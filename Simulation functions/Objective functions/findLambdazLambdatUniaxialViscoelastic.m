function [error] = findLambdazLambdatUniaxialViscoelastic(newLambda,lambdatV_past,lambdazV_past,timeV,...
    parameterV_elastic,parameterV_viscous,parameterV_active,elasticModel,viscousModel,activeModel,...
    targetPressure,refWallThickness,refMidWallRadius,activeActivator)
% [error] = findLambdazLambdatViscoelastic(newLambda,lambdatV_past,...
% lambdazV_past,timeV,parameterV_elastic,parameterV_viscous,parameterV_active,...
% elasticModel,viscousModel,activeModel,targetPressure,refWallThickness,...
% refMidWallRadius,activeActivator)
%
% This function is the cost function used for the iterative estimation of
% the circumferential and axial stretches in a simulated uniaxial tensile
% test in the circumferential direction, given the artery viscoelastic model 
% (which accounts for the summed contribution of elastin, collagen and, 
% possibly, VSMC). The VSMC contribution may be switched off by setting
% activatorActive to 0.
%
% The cost is defined in terms of the difference between the modelled
% circumferential stress and the pressure-equivalent stress (given the
% target pressure "targetPressure") and the difference between the modelled
% axial stress and 0 (given that the axial stress must be null in a
% uniaxial tensile test in the circumferential direction).
%
% The model formulation ("elasticModel", "viscousModel" "activeModel"), 
% the model parameters ("parameterV_elastic", "parameterV_viscous and
% "parameterV_active"), the reference geometry of the vessel ("refthickness" 
% and "refMidWallRadius") are required inputs. 

    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    
    lambdatV = [lambdatV_past; newLambda(1)]; % adding guessed stretches to
    lambdazV = [lambdazV_past; newLambda(2)]; % deformation history.
    
    lambdarV = 1./lambdatV./lambdazV;
    
    tV_inverted = ones(1,length(timeV))*max(timeV)-timeV;
    tM_temp = ones(length(timeV),length(timeV)).*tV_inverted-tV_inverted';
    tM = (tM_temp(2:end,1:end-1)+tM_temp(1:end-1,1:end-1))/2;
    locator = tM<0;
    tM(locator) = inf;
            
    sigmatt_veV=zeros(length(timeV),1);
    sigmazz_veV=zeros(length(timeV),1);
        
    for iConstituent=1:elasticModel.nConstituents
        activator=zeros(1,2);
        activator(iConstituent)=1;
    
        % Elastic 2ndPK extra stress in the three principal directions
        [sigmatt_eV,sigmazz_eV,~,~,~]=elasticModel.fun(parameterV_elastic,lambdatV,lambdazV,lambdarV,activator); 
        
        % Convert Cauchy stresses to 2ndPK stresses
        Stt_eV = sigmatt_eV./lambdatV.^2; 
        Szz_eV = sigmazz_eV./lambdazV.^2; 
        
        delta_Stt_eV = Stt_eV(2:end)-Stt_eV(1:end-1); % Circ. elastic 2ndPK extra stress increments
        delta_Szz_eV = Szz_eV(2:end)-Szz_eV(1:end-1); % Axial elastic 2ndPK extra stress increments
        
        % Relaxation function calculated at each instant of the time interval matrix
        G_fun = viscousModel.fun(parameterV_viscous((iConstituent-1)*viscousModel.nParam+1:(iConstituent-1)*viscousModel.nParam+3),tM,0);
        locator = isnan(G_fun);
        G_fun(locator) = 0;
        
        % Viscoelastic 2ndPK stresses
        Stt_veV = G_fun*delta_Stt_eV+Stt_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
        Szz_veV = G_fun*delta_Szz_eV+Szz_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
        
        % Viscoelastic Cauchy circumferential extra stress
        sigmatt_veV(1) = sigmatt_veV(1)+Stt_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdatV(1)^2;
        sigmazz_veV(1) = sigmazz_veV(1)+Szz_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdazV(1)^2;
    
        sigmatt_veV(2:end)=sigmatt_veV(2:end)+Stt_veV.*lambdatV(2:end).^2; % Circumferential
        sigmazz_veV(2:end)=sigmazz_veV(2:end)+Szz_veV.*lambdazV(2:end).^2; % Axial
    end

    if(activeModel.activator ==1 && activeActivator == 1)
        [sigmatt_eV,sigmazz_eV,~,~,~]=activeModel.fun(parameterV_active,lambdatV,lambdazV,lambdarV); 
       
        Stt_eV = sigmatt_eV./lambdatV.^2;
        Szz_eV = sigmazz_eV./lambdazV.^2;
        
        delta_Stt_eV = Stt_eV(2:end)-Stt_eV(1:end-1); % Circ. elastic 2ndPK extra stress increments
        delta_Szz_eV = Szz_eV(2:end)-Szz_eV(1:end-1); % Radial elastic 2ndPK extra stress increments
        
        G_fun = viscousModel.fun(parameterV_viscous((iConstituent)*viscousModel.nParam+1:(iConstituent)*viscousModel.nParam+3),tM,0);
        locator = isnan(G_fun);
        G_fun(locator) = 0;
        
        % Stt_veV = sum(G_fun*delta_Stt_eV)+Stt_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2))); % Circumferential 2nd Piola-Kirchhoff viscous stress
        % Szz_veV = sum(G_fun*delta_Szz_eV)+Szz_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2))); % Circumferential 2nd Piola-Kirchhoff viscous stress
        
        Stt_veV = G_fun*delta_Stt_eV+Stt_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)));
        Szz_veV = G_fun*delta_Szz_eV+Szz_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)));
        
        % sigmatt_veV=sigmatt_veV+Stt_veV*lambdatV(end)^2; % Circumferential Cauchy viscous stress
        % sigmazz_veV=sigmazz_veV+Szz_veV*lambdazV(end)^2; % Circumferential Cauchy viscous stress

        sigmatt_veV(1) = sigmatt_veV(1)+Stt_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)))*lambdatV(1)^2;
        sigmazz_veV(1) = sigmazz_veV(1)+Szz_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)))*lambdazV(1)^2;

        sigmatt_veV(2:end)=sigmatt_veV(2:end)+Stt_veV.*lambdatV(2:end).^2; % Circumferential
        sigmazz_veV(2:end)=sigmazz_veV(2:end)+Szz_veV.*lambdazV(2:end).^2; % Axial
    end
    
    % Vessel deformed geometry
    rmV = refMidWallRadius*lambdatV; % mid-wall radius
    riV = sqrt(rmV.^2-(refMidWallRadius^2-(refMidWallRadius-refWallThickness/2)^2)./lambdazV); % inner radius
    roV = sqrt(rmV.^2+((refMidWallRadius+refWallThickness/2)^2-refMidWallRadius^2)./lambdazV); % outer radius
    hV = roV-riV; % wall thickness
    
    % Convert target pressure into stress
    sigmatt_target = targetPressure*mtN*(rmV(end))/hV(end);
   
    error = [sigmatt_target-sigmatt_veV(end), -sigmazz_veV(end)];
end