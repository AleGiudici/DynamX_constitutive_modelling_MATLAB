function [outputMat] = simulateViscoelasticBehaviour(timeV,lambdatV,lambdazV,modelFormulation,modelParameters,...
    reference_configuration,activeActivator,trigger)

    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2

    % Importing model formulation
    elasticModel = modelFormulation.elasticModel;
    viscousModel = modelFormulation.viscousModel;
    
    % Importing model parameters
    parameterV_elastic = modelParameters.staticSEFparameters.("Final QLV values");
    if(isfield(modelFormulation,'activeModel'))
        activeModel = modelFormulation.activeModel;
        if(activeModel.activator == 1)
            parameterV_active = modelParameters.activeSEFparameters.("Final values");
        end
    else
        activeModel.activator = 0;
    end

    parameterV_viscous = modelParameters.viscousParameters.("Values");
    
    % Reference vessel geometry
    refMidWallRadius = (reference_configuration.(1)+reference_configuration.(2))/2;
    refOuterRadius = reference_configuration.(2);
    refInnerRadius = reference_configuration.(1);
    
    lambdarV = 1./lambdatV./lambdazV;
    
    rmV = refMidWallRadius*lambdatV;
    riV = sqrt(rmV.^2-(refMidWallRadius^2-refInnerRadius^2)./lambdazV);
    roV = sqrt(rmV.^2+(refOuterRadius^2-refMidWallRadius^2)./lambdazV);
    hV = roV-riV;
    
    tV_inverted=ones(1,length(timeV))*max(timeV)-timeV';
    tM_temp=ones(length(timeV),length(timeV)).*tV_inverted-tV_inverted';
    tM = (tM_temp(2:end,1:end-1)+tM_temp(1:end-1,1:end-1))/2;
    locator=tM<0;
    tM(locator) = inf;
            
    sigmatt_ve_compV=zeros(length(timeV),1);
    sigmazz_ve_compV=zeros(length(timeV),1);
    sigmarr_ve_compV=zeros(length(timeV),1);
        
    for iConstituent=1:elasticModel.nConstituents
        activator=zeros(1,2);
        activator(iConstituent)=1;
    
        % Elastic 2ndPK extra stress in the three principal directions
        [sigmatt_rr_eV,sigmazz_rr_eV,~,~,~]=elasticModel.fun(parameterV_elastic,lambdatV,lambdazV,lambdarV,activator); 
        [P_eV,~] = stress2pressure_force(sigmatt_rr_eV,sigmazz_rr_eV,rmV,hV);
        
        sigmatt_eV = sigmatt_rr_eV-P_eV/2*mtN; % Lagrange multiplier to get actual circ. stress
        sigmazz_eV = sigmazz_rr_eV-P_eV/2*mtN; % Lagrange multiplier to get actual axial stress
        sigmarr_eV = -P_eV/2*mtN; % Lagrange multiplier to get actual radial stress
    
        Stt_eV = sigmatt_eV./lambdatV.^2;
        Srr_eV = sigmarr_eV./lambdarV.^2;
        Szz_eV = sigmazz_eV./lambdazV.^2;
        
        delta_Stt_eV = Stt_eV(2:end)-Stt_eV(1:end-1); % Circ. elastic 2ndPK extra stress increments
        delta_Srr_eV = Srr_eV(2:end)-Srr_eV(1:end-1); % Radial elastic 2ndPK extra stress increments
        delta_Szz_eV = Szz_eV(2:end)-Szz_eV(1:end-1); % Axial elastic 2ndPK extra stress increments
        
        % Relaxation function calculated at each instant of the time interval matrix
        G_fun = viscousModel.fun(parameterV_viscous((iConstituent-1)*viscousModel.nParam+1:(iConstituent-1)*viscousModel.nParam+3),tM,0);
        locator = isnan(G_fun);
        G_fun(locator) = 0;
        
        % Viscoelastic 2ndPK circumferential extra stress
        Stt_veV = G_fun*delta_Stt_eV+Stt_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
        
        % Viscoelastic 2ndPK radial extra stress
        Srr_veV = G_fun*delta_Srr_eV+Srr_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
    
        % Viscoelastic 2ndPK axial extra stress
        Szz_veV = G_fun*delta_Szz_eV+Szz_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
        
        % Viscoelastic Cauchy circumferential extra stress
        sigmatt_ve_compV(1) = sigmatt_ve_compV(1)+Stt_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdatV(1)^2;
        sigmarr_ve_compV(1) = sigmarr_ve_compV(1)+Srr_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdarV(1)^2;
        sigmazz_ve_compV(1) = sigmazz_ve_compV(1)+Szz_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdazV(1)^2;
    
        sigmatt_ve_compV(2:end)=sigmatt_ve_compV(2:end)+Stt_veV.*lambdatV(2:end).^2; % Circumferential
        sigmarr_ve_compV(2:end)=sigmarr_ve_compV(2:end)+Srr_veV.*lambdarV(2:end).^2; % Radial
        sigmazz_ve_compV(2:end)=sigmazz_ve_compV(2:end)+Szz_veV.*lambdazV(2:end).^2; % Axial

        if(iConstituent == 1)
            sigmatt_elastin = sigmatt_ve_compV;
            sigmazz_elastin = sigmazz_ve_compV;
        elseif(iConstituent == 2)
            sigmatt_collagen = sigmatt_ve_compV-sigmatt_elastin;
            sigmazz_collagen = sigmazz_ve_compV-sigmazz_elastin;
        end
    end
    
    if(activeModel.activator ==1 && activeActivator == 1)
        [sigmatt_rr_eV,sigmazz_rr_eV,~,~,~]=activeModel.fun(parameterV_active,lambdatV,lambdazV,lambdarV); 
        [P_eV,~] = stress2pressure_force(sigmatt_rr_eV,sigmazz_rr_eV,rmV,hV);
    
        sigmatt_eV = sigmatt_rr_eV-P_eV/2*mtN; % Lagrange multiplier to get actual circ. stress
        sigmazz_eV = sigmazz_rr_eV-P_eV/2*mtN; % Lagrange multiplier to get actual axial stress
        sigmarr_eV = -P_eV/2*mtN; % Lagrange multiplier to get actual radial stress
    
        Stt_eV = sigmatt_eV./lambdatV.^2;
        Srr_eV = sigmarr_eV./lambdarV.^2;
        Szz_eV = sigmazz_eV./lambdazV.^2;
        
        delta_Stt_eV = Stt_eV(2:end)-Stt_eV(1:end-1); % Circ. elastic 2ndPK extra stress increments
        delta_Srr_eV = Srr_eV(2:end)-Srr_eV(1:end-1); % Radial elastic 2ndPK extra stress increments
        delta_Szz_eV = Szz_eV(2:end)-Szz_eV(1:end-1); % Axial elastic 2ndPK extra stress increments
        
        % Relaxation function calculated at each instant of the time interval matrix
        G_fun = viscousModel.fun(parameterV_viscous((iConstituent)*viscousModel.nParam+1:(iConstituent)*viscousModel.nParam+3),tM,0);
        locator = isnan(G_fun);
        G_fun(locator) = 0;
        
        % Viscoelastic 2ndPK circumferential extra stress
        Stt_veV = G_fun*delta_Stt_eV+Stt_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)));
        
        % Viscoelastic 2ndPK radial extra stress
        Srr_veV = G_fun*delta_Srr_eV+Srr_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)));
    
        % Viscoelastic 2ndPK axial extra stress
        Szz_veV = G_fun*delta_Szz_eV+Szz_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)));
        
        % Viscoelastic Cauchy circumferential extra stress
        sigmatt_ve_compV(1) = sigmatt_ve_compV(1)+Stt_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)))*lambdatV(1)^2;
        sigmarr_ve_compV(1) = sigmarr_ve_compV(1)+Srr_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)))*lambdarV(1)^2;
        sigmazz_ve_compV(1) = sigmazz_ve_compV(1)+Szz_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2)))*lambdazV(1)^2;
    
        sigmatt_ve_compV(2:end)=sigmatt_ve_compV(2:end)+Stt_veV.*lambdatV(2:end).^2; % Circumferential
        sigmarr_ve_compV(2:end)=sigmarr_ve_compV(2:end)+Srr_veV.*lambdarV(2:end).^2; % Radial
        sigmazz_ve_compV(2:end)=sigmazz_ve_compV(2:end)+Szz_veV.*lambdazV(2:end).^2; % Axial

        sigmatt_VSMC = sigmatt_ve_compV-sigmatt_elastin-sigmatt_collagen;
        sigmazz_VSMC = sigmazz_ve_compV-sigmazz_elastin-sigmazz_collagen;
    else
        sigmatt_VSMC = zeros(size(sigmatt_elastin,1),size(sigmatt_elastin,2));
        sigmazz_VSMC = zeros(size(sigmazz_elastin,1),size(sigmazz_elastin,2));
    end
        
    % Difference between viscoelastic Cauchy stress in the circ. and radial direction
    sigmatt_rr_ve_compV = sigmatt_ve_compV-sigmarr_ve_compV;
    
    % Difference between viscoelastic Cauchy stress in the axial and radial direction
    sigmazz_rr_ve_compV = sigmazz_ve_compV-sigmarr_ve_compV;
    
    [PV,fTV] = stress2pressure_force(sigmatt_rr_ve_compV,sigmazz_rr_ve_compV,rmV,hV);

    locator = find(trigger == 1);
       
    variableNames = ["Time [s]", "Pressure [mmHg]", "Inner diameter [mm]", "Outer diameter [mm]", "Transducer force [g]",...
                        "Circ. stretch [-]", "Circ. stress [MPa]", "Axial stretch [-]","Axial stress [MPa]","Elastin circ. load bearing [%]",...
                        "Elastin axial load bearing [%]","Collagen circ. load bearing [%]","Collagen axial load bearing [%]",...
                        "VSMC circ. load bearing [%]", "VSMC axial load bearing [%]"];
    
    outputMat = array2table([timeV(locator:end)-timeV(locator-1) PV(locator:end), 2*riV(locator:end), 2*roV(locator:end), fTV(locator:end), lambdatV(locator:end),...
        sigmatt_ve_compV(locator:end), lambdazV(locator:end), sigmazz_ve_compV(locator:end), sigmatt_elastin(locator:end)./sigmatt_ve_compV(locator:end)*100,...
        sigmazz_elastin(locator:end)./sigmazz_ve_compV(locator:end)*100,sigmatt_collagen(locator:end)./sigmatt_ve_compV(locator:end)*100,...
        sigmazz_collagen(locator:end)./sigmazz_ve_compV(locator:end)*100,sigmatt_VSMC(locator:end)./sigmatt_ve_compV(locator:end)*100,...
        sigmazz_VSMC(locator:end)./sigmazz_ve_compV(locator:end)*100],'VariableNames',variableNames);
end