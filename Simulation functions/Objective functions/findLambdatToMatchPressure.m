function [error] = findLambdatToMatchPressure(lambdat,lambdatV,lambdazV,timeV,modelFormulation,modelParameters,...
    targetPressure,referenceConfiguration,iTime,deltaTime,activeActivator)
% [error] = findLambdatToMatchPressure(lambdat,lambdatV,lambdazV,timeV,...
% modelFormulation,modelParameters,targetPressure,referenceConfiguration,...
% iTime,deltaTime,activeActivator)
% 
% This cost function is used to iteratively estimate the circumferential
% stretch ("lambdat") that, given the artery viscoelastic model (which
% accounts for the summed contribution of elastin, collagen and, possibly,
% VSMC), the circumferential ("lambdatV) and axial stretch history
% ("lambdazV") and the time vector ("timeV"), yields the desired target 
% pressure ("targetPressure").
%
% The model formulation ("modelFormulation" with substructures 
% "elasticModel", "viscousModel", and "activeModel"),
% the model parameters ("modelParameters"), the reference geometry of the 
% vessel ("referenceConfiguration") are also required inputs. 
% Note that parameterV_active may be unused if the field "activator" in
% "activeModel" is set to 0.
%
% "iTime" is the current time stamp, "deltaTime" is the sampling interval,
% and "activeActivator" determines whether the VSMC part of the model is
% activated or not.

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

% Reference configuration geometry
refMidWallRadius = (referenceConfiguration.("Inner radius [mm]")+referenceConfiguration.("Outer radius [mm]"))/2;
refThickness = referenceConfiguration.("Outer radius [mm]")-referenceConfiguration.("Inner radius [mm]");

% lambdattV = [lambdattV,lambdatt];
lambdatV(iTime+1) = lambdat;
lambdarV = 1./lambdatV./lambdazV;

% tV = [tV, 500+delta_t*(i-1)];
% tV(i+1) = 500+delta_t*(i-1);
timeV(iTime+1) = timeV(iTime)+deltaTime;

tM = timeV(end)-(timeV(2:end)+timeV(1:end-1))/2;
 
sigmatt_ve_compV = 0;
sigmarr_ve_compV = 0;

rmV = lambdatV*refMidWallRadius;

riV = sqrt(rmV.^2-(refMidWallRadius^2-(refMidWallRadius-refThickness/2)^2)./lambdazV);
roV = sqrt(rmV.^2+((refMidWallRadius+refThickness/2)^2-refMidWallRadius^2)./lambdazV);

hV = roV-riV;

for iConstituent=1:elasticModel.nConstituents
    activator=zeros(1,2);
    activator(iConstituent)=1;
    
    % First coarse time grid over the entire experiment
    [sigmatt_rr_eV,sigmazz_rr_eV,~,~,~]=elasticModel.fun(parameterV_elastic,lambdatV,lambdazV,lambdarV,activator); 
    [P_eV,~] = stress2pressure_force(sigmatt_rr_eV,sigmazz_rr_eV,rmV,hV);
        
    sigmatt_eV = sigmatt_rr_eV-P_eV/2*mtN; % Lagrange multiplier to get actual circ. stress
    sigmarr_eV = -P_eV/2*mtN; % Lagrange multiplier to get actual radial stress

    Stt_eV = sigmatt_eV./lambdatV.^2;
    Srr_eV = sigmarr_eV./lambdarV.^2;
    
    delta_Stt_eV = Stt_eV(2:end)-Stt_eV(1:end-1); % Circ. elastic 2ndPK extra stress increments
    delta_Srr_eV = Srr_eV(2:end)-Srr_eV(1:end-1); % Radial elastic 2ndPK extra stress increments
    
    G_fun = viscousModel.fun(parameterV_viscous((iConstituent-1)*viscousModel.nParam+1:(iConstituent-1)*viscousModel.nParam+3),tM,0);
    locator = isnan(G_fun);
    G_fun(locator) = 0;
    
    Stt_veV = sum(G_fun.*delta_Stt_eV)+Stt_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2))); % Circumferential 2nd Piola-Kirchhoff viscous stress
    Srr_veV = sum(G_fun.*delta_Srr_eV)+Srr_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2))); % Circumferential 2nd Piola-Kirchhoff viscous stress

    sigmatt_ve_compV=sigmatt_ve_compV+Stt_veV*lambdatV(end)^2; % Circumferential Cauchy viscous stress
    sigmarr_ve_compV=sigmarr_ve_compV+Srr_veV*lambdarV(end)^2; % Circumferential Cauchy viscous stress
end

if(activeModel.activator == 1 && activeActivator == 1)
    [sigmatt_rr_eV,sigmazz_rr_eV,~,~,~]=activeModel.fun(parameterV_active,lambdatV,lambdazV,lambdarV); 
    [P_eV,~] = stress2pressure_force(sigmatt_rr_eV,sigmazz_rr_eV,rmV,hV);

    sigmatt_eV = sigmatt_rr_eV-P_eV/2*mtN; % Lagrange multiplier to get actual circ. stress
    sigmarr_eV = -P_eV/2*mtN; % Lagrange multiplier to get actual radial stress

    Stt_eV = sigmatt_eV./lambdatV.^2;
    Srr_eV = sigmarr_eV./lambdarV.^2;
    
    delta_Stt_eV = Stt_eV(2:end)-Stt_eV(1:end-1); % Circ. elastic 2ndPK extra stress increments
    delta_Srr_eV = Srr_eV(2:end)-Srr_eV(1:end-1); % Radial elastic 2ndPK extra stress increments
    
    G_fun = viscousModel.fun(parameterV_viscous((iConstituent)*viscousModel.nParam+1:(iConstituent)*viscousModel.nParam+3),tM,0);
    locator = isnan(G_fun);
    G_fun(locator) = 0;
    
    Stt_veV = sum(G_fun.*delta_Stt_eV)+Stt_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2))); % Circumferential 2nd Piola-Kirchhoff viscous stress
    Srr_veV = sum(G_fun.*delta_Srr_eV)+Srr_eV(1)/(1+parameterV_viscous((iConstituent)*3+1)*log(parameterV_viscous((iConstituent)*3+3)/parameterV_viscous((iConstituent)*3+2))); % Circumferential 2nd Piola-Kirchhoff viscous stress

    sigmatt_ve_compV=sigmatt_ve_compV+Stt_veV*lambdatV(end)^2; % Circumferential Cauchy viscous stress
    sigmarr_ve_compV=sigmarr_ve_compV+Srr_veV*lambdarV(end)^2; % Circumferential Cauchy viscous stress
end

sigmatt_rr_ve_compV = sigmatt_ve_compV-sigmarr_ve_compV;

PV_compV = (sigmatt_rr_ve_compV.*hV(end)./rmV(end))/mtN;

error = PV_compV-targetPressure;
