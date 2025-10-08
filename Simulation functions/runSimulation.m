function [exportMat] = runSimulation(fit_data,targetPressureV,timeV,app,triggerV)
% [exportMat] = runSimulation(fit_data,targetPressureV,timeV,app,triggerV)
% 
% This function is used to run the iterative steps that yield the simulated
% curves/experiments ("exportMat". These simulations are run by iteratively
% estimating the circumferential deformation that gives the desired pressure
%
% The required input parameters are the modelling structure "fit_data", the 
% target pressure vector "targetPressureV", the time vector ("timeV"), the
% app structure ("app"), and the trigger vector ("triggerV"). The trigger
% vector is used to indicate at what time instant the actual simulation
% starts (i.e., removing the time instants that are required to obtain a
% steady state solution, when triggerV = 1).

%% Importing reference and unloaded configurations
    ri_ref = fit_data.vesselConfigurations.reference_configuration.("Inner radius [mm]");
    ro_ref = fit_data.vesselConfigurations.reference_configuration.("Outer radius [mm]");
    l_ref = fit_data.vesselConfigurations.reference_configuration.("Axial length [mm]");
    Ri = fit_data.vesselConfigurations.unloaded_configuration.("Inner radius [mm]");
    Ro = fit_data.vesselConfigurations.unloaded_configuration.("Outer radius [mm]");
    L = fit_data.vesselConfigurations.unloaded_configuration.("Axial length [mm]");
    
    lambdaz_to_use = str2double(app.Axial_stretch_sim.Value)/100*fit_data.in_vivo_axial_stretch;

%% Activating VSMC according to app selection
    if(strcmp(app.VSMCcontractionSwitch.Value,'On') && fit_data.modelFormulation.activeModel.activator == 1)
        activeActivator = 1;
    else
        activeActivator = 0;
    end    
    
%% Iterative estimation of the circumferential deformation vector
    timeV = [0, 500+timeV];
    lambdazV = zeros(1,length(timeV));
    lambdatV = zeros(1,length(timeV));
    
    lambdazV(1) = L/l_ref;
    lambdatV(1) = (Ro+Ri)/(ri_ref+ro_ref);
    
    for iTime = 2:length(timeV)
        lambdazV(iTime) = lambdaz_to_use;
    
        options = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 8000,'OptimalityTolerance',10^-10,'FunctionTolerance',10^-10); % setting max. no. of iterations to 8000
            
        lambdat_new = lsqnonlin(@(lambdat_new) findLambdatToMatchPressure(lambdat_new,lambdatV(1:iTime),lambdazV(1:iTime),timeV(1:iTime),fit_data.modelFormulation,...
            fit_data.modelParameters,targetPressureV(iTime-1),fit_data.vesselConfigurations.reference_configuration,iTime-1,timeV(iTime)-timeV(iTime-1),activeActivator),1,0,10,options);
        lambdatV(iTime) = lambdat_new;
    end

%% Output data
    exportMat = simulateViscoelasticBehaviour(timeV',lambdatV',lambdazV',fit_data.modelFormulation,fit_data.modelParameters,...
        fit_data.vesselConfigurations.reference_configuration,activeActivator,triggerV);
end