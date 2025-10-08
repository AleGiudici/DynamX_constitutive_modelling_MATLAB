function fit_data = pressureDrivenFittedCurves(fit_data)
% This function re-exports all the fitting results so that the modelled
% pressure matches the experimental pressures. Because the model works with
% deformation as input, this new curves are found interatively (i.e., for
% each input pressure, the code looks for the circumferential deformation
% that yields that pressure given the axial deformation as input.
% The function also compute the dynamic-to-quasi-static stiffening ratio.

    %% Retriving reference and unloaded configurations
    ri_ref = fit_data.vesselConfigurations.reference_configuration.("Inner radius [mm]");
    ro_ref = fit_data.vesselConfigurations.reference_configuration.("Outer radius [mm]");
    rm_ref = (ri_ref+ro_ref)/2;
    h_ref = fit_data.vesselConfigurations.reference_configuration.("Thickness [mm]");
    l_ref = fit_data.vesselConfigurations.reference_configuration.("Axial length [mm]");
    
    Ri = fit_data.vesselConfigurations.unloaded_configuration.("Inner radius [mm]");
    Ro = fit_data.vesselConfigurations.unloaded_configuration.("Outer radius [mm]");
    H = fit_data.vesselConfigurations.unloaded_configuration.("Thickness [mm]");
    L = fit_data.vesselConfigurations.unloaded_configuration.("Axial length [mm]");
    
    if(fit_data.modelFormulation.activeModel.activator == 1)
        ri_ref_active = fit_data.vesselConfigurations.reference_configuration_active.("Inner radius [mm]");
        ro_ref_active = fit_data.vesselConfigurations.reference_configuration_active.("Outer radius [mm]");
        rm_ref_active = (ri_ref_active+ro_ref_active)/2;
        h_ref_active = fit_data.vesselConfigurations.reference_configuration_active.("Thickness [mm]");
        l_ref_active = fit_data.vesselConfigurations.reference_configuration_active.("Axial length [mm]");
    end
    
    %% Retriving Model
    % elasticModel = fit_data.modelFormulation.elasticModel;
    % parameterV_elastic = fit_data.modelParameters.staticSEFparameters.("Final QLV values");
    % 
    % viscousModel = fit_data.modelFormulation.viscousModel;
    % parameterV_viscous = fit_data.modelParameters.viscousParameters.Values;
    
    activeModel = fit_data.modelFormulation.activeModel;
    % parameterV_active = fit_data.modelParameters.activeSEFparameters.("Final values");

    %% Quasi-static passive pressure sweeps
    PD = fit_data.PD_static;
        
    for i = 1:length(PD)
        lambdaz = fit_data.fittedCurves.quasiStaticFit.(PD{i}).('Axial stretch [-]');
        lambdazV = [L/l_ref; lambdaz];
        lambdatV = zeros(length(lambdazV),1);
        lambdatV(1) = (Ro+Ri)/(ri_ref+ro_ref);
    
        time = fit_data.fittedCurves.quasiStaticFit.(PD{i}).('Time [s]');
        timeV = zeros(length(lambdazV),1);
        timeV(1) = 0;
        timeV(2:end) = time;
    
        targetPressure = fit_data.fittedCurves.quasiStaticFit.(PD{i}).('Exp. pressure [mmHg]');
    
        for j = 1:length(targetPressure)
            lambdatt_new = lsqnonlin(@(lambdat) findLambdatToMatchPressure(lambdat,lambdatV(1:j+1),lambdazV(1:j+1),timeV(1:j+1),...
                fit_data.modelFormulation,fit_data.modelParameters,targetPressure(j),fit_data.vesselConfigurations.reference_configuration,...
                j,timeV(j+1)-timeV(j),0),1);
            
            lambdatV(j+1) = lambdatt_new;
        end
        
        trigger = zeros(1,length(timeV));
        trigger(2) = 1;
    
        pressureDrivenFitting.passive.(PD{i}) = simulateViscoelasticBehaviour(timeV,lambdatV,lambdazV,fit_data.modelFormulation,...
            fit_data.modelParameters,fit_data.vesselConfigurations.reference_configuration,0,trigger);
    end
    
    %% Quasi-static passive force sweeps
    FL = fit_data.FL;
        
    for i = 1:length(FL)
        lambdaz = fit_data.fittedCurves.quasiStaticFit.(FL{i}).('Axial stretch [-]');
        lambdazV = [L/l_ref; lambdaz];
        lambdatV = zeros(length(lambdazV),1);
        lambdatV(1) = (Ro+Ri)/(ri_ref+ro_ref);
    
        time = fit_data.fittedCurves.quasiStaticFit.(FL{i}).('Time [s]');
        timeV = zeros(length(lambdazV),1);
        timeV(1) = 0;
        timeV(2:end) = time;
    
        targetPressure = fit_data.fittedCurves.quasiStaticFit.(FL{i}).('Exp. pressure [mmHg]');
    
        for j = 1:length(targetPressure)
            lambdatt_new = lsqnonlin(@(lambdat) findLambdatToMatchPressure(lambdat,lambdatV(1:j+1),lambdazV(1:j+1),timeV(1:j+1),...
                fit_data.modelFormulation,fit_data.modelParameters,targetPressure(j),fit_data.vesselConfigurations.reference_configuration,...
                j,timeV(j+1)-timeV(j),0),0.5);
    
            lambdatV(j+1) = lambdatt_new;
        end
    

        

        trigger = zeros(1,length(timeV));
        trigger(2) = 1;
    
        pressureDrivenFitting.passive.(FL{i}) = simulateViscoelasticBehaviour(timeV,lambdatV,lambdazV,fit_data.modelFormulation,...
            fit_data.modelParameters,fit_data.vesselConfigurations.reference_configuration,0,trigger);
    end
    
    %% Dynamic passive pressure loops
    PD = fit_data.PD_dynamic;
        
    for i = 1:length(PD)
        lambdaz = fit_data.fittedCurves.dynamicFit.(PD{i}).('Axial stretch [-]');
        lambdazV = [L/l_ref; lambdaz; lambdaz; lambdaz; lambdaz];
        lambdatV = zeros(length(lambdazV),1);
        lambdatV(1) = (Ro+Ri)/(ri_ref+ro_ref);
    
        time = fit_data.fittedCurves.dynamicFit.(PD{i}).('Time [s]');
        time = time-2*time(1)+time(2);
        timeV = zeros(length(lambdazV),1);
        timeV(1) = 0;
        timeV(2:1+length(time)) = 500+time;
        timeV(2+length(time):1+2*length(time)) = 500+time+time(end);
        timeV(2+2*length(time):1+3*length(time)) = 500+time+2*time(end);
        timeV(2+3*length(time):1+4*length(time)) = 500+time+3*time(end);
    
        templatePressure = fit_data.fittedCurves.dynamicFit.(PD{i}).('Exp. pressure [mmHg]');
        targetPressure = [templatePressure; templatePressure; templatePressure; templatePressure];
    
        for j = 1:length(targetPressure)
            lambdatt_new = lsqnonlin(@(lambdat) findLambdatToMatchPressure(lambdat,lambdatV(1:j+1),lambdazV(1:j+1),timeV(1:j+1),...
                fit_data.modelFormulation,fit_data.modelParameters,targetPressure(j),fit_data.vesselConfigurations.reference_configuration,...
                j,timeV(j+1)-timeV(j),0),0.5);
            
            lambdatV(j+1) = lambdatt_new;
        end
   
        trigger = zeros(1,length(timeV));
        trigger(end-length(time)) = 1;
    
        pressureDrivenFitting.passive.(PD{i}) = simulateViscoelasticBehaviour(timeV,lambdatV,lambdazV,fit_data.modelFormulation,...
            fit_data.modelParameters,fit_data.vesselConfigurations.reference_configuration,0,trigger);
    end
    
    %% Stiffness ratio passive tests
    ref_P = pressureDrivenFitting.passive.PD_100.("Pressure [mmHg]");
    ref_lambdat = pressureDrivenFitting.passive.PD_100.("Circ. stretch [-]");
    ref_sigmatt = pressureDrivenFitting.passive.PD_100.("Circ. stress [MPa]");
    
    [~,locator_max] = max(ref_P);
    ref_P_loading = ref_P(1:locator_max);
    ref_P_unloading = ref_P(locator_max:end);
    
    ref_lambdat_loading = ref_lambdat(1:locator_max);
    ref_lambdat_unloading = ref_lambdat(locator_max:end);
    
    ref_sigmatt_loading = ref_sigmatt(1:locator_max);
    ref_sigmatt_unloading = ref_sigmatt(locator_max:end);

    ratio_mod_passive = zeros(3,size(fit_data.PD_dynamic,2)/3);
    ratio_exp_passive = zeros(3,size(fit_data.PD_dynamic,2)/3);
    
    for iTest = 1:length(PD) 
        P_loop = pressureDrivenFitting.passive.(PD{iTest}).("Pressure [mmHg]");
        lambdat_loop = pressureDrivenFitting.passive.(PD{iTest}).("Circ. stretch [-]");
        sigmatt_loop = pressureDrivenFitting.passive.(PD{iTest}).("Circ. stress [MPa]");
    
        [~,locator_max] = max(P_loop);
        P_loop_loading = P_loop(1:locator_max);
        P_loop_unloading = P_loop(locator_max:end);
    
        lambdat_loop_loading = lambdat_loop(1:locator_max);
        lambdat_loop_unloading = lambdat_loop(locator_max:end);
    
        sigmatt_loop_loading = sigmatt_loop(1:locator_max);
        sigmatt_loop_unloading = sigmatt_loop(locator_max:end);
    
        regression_line = polyfit(lambdat_loop_loading,sigmatt_loop_loading,1);
        dynamic_stiffness_loading = regression_line(1);
        regression_line = polyfit(lambdat_loop_unloading,sigmatt_loop_unloading,1);
        dynamic_stiffness_unloading = regression_line(1);
    
        sigmatt_qs_loading = interp1(ref_P_loading,ref_sigmatt_loading,P_loop_loading,"spline","extrap");
        lambdat_qs_loading = interp1(ref_P_loading,ref_lambdat_loading,P_loop_loading,"spline","extrap");
        regression_line = polyfit(lambdat_qs_loading,sigmatt_qs_loading,1);
        qs_stiffness_loading = regression_line(1);
    
        sigmatt_qs_unloading = interp1(ref_P_unloading,ref_sigmatt_unloading,P_loop_unloading,"spline","extrap");
        lambdat_qs_unloading = interp1(ref_P_unloading,ref_lambdat_unloading,P_loop_unloading,"spline","extrap");
        regression_line = polyfit(lambdat_qs_unloading,sigmatt_qs_unloading,1);
        qs_stiffness_unloading = regression_line(1);
        if(iTest<size(fit_data.PD_dynamic,2)/3+1)
            ratio_mod_passive(1,iTest) = (dynamic_stiffness_loading+dynamic_stiffness_unloading)/(qs_stiffness_loading+qs_stiffness_unloading);
        elseif(iTest>size(fit_data.PD_dynamic,2)/3 && iTest<2*size(fit_data.PD_dynamic,2)/3+1)
            ratio_mod_passive(2,iTest-size(fit_data.PD_dynamic,2)/3) = (dynamic_stiffness_loading+dynamic_stiffness_unloading)/(qs_stiffness_loading+qs_stiffness_unloading);
        else
            ratio_mod_passive(3,iTest-2*size(fit_data.PD_dynamic,2)/3) = (dynamic_stiffness_loading+dynamic_stiffness_unloading)/(qs_stiffness_loading+qs_stiffness_unloading);
        end
    end
    
    ref_P = fit_data.fittedCurves.quasiStaticFit.PD_100.("Exp. pressure [mmHg]");
    ref_lambdat = fit_data.fittedCurves.quasiStaticFit.PD_100.("Circ. stretch [-]");
    ref_sigmatt = fit_data.fittedCurves.quasiStaticFit.PD_100.("Exp. circ. stress [MPa]");
    
    [~,locator_max] = max(ref_P);
    ref_P_loading = ref_P(1:locator_max);
    ref_P_unloading = ref_P(locator_max:end);
    
    ref_lambdat_loading = ref_lambdat(1:locator_max);
    ref_lambdat_unloading = ref_lambdat(locator_max:end);
    
    ref_sigmatt_loading = ref_sigmatt(1:locator_max);
    ref_sigmatt_unloading = ref_sigmatt(locator_max:end);
    
    for iTest = 1:length(PD) 
        P_loop = fit_data.fittedCurves.dynamicFit.(PD{iTest}).("Exp. pressure [mmHg]");
        lambdat_loop = fit_data.fittedCurves.dynamicFit.(PD{iTest}).("Circ. stretch [-]");
        sigmatt_loop = fit_data.fittedCurves.dynamicFit.(PD{iTest}).("Exp. circ. stress [MPa]");
    
        [~,locator_max] = max(P_loop);
        P_loop_loading = P_loop(1:locator_max);
        P_loop_unloading = P_loop(locator_max:end);
    
        lambdat_loop_loading = lambdat_loop(1:locator_max);
        lambdat_loop_unloading = lambdat_loop(locator_max:end);
    
        sigmatt_loop_loading = sigmatt_loop(1:locator_max);
        sigmatt_loop_unloading = sigmatt_loop(locator_max:end);
    
        regression_line = polyfit(lambdat_loop_loading,sigmatt_loop_loading,1);
        dynamic_stiffness_loading = regression_line(1);
        regression_line = polyfit(lambdat_loop_unloading,sigmatt_loop_unloading,1);
        dynamic_stiffness_unloading = regression_line(1);
    
        sigmatt_qs_loading = interp1(ref_P_loading,ref_sigmatt_loading,P_loop_loading,"spline","extrap");
        lambdat_qs_loading = interp1(ref_P_loading,ref_lambdat_loading,P_loop_loading,"spline","extrap");
        regression_line = polyfit(lambdat_qs_loading,sigmatt_qs_loading,1);
        qs_stiffness_loading = regression_line(1);
    
        sigmatt_qs_unloading = interp1(ref_P_unloading,ref_sigmatt_unloading,P_loop_unloading,"spline","extrap");
        lambdat_qs_unloading = interp1(ref_P_unloading,ref_lambdat_unloading,P_loop_unloading,"spline","extrap");
        regression_line = polyfit(lambdat_qs_unloading,sigmatt_qs_unloading,1);
        qs_stiffness_unloading = regression_line(1);
    
        if(iTest<size(fit_data.PD_dynamic,2)/3+1)
            ratio_exp_passive(1,iTest) = (dynamic_stiffness_loading+dynamic_stiffness_unloading)/(qs_stiffness_loading+qs_stiffness_unloading);
        elseif(iTest>size(fit_data.PD_dynamic,2)/3 && iTest<2*size(fit_data.PD_dynamic,2)/3+1)
            ratio_exp_passive(2,iTest-size(fit_data.PD_dynamic,2)/3) = (dynamic_stiffness_loading+dynamic_stiffness_unloading)/(qs_stiffness_loading+qs_stiffness_unloading);
        else
            ratio_exp_passive(3,iTest-2*size(fit_data.PD_dynamic,2)/3) = (dynamic_stiffness_loading+dynamic_stiffness_unloading)/(qs_stiffness_loading+qs_stiffness_unloading);
        end
    end
        
    if(activeModel.activator == 1)
        %% Quasi-static active pressure sweep
        PD = fieldnames(fit_data.fittedCurves.activeFit.Lname);
        
        lambdaz = fit_data.fittedCurves.activeFit.Lname.(PD{1}).('Axial stretch [-]');
        lambdazV = [L/l_ref; lambdaz];
        lambdatV = zeros(length(lambdazV),1);
        lambdatV(1) = (Ro+Ri)/(ri_ref+ro_ref);
        
        time = fit_data.fittedCurves.activeFit.Lname.(PD{1}).('Time [s]');
        time = time-2*time(1)+time(2);
        timeV = zeros(length(lambdazV),1);
        timeV(1) = 0;
        timeV(2:1+length(time)) = 1000+time;
        
        targetPressure = fit_data.fittedCurves.activeFit.Lname.(PD{1}).('Exp. pressure [mmHg]');
        
        for j = 1:length(targetPressure)
            lambdatt_new = lsqnonlin(@(lambdat) findLambdatToMatchPressure(lambdat,lambdatV(1:j+1),lambdazV(1:j+1),timeV(1:j+1),...
                fit_data.modelFormulation,fit_data.modelParameters,targetPressure(j),fit_data.vesselConfigurations.reference_configuration,...
                j,timeV(j+1)-timeV(j),1),0.3);
            
            lambdatV(j+1) = lambdatt_new;
        end
        
        trigger = zeros(1,length(timeV));
        trigger(2) = 1;
        
        pressureDrivenFitting.active.(PD{1}) = simulateViscoelasticBehaviour(timeV,lambdatV,lambdazV,fit_data.modelFormulation,...
            fit_data.modelParameters,fit_data.vesselConfigurations.reference_configuration,1,trigger);
        
        %% Dynamic active pressure loops
        for i = 2:length(PD)
            lambdaz = fit_data.fittedCurves.activeFit.Lname.(PD{i}).('Axial stretch [-]');
            lambdazV = [L/l_ref; lambdaz; lambdaz; lambdaz; lambdaz; lambdaz(1)];
            lambdatV = zeros(length(lambdazV),1);
            lambdatV(1) = (Ro+Ri)/(ri_ref+ro_ref);
        
            time = fit_data.fittedCurves.activeFit.Lname.(PD{i}).('Time [s]');
            time = time-2*time(1)+time(2);
            timeV = zeros(length(lambdazV),1);
            timeV(1) = 0;
            timeV(2:1+length(time)) = 1000+time;
            timeV(2+length(time):1+2*length(time)) = 1000+time+time(end);
            timeV(2+2*length(time):1+3*length(time)) = 1000+time+2*time(end);
            timeV(2+3*length(time):1+4*length(time)) = 1000+time+3*time(end);
            timeV(2+4*length(time)) = 1000+time(1)+4*time(end);
        
            templatePressure = fit_data.fittedCurves.activeFit.Lname.(PD{i}).('Exp. pressure [mmHg]');
            targetPressure = [templatePressure; templatePressure; templatePressure; templatePressure; templatePressure(1)];
        
            for j = 1:length(targetPressure)
                lambdatt_new = lsqnonlin(@(lambdat) findLambdatToMatchPressure(lambdat,lambdatV(1:j+1),lambdazV(1:j+1),timeV(1:j+1),...
                    fit_data.modelFormulation,fit_data.modelParameters,targetPressure(j),fit_data.vesselConfigurations.reference_configuration,...
                    j,timeV(j+1)-timeV(j),1),0.5);
                
                lambdatV(j+1) = lambdatt_new;
            end
        
            trigger = zeros(1,length(timeV));
            trigger(end-length(time)) = 1;
        
            pressureDrivenFitting.active.(PD{i}) = simulateViscoelasticBehaviour(timeV,lambdatV,lambdazV,fit_data.modelFormulation,...
                fit_data.modelParameters,fit_data.vesselConfigurations.reference_configuration,1,trigger);
        end
                
        %% Stiffness ratio active tests
        ref_P = pressureDrivenFitting.active.PD_100.("Pressure [mmHg]");
        ref_lambdat = pressureDrivenFitting.active.PD_100.("Circ. stretch [-]");
        ref_sigmatt = pressureDrivenFitting.active.PD_100.("Circ. stress [MPa]");
        
        [~,locator_max] = max(ref_P);
        ref_P_loading = ref_P(1:locator_max);
        ref_P_unloading = ref_P(locator_max:end);
        
        ref_lambdat_loading = ref_lambdat(1:locator_max);
        ref_lambdat_unloading = ref_lambdat(locator_max:end);
        
        ref_sigmatt_loading = ref_sigmatt(1:locator_max);
        ref_sigmatt_unloading = ref_sigmatt(locator_max:end);
    
        ratio_mod_active = zeros(6,1);
        ratio_exp_active = zeros(6,1);
        
        for iTest = 2:length(PD) 
            P_loop = pressureDrivenFitting.active.(PD{iTest}).("Pressure [mmHg]");
            lambdat_loop = pressureDrivenFitting.active.(PD{iTest}).("Circ. stretch [-]");
            sigmatt_loop = pressureDrivenFitting.active.(PD{iTest}).("Circ. stress [MPa]");
        
            [~,locator_max] = max(P_loop);
            P_loop_loading = P_loop(1:locator_max);
            P_loop_unloading = P_loop(locator_max:end);
        
            lambdat_loop_loading = lambdat_loop(1:locator_max);
            lambdat_loop_unloading = lambdat_loop(locator_max:end);
        
            sigmatt_loop_loading = sigmatt_loop(1:locator_max);
            sigmatt_loop_unloading = sigmatt_loop(locator_max:end);
        
            regression_line = polyfit(lambdat_loop_loading,sigmatt_loop_loading,1);
            dynamic_stiffness_loading = regression_line(1);
            regression_line = polyfit(lambdat_loop_unloading,sigmatt_loop_unloading,1);
            dynamic_stiffness_unloading = regression_line(1);
        
            sigmatt_qs_loading = interp1(ref_P_loading,ref_sigmatt_loading,P_loop_loading,"spline","extrap");
            lambdat_qs_loading = interp1(ref_P_loading,ref_lambdat_loading,P_loop_loading,"spline","extrap");
            regression_line = polyfit(lambdat_qs_loading,sigmatt_qs_loading,1);
            qs_stiffness_loading = regression_line(1);
        
            sigmatt_qs_unloading = interp1(ref_P_unloading,ref_sigmatt_unloading,P_loop_unloading,"spline","extrap");
            lambdat_qs_unloading = interp1(ref_P_unloading,ref_lambdat_unloading,P_loop_unloading,"spline","extrap");
            regression_line = polyfit(lambdat_qs_unloading,sigmatt_qs_unloading,1);
            qs_stiffness_unloading = regression_line(1);
        
            ratio_mod_active(iTest-1,1) = (dynamic_stiffness_loading+dynamic_stiffness_unloading)/(qs_stiffness_loading+qs_stiffness_unloading);
        end
        
        ref_P = fit_data.fittedCurves.activeFit.Lname.PD_100.("Exp. pressure [mmHg]");
        ref_lambdat = fit_data.fittedCurves.activeFit.Lname.PD_100.("Circ. stretch [-]");
        ref_sigmatt = fit_data.fittedCurves.activeFit.Lname.PD_100.("Exp. circ. stress [MPa]");
        
        [~,locator_max] = max(ref_P);
        ref_P_loading = ref_P(1:locator_max);
        ref_P_unloading = ref_P(locator_max:end);
        
        ref_lambdat_loading = ref_lambdat(1:locator_max);
        ref_lambdat_unloading = ref_lambdat(locator_max:end);
        
        ref_sigmatt_loading = ref_sigmatt(1:locator_max);
        ref_sigmatt_unloading = ref_sigmatt(locator_max:end);
        
        for iTest = 2:length(PD) 
            P_loop = fit_data.fittedCurves.activeFit.Lname.(PD{iTest}).("Exp. pressure [mmHg]");
            lambdat_loop = fit_data.fittedCurves.activeFit.Lname.(PD{iTest}).("Circ. stretch [-]");
            sigmatt_loop = fit_data.fittedCurves.activeFit.Lname.(PD{iTest}).("Exp. circ. stress [MPa]");
        
            [~,locator_max] = max(P_loop);
            P_loop_loading = P_loop(1:locator_max);
            P_loop_unloading = P_loop(locator_max:end);
        
            lambdat_loop_loading = lambdat_loop(1:locator_max);
            lambdat_loop_unloading = lambdat_loop(locator_max:end);
        
            sigmatt_loop_loading = sigmatt_loop(1:locator_max);
            sigmatt_loop_unloading = sigmatt_loop(locator_max:end);
        
            regression_line = polyfit(lambdat_loop_loading,sigmatt_loop_loading,1);
            dynamic_stiffness_loading = regression_line(1);
            regression_line = polyfit(lambdat_loop_unloading,sigmatt_loop_unloading,1);
            dynamic_stiffness_unloading = regression_line(1);
        
            sigmatt_qs_loading = interp1(ref_P_loading,ref_sigmatt_loading,P_loop_loading,"spline","extrap");
            lambdat_qs_loading = interp1(ref_P_loading,ref_lambdat_loading,P_loop_loading,"spline","extrap");
            regression_line = polyfit(lambdat_qs_loading,sigmatt_qs_loading,1);
            qs_stiffness_loading = regression_line(1);
        
            sigmatt_qs_unloading = interp1(ref_P_unloading,ref_sigmatt_unloading,P_loop_unloading,"spline","extrap");
            lambdat_qs_unloading = interp1(ref_P_unloading,ref_lambdat_unloading,P_loop_unloading,"spline","extrap");
            regression_line = polyfit(lambdat_qs_unloading,sigmatt_qs_unloading,1);
            qs_stiffness_unloading = regression_line(1);
        
            ratio_exp_active(iTest-1,1) = (dynamic_stiffness_loading+dynamic_stiffness_unloading)/(qs_stiffness_loading+qs_stiffness_unloading);
        end
    end

    %% Saving
    fit_data.fittedCurves.pressureDrivenFit = pressureDrivenFitting;
    freq = cell(1,size(fit_data.PD_dynamic,2)/3);
    for iTest = 1:size(fit_data.PD_dynamic,2)/3
        loc_underscore = find(fit_data.PD_dynamic{iTest} == '_');
        if(loc_underscore(end)-loc_underscore(end-1)-1 == 4)
            freq{iTest} = [fit_data.PD_dynamic{iTest}(loc_underscore(end-1)+1) '.' fit_data.PD_dynamic{iTest}(loc_underscore(end-1)+2:loc_underscore(end)-1) ' Hz'];
        elseif(loc_underscore(end)-loc_underscore(end-1)-1 == 5)
            freq{iTest} = [fit_data.PD_dynamic{iTest}(loc_underscore(end-1)+1:loc_underscore(end-1)+2) '.' fit_data.PD_dynamic{iTest}(loc_underscore(end-1)+3:loc_underscore(end)-1) ' Hz'];    
        else
           freq{iTest} = ['0.' fit_data.PD_dynamic{iTest}(loc_underscore(end-1)+1:loc_underscore(end)-1) ' Hz'];
        end
    end
    fit_data.D_to_QS_stiffeningRatio.passive.Experimental = array2table(ratio_exp_passive,"RowNames",{'40-80 mmHg' '80-120 mmHg' '120-160 mmHg'},...
        "VariableNames",freq);
    fit_data.D_to_QS_stiffeningRatio.passive.Modelled = array2table(ratio_mod_passive,"RowNames",{'40-80 mmHg' '80-120 mmHg' '120-160 mmHg'},...
        "VariableNames",freq);
    if(activeModel.activator == 1)
        fit_data.D_to_QS_stiffeningRatio.active.Experimental = array2table(ratio_exp_active,"RowNames",{'50-90 mmHg' '70-110 mmHg' '90-130 mmHg' '110-150 mmHg'...
            '130-170 mmHg' '150-190 mmHg'},"VariableNames",{'10 Hz'});
        fit_data.D_to_QS_stiffeningRatio.active.Modelled = array2table(ratio_mod_active,"RowNames",{'50-90 mmHg' '70-110 mmHg' '90-130 mmHg' '110-150 mmHg'...
            '130-170 mmHg' '150-190 mmHg'},"VariableNames",{'10 Hz'});
    end
    
    % save(fileName,'fit_data');
end