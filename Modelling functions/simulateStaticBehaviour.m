function fit_data = simulateStaticBehaviour(fit_data)
    PD = {'PD_95' 'PD_100' 'PD_105'};

    % Scaling factors
    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2

    % Importing model formulation
    elasticModel = fit_data.modelFormulation.elasticModel;
    
    % Reference configuration geometry
    refInnerRadius = fit_data.vesselConfigurations.reference_configuration.("Inner radius [mm]");
    refOuterRadius = fit_data.vesselConfigurations.reference_configuration.("Outer radius [mm]");
    refThickness = fit_data.vesselConfigurations.reference_configuration.("Thickness [mm]");
    refAxialLength = fit_data.vesselConfigurations.reference_configuration.("Axial length [mm]");
    vessel_volume = pi*(refOuterRadius^2-refInnerRadius^2)*refAxialLength;

    overallAxStretch = fit_data.in_vivo_axial_stretch;
    
    targetPressureV = (10:1:200)';
    lambdaZPlotV = [0.95*overallAxStretch overallAxStretch 1.05*overallAxStretch]'; % list of axial stretches to be used for the simulation of inflation
    
    if(fit_data.modelFormulation.viscousModel.activator==1)
        [parameterV_elastic_static,~] = viscoElastic2ElasticModel(fit_data.modelFormulation,fit_data.modelParameters);
    else
        parameterV_elastic_static = fit_data.modelParameters.staticSEFparameters.("First QS values");
    end
    
    for iStretch = 1:length(lambdaZPlotV)
        lambdaZPlot = lambdaZPlotV(iStretch);
    
        riFoundV = zeros(length(targetPressureV),1);
        
        %Inner radius calculation by sweeping pressure for fixed axial stretch
        for iP = 1:length(targetPressureV)        
            targetPressure = targetPressureV(iP);    
            riFoundV(iP) = fzero(@(ri) (getPressure(parameterV_elastic_static, ri,lambdaZPlot,refThickness,refOuterRadius,refAxialLength,elasticModel.fun) - targetPressure), refOuterRadius); % estimated internal radius
        end
        
        roV = radiusFromIincompressibility(riFoundV,lambdaZPlot*refAxialLength,vessel_volume,2); % outer radius - conservation of volume
        hV = roV-riFoundV; % wall thickness
                   
        lambdattV = (riFoundV+hV/2)/(refInnerRadius+refThickness/2);
        lambdarrV = 1./(lambdaZPlot*lambdattV); 
    
        [sigmatt_rrV,sigmazz_rrV,strain_energyV,CttttV,CzzzzV] = elasticModel.fun(parameterV_elastic_static,lambdattV,lambdaZPlot,lambdarrV);
        [sigmatt_rrV_elastin,sigmazz_rrV_elastin,strain_energyV_elastin,~,~] = elasticModel.fun(parameterV_elastic_static,lambdattV,lambdaZPlot,lambdarrV,[1,0]);
        
        [PV,fTV] = stress2pressure_force(sigmatt_rrV,sigmazz_rrV,riFoundV+hV/2,hV);
    
        sigmattV=sigmatt_rrV-PV/2*mtN; % Lagrange multiplier to get the actual stress
        sigmazzV=sigmazz_rrV-PV/2*mtN; % Lagrange multiplier to get the actual stress
    
        elastin_circ_load_bearing = sigmatt_rrV_elastin./sigmatt_rrV*100;
        collagen_circ_load_bearing = (1-sigmatt_rrV_elastin./sigmatt_rrV)*100;
    
        elastin_axial_load_bearing = sigmazz_rrV_elastin./sigmazz_rrV*100;
        collagen_axial_load_bearing = (1-sigmazz_rrV_elastin./sigmazz_rrV)*100;
    
        elastin_energy_bearing = strain_energyV_elastin./strain_energyV*100;
        collagen_energy_bearing = (1-strain_energyV_elastin./strain_energyV)*100;
    
        temp_data = [PV (roV-hV) roV fTV sigmattV...
        lambdattV sigmazzV ones(length(fTV),1)*lambdaZPlot CttttV CzzzzV strain_energyV...
        elastin_circ_load_bearing collagen_circ_load_bearing elastin_axial_load_bearing...
        collagen_axial_load_bearing elastin_energy_bearing collagen_energy_bearing];
    
        headings = ["Pressure [mmHg]","Inner radius [mm]","Outer radius [mm]",...
            "Axial force [g]","Circumferential stress [MPa]","Circumferential stretch [-]",...
            "Axial stress [MPa]","Axial stretch [-]","Circumferential stiffness [MPa]",...
            "Axial stiffness [MPa]", "Elastic energy [MPa]", "Elastin circ. load bearing [%]",...
            "Collagen circ. load bearing [%]", "Elastin axial load bearing [%]",...
            "Collagen axial load bearing [%]", "Elastin elastic energy [%]", "Collagen elastic energy [%]"];
        
        fit_data.modBehaviour.(PD{iStretch}) = array2table(temp_data,"VariableNames",headings);
    end
end