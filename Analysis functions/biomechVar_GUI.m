function fit_data = biomechVar_GUI(targetP, SBP, DBP, fit_data, app)

    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    rho_blood = 1050;

%% Importing data from simulated static behaviour at in vivo stretch
    lambdattV = fit_data.modBehaviour.PD_100.("Circumferential stretch [-]");
    lambdazzV = fit_data.modBehaviour.PD_100.("Axial stretch [-]");
    lambdarrV = 1./lambdattV./lambdazzV;
    pressureV = fit_data.modBehaviour.PD_100.("Pressure [mmHg]");

%% Importing model formulation and parameters
    % Importing model formulation
    elasticModel = fit_data.modelFormulation.elasticModel;
    
    if(fit_data.modelFormulation.viscousModel.activator == 1)
        [parameterV_elastic_static,~] = viscoElastic2ElasticModel(fit_data.modelFormulation,fit_data.modelParameters);
    else
        parameterV_elastic_static = fit_data.modelParameters.staticSEFparameters.("First QS values");
    end

%% Checking for previously calculated variables
    if(isfield(fit_data, 'bioMechVar'))
        k = size(fit_data.bioMechVar,2);  
        bioMechVarMat = zeros(17,k+1);
        bioMechVarMat(:,1:end-length(targetP)) = table2array(fit_data.bioMechVar);
        temp = fit_data.bioMechVar;
        column_headings = strings(1,k+1);
        for i = 1:k
            column_headings{i} = temp.Properties.VariableNames{i};
        end
    else
        k = 0;
    end
    
    for j = 1:size(targetP,2)       
        [~,locator_P] = min(abs(pressureV-targetP(1,j)));
            
        [sigmatt_rr_P,sigmazz_rr_P,~,Ctttt_P,Czzzz_P] = elasticModel.fun(parameterV_elastic_static,lambdattV(locator_P),lambdazzV(locator_P),lambdarrV(locator_P));
        [sigmatt_rr_P_elastin,sigmazz_rr_P_elastin,~,~,~] = elasticModel.fun(parameterV_elastic_static,lambdattV(locator_P),lambdazzV(locator_P),lambdarrV(locator_P),[1,0]);
        [sigmatt_rr_P_collagen,sigmazz_rr_P_collagen,~,~,~] = elasticModel.fun(parameterV_elastic_static,lambdattV(locator_P),lambdazzV(locator_P),lambdarrV(locator_P),[0,1]);
        
        if(~isnan(SBP) && ~isnan(DBP))
            % Elastic energy at diastole
            [~,locator_P_lb]=min(abs(pressureV-DBP));
            [~,~,W_lb,~,~]=elasticModel.fun(parameterV_elastic_static,lambdattV(locator_P_lb),lambdazzV(locator_P),lambdarrV(locator_P_lb));
            [~,~,W_lb_elastin,~,~]=elasticModel.fun(parameterV_elastic_static,lambdattV(locator_P_lb),lambdazzV(locator_P),lambdarrV(locator_P_lb),[1,0]);
            [~,~,W_lb_collagen,~,~]=elasticModel.fun(parameterV_elastic_static,lambdattV(locator_P_lb),lambdazzV(locator_P),lambdarrV(locator_P_lb),[0,1]);
            
            % Elastic energy at systole
            [~,locator_P_ub]=min(abs(pressureV-SBP));
            [~,~,W_ub,~,~]=elasticModel.fun(parameterV_elastic_static,lambdattV(locator_P_ub),lambdazzV(locator_P),lambdarrV(locator_P_ub));
            [~,~,W_ub_elastin,~,~]=elasticModel.fun(parameterV_elastic_static,lambdattV(locator_P_ub),lambdazzV(locator_P),lambdarrV(locator_P_ub),[1,0]);
            [~,~,W_ub_collagen,~,~]=elasticModel.fun(parameterV_elastic_static,lambdattV(locator_P_ub),lambdazzV(locator_P),lambdarrV(locator_P_ub),[0,1]);
    
            % PWV, distensibility and compliance (systolic-to-diastolic
            ri_dia = fit_data.modBehaviour.PD_100.("Inner radius [mm]")(locator_P_lb);
            ri_sys = fit_data.modBehaviour.PD_100.("Inner radius [mm]")(locator_P_ub);
                
            PWV = ri_dia*sqrt((SBP-DBP)*133.32/rho_blood/(ri_sys^2-ri_dia^2));
    
            distensibility = 1/rho_blood/PWV^2*10^6;
    
            compliance = pi*ri_dia^2*distensibility;
            
            % Plotting on app graphs
            if(strcmp(app.graph_choice.Value,'Pressure and force'))
                hold(app.axes1, 'on')
                plot(app.axes1,2*fit_data.modBehaviour.PD_100.("Inner radius [mm]")(locator_P_lb:locator_P_ub), fit_data.modBehaviour.PD_100.("Pressure [mmHg]")(locator_P_lb:locator_P_ub),'c','LineWidth',4)
                plot(app.axes1,2*fit_data.modBehaviour.PD_100.("Inner radius [mm]")(locator_P), fit_data.modBehaviour.PD_100.("Pressure [mmHg]")(locator_P),'r-*','LineWidth',10)
                hold(app.axes1, 'off')
    
                hold(app.axes2, 'on')
                plot(app.axes2,fit_data.modBehaviour.PD_100.("Pressure [mmHg]")(locator_P_lb:locator_P_ub), fit_data.modBehaviour.PD_100.("Axial force [g]")(locator_P_lb:locator_P_ub),'c','LineWidth',4)
                plot(app.axes2,fit_data.modBehaviour.PD_100.("Pressure [mmHg]")(locator_P), fit_data.modBehaviour.PD_100.("Axial force [g]")(locator_P),'r-*','LineWidth',10)
                hold(app.axes2, 'off')
    
            elseif(strcmp(app.graph_choice.Value,'Stresses'))
                hold(app.axes1, 'on')
                plot(app.axes1,fit_data.modBehaviour.PD_100.("Circumferential stretch [-]")(locator_P_lb:locator_P_ub), fit_data.modBehaviour.PD_100.("Circumferential stress [MPa]")(locator_P_lb:locator_P_ub),'c','LineWidth',4)
                plot(app.axes1,fit_data.modBehaviour.PD_100.("Circumferential stretch [-]")(locator_P), fit_data.modBehaviour.PD_100.("Circumferential stress [MPa]")(locator_P),'r-*','LineWidth',10)
                hold(app.axes1, 'off')
    
                hold(app.axes2, 'on')
                plot(app.axes2,fit_data.modBehaviour.PD_100.("Circumferential stretch [-]")(locator_P_lb:locator_P_ub), fit_data.modBehaviour.PD_100.("Axial stress [MPa]")(locator_P_lb:locator_P_ub),'c','LineWidth',4)
                plot(app.axes2,fit_data.modBehaviour.PD_100.("Circumferential stretch [-]")(locator_P), fit_data.modBehaviour.PD_100.("Axial stress [MPa]")(locator_P),'r-*','LineWidth',10)
                hold(app.axes2, 'off')
            end
        else
            W_ub = NaN;
            W_lb = NaN;
            PWV = NaN;
            distensibility = NaN;
            compliance = NaN;
        end
    end
    
    bioMechVarMat(1,k+1:k+j) = targetP;
    bioMechVarMat(2,k+1:k+j) = sigmatt_rr_P-targetP(1,1:j)/2*mtN; % circumferential stress
    bioMechVarMat(3,k+1:k+j) = sigmazz_rr_P-targetP(1,1:j)/2*mtN; % axail stress
    bioMechVarMat(4,k+1:k+j) = Ctttt_P(1,1:j); % circumferential stiffness
    bioMechVarMat(5,k+1:k+j) = Czzzz_P(1,1:j); % axial stiffness
    bioMechVarMat(6,k+1:k+j) = SBP(1,1:j);
    bioMechVarMat(7,k+1:k+j) = DBP(1,1:j);
    bioMechVarMat(8,k+1:k+j) = (W_ub(1,1:j)-W_lb(1,1:j))*1000; % systolic-diastolic elastic energy
    bioMechVarMat(9,k+1:k+j) = PWV(1,1:j); 
    bioMechVarMat(10,k+1:k+j) = distensibility(1,1:j); 
    bioMechVarMat(11,k+1:k+j) = compliance(1,1:j);
    bioMechVarMat(12,k+1:k+j) = sigmatt_rr_P_elastin/sigmatt_rr_P*100; % elastin's circumferential load bearing
    bioMechVarMat(13,k+1:k+j) = sigmatt_rr_P_collagen/sigmatt_rr_P*100; % collagen's circumferential load bearing
    bioMechVarMat(14,k+1:k+j) = sigmazz_rr_P_elastin/sigmazz_rr_P*100; % elastin's axial load bearing
    bioMechVarMat(15,k+1:k+j) = sigmazz_rr_P_collagen/sigmazz_rr_P*100; % collagen's axial load bearing
    if(~isnan(SBP) && ~isnan(DBP))
        bioMechVarMat(16,k+1:k+j) = (W_ub_elastin(1,1:j)-W_lb_elastin(1,1:j))/(W_ub(1,1:j)-W_lb(1,1:j))*100; % elastin's elastic energy
        bioMechVarMat(17,k+1:k+j) = (W_ub_collagen(1,1:j)-W_lb_collagen(1,1:j))/(W_ub(1,1:j)-W_lb(1,1:j))*100; % collagen's elastic energy
    else
        bioMechVarMat(16,k+1:k+j) = NaN;
        bioMechVarMat(17,k+1:k+j) = NaN;
    end
    
    row_headings = ["Target pressure [mmHg]","Circumferential stress [MPa]","Axial stress [MPa]",...
        "Circumferential stiffness [MPa]","Axial stiffness [MPa]","Systolic pressure [mmHg]",...
        "Diastolic pressure [mmHg]","Stored elastic energy [kPa]","PWV [m/s]",...
        "Distensibility [1/MPa]","Compliance [mm^4/N]","Elastin circ. load-bearing [%]",...
        "Collagen circ.load-bearing [%]","Elastin axial load-bearing [%]","Collagen axial load-bearing [%]",...
        "Elastin stored energy [%]","Collagen stored energy [%]"];
    column_headings{k+1} = datestr(datetime);
    
    fit_data.bioMechVar = array2table(bioMechVarMat,"RowNames",row_headings,"VariableNames",column_headings);
end