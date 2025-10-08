function fit_data = biomechVar_script(targetP, lambdattV, lambdazzV, lambdarrV, P, elastic, fit_data, SBP, DBP)

check_PP = (~isnan(SBP))*(~isnan(DBP));

mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
rho_blood = 1050;

if(isfield(fit_data, 'bioMechVar'))
    k = size(fit_data.bioMechVar,2);  
    bioMechVarMat = zeros(18,k+1);
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
    
    i = j+k;

    [SEF_param_final,~] = viscoElastic2ElasticModel(fit_data.modelFormulation,fit_data.modelParameters);
    
    [~,locator_P]=min(abs(P-targetP(1,j)));

    % SEF_param = table2array(fit_data.modelParameters.staticSEFparameters);
    % if(fit_data.modelFormulation.viscousModel.activator == 1)
    %     viscous_param = table2array(fit_data.modelParameters.viscousParameters);
    % end
    % SEF_param_final = SEF_param(:,1);
    % 
    % if(fit_data.modelFormulation.viscousModel.activator == 1)
    %     if(strcmp(fit_data.modelFormulation.elasticModel.modelName,'4-fiber families 2') || strcmp(fit_data.modelFormulation.elasticModel.modelName,'4-fiber familiy model') || strcmp(fit_data.modelFormulation.elasticModel.modelName,'4-fiber families (not neo-Hookean)'))
    %         SEF_param_final(1) = SEF_param_final(1)/(1+viscous_param(1,1)*log(viscous_param(3,1)/viscous_param(2,1)));
    %         SEF_param_final(2) = SEF_param_final(2)/(1+viscous_param(4,1)*log(viscous_param(6,1)/viscous_param(5,1)));
    %         SEF_param_final(4) = SEF_param_final(4)/(1+viscous_param(4,1)*log(viscous_param(6,1)/viscous_param(5,1)));
    %         SEF_param_final(7) = SEF_param_final(7)/(1+viscous_param(4,1)*log(viscous_param(6,1)/viscous_param(5,1)));
    % 
    %     elseif(strcmp(fit_data.modelFormulation.elasticModel.modelName,'2-fibre family model') || strcmp(fit_data.modelFormulation.elasticModel.modelName,'2-fibre family model with dispersion'))
    %         SEF_param_final(1) = SEF_param_final(1)/(1+viscous_param(1,1)*log(viscous_param(3,1)/viscous_param(2,1)));
    %         SEF_param_final(2) = SEF_param_final(2)/(1+viscous_param(4,1)*log(viscous_param(6,1)/viscous_param(5,1)));
    % 
    %     elseif(strcmp(fit_data.modelFormulation.elasticModel.modelName,'Zulliger´s model'))
    %         SEF_param_final(1) = SEF_param_final(1)/(1+viscous_param(1,1)*log(viscous_param(3,1)/viscous_param(2,1)));
    %         SEF_param_final(2) = SEF_param_final(2)/(1+viscous_param(4,1)*log(viscous_param(6,1)/viscous_param(5,1)));
    %         SEF_param_final(3) = SEF_param_final(3)/(1+viscous_param(4,1)*log(viscous_param(6,1)/viscous_param(5,1)));
    %     end
    % end

    [sigmatt_rr_P,sigmazz_rr_P,~,Ctttt_P,Czzzz_P,Czztt_P]=elastic.fun(SEF_param_final,lambdattV(locator_P),lambdazzV,lambdarrV(locator_P));
    [sigmatt_rr_P_elastin,sigmazz_rr_P_elastin,~,~,~]=elastic.fun(SEF_param_final,lambdattV(locator_P),lambdazzV,lambdarrV(locator_P),[1,0]);
    [sigmatt_rr_P_collagen,sigmazz_rr_P_collagen,~,~,~]=elastic.fun(SEF_param_final,lambdattV(locator_P),lambdazzV,lambdarrV(locator_P),[0,1]);
        
    if(check_PP)
        [~,locator_P_lb]=min(abs(P-DBP));
        [~,~,W_lb,~,~]=elastic.fun(SEF_param_final,lambdattV(locator_P_lb),lambdazzV,lambdarrV(locator_P_lb));
        [~,~,W_lb_elastin,~,~]=elastic.fun(SEF_param_final,lambdattV(locator_P_lb),lambdazzV,lambdarrV(locator_P_lb),[1,0]);
        [~,~,W_lb_collagen,~,~]=elastic.fun(SEF_param_final,lambdattV(locator_P_lb),lambdazzV,lambdarrV(locator_P_lb),[0,1]);

        [~,locator_P_ub]=min(abs(P-SBP));
        [~,~,W_ub,~,~]=elastic.fun(SEF_param_final,lambdattV(locator_P_ub),lambdazzV,lambdarrV(locator_P_ub));
        [~,~,W_ub_elastin,~,~]=elastic.fun(SEF_param_final,lambdattV(locator_P_ub),lambdazzV,lambdarrV(locator_P_ub),[1,0]);
        [~,~,W_ub_collagen,~,~]=elastic.fun(SEF_param_final,lambdattV(locator_P_ub),lambdazzV,lambdarrV(locator_P_ub),[0,1]);

        % PWV, distensibility and compliance
        
        ri_dia = fit_data.modBehaviour.PD_100.("Inner radius [mm]")(locator_P_lb);
        ri_sys = fit_data.modBehaviour.PD_100.("Inner radius [mm]")(locator_P_ub);
            
        PWV = ri_dia*sqrt((SBP-DBP)*133.32/rho_blood/(ri_sys^2-ri_dia^2));

        distensibility = 1/rho_blood/PWV^2*10^6;

        compliance = pi*ri_dia^2*distensibility;
    else
        W_ub = NaN;
        W_lb = NaN;
        PWV = NaN;
        distensibility = NaN;
        compliance = NaN;
    end
end

bioMechVarMat(1,k+1:k+j) = targetP;
bioMechVarMat(2,k+1:k+j) = sigmatt_rr_P-targetP(1,1:j)/2*mtN;
bioMechVarMat(3,k+1:k+j) = sigmazz_rr_P-targetP(1,1:j)/2*mtN;
bioMechVarMat(4,k+1:k+j) = Ctttt_P(1,1:j);
bioMechVarMat(5,k+1:k+j) = Czzzz_P(1,1:j);
bioMechVarMat(6,k+1:k+j) = Czztt_P(1,1:j);
bioMechVarMat(7,k+1:k+j) = SBP(1,1:j);
bioMechVarMat(8,k+1:k+j) = DBP(1,1:j);
bioMechVarMat(9,k+1:k+j) = (W_ub(1,1:j)-W_lb(1,1:j))*1000;
bioMechVarMat(10,k+1:k+j) = PWV(1,1:j);
bioMechVarMat(11,k+1:k+j) = distensibility(1,1:j);
bioMechVarMat(12,k+1:k+j) = compliance(1,1:j);
bioMechVarMat(13,k+1:k+j) = sigmatt_rr_P_elastin/sigmatt_rr_P*100;
bioMechVarMat(14,k+1:k+j) = sigmatt_rr_P_collagen/sigmatt_rr_P*100;
bioMechVarMat(15,k+1:k+j) = sigmazz_rr_P_elastin/sigmazz_rr_P*100;
bioMechVarMat(16,k+1:k+j) = sigmazz_rr_P_collagen/sigmazz_rr_P*100;
bioMechVarMat(17,k+1:k+j) = (W_ub_elastin(1,1:j)-W_lb_elastin(1,1:j))/(W_ub(1,1:j)-W_lb(1,1:j))*100;
bioMechVarMat(18,k+1:k+j) = (W_ub_collagen(1,1:j)-W_lb_collagen(1,1:j))/(W_ub(1,1:j)-W_lb(1,1:j))*100;

row_headings = ["Target pressure [mmHg]","Circumferential stress [MPa]","Axial stress [MPa]",...
    "Circumferential stiffness [MPa]","Axial stiffness [MPa]","Axial-circumferential stiffness [MPa]","Systolic pressure [mmHg]",...
    "Diastolic pressure [mmHg]","Stored elastic energy [kPa]","PWV [m/s]",...
    "Distensibility [1/MPa]","Compliance [mm^4/N]","Elastin circ. load-bearing [%]",...
    "Collagen circ.load-bearing [%]","Elastin axial load-bearing [%]","Collagen axial load-bearing [%]",...
    "Elastin stored energy [%]","Collagen stored energy [%]"];
column_headings{k+1} = datestr(datetime);

fit_data.bioMechVar = array2table(bioMechVarMat,"RowNames",row_headings,"VariableNames",column_headings);
