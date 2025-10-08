function [costV] = objectiveStatic(parameterV_normalised,inputDataMatrixQS,elasticModel)
% Objective function for the optimisation of the model parameters of a
% fully elastic passive wall model. The error is calculated using the
% normalised differences between modelled and experimental pressure and
% modelled and experimental transducer axial force.

    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    gtN = 0.0098; %scaling factor to convert g to N

    parameterV_temp=parameterV_normalised.*elasticModel.normV;
    
    j=0;
    parameterV_elastic=zeros(1,elasticModel.nParam_baseline);
    for i=1:elasticModel.nParam_baseline
        if(isnan(elasticModel.imposedParam(i)))
            j=j+1;
            parameterV_elastic(i)=parameterV_temp(j);
        elseif(isinf(elasticModel.imposedParam(i)))
            if(strcmp(elasticModel.modelName,'4-fiber families (not neo-Hookean)') || strcmp(elasticModel.modelName,'4-fiber families'))
                if(i==4 || i==7)
                    parameterV_elastic(i)=parameterV_temp(2);
                elseif(i==6 || i==8)
                    parameterV_elastic(i)=parameterV_temp(3);
                end
            end
        else
            parameterV_elastic(i)=elasticModel.imposedParam(i);
        end
    end
    
    lambdatV = inputDataMatrixQS(:,7);
    lambdazV = inputDataMatrixQS(:,8);
    lambdarV = inputDataMatrixQS(:,9);
    
    riV = inputDataMatrixQS(:,3);
    hV = inputDataMatrixQS(:,4);
    
    P_expV = inputDataMatrixQS(:,1);
    fT_expV = inputDataMatrixQS(:,5);
    
    weights_qs = inputDataMatrixQS(:,end-1);
       
    [sigmatt_rr_compV,sigmazz_rr_compV,~,~,~]=elasticModel.fun(parameterV_elastic,lambdatV,lambdazV,lambdarV);
    [P_compV,fT_compV] = stress2pressure_force(sigmatt_rr_compV,sigmazz_rr_compV,riV+hV/2,hV);
    
    cost_tt = (P_compV*mtN - P_expV)./mean(P_expV).*weights_qs;
    
    cost_zz = (fT_compV*gtN - fT_expV)./mean(fT_expV).*weights_qs;
    
    costV = cat(1,cost_tt,cost_zz); %cost function for optimization
    
    locator=isinf(costV); % looks for inf values in the cost vector
    costV(locator)=10^6;
    
    locator_NaN = isnan(costV);
    costV(locator_NaN)=10^10;
    
end