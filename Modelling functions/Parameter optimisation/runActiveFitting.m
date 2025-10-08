function [parameterV_active,parameterV_viscous,R2_active,RMSE_active,fit_data] = runActiveFitting(inputMatActiveFitting,elasticModel,viscousModel,...
    activeModel,parameterV_elastic_2ndQLV,parameterV_viscous,reference_configuration_for_active,unloaded_configuration,n_dataPoints,weight_vector,fit_data)

    elasticModel_for_active = elasticModel;
    elasticModel_for_active.imposedParam = NaN(1,elasticModel_for_active.nParam_baseline);
    elasticModel_for_active.imposedParam([2,3,7,8,9,10,11]) = [0 0 0 0 0.000001 parameterV_elastic_2ndQLV(10) reference_configuration_for_active(4)/unloaded_configuration(4)];
    elasticModel_for_active.nParam = 5;
    elasticModel_for_active.lb = [0 0 0 0 0.5];
    elasticModel_for_active.ub = [1 1 1 1 1.5];
    elasticModel_for_active.normV = [0.2 0.2 pi/2 50, 1];

    %% Optomisation steps 1
    lb = zeros(1,elasticModel_for_active.nParam+elasticModel_for_active.nConstituents*viscousModel.nParam);
    ub = zeros(1,elasticModel_for_active.nParam+elasticModel_for_active.nConstituents*viscousModel.nParam);

    lb(1:elasticModel_for_active.nParam) = elasticModel_for_active.lb;
    ub(1:elasticModel_for_active.nParam) = elasticModel_for_active.ub;

    for i = 1:elasticModel_for_active.nConstituents
        if(parameterV_viscous(1+(i-1)*3)<0.00001)
            lb(elasticModel_for_active.nParam+1+(i-1)*3:elasticModel_for_active.nParam+i*3) = [0 0.001 100]./viscousModel.normV;
            ub(elasticModel_for_active.nParam+1+(i-1)*3:elasticModel_for_active.nParam+i*3) = [0 0.001 100]./viscousModel.normV;
        else
            lb(elasticModel_for_active.nParam+1+(i-1)*3:elasticModel_for_active.nParam+i*3) = viscousModel.lb;
            ub(elasticModel_for_active.nParam+1+(i-1)*3:elasticModel_for_active.nParam+i*3) = viscousModel.ub;
            lb(elasticModel_for_active.nParam+i*3) = parameterV_viscous(i*3)/viscousModel.normV(3);
            ub(elasticModel_for_active.nParam+i*3) = parameterV_viscous(i*3)/viscousModel.normV(3);
        end
    end

    x0 = (lb+ub)/10;
    test_list = fieldnames(inputMatActiveFitting.passive_reference);

    options = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 8000); % setting max. no. of iterations to 8000
    
    problem = createOptimProblem('lsqnonlin','objective',@(xV) objectiveViscoelastic(xV,inputMatActiveFitting.passive_reference,elasticModel_for_active,viscousModel,149,ones(1,length(test_list)),test_list),'x0',x0,'lb',lb,'ub',ub,'options',options);
    ms = MultiStart('UseParallel',true);
    
    parpool;
    [parameterV_norm_1] = run(ms,problem,50);
    delete(gcp);

    %%
    lb = zeros(1,elasticModel_for_active.nParam+activeModel.nParam_baseline+...
        (elasticModel_for_active.nConstituents+activeModel.nConstituents)*viscousModel.nParam);
    ub = zeros(1,elasticModel_for_active.nParam+activeModel.nParam_baseline+...
        (elasticModel_for_active.nConstituents+activeModel.nConstituents)*viscousModel.nParam);

    lb(1:elasticModel_for_active.nParam) = parameterV_norm_1(1:elasticModel_for_active.nParam);
    ub(1:elasticModel_for_active.nParam) = parameterV_norm_1(1:elasticModel_for_active.nParam);

    lb(elasticModel_for_active.nParam+1:elasticModel_for_active.nParam+activeModel.nParam_baseline) = activeModel.lb_baseline;
    ub(elasticModel_for_active.nParam+1:elasticModel_for_active.nParam+activeModel.nParam_baseline) = activeModel.ub_baseline;

    for i = 1:elasticModel_for_active.nConstituents+activeModel.nConstituents
        if(i<=2)
            lb(elasticModel_for_active.nParam+activeModel.nParam_baseline+1+(i-1)*3:elasticModel_for_active.nParam+activeModel.nParam_baseline+i*3) = parameterV_norm_1(elasticModel_for_active.nParam+1+(i-1)*3:elasticModel_for_active.nParam+i*3);
            ub(elasticModel_for_active.nParam+activeModel.nParam_baseline+1+(i-1)*3:elasticModel_for_active.nParam+activeModel.nParam_baseline+i*3) = parameterV_norm_1(elasticModel_for_active.nParam+1+(i-1)*3:elasticModel_for_active.nParam+i*3);
        else
            lb(elasticModel_for_active.nParam+activeModel.nParam_baseline+1+(i-1)*3:elasticModel_for_active.nParam+activeModel.nParam_baseline+i*3) = viscousModel.lb;
            ub(elasticModel_for_active.nParam+activeModel.nParam_baseline+1+(i-1)*3:elasticModel_for_active.nParam+activeModel.nParam_baseline+i*3) = viscousModel.ub;
        end
    end

    ub(end) = 10;

    x0 = (lb+ub)/10;

    options = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 8000); % setting max. no. of iterations to 8000
    
    problem = createOptimProblem('lsqnonlin','objective',@(xV) objectiveActive(xV,inputMatActiveFitting,elasticModel_for_active,viscousModel,activeModel,n_dataPoints/2,weight_vector),'x0',x0,'lb',lb,'ub',ub,'options',options);
    ms = MultiStart('UseParallel',true);
    
    parpool;
    [parameterV_norm] = run(ms,problem,50);
    delete(gcp);

    parameterV_temp=parameterV_norm(1:elasticModel_for_active.nParam).*elasticModel_for_active.normV; 
    
    j=0;
    parameterV_elastic=zeros(1,elasticModel_for_active.nParam_baseline);
    for i=1:elasticModel_for_active.nParam_baseline
        if(isnan(elasticModel_for_active.imposedParam(i)))
            j=j+1;
            parameterV_elastic(i)=parameterV_temp(j);
        elseif(isinf(elasticModel_for_active.imposedParam(i)))
            if(strcmp(elasticModel_for_active.modelName,'4-fiber families (not neo-Hookean)') || strcmp(elasticModel_for_active.modelName,'4-fiber families'))
                if(i==4 || i==7)
                    parameterV_elastic(i)=parameterV_temp(2);
                elseif(i==6 || i==8)
                    parameterV_elastic(i)=parameterV_temp(3);
                end
            end
        else
            parameterV_elastic(i)=elasticModel_for_active.imposedParam(i);
        end
    end

    % De-normalise active parameter vector
    parameterV_active = parameterV_norm(elasticModel_for_active.nParam+1:elasticModel_for_active.nParam+activeModel.nParam_baseline).*activeModel.normV_baseline;

    % De-normalise viscous parameter vector
    visco_normV = zeros(1,(elasticModel_for_active.nConstituents+activeModel.nConstituents)*3);
    for i = 1:elasticModel_for_active.nConstituents+activeModel.nConstituents
        visco_normV(1+(i-1)*3:i*3) = viscousModel.normV;
    end
    parameterV_viscous = parameterV_norm(elasticModel_for_active.nParam+activeModel.nParam_baseline+1:end).*visco_normV;
    
    [fit_data,P_final,fT_final] = getActiveBehaviour(parameterV_elastic,parameterV_active,parameterV_viscous,inputMatActiveFitting,...
        elasticModel,activeModel,viscousModel,n_dataPoints,fit_data);

    % determine squared residual error
    e = objectiveActive(parameterV_norm,inputMatActiveFitting,elasticModel_for_active,viscousModel,activeModel,n_dataPoints/2,weight_vector);
    
    % Calculate the R2
    SS_tot=sum([(P_final(round(end/2):end)-mean(P_final(round(end/2):end)))/mean(P_final(round(end/2):end));(fT_final(round(end/2):end)-mean(fT_final(round(end/2):end)))/mean(fT_final(round(end/2):end))].^2);
    R2_active=1-sum(e.^2)/SS_tot; % coefficient of determination of the static fitting
    RMSE_active = sqrt(sum(e.^2)/length(e));
end
