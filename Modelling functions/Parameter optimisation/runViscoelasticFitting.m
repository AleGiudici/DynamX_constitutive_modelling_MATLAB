function [parameterV_elastic,parameterV_viscous,RMSE,R2,fit_data] = runViscoelasticFitting(inputMatViscoelasticFitting,elastic,visco,weightsV,fit_data,optimisation_step)
% Function that interatively estimates the elastic parameters that fit the
% measured pseudoelastic arterial behaviour. 
% The input parameters are the structure with the experimental data
% "inputMatElasticFitting", the chosen elastic and viscous models "elastic" 
% and "viscous", the vector of weights "weightsV" (used to differentially 
% weight protocol experiments), the output matrix "fit_data" and the 
% optimisation step number "optimisation_step" (either 1 or 2).

%% Define the total number of experimental data points in the cost function
    n_dataPoints = 0; 
    test_list = fieldnames(inputMatViscoelasticFitting);
    for i = 1:length(test_list)
        locator = find(inputMatViscoelasticFitting.(test_list{i}).processed_data(:,end)==1);
        n_dataPoints = n_dataPoints+length(inputMatViscoelasticFitting.(test_list{i}).processed_data(locator:end,end));
    end
    
    if(optimisation_step == 1)
        for i = 1:elastic.nConstituents
            visco_lb((i-1)*3+1:i*3) = visco.lb;
            visco_ub((i-1)*3+1:i*3) = visco.ub;
        end
    else
        visco_lb = visco.lb;
        visco_ub = visco.ub;
    end
    
%% Define parameter space based on modelling choices
    lb=[elastic.lb visco_lb]; % lower bound
    ub=[elastic.ub visco_ub]; % normalised upper bound
    x0=[(elastic.lb+elastic.ub)/10, (visco_lb+visco_ub)/2];

%% Run viscoelastic fitting
    options = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 8000); % setting max. no. of iterations to 8000

    problem = createOptimProblem('lsqnonlin','objective',@(xV) objectiveViscoelastic(xV,inputMatViscoelasticFitting,elastic,visco,n_dataPoints,weightsV,test_list),'x0',x0,'lb',lb,'ub',ub,'options',options);
    ms = MultiStart('UseParallel',true);
    
    parpool;
    [parameterV_normalised] = run(ms,problem,50);
    delete(gcp);

%% Split model parameters in elastic and viscous parameters
    parameterV_elastic_normalised = parameterV_normalised(1:elastic.nParam);
    parameterV_viscous_normalised = parameterV_normalised(elastic.nParam+1:end);
        
%     % NEEDS TO BE DELETED
%     parameterV_viscous_temp(4:6) = parameterV_viscous_temp(1:3);
%     % NEEDS TO BE DELETED

%% Finilising the parameter vector: 1) de-normalise and 2) add imposed parameters
    [parameterV_elastic] = finaliseParameterVector(parameterV_elastic_normalised,elastic);
    [parameterV_viscous] = finaliseViscousParameterVector(parameterV_viscous_normalised,elastic,visco);
            
    [fit_data,P_final,fT_final] = getViscoelasticBehaviour(parameterV_elastic,parameterV_viscous,inputMatViscoelasticFitting,elastic,visco,n_dataPoints,fit_data,optimisation_step);
        
    % determine squared residual error
    e = objectiveViscoelastic(parameterV_normalised,inputMatViscoelasticFitting,elastic,visco,n_dataPoints,weightsV,test_list);
    
    % Calculate the R2
    SS_tot=sum([(P_final-mean(P_final))/mean(P_final);(fT_final-mean(fT_final))/mean(fT_final)].^2);
    R2=1-sum(e.^2)/SS_tot; % coefficient of determination of the static fitting
    RMSE = sqrt(sum(e.^2)/length(e));
end

