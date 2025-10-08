function [parameterV,R2,RMSE] = runElasticFitting(inputMatElasticFitting,elastic)
% Function that interatively estimates the elastic parameters that fit the
% measured pseudoelastic arterial behaviour.

    options = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 8000); % setting max. no. of iterations to 8000

%% Running fitting iterations
    problem = createOptimProblem('lsqnonlin','objective',@(xV) objectiveStatic(xV,inputMatElasticFitting,elastic),'x0',(elastic.lb+elastic.ub)/10,'lb',elastic.lb,'ub',elastic.ub,'options',options);
    ms = MultiStart('UseParallel',true);
    parpool;
    [parameterV_normalised] = run(ms,problem,50);
    delete(gcp);

%% Determining fitting quality
    error = objectiveStatic(parameterV_normalised,inputMatElasticFitting,elastic);
    SS_tot=sum([(inputMatElasticFitting(:,1)-mean(inputMatElasticFitting(:,1)))/mean(inputMatElasticFitting(:,1));...
        (inputMatElasticFitting(:,5)-mean(inputMatElasticFitting(:,5)))/mean(inputMatElasticFitting(:,5))].^2);
    R2=1-sum(error.^2)/SS_tot; % coefficient of determination of the static fitting
    RMSE = sqrt(sum(error.^2)/length(error)); % root mean square error

%% Finilising the parameter vector: 1) de-normalise and 2) add imposed parameters
    [parameterV] = finaliseParameterVector(parameterV_normalised,elastic);
end