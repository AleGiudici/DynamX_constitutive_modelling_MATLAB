function [parameterV_elastic_static,parameterV_active_static] = viscoElastic2ElasticModel(modelFormulation,modelParameters)
%parameterV_elastic,parameterV_viscous,elasticModel)
% This function converts all the stiffness-like parameters of the passive
% strain energy function of a viscoelastic into what they would be if the 
% model was purely elastic. This conversion is done by dividing the
% stiffness-like parameter by 
% (1 + viscoelastic_gain * log(viscoelastic_time_constant_2 / viscoelastic_time_constant_1)

%% Importing model formulation and parameters
    % Importing model formulation
    elasticModel = modelFormulation.elasticModel;
    activeModel = modelFormulation.activeModel;
    
    % Importing model parameters
    parameterV_elastic = modelParameters.staticSEFparameters.("Final QLV values"); 
    parameterV_viscous = modelParameters.viscousParameters.("Values");

    parameterV_elastic_static = parameterV_elastic;

    if(strcmp(elasticModel.modelName,'4-fiber families 2') || strcmp(elasticModel.modelName,'4-fiber familiy model') || strcmp(elasticModel.modelName,'4-fiber families (not neo-Hookean)'))
        parameterV_elastic_static(1) = parameterV_elastic(1)/(1+parameterV_viscous(1)*log(parameterV_viscous(3)/parameterV_viscous(2)));
        parameterV_elastic_static(2) = parameterV_elastic(2)/(1+parameterV_viscous(4)*log(parameterV_viscous(6)/parameterV_viscous(5)));
        parameterV_elastic_static(4) = parameterV_elastic(4)/(1+parameterV_viscous(4)*log(parameterV_viscous(6)/parameterV_viscous(5)));
        parameterV_elastic_static(7) = parameterV_elastic(7)/(1+parameterV_viscous(4)*log(parameterV_viscous(6)/parameterV_viscous(5)));
    
    elseif(strcmp(elasticModel.modelName,'2-fibre family model') || strcmp(elasticModel.modelName,'2-fibre family model with dispersion'))
        parameterV_elastic_static(1) = parameterV_elastic(1)/(1+parameterV_viscous(1)*log(parameterV_viscous(3,1)/parameterV_viscous(2)));
        parameterV_elastic_static(2) = parameterV_elastic(2)/(1+parameterV_viscous(4)*log(parameterV_viscous(6,1)/parameterV_viscous(5)));
    
    elseif(strcmp(elasticModel.modelName,'Zulliger´s model'))
        parameterV_elastic_static(1) = parameterV_elastic(1)/(1+parameterV_viscous(1)*log(parameterV_viscous(3)/parameterV_viscous(2)));
        parameterV_elastic_static(2) = parameterV_elastic(2)/(1+parameterV_viscous(4)*log(parameterV_viscous(6)/parameterV_viscous(5)));
        parameterV_elastic_static(3) = parameterV_elastic(3)/(1+parameterV_viscous(4)*log(parameterV_viscous(6)/parameterV_viscous(5)));
    end
    
    if(activeModel.activator == 1)
        parameterV_active_static = modelParameters.activeSEFparameters.("Final values");
        parameterV_active_static(1) = parameterV_active_static(1)/(1+parameterV_viscous(7)*log(parameterV_viscous(9)/parameterV_viscous(8)));
    else
        parameterV_active_static = NaN;
    end
end
