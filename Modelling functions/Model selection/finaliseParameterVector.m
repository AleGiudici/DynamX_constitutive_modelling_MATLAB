function [parameterV] = finaliseParameterVector(parameterV_normalised,elastic)
% finaliseParameterVector(parameterV_normalised,elastic) builds the final SEF
% parameter vector by 1) de-normalising the parameter values, and 2) adding
% imposed parameters. Function inputs are normalised parameter vecotor
% "parameterV_normalised" and the structure "elastic" which contains all
% the information on the chosen elastic model.

%% Denormalising the parameter vector
    parameterV_normalised=parameterV_normalised.*elastic.normV;

%% Finalising parameter vector with imposed parameters
    j=0;
    parameterV=zeros(1,elastic.nParam_baseline);
    for i=1:elastic.nParam_baseline
        if(isnan(elastic.imposedParam(i)))
            j=j+1;
            parameterV(i)=parameterV_normalised(j);
        elseif(isinf(elastic.imposedParam(i)))
            if(strcmp(elastic.modelName,'4-fiber families (not neo-Hookean)') || strcmp(elastic.modelName,'4-fiber families'))
                if(i==4 || i==7)
                    parameterV(i)=parameterV_normalised(2);
                elseif(i==6 || i==8)
                    parameterV(i)=parameterV_normalised(3);
                end
            end
        else
            parameterV(i)=elastic.imposedParam(i);
        end
    end
end