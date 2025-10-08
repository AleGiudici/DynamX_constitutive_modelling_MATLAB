function [parameterV] = finaliseViscousParameterVector(parameterV_normalised,elasticModel,viscousModel)
% finaliseViscousParameterVector(parameterV_normalised,elastic) builds the 
% final (de-normalised) viscous parameter vector. Function inputs are the 
% normalised viscous parameter vecotor "parameterV_normalised" and the
% structures "elasticModel"  and "viscousModel" which contains the  
% formulation of the chosen elastic and viscous models.


%% Unnormalising the viscous parameter
    viscous_normV = zeros(1,3*elasticModel.nConstituents);
    for i = 1:elasticModel.nConstituents
        viscous_normV(1+(i-1)*3:i*3) = viscousModel.normV;
    end
    
    parameterV = parameterV_normalised.*viscous_normV;
end