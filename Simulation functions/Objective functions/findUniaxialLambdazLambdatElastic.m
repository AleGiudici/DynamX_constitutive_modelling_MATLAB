function [error] = findUniaxialLambdazLambdatElastic(newLambdas,parameterV_elastic,parameterV_active,elasticModel,activeModel,...
    targetPressure,refWallThickness,refMidWallRadius,activatorActive)
% [error] = findUniaxialLambdazLambdatElastic(newLambdas,...
% parameterV_elastic,parameterV_active,elasticModel,activeModel,...
% targetPressure,refWallThickness,refMidWallRadius,activatorActive)
%
% This function is the cost function used for the iterative estimation of
% the circumferential and axial stretches in a simulated uniaxial tensile
% test in the circumferential direction, given the artery elastic model 
% (which accounts for the summed contribution of elastin, collagen and, 
% possibly, VSMC). The VSMC contribution may be switched off by setting
% activatorActive to 0.
%
% The cost is defined in terms of the difference between the modelled
% circumferential stress and the pressure-equivalent stress (given the
% target pressure "targetPressure") and the difference between the modelled
% axial stress and 0 (given that the axial stress must be null in a
% uniaxial tensile test in the circumferential direction).
%
% Note that parameterV_elastic and parameterV_active need to reflect the
% behaviour of a purerly elastic artery. If these are taken from a
% viscoelastic model, some adjustements are necessary (see the function
% "viscoElastic2ElasticModel" in the "Modelling function" folder).
%
% The model formulation ("elasticModel" and "activeModel"), the model 
% parameters ("parameterV_elastic" and "parameterV_active"), the reference 
% geometry of the vessel ("refthickness" and "refMidWallRadius") are 
% required inputs. 
    
    lambdat = newLambdas(1);
    lambdaz = newLambdas(2);
    
    lambdar = 1/lambdat/lambdaz;
           
    [sigmatt,sigmazz] = elasticModel.fun(parameterV_elastic,lambdat,lambdaz,lambdar);
    
    if(activatorActive)
        [sigmatt_act,sigmazz_act] = activeModel.fun(parameterV_active,lambdat,lambdaz,lambdar);
    
        sigmatt = sigmatt+sigmatt_act;
        sigmazz = sigmazz+sigmazz_act;
    end
    
    rm = refMidWallRadius*lambdat;
    ri = sqrt(rm^2-(refMidWallRadius^2-(refMidWallRadius-refWallThickness/2)^2)/lambdaz);
    ro = sqrt(rm^2+((refMidWallRadius+refWallThickness/2)^2-refMidWallRadius^2)/lambdaz);
    h = ro-ri;
    
    sigmatt_target = targetPressure*133.32/10^6*rm/h;
        
    error = [sigmatt_target-sigmatt, -sigmazz];    
end