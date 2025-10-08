function weight = FungMemoryFun(parameterV_viscous,timeV,derivative)
% Fung's reduced order relaxation function (for detailed description see
% Y.C. Fung, Biomechanics: Mechanical properties of living tissues 
% (2nd edition), 1993 - Chapter 7.6)
%
% "derivative controls" the output:
%   if derivative = 0 --> return the standard relaxation function
%   if derivative ~= 0 --> return the first derivative of the relaxation
%   function (for integration by parts)

if(~exist('derivative'))
    derivative = 1;
end

if(derivative)
    weight = (parameterV_viscous(1)./timeV.*(exp(-timeV/parameterV_viscous(2))-exp(-timeV/parameterV_viscous(3)))/...
        (1+parameterV_viscous(1)*log(parameterV_viscous(3)/parameterV_viscous(2))));
else
    weight = (1+parameterV_viscous(1)*(expint(timeV/parameterV_viscous(3))-expint(timeV/parameterV_viscous(2))))/...
        (1+parameterV_viscous(1)*log(parameterV_viscous(3)/parameterV_viscous(2)));
end

end
