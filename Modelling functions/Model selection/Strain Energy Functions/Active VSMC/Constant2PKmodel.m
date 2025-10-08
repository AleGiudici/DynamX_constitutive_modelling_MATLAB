function [sigmatt_rr,sigmazz_rr, W, Ctttt, Czzzz]=Constant2PKmodel(parameterV,lambdattV,lambdazzV,lambdarrV)
% Active stress  model that models the VSMC's stress as a constat 2nd
% Piola-Kirchhoff stress.
%
% parameterV is the 2-element vector of model parameters where:
%   position 1: VSMC stiffness like parameter [MPa]
%   position 2: VSMC circumferential deposition stretch [-]

sigmatt_rr = parameterV(1)*(lambdattV*parameterV(2)).^2;
sigmazz_rr = 0;

W = parameterV(1)/2*(lambdattV*parameterV(2)).^2;

Ctttt = 2*parameterV(1)*(lambdattV*parameterV(2)).^2;
Czzzz = 0;

end