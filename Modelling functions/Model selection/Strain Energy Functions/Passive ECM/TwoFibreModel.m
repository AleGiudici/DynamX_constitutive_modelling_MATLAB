function [sigmatt_rr,sigmazz_rr,W,Ctttt,Czzzz]=TwoFibreModel(parameterV,lambdatV,lambdazV,lambdarV,activator)
% The function returns the circumferential and axial stresses, the stored
% elastic energy, and the small-on-large circumferential and axial material
% stiffnesses for four-fibre family constitutive model given a parameter
% vector xV, the deformation in the three principal directions, and an activator 
% that activates/deactivates the contribution from each wal constituent.
%
% For a detailed desription of the constitutive model refer to Gasser et al. 
% Hyperelastic modelling of arterial layers with distributed collagen fibre orientations.  
% Journal of the Royal Society Interface, 3, pp.15-35, 2006.
%
% parameterV is the 8-element vector of model parameters where:
%   position 1: elastin stiffness-like parameter
%   position 2: diagonally oriented collagen fibre stiffness-like parameter
%   position 3: diagonally oriented collagen fibre non-linearity parameter
%   position 4: diagonally oriented collagen fibre orientation parameter
%   position 5: diagonally oriented collagen fibre non-linearity parameter for compression
%   position 6: elastin circumferential deposition stretch
%   position 7: elastin axial deposition stretch
%   position 8: collagen deposition stretch along the fibre direction
%
% For a detailed description of the small-on-large formulation refer to
% Baek et al. Theory of small on large: Potential utility in computations
% of fluid-solid interactions in arteries. Computer methods in applied
% mechanics nad engineering, 196(31-32), pp.3070-3078.

if(~exist('activator'))
    activator=[1,1]; % selective activation of collagen and elastin --> if not given as input, is set to [1, 1] 
    % (i.e., both active).
end

Ic=(lambdatV*parameterV(6)).^2+(lambdazV*parameterV(7)).^2+(lambdarV/(parameterV(6)*parameterV(7))).^2;
Id=(lambdatV.^2*sin(parameterV(4))^2+lambdazV.^2*cos(parameterV(4))^2)*parameterV(8)^2;

sigmatt_rr=activator(1)*parameterV(1)*((lambdatV*parameterV(6)).^2-(lambdarV/(parameterV(6)*parameterV(7))).^2) +...
    activator(2)*4*parameterV(2)*sin(parameterV(4))^2*(lambdatV*parameterV(8)).^2.*(Id-1).*exp((parameterV(3)*(Id>=1)+parameterV(5)*(Id<1)).*(Id-1).^2);

sigmazz_rr=activator(1)*parameterV(1)*((lambdazV*parameterV(7)).^2-(lambdarV/(parameterV(6)*parameterV(7))).^2) +...
    activator(2)*4*parameterV(2)*cos(parameterV(4))^2*(lambdazV*parameterV(8)).^2.*(Id-1).*exp((parameterV(3)*(Id>=1)+parameterV(5)*(Id<1)).*(Id-1).^2);

W=activator(1)*parameterV(1)*(Ic-3)+activator(2)*parameterV(2)/parameterV(3)*(exp((parameterV(3)*(Id>1)+parameterV(5)*(Id<1)).*(Id-1).^2)-1);

Ctttt=2*(sigmatt_rr+activator(1)*parameterV(1)*(lambdarV/(parameterV(6)*parameterV(7))).^2) + activator(2)*8*...
    (lambdatV*parameterV(8)).^4.*(1+2*(parameterV(3)*(Id>=1)+parameterV(5)*(Id<1)).*(Id-1).^2)*...
    parameterV(2)*sin(parameterV(4))^4.*exp((parameterV(3)*(Id>=1)+parameterV(5)*(Id<1)).*(Id-1).^2);
Czzzz=2*(sigmazz_rr+activator(1)*parameterV(1)*(lambdarV/(parameterV(6)*parameterV(7))).^2) + activator(2)*8*...
    (lambdazV*parameterV(8)).^4.*(1+2*(parameterV(3)*(Id>=1)+parameterV(5)*(Id<1)).*(Id-1).^2)*...
    parameterV(2)*cos(parameterV(4))^4.*exp((parameterV(3)*(Id>=1)+parameterV(5)*(Id<1)).*(Id-1).^2);

end