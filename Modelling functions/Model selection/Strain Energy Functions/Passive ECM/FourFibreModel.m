function [sigmatt_rr,sigmazz_rr, W, Ctttt, Czzzz]=FourFibreModel(parameterV,lambdatV,lambdazV,lambdarV,activator)
% The function returns the circumferential and axial stresses, the stored
% elastic energy, and the small-on-large circumferential and axial material
% stiffnesses for four-fibre family constitutive model given a parameter
% vector xV, the deformation in the three principal directions, and an activator 
% that activates/deactivates the contribution from each wal constituent.
%
% For a detailed desription of the constitutive model refer to Ferruzzi et al. 
% On constitutive descriptors of the biaxial mechanical behaviour of human abdominal  
% aorta and aneurysms. Journal of the Royal Society Interface, 8(56), pp.435-450, 2011.
%
% parameterV is the 12-element vector of model parameters where:
%   position 1: elastin stiffness-like parameter
%   position 2: circumferentially oriented collagen fibre stiffness-like parameter
%   position 3: circumferentially oriented collagen fibre non-linearity parameter
%   position 4: diagonally oriented collagen fibre stiffness-like parameter
%   position 5: diagonally oriented collagen fibre orientation parameter
%   position 6: diagonally oriented collagen fibre non-linearity parameter
%   position 7: axially oriented collagen fibre stiffness-like parameter
%   position 8: axially oriented collagen fibre non-linearity parameter
%   position 9: collagen non-linearity parameter for compression
%   position 10: elastin circumferential deposition stretch
%   position 11: elastin axial deposition stretch
%   position 12: collagen deposition stretch along the fibre direction
%   
%
% For a detailed description of the small-on-large formulation refer to
% Baek et al. Theory of small on large: Potential utility in computations
% of fluid-solid interactions in arteries. Computer methods in applied
% mechanics nad engineering, 196(31-32), pp.3070-3078.

if(~exist('activator'))
    activator=[1,1]; % selective activation of collagen and elastin --> if not given as input, is set to [1, 1]
end

I1 = (lambdatV*parameterV(10)).^2+(lambdazV*parameterV(11)).^2+(lambdarV/(parameterV(10)*parameterV(11))).^2; % First invariant

Ld = (lambdazV.^2*(cos(parameterV(5)))^2+lambdatV.^2*(sin(parameterV(5)))^2)*parameterV(12)^2; % squared stretch of the diagonal fibre
Lt = lambdatV.^2*parameterV(12)^2; % squared stretch of the circumferential fibre
Lz = lambdazV.^2*parameterV(12)^2; % squared stretch of the axial fibre

% sigmatt_rr = xV(1)*((lambdattV*xV(9)).^2-(lambdarrV/(xV(9)*xV(10))).^2)+ xV(2)*((lambdattV*xV(11)).^2-1).*exp(xV(3)*(((lambdattV*xV(11)).^2-1).^2)).*(lambdattV*xV(11)).^2+...
%     + 2*xV(4)*(Id-1).*exp(xV(6)*((Id-1).^2)).*(lambdattV*xV(12)).^2*sin(xV(5))^2; % circumeferential Cauchy stress N/mm2
% 
% sigmazz_rr = xV(1)*((lambdazzV*xV(10)).^2-(lambdarrV/(xV(9)*xV(10))).^2)+ xV(7)*((lambdazzV*xV(13)).^2-1).*exp(xV(8)*(((lambdazzV*xV(13)).^2-1).^2)).*(lambdazzV*xV(13)).^2+...
%     + 2*xV(4)*(Id-1).*exp(xV(6)*((Id-1).^2)).*(lambdazzV*xV(12)).^2*cos(xV(5))^2; % axial Cauchy stress N/mm2

k2_circ = parameterV(3)*(Lt>=1)+parameterV(9)*(Lt<1);
k2_diag = (parameterV(6)*(Ld>=1)+parameterV(9)*(Ld<1));
k2_axial = (parameterV(8)*(Lz>=1)+parameterV(9)*(Lz<1));

sigmatt_rr = activator(1)*parameterV(1)*((lambdatV*parameterV(10)).^2-(lambdarV/(parameterV(10)*parameterV(11))).^2) +...
    activator(2)*(parameterV(2)*(Lt-1).*exp(k2_circ.*((Lt-1).^2)).*Lt +...
    2*parameterV(4)*(Ld-1).*exp(k2_diag.*((Ld-1).^2)).*Lt*sin(parameterV(5))^2); % circumeferential Cauchy stress N/mm2...
     
sigmazz_rr = activator(1)*parameterV(1)*((lambdazV*parameterV(11)).^2-(lambdarV/(parameterV(10)*parameterV(11))).^2) +...
    activator(2)*(parameterV(7)*(Lz-1).*exp(k2_axial.*((Lz-1).^2)).*Lz +...
    2*parameterV(4)*(Ld-1).*exp(k2_diag.*((Ld-1).^2)).*Lz*cos(parameterV(5))^2); % axial Cauchy stress N/mm2

W = activator(1)*(parameterV(1)/2)*(I1-3) +...
    activator(2)*(parameterV(2)./(4*k2_circ).*(exp(k2_circ.*(Lt-1).^2)-1) +... 
    parameterV(7)./(4*k2_axial).*(exp(k2_axial.*(Lz-1).^2)-1) +...
    parameterV(4)./(2*k2_diag).*(exp(k2_diag.*(Ld-1).^2)-1)); %Nmm2

Ctttt = 2*(sigmatt_rr+activator(1)*parameterV(1)*(lambdarV/(parameterV(10)*parameterV(11))).^2) +...
    activator(2)*(2*parameterV(2)*Lt.^2.*(1+2*k2_circ.*(Lt-1).^2).*exp(k2_circ.*(Lt-1).^2) +...
    4*parameterV(4)*sin(parameterV(5))^4.*Lt.^2.*(1+2*k2_diag.*(Ld-1).^2).*exp(k2_diag.*(Ld-1).^2));

Czzzz = 2*(sigmazz_rr+activator(1)*parameterV(1)*(lambdarV/(parameterV(10)*parameterV(11))).^2) +...
    activator(2)*(2*parameterV(7)*Lz.^2.*(1+2*k2_axial.*(Lz-1).^2).*exp(k2_axial.*(Lz-1).^2) +...
    4*parameterV(4)*cos(parameterV(5))^4.*Lz.^2.*(1+2*k2_diag.*(Ld-1).^2).*exp(k2_diag.*(Ld-1).^2));

Czztt = activator(2)*4*parameterV(4)*lambdazV.^2*cos(parameterV(5))^2.*lambdatV.^2*sin(parameterV(5))^2.*...
    exp((parameterV(6)*(Ld>=1)+parameterV(9)*(Ld<1)).*(Ld-1).^2).*(1+2*(parameterV(6)*(Ld>=1)+parameterV(9)*(Ld<1)).*(Ld-1).^2);
end