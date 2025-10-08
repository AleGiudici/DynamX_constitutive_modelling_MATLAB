function [sigmatt_rr,sigmazz_rr, W, Ctttt, Czzzz, Czztt]=UFDmodelNNH(parameterV,lambdatV,lambdazV,lambdarV,activator,not_nH_coeff)
% The function returns the circumferential and axial stresses, the stored
% elastic energy, and the small-on-large circumferential and axial material
% stiffnesses for four-fibre family constitutive model given a parameter
% vector "parameterV", the deformation in the three principal directions, 
% an "activator" vector that activates/deactivates the contribution from
% each wall constituent and non-Neo-Hookean coefficient "not_nH_coeff" 
% which determines the shape of elastin response (not_nH_coeff = 0 yields
% standard four fiber model).
%
% For a detailed desription of the constitutive model refer to Giudici et al. 
% Constituent-based quasi-linear viscoelasticity: A revised quasi-linear 
% modelling framework to capture non-linear viscoelasticity in arteries. 
% Biomechanics and Modellining in Mechanobiology, 2023 [in press].
%
% parameterV is the 12-element vector of model parameters where:
%   position 1: elastin stiffness-like parameter
%   position 2: collagen fibre stiffness-like parameter
%   position 3: collagen fibre non-linearity parameter
%   position 4: fibre distribution
%   position 5: collagen non-linearity parameter for compression
%   position 6: elastin circumferential deposition stretch
%   position 7: elastin axial deposition stretch
%   position 8: collagen deposition stretch along the fibre direction
%   
%
% For a detailed description of the small-on-large formulation refer to
% Baek et al. Theory of small on large: Potential utility in computations
% of fluid-solid interactions in arteries. Computer methods in applied
% mechanics nad engineering, 196(31-32), pp.3070-3078.

if(~exist('activator'))
    activator=[1,1]; % selective activation of collagen and elastin --> if not given as input, is set to [1, 1]
end

if(~exist('not_nH_coeff'))
    not_nH_coeff=0.15;
end

I1 = (lambdatV*parameterV(6)).^2+(lambdazV*parameterV(7)).^2+(lambdarV/(parameterV(6)*parameterV(7))).^2; % First invariant

Lt = lambdatV.^2*parameterV(8)^2; % squared stretch of the circumferential fibre
Lz = lambdazV.^2*parameterV(8)^2; % squared stretch of the axial fibre

exponent = (parameterV(3)*(Lt>=1)+parameterV(5)*(Lt<1))...
    *parameterV(4)^2.*(Lt-1).^2+(parameterV(3)*(Lz>=1)+parameterV(5)*(Lz<1))*(1-parameterV(4)^2).*(Lz-1).^2;

sigmatt_rr = activator(1)*parameterV(1)*(1+not_nH_coeff)*((lambdatV*parameterV(6)).^2-(lambdarV/(parameterV(6)*parameterV(7))).^2).*(I1-3).^not_nH_coeff +...
    activator(2)*(2*parameterV(2)*parameterV(4)^2*(Lt-1).*Lz.*exp(exponent)); % circumferential Cauchy stress N/mm2
     
sigmazz_rr = activator(1)*parameterV(1)*(1+not_nH_coeff)*((lambdazV*parameterV(7)).^2-(lambdarV/(parameterV(6)*parameterV(7))).^2).*(I1-3).^not_nH_coeff +...
    activator(2)*(2*parameterV(2)*(1-parameterV(4)^2)*(Lz-1).*Lz.*exp(exponent)); % axial Cauchy stress N/mm2
    
W = activator(1)*(parameterV(1)/2)*(I1-3).^(1+not_nH_coeff) +...
    activator(2)*(parameterV(2)/(2*parameterV(4))*exp(exponent));

Ctttt = 2*(sigmatt_rr+activator(1)*parameterV(1)*(1+not_nH_coeff)*(I1-3).^not_nH_coeff.*(lambdarV/parameterV(6)*parameterV(7)).^2) +...
    activator(1)*2*(parameterV(6)*lambdatV).^4*parameterV(1)*(1+not_nH_coeff)*not_nH_coeff.*(I1-3).^(not_nH_coeff-1) +...
    activator(2)*(parameterV(2)*parameterV(4)^2*Lt.^4.*(1+2*(parameterV(3)*(Lt>=1)+parameterV(5)*(Lt<1))*parameterV(4)^2.*...
    (Lt-1).^2).*exp(exponent));

Czzzz = 2*(sigmazz_rr+activator(1)*activator(1)*parameterV(1)*(1+not_nH_coeff)*(I1-3).^not_nH_coeff.*(lambdarV/parameterV(6)*parameterV(7)).^2) +...
    activator(1)*2*(parameterV(7)*lambdazV).^4*parameterV(1)*(1+not_nH_coeff)*not_nH_coeff.*(I1-3).^(not_nH_coeff-1) +...
    activator(2)*(parameterV(2)*(1-parameterV(4))^2*Lz.^4.*(1+2*(parameterV(3)*(Lz>=1)+parameterV(5)*(Lz<1))*(1-parameterV(4))^2.*...
    (Lz-1).^2).*exp(exponent));

Czztt = 0; %2*activator(1)*(lambdatV*parameterV(10).*lambdazV*parameterV(11)).^2*parameterV(1)*(1+not_nH_coeff)*not_nH_coeff.*(I1-3).^(not_nH_coeff-1) +...
    % activator(2)*4*parameterV(4)*Lz*cos(parameterV(5))^2.*Lt*sin(parameterV(5))^2.*exp((parameterV(6)*(Ld>=1)+parameterV(9)*(Ld<1)).*(Ld-1).^2).*...
    % (1+2*(parameterV(6)*(Ld>=1)+parameterV(9)*(Ld<1)).*(Ld-1).^2);
end