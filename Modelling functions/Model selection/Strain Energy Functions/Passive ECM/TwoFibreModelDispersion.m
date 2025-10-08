function [sigmatt_rr,sigmazz_rr,W,Ctttt,Czzzz]=TwoFibreModelDispersion(parameterV,lambdatV,lambdazV,lambdarV,activator)
% The function returns the circumferential and axial stresses, the stored
% elastic energy, and the small-on-large circumferential and axial material
% stiffnesses for four-fibre family constitutive model given a parameter
% vector parameterV, the deformation in the three principal directions, and an activator 
% that activates/deactivates the contribution from each wal constituent.
%
% For a detailed desription of the constitutive model refer to Gasser et al. 
% Hyperelastic modelling of arterial layers with distributed collagen fibre orientations.  
% Journal of the Royal Society Interface, 3, pp.15-35, 2006.
%
% parameterV is the 5-element vector of model parameters where:
%   position 1: elastin stiffness-like parameter
%   position 2: diagonally oriented collagen fibre stiffness-like parameter
%   position 3: diagonally oriented collagen fibre non-linearity parameter
%   position 4: diagonally oriented collagen fibre orientation parameter
%   position 5: collagen fibre dispersion parameter
%   position 6: diagonally oriented collagen fibre non-linearity parameter for compression
%   position 7: elastin circumferential deposition stretch
%   position 8: elastin axial deposition stretch
%   position 9: collagen deposition stretch along the fibre direction
%
% For a detailed description of the small-on-large formulation refer to
% Baek et al. Theory of small on large: Potential utility in computations
% of fluid-solid interactions in arteries. Computer methods in applied
% mechanics nad engineering, 196(31-32), pp.3070-3078.

    if(~exist('activator'))
        activator=[1,1]; % selective activation of collagen and elastin --> if not given as input, is set to [1, 1]
    end
    
    Ic = (lambdatV*parameterV(7)).^2 + (lambdazV*parameterV(8)).^2 + (lambdarV/(parameterV(7)*parameterV(8))).^2; % First invariant of C for elastin
    Ic_f = lambdatV.^2 + lambdazV.^2 + lambdarV.^2;
    Id = lambdatV.^2*sin(parameterV(4))^2 + lambdazV.^2*cos(parameterV(4))^2; % Fourth invariant for collagen
    
    Ld = (parameterV(5)*Ic_f + (1-3*parameterV(5))*Id) * parameterV(9)^2; % Fibre stretch
    
    sigmatt_rr = activator(1)*2*parameterV(1)*((lambdatV*parameterV(7)).^2 - (lambdarV/(parameterV(7)*parameterV(8))).^2) +...
        activator(2)*4*parameterV(2)*(parameterV(5)*(lambdatV.^2-lambdarV.^2) + (1-3*parameterV(5))*(sin(parameterV(4))*lambdatV).^2)*parameterV(9)^2.*...
        (Ld-1).*exp((parameterV(3)*(Ld>=1) + parameterV(6)*(Ld<1)).*(Ld-1).^2);
        
    sigmazz_rr = activator(1)*2*parameterV(1)*((lambdazV*parameterV(8)).^2 - (lambdarV/(parameterV(7)*parameterV(8))).^2) +...
        activator(2)*4*parameterV(2)*(parameterV(5)*(lambdazV.^2-lambdarV.^2) + (1-3*parameterV(5))*(cos(parameterV(4))*lambdazV).^2)*parameterV(9)^2.*...
        (Ld-1).*exp((parameterV(3)*(Ld>=1) + parameterV(6)*(Ld<1)).*(Ld-1).^2);
    
    W = activator(1) * parameterV(1)*(Ic-3)+activator(2)*parameterV(2)/parameterV(3)*(exp((parameterV(3)*(Ld>=1) + parameterV(6)*(Ld<1)).*(Ld-1).^2)-1);
    
    Ctttt = 2*(sigmatt_rr + activator(1)*2*parameterV(1)*(lambdarV/(parameterV(7)*parameterV(8))).^2) +...
        activator(2)*8*(lambdatV*parameterV(9)).^4.*(1+2*(parameterV(3)*(Ld>=1)+parameterV(6)*(Ld<1)).*(Ld-1).^2)*parameterV(2)*...
        (parameterV(5)+(1-3*parameterV(5))*sin(parameterV(4))^2)^2.*exp((parameterV(3)*(Ld>=1)+parameterV(6)*(Ld<1)).*(Ld-1).^2);
    
    Czzzz = 2*(sigmazz_rr + activator(1)*2*parameterV(1)*(lambdarV/(parameterV(7)*parameterV(8))).^2) +...
        activator(2)*8*(lambdazV*parameterV(9)).^4.*(1+2*(parameterV(3)*(Ld>=1)+parameterV(6)*(Ld<1)).*(Ld-1).^2)*parameterV(2)*...
        (parameterV(5)+(1-3*parameterV(5))*cos(parameterV(4))^2)^2.*exp((parameterV(3)*(Ld>=1)+parameterV(6)*(Ld<1)).*(Ld-1).^2);

end