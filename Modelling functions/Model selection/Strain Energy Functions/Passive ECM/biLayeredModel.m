function [sigmatt_rr,sigmazz_rr, W, Ctttt, Czzzz]=biLayeredModel(parameterV,lambdatV,lambdazV,lambdarV,activator,not_nH_coeff)

if(~exist('activator'))
    activator=[1,1]; % selective activation of collagen and elastin --> if not given as input, is set to [1, 1]
end

if(~exist('not_nH_coeff'))
    not_nH_coeff=0.15;
end

% parameterV is the 12-element vector of model parameters where:
%   position 1: elastin stiffness-like parameter
%   position 2: medial collagen fibre stiffness-like parameter
%   position 3: medial collagen fibre non-linearity parameter
%   position 4: medial collagen fibre angle
%   position 5: medial collagen fibre distribution
%   position 6: adventitial collagen fibre stiffness-like parameter
%   position 7: adventitial collagen fibre non-linearity parameter multiplier
%   position 8: adventitial collagen fibre angle multiplier
%   position 9: adventitial collagen fibre distribution
%   position 10: collagen compression non-linearity parameter
%   position 11: elastin circumferential deposition stretch
%   position 12: elastin axial deposition stretch
%   position 13: medial collagen deposition stretch along the fibre direction
%   position 14: adventitial collagen deposition stretch along the fibre direction

% parameterV(7) = parameterV(7)*parameterV(3);
% parameterV(8) = parameterV(4)*parameterV(8);

I1e = (lambdatV*parameterV(11)).^2+(lambdazV*parameterV(12)).^2+(lambdarV/(parameterV(11)*parameterV(12))).^2; % First invariant
I1c = lambdatV.^2+lambdazV.^2+lambdarV.^2; % First invariant
I4m = (lambdatV.^2*sin(parameterV(4))^2+lambdazV.^2*cos(parameterV(4))^2)*parameterV(13)^2;
I4a = (lambdatV.^2*sin(parameterV(8))^2+lambdazV.^2*cos(parameterV(8))^2)*parameterV(14)^2;

Ldm = parameterV(5)*(I1c-3)*parameterV(13)^2+(1-3*parameterV(5))*I4m;
Lda = parameterV(9)*(I1c-3)*parameterV(14)^2+(1-3*parameterV(9))*I4a;

sigmatt_rr = activator(1)*(parameterV(1)*(1+not_nH_coeff)*(lambdatV.^2*parameterV(11)^2-lambdarV.^2/(parameterV(11)^2*parameterV(12)^2)).*(I1e-3).^not_nH_coeff) +...
    activator(2)*(2*parameterV(2)*(Ldm-1).*(parameterV(5)*(lambdatV.^2-lambdarV.^2)*parameterV(13)^2+(1-3*parameterV(5))*lambdatV.^2*parameterV(13)^2*sin(parameterV(4))^2)...
    .*exp((parameterV(3)*(Ldm>=1)+parameterV(10)*(Ldm<1)).*(Ldm-1).^2)) +...
    activator(2)*(2*parameterV(6)*(Lda-1).*(parameterV(9)*(lambdatV.^2-lambdarV.^2)*parameterV(14)^2+(1-3*parameterV(9))*lambdatV.^2*parameterV(14)^2*sin(parameterV(8))^2)...
    .*exp((parameterV(7)*(Lda>=1)+parameterV(10)*(Lda<1)).*(Lda-1).^2));

sigmazz_rr = activator(1)*(parameterV(1)*(1+not_nH_coeff)*(lambdazV.^2*parameterV(11)^2-lambdarV.^2/(parameterV(11)^2*parameterV(12)^2)).*(I1e-3).^not_nH_coeff) +...
    activator(2)*(2*parameterV(2)*(Ldm-1).*(parameterV(5)*(lambdazV.^2-lambdarV.^2)*parameterV(13)^2+(1-3*parameterV(5))*lambdazV.^2*parameterV(13)^2*cos(parameterV(4))^2)...
    .*exp((parameterV(3)*(Ldm>=1)+parameterV(10)*(Ldm<1)).*(Ldm-1).^2)) +...
    activator(2)*(2*parameterV(6)*(Lda-1).*(parameterV(9)*(lambdazV.^2-lambdarV.^2)*parameterV(14)^2+(1-3*parameterV(9))*lambdazV.^2*parameterV(14)^2*cos(parameterV(8))^2)...
    .*exp((parameterV(7)*(Lda>=1)+parameterV(10)*(Lda<1)).*(Lda-1).^2));

Ctttt = 0;
Czzzz = 0;
W = 0;

end


