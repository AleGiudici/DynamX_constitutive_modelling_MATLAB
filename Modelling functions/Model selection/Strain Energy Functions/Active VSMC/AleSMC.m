function [sigmatt_rr,sigmazz_rr,W,Ctttt,Czzzz,Cttzz] = AleSMC(parameterV,lambdat,lambdaz,lambdar)
% This strain energy function models the smooth muscle contraction as a
% Gaussian-shaped function with principal contribution in the
% circumferential direction and in-plane (circumferential-axial) dispersion.
%
% parameterV is the 5-element vector of model parameters where:
%   position 1: VSMC peak stress parameter
%   position 2: VSMC dispersion coefficient (between 0 and 0.5)
%   position 3: width of the Gaussian function
%   position 4: VSMC optimal stretch in the circumferential direction
%   position 5: VSMC optimal stretch in the axial direction
%
% For a detailed description of the small-on-large formulation refer to
% Baek et al. Theory of small on large: Potential utility in computations
% of fluid-solid interactions in arteries. Computer methods in applied
% mechanics nad engineering, 196(31-32), pp.3070-3078.

    
    K = (1-2*parameterV(2))*(lambdat-parameterV(4))+parameterV(2)*(lambdat+lambdaz-parameterV(4)-parameterV(5));

    sigmatt_rr = parameterV(1)*lambdat*(1-parameterV(2)).*exp(-K.^2/(2*parameterV(3)^2));
    sigmazz_rr = parameterV(1)*lambdaz*parameterV(2).*exp(-K.^2/(2*parameterV(3)^2));

    W = parameterV(1)*sqrt(pi)/2/sqrt(-1/(2*parameterV(3)^2))*erfi(sqrt(-1/(2*parameterV(3)^2))*K);

    Ctttt = 2*parameterV(1)*lambdat*(1-parameterV(2)).*exp(-K.^2/(2*parameterV(3)^2)).*...
        (1-(1-parameterV(2))*K.*lambdat/(2*parameterV(3)^2));

    Czzzz = 2*parameterV(1)*lambdaz*parameterV(2).*exp(-K.^2/(2*parameterV(3)^2)).*...
        (1-parameterV(2)*K.*lambdaz/(2*parameterV(3)^2));

    Cttzz = -2*parameterV(1)*(1-parameterV(2))*parameterV(2)*lambdat.*lambdaz.*K/(2*parameterV(3)^2).*...
        exp(-K.^2/(2*parameterV(3)^2));
end  