function [sigmatt_rr,sigmazz_rr,W,Ctttt,Czzzz] = rachevSMC(parameterV,lambdat,lambdaz,lambdar)
% This strain energy function implements Rachev's parabolic VSMC model.
% VSMCs are assumed to contribute to the circumferential response only.
%
% parameterV is the 3-element vector of model parameters where:
%   position 1: VSMC peak stress parameter
%   position 2: VSMC optimal stretch in the circumferential direction
%   position 3: controls the width of the parabola
%
% For a detailed description of the model, please refer to Rachev et al.
% ﻿Theoretical study of the effects of vascular smooth muscle contraction
% on strain and stress sistributions in arteries. Annals of Biomedical
% Engineering 1999, 27(4): 459-468.
%
% For a detailed description of the small-on-large formulation refer to
% Baek et al. Theory of small on large: Potential utility in computations
% of fluid-solid interactions in arteries. Computer methods in applied
% mechanics nad engineering, 196(31-32), pp.3070-3078.

    sigmatt_rr = parameterV(1)*lambdat.*(1-((parameterV(2)-lambdat)/(parameterV(2)-parameterV(3)).^2)).*...
        (lambdat >= parameterV(2)-parameterV(3) & lambdat <= parameterV(2)+parameterV(3));
    sigmazz_rr = 0;

    W = 0;
    Ctttt = 2*sigmatt_rr+2*parameterV(1)*lambdat.^2.*(parameterV(2)-lambdat)/(parameterV(2)-parameterV(3)).*...
        (lambdat >= parameterV(2)-parameterV(3) & lambdat <= parameterV(2)+parameterV(3));
    Czzzz = 0;
end  