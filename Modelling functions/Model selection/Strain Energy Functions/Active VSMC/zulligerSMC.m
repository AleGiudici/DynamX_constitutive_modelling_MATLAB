function [sigmatt_rr,sigmazz_rr,W,Ctttt,Czzzz] = zulligerSMC(parameterV,lambdat,lambdaz,lambdar)
% This strain energy function implements Zulliger's VSMC model.
% VSMCs are assumed to contribute to the circumferential response only.
%
% parameterV is the 2-element vector of model parameters where:
%   position 1: VSMC stiffness-like parameter
%   position 2: VSMC deposition stretch in the circumferential direction
%
% For a detailed description of the model, please refer to Zulliger et al.
% ﻿﻿A constitutive formulation of arterial mechanics including vascular
% smooth muscle tone . American Journal of Physiology - Heart and
% Circulatory Physiology 2004, 287(2 56-3): 1335-1343.
%
% For a detailed description of the small-on-large formulation refer to
% Baek et al. Theory of small on large: Potential utility in computations
% of fluid-solid interactions in arteries. Computer methods in applied
% mechanics nad engineering, 196(31-32), pp.3070-3078.

    sigmatt_rr = parameterV(1)*(lambdat*parameterV(2)-1);
    sigmazz_rr = 0;
    
    W = parameterV(1)*(lambdat*parameterV(2)-log(lambdat*parameterV(2))-1);
    
    Ctttt = 2*sigmatt_rr+parameterV(1);
    Czzzz = 0;
end  