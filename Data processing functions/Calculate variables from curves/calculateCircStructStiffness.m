function [haemodynamicVariables] = calculateCircStructStiffness(dataMat,static_P_ref_res,static_ri_ref_res)
% This function calculates circumferential structural stiffness parameter
% (i.e., PWV, area distensibility and area compliance) from the
% experimental data. Calculations are performed on both the dynamic pressure-
% diameter loops and in matching pressure ranges of the quasi-static
% pressure-diameter curves. PWV is [m/s], distensibility is [1/MPa] and
% compliance is [mm^4/Pa].

    Mtm = 1000;
    bloodDensity = 1050;

    locator = find(~isnan(dataMat(:,3)));

    [max_P,loc_max_P] = max(dataMat(locator,2));
    [min_P,loc_min_P] = min(dataMat(locator,2));

    max_ri = dataMat(locator(1)+loc_max_P-1,3)/Mtm/2;
    min_ri = dataMat(locator(1)+loc_min_P-1,3)/Mtm/2;

    PWV = min_ri*sqrt((max_P-min_P)*133.32/(max_ri^2-min_ri^2)/bloodDensity);

    Distensibility = 1/(bloodDensity*PWV^2)*10^6;
    Compliance = (pi*min_ri^2)*Distensibility;

    [~, loc_max_P_static] = min(abs(max_P-static_P_ref_res));
    [~, loc_min_P_static] = min(abs(min_P-static_P_ref_res));
    
    max_P_static = static_P_ref_res(loc_max_P_static);
    min_P_static = static_P_ref_res(loc_min_P_static);
    max_ri_static = static_ri_ref_res(loc_max_P_static)/Mtm;
    min_ri_static = static_ri_ref_res(loc_min_P_static)/Mtm;

    PWV_static = min_ri_static*sqrt((max_P_static-min_P_static)*133.32/(max_ri_static^2-min_ri_static^2)/bloodDensity);
    Distensibility_static = 1/(bloodDensity*PWV_static^2)*10^6;
    Compliance_static = (pi*min_ri_static^2)*Distensibility_static;

    haemodynamicVariables(1:8,1) = [max_P; min_P; PWV; PWV_static; Distensibility;...
        Distensibility_static; Compliance; Compliance_static];
end