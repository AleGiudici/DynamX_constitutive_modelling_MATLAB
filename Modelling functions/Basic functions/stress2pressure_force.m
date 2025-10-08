function [pressureV,transducerForceV] = stress2pressure_force(sigmatt_rr,sigmazz_rr,midWallRadiusV,thicknessV)
% Given the circumferential and axial extra stresses, the deformed mid-wall 
% radius and wall thickness, this function returns the transmural pressure 
% in mmHg and the transducer axial force in grams.

    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    gtN = 0.0098; %scaling factor to convert g to N

    pressureV = sigmatt_rr.*thicknessV./midWallRadiusV/mtN;
    transducerForceV = (2*sigmazz_rr-sigmatt_rr)*pi.*thicknessV.*midWallRadiusV/gtN;

end