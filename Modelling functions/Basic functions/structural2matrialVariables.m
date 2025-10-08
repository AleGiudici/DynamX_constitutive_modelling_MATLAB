function [materialMat] = structural2matrialVariables(pressureV,transducerForceV,midWallRadiusV,thicknessV,axialLengthV,refMidWallRadius,refAxialLength)
% This function takes structural vessel variables (i.e., pressure,
% transducer axial force, mid-wall radius and wall thickness) and
% geometrical features in the reference configuration (i.e., mid-wall
% radius and axial length), and returns a matrix of material variables 
% (i.e., circumferential and axial stresses and circumferential, axial and
% radial stretches).

    sigmatt_rr = pressureV.*midWallRadiusV./thicknessV; % experimental circumferential extra stress (Laplace's law) (N/mm^2)
    sigmazz_rr = (transducerForceV./(pi*midWallRadiusV.*thicknessV)+sigmatt_rr)/2; % experimental axial extra stress (Laplace's law) (N/mm^2)
    
    sigmatt = sigmatt_rr-pressureV/2;
    sigmazz = sigmazz_rr-pressureV/2;

    lambdatt_qs = midWallRadiusV/refMidWallRadius; % circumferential stretch
    lambdazz_qs = axialLengthV/refAxialLength; % Axial stretch
    lambdarr_qs = 1./(lambdatt_qs.*lambdazz_qs); % radial stretch (to be determined from incompressibility condition)
    
    materialMat = [sigmatt,sigmazz,lambdatt_qs,lambdazz_qs,lambdarr_qs];
end