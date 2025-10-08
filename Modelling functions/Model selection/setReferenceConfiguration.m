function [l_ref,ri_ref,h_ref,ro_ref,lambdat_tf,lambdaz_tf,lambdar_tf] = setReferenceConfiguration(passive,in_vivo_ref)

% setReferenceConfiguration builds the reference configuration (i.e., the
% reference axial length (l_ref), inner radius (ri_ref), wall thickness
% (h_ref), and outer radius (ro_ref)) based on the users's modelling choices. 
% Additionally,setReferenceConfiguration calculates the stretches in the
% three principal directions for the traction-free (tf) configuration
% (lambdat_tf, lambdaz_tf and lambdar_tf correspond to the circumferential,
% axial and radial directions, respectively.
%
% Function imputs are the structure "passive" with the experimentally
% measured passive arterial behaviour and "in_vivo_ref" which is a logical
% variable determining whether the reference configuration corresponds to
% the unloaded configuration (=false) or to an in vivo homeostatic reference 
% (i.e., when pressure is 100 mmHg and the axial stretch is that in vivo).

%% Scaling factors
    mtm = 0.001; %scaling factor to convert mircometer to mm
    
%% Geometrical parameters from unloaded intact configuration
    Ro = ((passive.OD)/2)*mtm;%outer radius (mm)
    H = (passive.H)*mtm; %wall thickness (mm)
    Ri = Ro-H; %inner radius (mm)
    L = passive.L; %unloaded intact length (mm)
    
%% Setting reference configuration

    % Geometry
    if(in_vivo_ref)
        dM=(table2array(passive.static_data.loading.biaxial_Pd.PD_100)+flip(table2array(passive.static_data.unloading.biaxial_Pd.PD_100)))/2;
        [~,loc_ref]=min(abs(dM(:,2)-100));
        l_ref = dM(loc_ref,6);
        ri_ref = dM(loc_ref,7)/2*mtm; 
        h_ref = (dM(loc_ref,3)-dM(loc_ref,7))/2*mtm;
        ro_ref = dM(loc_ref,3)/2*mtm;
    else
        l_ref = L;
        ri_ref = Ri;
        h_ref = H;
        ro_ref = Ro;
    end
    
    % Define principal stretches in the unloaded configuration
    if(in_vivo_ref)
        lambdat_tf = (Ro-H/2)/(ri_ref+(h_ref/2));
        lambdaz_tf = (L)/l_ref;
        lambdar_tf = 1/lambdat_tf/lambdaz_tf;
    else
        % Set to 1 if the unloaded configuration is the reference configuration
        lambdat_tf = 1;
        lambdaz_tf = 1;
        lambdar_tf = 1;
    end
end