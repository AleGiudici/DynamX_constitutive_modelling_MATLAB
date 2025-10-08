function [inputMatActiveFitting,n_dataPoints,reference_configuration] = createInputStructureActiveFitting(experimental_data,reference_configuration)
% createInputStructureActiveFitting creates the input structure that is
% active VSMC parameter identification step. The structure is divided in
% a number of substructures corresponding to the number of contraction
% states. At least two contraction states are needed (i.e., a reference
% uncontracted state and a contracted state).
% Each substrcuture contains a list of experimental test matrices, in which:
%       Column 1: time [s]
%       Column 2: pressure [MPa]
%       Column 3: outer radius [mm]
%       Column 4: inner radius [mm]
%       Column 5: wall thickness [mm]
%       Column 6: transducer force [N]
%       Column 7: axial length [mm]
%       Column 8: circumferential stretch [-]
%       Column 9: axial stretch [-]
%       Column 10: radial stretch [–]
%       Column 11: circumferential stress [MPa]
%       Column 12: axial stress [MPa]
%       Column 13: trigger (to separate the concatenated experiments) [-]
%
% The trigger column is made of zeros except for ones in correspondence
% of the first datapoint of each experiment. In harmonic loading
% experiments, each cycle is identified by a "2" in the trigger column and
% the last cycle is indicated by a "1".
%
% Function inputs are the structure "experimental data" with the experimentally
% measured passive and active arterial behaviours and "reference
% configuration" which contains information on the vessel geometry in the
% reference configuration (i.e., inner radius, outer radius, wall thickness
% and axial length).

%% Scaling factors
    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    mtm = 0.001; %scaling factor to convert mircometer to mm
    gtN = 0.0098; %scaling factor to convert g to N
    conversion_vector = [1; mtN; mtm/2; gtN; 1; 1; mtm/2];

%% Unloaded geometry
    Ro = experimental_data.passive.OD/2*mtm;
    H = experimental_data.passive.H*mtm;
    L = experimental_data.passive.L;
    vesselVolume = pi*(Ro^2-(Ro-H)^2)*L;

%% Reference configuration
    contraction_states = fieldnames(experimental_data.active);
    
    n_dataPoints = 0;

    for iState = 1:length(contraction_states)
        test_list = fieldnames(experimental_data.active.(contraction_states{iState,1}));
        
        dataMat = zeros(size(table2array(experimental_data.active.(contraction_states{iState,1}).(test_list{1,1})),1)+1,...
            size(table2array(experimental_data.active.(contraction_states{iState,1}).(test_list{1,1})),2));

        dataMat(2:end,:) = table2array(experimental_data.active.(contraction_states{iState,1}).PD_100).*conversion_vector';

        if(iState == 1)
            [~,loc_max] = max(dataMat(2:end,2));
            [~,loc_ref_loading] = min(abs(dataMat(2:loc_max,2)-100));
            [~,loc_ref_unloading] = min(abs(dataMat(loc_max:end,2)-100));

            reference_configuration(2) = (dataMat(loc_ref_loading,3)+dataMat(loc_ref_unloading+loc_max-1,3))/2;
            reference_configuration(4) = (dataMat(loc_ref_loading,6)+dataMat(loc_ref_unloading+loc_max-1,6))/2;

            reference_configuration(1) = radiusFromIincompressibility(reference_configuration(2),reference_configuration(4),vesselVolume,1);
            reference_configuration(3) = reference_configuration(2)-reference_configuration(1);

            refMidWallRadius = (reference_configuration(1)+reference_configuration(2))/2;
            refAxialLength = reference_configuration(4);
        end

        n_dataPoints = n_dataPoints+size(dataMat(2:end,:),1);

        locator_NaN = ~isnan(dataMat(:,3)); % locate video tracking issues
        dataMat = dataMat(locator_NaN,:);
               
        triggerV = zeros(size(dataMat,1),1);
        triggerV(2) = 1;

        dataMat(2:end,1) = dataMat(2:end,1)+5000-dataMat(2,1);
        
        dataMat(1,:) = [0 0 Ro 0 1 L Ro-H-H];
    
% Converting experimental data into fitting data
        ri = radiusFromIincompressibility(dataMat(:,3),dataMat(:,6),vesselVolume,1); % inner radii in loaded configurations (mm)
        h = dataMat(:,3)-ri; % thickness in loaded configurations (mm)
        rm = (ri+dataMat(:,3))/2;
    
% Converting structural variables (pressure,force,radius,length) into
% material variables (stresses and stretches)
        [materialMat] = structural2matrialVariables(dataMat(:,2),dataMat(:,4),rm,h,dataMat(:,6),refMidWallRadius,refAxialLength);
    
% Creating data structure
        inputMatActiveFitting.(contraction_states{iState,1}).(test_list{1,1}).processed_data = [dataMat(:,1), dataMat(:,2) dataMat(:,3) ri, h, dataMat(:,4), dataMat(:,6),...
            materialMat(:,3), materialMat(:,4) materialMat(:,5) materialMat(:,1) materialMat(:,2), triggerV];

        for iTest = 2:length(test_list)
            locator_name = find(test_list{iTest} == '_');
            HR = str2double(test_list{iTest}(locator_name(2)+1:end))/1000;
    
            dataMat = ensambleAveragingDynamicData(table2array(experimental_data.active.(contraction_states{iState}).(test_list{iTest})).*conversion_vector',HR);
            
            n_dataPoints = n_dataPoints+size(dataMat,1);

            triggerV = zeros(size(dataMat,1),1);
            triggerV(1) = 2;
    
            concatenationMat = [dataMat triggerV];
    
            concatenated_dataMat = [0 0 Ro 0 1 L Ro-H 0; concatenationMat; concatenationMat; concatenationMat; concatenationMat];
    
            concatenated_dataMat(2:end,1) = concatenated_dataMat(2:end,1)+5000;
            concatenated_dataMat(2:end,1) = find(concatenated_dataMat(2:end,1)>0).*(concatenated_dataMat(3,1)-concatenated_dataMat(2,1))+500;
       
    % Converting experimental data into fitting data
            ri = radiusFromIincompressibility(concatenated_dataMat(:,3),concatenated_dataMat(:,6),vesselVolume,1); % inner radii in loaded configurations (mm)
            h = concatenated_dataMat(:,3)-ri; % thickness in loaded configurations (mm)
            rm = (ri+concatenated_dataMat(:,3))/2;
        
    % Converting structural variables (pressure,force,radius,length) into
    % material variables (stresses and stretches)
            [materialMat] = structural2matrialVariables(concatenated_dataMat(:,2),concatenated_dataMat(:,4),rm,h,concatenated_dataMat(:,6),refMidWallRadius,refAxialLength);
        
            locator2inTrigger = find(concatenated_dataMat(:,end) == 2);
            concatenated_dataMat(locator2inTrigger(end),end) = 1;
    
    % Creating data structure
            inputMatActiveFitting.(contraction_states{iState,1}).(test_list{iTest}).processed_data = [concatenated_dataMat(:,1), concatenated_dataMat(:,2) concatenated_dataMat(:,3) ri,...
                h, concatenated_dataMat(:,4), concatenated_dataMat(:,6), materialMat(:,3), materialMat(:,4), materialMat(:,5), materialMat(:,1),...
                materialMat(:,2), concatenated_dataMat(:,end)]; 
        end
    end
end