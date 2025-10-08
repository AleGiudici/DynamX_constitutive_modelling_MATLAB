function [inputMatViscoelasticFitting1,inputMatViscoelasticFitting2] = createInputStructuresViscoelasticFitting(passive,pressure_sweep_list,axial_sweeps_check,reference_configuration)
% createInputStructuresViscoelasticFitting creates the 2 input structures that are
% used in the 2 viscoelastic parameter identification steps. Each structure
% contains a list of experimental tet matrices, in which:
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
%       Column 13:  trigger (to separate the concatenated experiments) [-]
%
% The trigger column is made of zeros except for ones in correspondence
% of the first datapoint of each experiment. In harmonic loading
% experiments, each cycle is identified by a "2" in the trigger column and
% the last cycle is indicated by a "1".
%
% Function imputs are the structure "passive" with the experimentally
% measured passive arterial behaviour, "pressure_sweep_list" containing the 
% list of pressure-sweeps experiment to be included in the parameter
% estimation (according to previous modelling choices), and 
% "axial_sweeps_check" which is a logical variable stating whether
% axial sweep experiments are included in the parameter estimation.
% "reference configuration" contains information on the vessel geometry in the
% reference configuration (i.e., inner radius, outer radius, wall thickness
% and axial length).

%% Scaling factors
    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    mtm = 0.001; %scaling factor to convert mircometer to mm
    gtN = 0.0098; %scaling factor to convert g to N
    conversion_vector = [1; mtN; mtm/2; gtN; 1; 1; mtm/2];

%% Unloaded geometry
    Ro = passive.OD/2*mtm;
    H = passive.H*mtm;
    L = passive.L;
    vesselVolume = pi*(Ro^2-(Ro-H)^2)*L;

%% Reference configuration
    refMidWallRadius = (reference_configuration(1)+reference_configuration(2))/2;
    refAxialLength = reference_configuration(4);
   
    for k=1:length(pressure_sweep_list) %for k=1:n; where n = number of PD tests

        dataMat = zeros(size(table2array(passive.static_data.loading.biaxial_Pd.(pressure_sweep_list{k})),1)+...
            size(table2array(passive.static_data.loading.biaxial_Pd.(pressure_sweep_list{k})),1)+1,...
            size(table2array(passive.static_data.loading.biaxial_Pd.(pressure_sweep_list{k})),2));

        dataMat(2:end,:) = [table2array(passive.static_data.loading.biaxial_Pd.(pressure_sweep_list{k}));...
            table2array(passive.static_data.unloading.biaxial_Pd.(pressure_sweep_list{k}))].*conversion_vector';

        locator_NaN = ~isnan(dataMat(:,3)); % locate video tracking issues
        dataMat = dataMat(locator_NaN,:);
               
        triggerV = zeros(size(dataMat,1),1);
        triggerV(2) = 1;

        dataMat(2:end,1) = dataMat(2:end,1)+500;
        dataMat(1,:) = [0 0 Ro 0 1 L Ro-H-H];
       
% Converting experimental data into fitting data
        ri = radiusFromIincompressibility(dataMat(:,3),dataMat(:,6),vesselVolume,1); % inner radii in loaded configurations (mm)
        h = dataMat(:,3)-ri; % thickness in loaded configurations (mm)
        rm = (ri+dataMat(:,3))/2;

% Converting structural variables (pressure,force,radius,length) into
% material variables (stresses and stretches)
        [materialMat] = structural2matrialVariables(dataMat(:,2),dataMat(:,4),rm,h,dataMat(:,6),refMidWallRadius,refAxialLength);

% Creating data structure
        inputMatViscoelasticFitting1.(pressure_sweep_list{k}).processed_data = [dataMat(:,1), dataMat(:,2) dataMat(:,3) ri, h, dataMat(:,4), dataMat(:,6),...
            materialMat(:,3), materialMat(:,4) materialMat(:,5) materialMat(:,1) materialMat(:,2), triggerV];

        if(strcmp(pressure_sweep_list{k},'PD_100'))
            inputMatViscoelasticFitting2.(pressure_sweep_list{k}).processed_data = [dataMat(:,1), dataMat(:,2) dataMat(:,3) ri, h, dataMat(:,4), dataMat(:,6),...
                materialMat(:,3), materialMat(:,4) materialMat(:,5) materialMat(:,1) materialMat(:,2), triggerV];
        end
    end

    if(axial_sweeps_check)
        FL = passive.static_data.FL;

        for k=1:length(FL) %for k=1:n; where n = number of PD tests

            dataMat = zeros(size(table2array(passive.static_data.loading.biaxial_Fl.(FL{k})),1)+...
                size(table2array(passive.static_data.loading.biaxial_Fl.(FL{k})),1)+1,...
                size(table2array(passive.static_data.loading.biaxial_Fl.(FL{k})),2));
    
            dataMat(2:end,:) = [table2array(passive.static_data.loading.biaxial_Fl.(FL{k}));...
                table2array(passive.static_data.unloading.biaxial_Fl.(FL{k}))].*conversion_vector';
    
            locator_NaN = ~isnan(dataMat(:,3)); % locate video tracking issues
            dataMat = dataMat(locator_NaN,:);
                   
            triggerV = zeros(size(dataMat,1),1);
            triggerV(2) = 1;
    
            dataMat(2:end,1) = dataMat(2:end,1)+500;
            dataMat(1,:) = [0 0 Ro 0 1 L Ro-H-H];
           
    % Converting experimental data into fitting data
            ri = radiusFromIincompressibility(dataMat(:,3),dataMat(:,6),vesselVolume,1); % inner radii in loaded configurations (mm)
            h = dataMat(:,3)-ri; % thickness in loaded configurations (mm)
            rm = (ri+dataMat(:,3))/2;
    
    % Converting structural variables (pressure,force,radius,length) into
    % material variables (stresses and stretches)
            [materialMat] = structural2matrialVariables(dataMat(:,2),dataMat(:,4),rm,h,dataMat(:,6),refMidWallRadius,refAxialLength);
    
    % Creating data structure
            inputMatViscoelasticFitting1.(FL{k}).processed_data = [dataMat(:,1), dataMat(:,2) dataMat(:,3) ri, h, dataMat(:,4), dataMat(:,6),...
                materialMat(:,3), materialMat(:,4) materialMat(:,5) materialMat(:,1) materialMat(:,2), triggerV];
        end
    end

    PD_dyna = passive.dynamic_data.PD(4:end);

    for k=1:length(PD_dyna) %for k=1:n; where n = number of PD tests
        locator_name = find(PD_dyna{k} == '_');
        HR = str2double(PD_dyna{k}(locator_name(end-1)+1:locator_name(end)-1))/1000;

        dataMat = ensambleAveragingDynamicData(table2array(passive.dynamic_data.dynamic_PD.(PD_dyna{k})).*conversion_vector',HR);

        triggerV = zeros(size(dataMat,1),1);
        triggerV(1) = 2;

        concatenationMat = [dataMat triggerV];

        concatenated_dataMat = [0 0 Ro 0 1 L Ro-H 0; concatenationMat; concatenationMat; concatenationMat; concatenationMat];

        concatenated_dataMat(2:end,1) = concatenated_dataMat(2:end,1)+500;
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
        inputMatViscoelasticFitting2.(PD_dyna{k}).processed_data = [concatenated_dataMat(:,1), concatenated_dataMat(:,2) concatenated_dataMat(:,3) ri,...
            h, concatenated_dataMat(:,4), concatenated_dataMat(:,6), materialMat(:,3), materialMat(:,4), materialMat(:,5), materialMat(:,1),...
            materialMat(:,2), concatenated_dataMat(:,end)];        
    end
end