function [inputMatElasticFitting,weights_vector] = createInputMatrixElasticFitting(passive,pressure_sweep_list,axial_sweeps_check,reference_configuration)
% createInputMatrixElasticFitting creates the input matrix for the elastic
% fitting with the concatenated experimental data, containing:
%       Column 1: pressure [MPa]
%       Column 2: outer radius [mm]
%       Column 3: inner radius [mm]
%       Column 4: wall thickness [mm]
%       Column 5: transducer force [N]
%       Column 6: axial length [mm]
%       Column 7: circumferential stretch [-]
%       Column 8: axial stretch [-]
%       Column 9: radial stretch [–]
%       Column 10: circumferential stress [MPa]
%       Column 11: axial stress [MPa]
%       Column 12: wieghts vector [-]
%       Column 13: trigger (to separate the concatenated experiments) [-]
%
% Note that loading and unloading curves of each experiment are averaged to
% yield a unique pseudoelastic relationship.
%
% The trigger column is made of zeros except for ones in correspondence
% of the first datapoint of each experiment. 
% The weights column states the weight of each datapoint in the parameter
% estimation cost function.
%
% Function imputs are the structure "passive" with the experimentally
% measured passive arterial behaviour, "pressure_sweep_list" containing the 
% list of pressure-sweeps experiment to be included in the parameter
% estimation (according to previous modelling choices), and 
% "axial_sweeps_check" which is a logical variable stating whether
% axial-sweep experiments are included in the parameter estimation.
% "refMidWallRadius" and "refAxialLength" which are the mid wall radius 
% and the axial length in the reference configuration.

%% Scaling factors
    mtm = 0.001; %scaling factor to convert mircometer to mm

%% Geometrical parameters from unloaded intact configuration
    Ro = ((passive.OD)/2)*mtm;%outer radius (mm)
    H = (passive.H)*mtm; %wall thickness (mm)
    Ri = Ro-H; %inner radius (mm)
    L = passive.L; %unloaded intact length (mm)
    vessel_volume = pi*((Ro^2-Ri^2)*L); %volume of segment (mm^3)

%% Geometrical parameters from reference configuration
    refMidWallRadius = (reference_configuration(1)+reference_configuration(2))/2;
    refAxialLength = reference_configuration(4);

    [dataMat] = importElasticData(passive,pressure_sweep_list,axial_sweeps_check);

%% Converting experimental data into fitting data
    ri_qs = radiusFromIincompressibility(dataMat(:,2),dataMat(:,4),vessel_volume,1); % inner radii in loaded configurations (mm)
    h_qs = dataMat(:,2)-ri_qs; % thickness in loaded configurations (mm)
    rm_qs = (ri_qs+dataMat(:,2))/2;

% Converting structural variables (pressure,force,radius,length) into
% material variables (stresses and stretches)
    [materialMat] = structural2matrialVariables(dataMat(:,1),dataMat(:,3),rm_qs,h_qs,dataMat(:,4),refMidWallRadius,refAxialLength);

    separator = find(dataMat(:,5) == 1);
    mult_factor = 1; % defines the cumulative weight of the F-L experiments

    width_test = zeros(1,length(separator));
    weights_vector = zeros(1,length(separator));
    
    for i = 1:length(separator)
        if(i == length(separator))
            lambdatt_test = materialMat(separator(i):end,3);
            lambdazz_test = materialMat(separator(i):end,4);
        else
            lambdatt_test = materialMat(separator(i):separator(i+1)-1,3);
            lambdazz_test = materialMat(separator(i):separator(i+1)-1,4);
        end
        
        width_test(i) = sum(sqrt(diff(lambdatt_test).^2+diff(lambdazz_test).^2));
    end
    for i = 1:length(separator)
        if(i<=3)
            dataMat(separator(i):end,6) = 1;
            weights_vector(i) = 1;
        elseif(i == length(separator))
            dataMat(separator(i):end,6) = width_test(i)/sum(width_test(4:end));
            weights_vector(i) = width_test(i)/sum(width_test(4:end));
        else
            dataMat(separator(i):separator(i+1)-1,6) = width_test(i)/sum(width_test(4:end))*mult_factor;
            weights_vector(i) = width_test(i)/sum(width_test(4:end))*mult_factor;
        end
    end
    
    inputMatElasticFitting = [dataMat(:,1),dataMat(:,2),ri_qs,h_qs,dataMat(:,3),dataMat(:,4),materialMat(:,3),materialMat(:,4),materialMat(:,5),materialMat(:,1),materialMat(:,2),dataMat(:,6),dataMat(:,5)];
end