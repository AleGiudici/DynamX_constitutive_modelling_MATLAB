function [outputMat] = importElasticData(passive,pressure_sweep_list,axial_sweeps_check)
% importElasticData concatenates the data from all experiments into a 
% matrix [pressure | outer radius | transducer axial force | axial length |
% trigger (experiment separator) | fitting weights]. 
%
% Note that loading and unloading curves of each experiment are averaged to
% yield a unique pseudoelastic relationship.
%
% The first four columns contain the actual experimental data. 
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
    
%% Scaling factors
    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    mtm = 0.001; %scaling factor to convert mircometer to mm
    gtN = 0.0098; %scaling factor to convert g to N

%% Initialising vectors
    pressureV = zeros(1000,1);
    outerRadiusV = zeros(1000,1);
    transducerForceV = zeros(1000,1);
    axialLengthV = zeros(1000,1);
    triggerV = zeros(1000,1);
    weightsV = zeros(1000,1);

    counter_1 = 1;
    counter_2 = 0;

    for k=1:length(pressure_sweep_list) %for k=1:n; where n = number of PD tests
    
        dataMat = (table2array(passive.static_data.unloading.biaxial_Pd.(pressure_sweep_list{k}))+flip(table2array(passive.static_data.loading.biaxial_Pd.(pressure_sweep_list{k}))))/2;
        
        locator_NaN = ~isnan(dataMat(:,3)); % locate video tracking issues
        dataMat = dataMat(locator_NaN,:);

        counter_2 = counter_2+size(dataMat,1); 
    
        trigger = zeros(size(dataMat,1),1);
        trigger(1) = 1;
        
        pressureV(counter_1:counter_2) = dataMat(:,2)*mtN;
        outerRadiusV(counter_1:counter_2) = dataMat(:,3)/2*mtm;
        transducerForceV(counter_1:counter_2) = dataMat(:,4)*gtN;
        axialLengthV(counter_1:counter_2) = dataMat(:,6);
        triggerV(counter_1:counter_2) = trigger;    
        weightsV(counter_1:counter_2) = ones(size(dataMat,1),1); % Pressure-diameter sweeps have weight 1

        counter_1 = counter_2+1;
    end
    
    if(axial_sweeps_check)
        j = 0;
        FL = passive.static_data.FL;
        
        for k=length(pressure_sweep_list)+1:length(pressure_sweep_list)+length(FL) %for k=1:n; where n = number of FL tests
            j = j+1;
    
            dataMat = (table2array(passive.static_data.unloading.biaxial_Fl.(FL{j}))+flip(table2array(passive.static_data.loading.biaxial_Fl.(FL{j}))))/2;
            
            locator_NaN = ~isnan(dataMat(:,3)); % locate video tracking issues
            dataMat = dataMat(locator_NaN,:);

            counter_2 = counter_2+size(dataMat,1); 
    
            trigger = zeros(size(dataMat,1),1);
            trigger(1) = 1;
        
            check = prod(isnan(dataMat(:,3))); % check if data was acquired with dynamice s1 (i.e., missing diameter recording
            % in FL-sweeps.
    
            if((~check) && (size(pressure_sweep_list,2)~=1))
                pressureV(counter_1:counter_2) = dataMat(:,2)*mtN;
                outerRadiusV(counter_1:counter_2) = dataMat(:,3)/2*mtm;
                transducerForceV(counter_1:counter_2) = dataMat(:,4)*gtN;
                axialLengthV(counter_1:counter_2) = dataMat(:,6);
                triggerV(counter_1:counter_2) = trigger;
                weightsV(counter_1:counter_2) = length(pressure_sweep_list)/length(FL)/5*ones(size(dataMat,1),1); % The cumulative weight of all force-length sweeps is equal 
                % to that of a single pressure-diameter sweep.
            end

            counter_1 = counter_2+1;
        end
    end
    
    outputMat = [pressureV,outerRadiusV,transducerForceV,axialLengthV,triggerV,weightsV];
    outputMat = outputMat(1:counter_2,:);
end