function [dataMat] = importDataDynamX2(app,directoryName,dataType,test_name,passive,forceLengthCheck,dynamicCheck)
    
    Mtm = 1000; % mm to micrometer
    mN2g = 1/9.81; % mN to g comversion
    vesselVolume = pi*(((passive.OD/2/Mtm)^2-((passive.OD-passive.H)/2/Mtm)^2)*passive.L); % Volume of unloaded geometry [mm^3]
    
%% Establish the file path
    path = findFilePathDynamX2(app,directoryName,forceLengthCheck,dynamicCheck,test_name);

%% Calibration coefficient
    if(strcmp(app.protocol_menu.Value,'Custom'))
        load(app.load_data_dat_proc.UserData.file_protocol)
        V2mN = protocolMat.V2mN;
        p2um = protocolMat.p2um;
    else
        % MOUSE CAROTID
        % V2mN = 29.5; % Volt to mN conversion for mouse carotids
        % p2um = 1.05; % pixel to um conversion for mouse carotids
        
        % % RAT AORTA
        % V2mN = 28.08; % Volt to mN conversion for rat aorta
        % p2um = 2.712; % pixel to um conversion for rat aorta
        
        % % MOUSE AORTA
        V2mN = 37.3; % Volt to mN conversion for rat aorta
        p2um = 1.638; % pixel to um conversion for rat aorta
    end

%% First attempt at uploading (internally) preprocessed data
    check = 0;
    clear exported_data extracted_data
    if(strcmp(dataType,'PD_static') || strcmp(dataType,'PD_dynamic'))
        if(isfile([fullfile(path,test_name) '.mat']))
            load([fullfile(path,test_name) '.mat'],"exported_data");
            extracted_data = table2array(exported_data);
            check = 1;
        end
    elseif(strcmp(dataType,'FL_static'))
        if(isfile([fullfile(path,test_name) '.mat']))
            load([fullfile(path,test_name) '.mat'],"exported_data");
            extracted_data = table2array(exported_data);
            check = 1;
        end
    end

    if(check == 0)
%% Second attempt at uploading (externally) preprocessed data
        if((ismac && isfile([path '/Data.mat'])) || (~ismac && isfile([path '\Data.mat'])))
            load([path '/Data.mat']);
            dataMat = matfile([path '/Data.mat']);
            name_fields = fieldnames(dataMat);

            timeDAQV = eval([name_fields{2,1} '.time_pressure']);
            pressure1V = eval([name_fields{2,1} '.p1']);
            pressure2V = eval([name_fields{2,1} '.p2']);
            transducerForceV = eval([name_fields{2,1} '.load']);
            stretchV = eval([name_fields{2,1} '.stretch']);
            interpolatedOuterDiameter = eval([name_fields{2,1} '.D']);
            timeMotorsV = eval([name_fields{2,1} '.time_stretch']);
            
            % interpolate motor data to daq data
            transducerForceV=interp1(timeMotorsV,transducerForceV,timeDAQV,'linear','extrap');
            interpolatedStretchV=interp1(timeMotorsV,stretchV,timeDAQV,'linear','extrap');
        else
%% Preprocessed data not found --> processing data
            [timeDAQV,pressure1V,pressure2V,interpolatedOuterDiameter,transducerForceV,interpolatedStretchV] = processDynamX2Data(path);
        end
        if(size(interpolatedOuterDiameter,2)==1)
            extracted_data = [timeDAQV, pressure1V, pressure2V, interpolatedOuterDiameter*p2um, transducerForceV*V2mN*mN2g, interpolatedStretchV];    
        else
            extracted_data = [timeDAQV, pressure1V, pressure2V, interpolatedOuterDiameter'*p2um, transducerForceV*V2mN*mN2g, interpolatedStretchV];    
        end
    end

%% Resampling experimental data
    if(strcmp(app.protocol_menu.Value,'Old'))
        [locatorPeak] = find(extracted_data(:,2)>200);
        [locatorBottom] = find(extracted_data(:,2)<=11);
    
        if((size(locatorPeak,1)~=0) && (size(locatorBottom,1)~=0))
            speedLocatorBottom = (locatorBottom(2:end)-locatorBottom(1:end-1));
            [~,jump] = max(speedLocatorBottom);
            
            loadingData = extracted_data(jump:locatorPeak(1),:);
            resampling_time = linspace(loadingData(1,1),loadingData(end,1),35);
            loading_resampled(:,1) = resampling_time';
            loading_resampled(:,2:6) = interp1(loadingData(:,1),loadingData(:,2:6),resampling_time');
                
            unloadingData = extracted_data(locatorPeak(end):locatorBottom(jump+1),:);
            resampling_time = linspace(unloadingData(1,1),unloadingData(end,1),35);
            unloading_resampled(:,1) = resampling_time';
            unloading_resampled(:,2:6) = interp1(unloadingData(:,1),unloadingData(:,2:6),resampling_time');
            
            extracted_data_resampled = [loading_resampled;unloading_resampled];
            extracted_data_resampled = double(extracted_data_resampled);
        else
            if(sum(abs(extracted_data(:,2)-mean(extracted_data(:,2)))<=2)==size(extracted_data,1))
                minLoading = min(extracted_data(round(end/2):round(3*end/4),5));
                minUnloading = min(extracted_data(round(3*end/4):end,5));
                
                [locatorBottom] = find(extracted_data(end/2:end,5)<=max([minLoading+minUnloading,0]));
                speedLocatorBottom = (locatorBottom(2:end)-locatorBottom(1:end-1));
                [~,jump] = max(speedLocatorBottom);
                extracted_data = extracted_data(end/2+jump:end/2+(locatorBottom(jump+1)),:);
            
                resampling_time = linspace(extracted_data(1,1),extracted_data(end,1),70);
                extracted_data_resampled(:,1) = resampling_time';
                extracted_data_resampled(:,2:6) = interp1(extracted_data(:,1),extracted_data(:,2:6),resampling_time');
                extracted_data_resampled = double(extracted_data_resampled);
            else
                locatorBottom = find(~isnan(extracted_data(:,4)));
                resampling_time = extracted_data(locatorBottom(1),1):0.001:extracted_data(locatorBottom(end),1);
                extracted_data_resampled(:,2:6) = interp1(extracted_data(locatorBottom,1),extracted_data(locatorBottom,2:6),resampling_time');
                extracted_data_resampled(:,1) = resampling_time;
                extracted_data_resampled = double(extracted_data_resampled);
            end
        end
    else
        if(strcmp(dataType,'PD_static'))
            [locatorPeak] = find(extracted_data(:,2)>0.99*max(extracted_data(:,2)));
            [locatorBottom] = find(extracted_data(:,2)<=11);
            
            % Find beginning and end of pressure sweep experiment
            speedLocatorBottom = (locatorBottom(2:end)-locatorBottom(1:end-1));
            [~,jump] = max(speedLocatorBottom);
            
            % Resampling loading curves to 35 data points
            loadingData = extracted_data(jump:locatorPeak(1),:);
            resampling_time = linspace(loadingData(1,1),loadingData(end,1),35);
            loading_resampled(:,1) = resampling_time';
            loading_resampled(:,2:6) = interp1(loadingData(:,1),loadingData(:,2:6),resampling_time');
            
            % Resampling unloading curves to 35 datapoints
            unloadingData = extracted_data(locatorPeak(end):locatorBottom(jump+1),:);
            resampling_time = linspace(unloadingData(1,1),unloadingData(end,1),35);
            unloading_resampled(:,1) = resampling_time';
            unloading_resampled(:,2:6) = interp1(unloadingData(:,1),unloadingData(:,2:6),resampling_time');
            
            extracted_data_resampled = double([loading_resampled;unloading_resampled]);

        elseif(strcmp(dataType,'FL_static'))
            locatorUnderscore = find(test_name == '_');
            FL_pressure = test_name(locatorUnderscore+1:end);

            % Find axial force at the beginning and end of the experiment
            minLoading = min(extracted_data(1:round(end/2),5));
            minUnloading = min(extracted_data(round(end/2):end,5));
            
            % Find beginning and end of the axial force experiment
            [locatorBottom] = find(extracted_data(1:end,5)<=max([minLoading+minUnloading,0]));
            speedLocatorBottom = (locatorBottom(2:end)-locatorBottom(1:end-1));
            [max_jump,jump] = max(speedLocatorBottom);
            if(max_jump>150)
                extracted_data = extracted_data(1+jump:1+(locatorBottom(jump+1)),:);
            else
                [locatorBottom] = find(extracted_data(round(end/2):end,5)<=max([minUnloading,0,0.1]));
                extracted_data = extracted_data(1:round(end/2)+(locatorBottom(1)),:);
            end
            
            % Resampling experiments to 70 datapoints
            resampling_time = linspace(extracted_data(1,1),extracted_data(end,1),70);
            extracted_data_resampled(:,1) = resampling_time';
            extracted_data_resampled(:,2:6) = interp1(extracted_data(:,1),extracted_data(:,2:6),resampling_time');
            extracted_data_resampled = double(extracted_data_resampled);
        else
            % Resampling dynamic experiment to 1 kHz
            locatorBottom = find(~isnan(extracted_data(:,4)));
            resampling_time = extracted_data(locatorBottom(1),1):0.001:extracted_data(locatorBottom(end),1);
            extracted_data_resampled(:,2:6) = interp1(extracted_data(locatorBottom,1),extracted_data(locatorBottom,2:6),resampling_time');
            extracted_data_resampled(:,1) = resampling_time;
            extracted_data_resampled = double(extracted_data_resampled);
        end
    end

    pressureV = extracted_data_resampled(:,3); % pressure data
    transducerForceV = extracted_data_resampled(:,5); % axial force data
    axialStretchV = extracted_data_resampled(:,6); % axial stretch data
    axialLengthV = axialStretchV*passive.L; % axial length data
    outerDiameterV = extracted_data_resampled(:,4); % outer diameter data
    timeV = extracted_data_resampled(:,1);

%% Resetting unloaded outer diameter according to force-length sweep at 10 mmHg
    if(forceLengthCheck ==1 && str2double(FL_pressure) == 10)
        Do_l = interp1(transducerForceV(1:5:35),outerDiameterV(1:5:35),0,'spline',"extrap");
        Do_ul = interp1(transducerForceV(35:5:70),outerDiameterV(35:5:70),0,'spline',"extrap");
        passive.OD = (Do_l+Do_ul)/2;
        vesselVolume = pi*(((passive.OD/2/Mtm)^2-((passive.OD-passive.H)/2/Mtm)^2)*passive.L); % Volume of unloaded geometry [mm^3]
    end

    roV = outerDiameterV/2/Mtm;
    riV = radiusFromIincompressibility(roV,axialLengthV,vesselVolume,1); % inner radius from incompressibility
    innerDiameterV = riV*2*Mtm; % inner diameter

    dataMat = [timeV,pressureV,innerDiameterV,outerDiameterV,transducerForceV,axialStretchV,axialLengthV];
end