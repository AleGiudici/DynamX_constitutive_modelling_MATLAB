% load DAQ data from files
% dont forget to change decimal seperator to points in DAQ data file

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

mN2g = 1/9.81; % mN to g comversion

% get folder location with DAQ, motors and CAM files.
% path=uigetdir(VEVODir_save);
% path=uigetdir('C:\Users\p70068218\Documents\BME PhD maastricht\Data\')
check = 0;
clear exported_data extracted_data
if(strcmp(data_type,'PD_static') || strcmp(data_type,'PD_dynamic'))
    if(isfile([fullfile(path,PD{k}) '.mat']))
        load([fullfile(path,PD{k}) '.mat']);
        extracted_data = table2array(exported_data);
        check = 1;
    end
elseif(strcmp(data_type,'FL_static'))
    if(isfile([fullfile(path,FL{k}) '.mat']))
        load([fullfile(path,FL{k}) '.mat']);
        extracted_data = table2array(exported_data);
        check = 1;
    end
end

if(check==0)
    if(ismac)
        if(~isfile([path '/Data.mat']))
            [A,B]=importdata([path,'/DAQ.txt'],'\t',5);
            if isempty(strfind(A{5},'.'))
                sep=',';
            else
                sep='.';
            end
            
            % load camera video
            mov=VideoReader([path '/CAM.avi']);
            frames=read(mov,[1,mov.NumFrames]);
            frames=squeeze(frames(:,:,1,:));
            clear mov 
            
            % find background intensity vallue to subtract from frames 
            int_count=hist(frames(:),0:255);
            frames=frames-find(int_count==max(int_count),1)-1;
            
            % set D tracking parameters
            cutoff=10;
            steep=0.5;
            struc_s=ones(21,21,1);
            struc_l=ones(1,301,1);
            
            % perform Diameter tracking
            
            %make psuedo binary grayscale image
            scaled=1./(1+steep.*exp(cutoff-single(frames)));
            clear frames
            % erode small strucutre
            scaled_er=imopen(scaled,struc_s);
            clear scaled
            %close dark center of vessel
            scaled_dl=imclose(scaled_er,struc_l);
            clear scaled_er
            
            % sum horizontal psuedo-binary pixels to determine diameter 
            d=sum(scaled_dl,2);
            clear scaled_dl
            D=squeeze(mean(d,1));
            clear d
            
            %load daq data
            DAQ=readmatrix([path '/DAQ.txt'],'NumHeaderLines',23,'DecimalSeparator',sep,'Delimiter','\t','ConsecutiveDelimitersRule', 'join');
            time_daq=DAQ(:,1);
            force=DAQ(:,2);
            p1=DAQ(:,3);
            p2=DAQ(:,4);
            flash=DAQ(:,9);
            
            % determine period/frequency of camera frames (DBA 2021 november)
            diff=(flash(2:end)-flash(1:end-1));
            diff=(flash(51:end)-flash(50:end-1));
            edge_p=find(diff>2);
            periods=edge_p(2:end)-edge_p(1:end-1);
            cam_period=mode(periods)*(time_daq(2)-time_daq(1));
            
            % determine starting time of camera frames (DBA 2021 november)
            bin_flash=flash.*0;
            bin_flash(flash>2)=1;
            a=mode(cumsum(bin_flash));
            time_cam_start=time_daq(find(cumsum(bin_flash)==a,1,'last')+1);
            
            % make time_cam from start and period
            time_cam=(0:length(D)-1).*cam_period+time_cam_start;
            
            % load motors data
            motors=readmatrix([path '/motors.txt'],'NumHeaderLines',13,'DecimalSeparator',sep,'Delimiter','\t');
            time_motors=motors(:,2);
            unloaded=motors(:,3);
            stretch=motors(:,4);
            m1=motors(:,5);
            m2=motors(:,6);
    
            % interpolatie camera diameter data
            if length(time_cam)<length(D)
               D=D(1:length(time_cam));
            else
                time_cam=time_cam(1:length(D));
            end
            int_D=interp1(time_cam,D,time_daq);
        else
           load([path '/Data.mat']);
           dataMat = matfile([path '/Data.mat']);
           name_fields = fieldnames(dataMat);
           time_daq = eval([name_fields{2,1} '.time_pressure']);
           p1 = eval([name_fields{2,1} '.p1']);
           p2 = eval([name_fields{2,1} '.p2']);
           force = eval([name_fields{2,1} '.load']);
           stretch = eval([name_fields{2,1} '.stretch']);
           D = eval([name_fields{2,1} '.D']);
           time_motors = eval([name_fields{2,1} '.time_stretch']);
    
           % interpolate motor data to daq data
           force=interp1(time_motors,force,time_daq,'linear','extrap');
        end
    else
        if(~isfile([path '\Data.mat']))
            [A,B]=importdata([path,'\DAQ.txt'],'\t',5);
            if isempty(strfind(A{5},'.'))
                sep=',';
            else
                sep='.';
            end
            
            % load camera video
            mov=VideoReader([path '\CAM.avi']);
            frames=read(mov,[1,mov.NumFrames]);
            frames=squeeze(frames(:,:,1,:));
            clear mov 
            
            % find background intensity vallue to subtract from frames 
            int_count=hist(frames(:),0:255);
            frames=frames-find(int_count==max(int_count),1)-1;
            
            % set D tracking parameters
            cutoff=10;
            steep=0.5;
            struc_s=ones(21,21,1);
            struc_l=ones(1,301,1);
            
            % perform Diameter tracking
            
            %make psuedo binary grayscale image
            scaled=1./(1+steep.*exp(cutoff-single(frames)));
            clear frames
            % erode small strucutre
            scaled_er=imopen(scaled,struc_s);
            clear scaled
            %close dark center of vessel
            scaled_dl=imclose(scaled_er,struc_l);
            clear scaled_er
            
            % sum horizontal psuedo-binary pixels to determine diameter 
            d=sum(scaled_dl,2);
            clear scaled_dl
            D=squeeze(mean(d,1));
            clear d
            
            %load daq data
            DAQ=readmatrix([path '\DAQ.txt'],'NumHeaderLines',23,'DecimalSeparator',sep,'Delimiter','\t','ConsecutiveDelimitersRule', 'join');
            time_daq=DAQ(:,1);
            force=DAQ(:,2);
            p1=DAQ(:,3);
            p2=DAQ(:,4);
            flash=DAQ(:,9);
            
            % determine period/frequency of camera frames (DBA 2021 november)
            diff=(flash(2:end)-flash(1:end-1));
            diff=(flash(51:end)-flash(50:end-1));
            edge_p=find(diff>2);
            periods=edge_p(2:end)-edge_p(1:end-1);
            cam_period=mode(periods)*(time_daq(2)-time_daq(1));
            
            % determine starting time of camera frames (DBA 2021 november)
            bin_flash=flash.*0;
            bin_flash(flash>2)=1;
            a=mode(cumsum(bin_flash));
            time_cam_start=time_daq(find(cumsum(bin_flash)==a,1,'last')+1);
            
            % make time_cam from start and period
            time_cam=(0:length(D)-1).*cam_period+time_cam_start;
            
            % load motors data
            motors=readmatrix([path '\motors.txt'],'NumHeaderLines',13,'DecimalSeparator',sep,'Delimiter','\t');
            time_motors=motors(:,2);
            unloaded=motors(:,3);
            stretch=motors(:,4);
            m1=motors(:,5);
            m2=motors(:,6);
    
            % interpolatie camera diameter data
            if length(time_cam)<length(D)
               D=D(1:length(time_cam));
            else
                time_cam=time_cam(1:length(D));
            end
            int_D=interp1(time_cam,D,time_daq);
        else
           load([path '\Data.mat']);
           dataMat = matfile([path '\Data.mat']);
           name_fields = fieldnames(dataMat);
           time_daq = eval([name_fields{2,1} '.time_pressure']);
           p1 = eval([name_fields{2,1} '.p1']);
           p2 = eval([name_fields{2,1} '.p2']);
           force = eval([name_fields{2,1} '.load']);
           stretch = eval([name_fields{2,1} '.stretch']);
           D = eval([name_fields{2,1} '.D']);
           time_motors = eval([name_fields{2,1} '.time_stretch']);
    
           % interpolate motor data to daq data
           force=interp1(time_motors,force,time_daq,'linear','extrap');
        end
    end
    
    % interpolate motor data to daq data
    int_stretch=interp1(time_motors,stretch,time_daq,'linear','extrap');

    if(isnan(int_stretch(1)))
        int_stretch = ones(length(int_stretch),1)*1.6138;
    end
    
    extracted_data = [time_daq, p1, p2, D*p2um, force*V2mN*mN2g, int_stretch];

%     locator_NaN = ~isnan(int_D);
%     extracted_data = extracted_data(locator_NaN,:);
%     extracted_data(:,6)=1.5908430;

end

clear extracted_data_resampled;
if(strcmp(protocol_type,'Old'))
    [locator] = find(extracted_data(:,2)>200);
    [locator_2] = find(extracted_data(:,2)<=11);

    if((size(locator,1)~=0)&(size(locator_2,1)~=0))
        speed_locator_2 = (locator_2(2:end)-locator_2(1:end-1));
        [~,jump] = max(speed_locator_2);
        
        loading = extracted_data(jump:locator(1),:);
        resampling_time = linspace(loading(1,1),loading(end,1),35);
        loading_resampled(:,1) = resampling_time';
        loading_resampled(:,2:6) = interp1(loading(:,1),loading(:,2:6),resampling_time');
            
        unloading = extracted_data(locator(end):locator_2(jump+1),:);
        resampling_time = linspace(unloading(1,1),unloading(end,1),35);
        unloading_resampled(:,1) = resampling_time';
        unloading_resampled(:,2:6) = interp1(unloading(:,1),unloading(:,2:6),resampling_time');
        
        extracted_data_resampled = [loading_resampled;unloading_resampled];
        extracted_data_resampled = double(extracted_data_resampled);
    else
        if(sum(abs(extracted_data(:,2)-mean(extracted_data(:,2)))<=2)==size(extracted_data,1))
            min_1 = min(extracted_data(round(end/2):round(3*end/4),5));
            min_2 = min(extracted_data(round(3*end/4):end,5));
            
            [locator_2] = find(extracted_data(end/2:end,5)<=max([min_1+min_2,0]));
            speed_locator = (locator_2(2:end)-locator_2(1:end-1));
            [~,jump] = max(speed_locator);
            extracted_data = extracted_data(end/2+jump:end/2+(locator_2(jump+1)),:);
        
            resampling_time = linspace(extracted_data(1,1),extracted_data(end,1),70);
            extracted_data_resampled(:,1) = resampling_time';
            extracted_data_resampled(:,2:6) = interp1(extracted_data(:,1),extracted_data(:,2:6),resampling_time');
            extracted_data_resampled = double(extracted_data_resampled);
        else
            locator_2 = find(~isnan(extracted_data(:,4)));
            resampling_time = extracted_data(locator_2(1),1):0.001:extracted_data(locator_2(end),1);
    %         resampling_time = linspace(extracted_data(locator(1),1),extracted_data(locator(end),1),2500);
            extracted_data_resampled(:,2:6) = interp1(extracted_data(locator_2,1),extracted_data(locator_2,2:6),resampling_time');
            extracted_data_resampled(:,1) = resampling_time;
            extracted_data_resampled = double(extracted_data_resampled);
        end
    end
else
    if(strcmp(data_type,'PD_static'))
        [locator] = find(extracted_data(:,2)>0.99*max(extracted_data(:,2)));
        [locator_2] = find(extracted_data(:,2)<=11);

        speed_locator_2 = (locator_2(2:end)-locator_2(1:end-1));
        [~,jump] = max(speed_locator_2);
        
        loading = extracted_data(jump:locator(1),:);
        resampling_time = linspace(loading(1,1),loading(end,1),35);
        loading_resampled(:,1) = resampling_time';
        loading_resampled(:,2:6) = interp1(loading(:,1),loading(:,2:6),resampling_time');
            
        unloading = extracted_data(locator(end):locator_2(jump+1),:);
        resampling_time = linspace(unloading(1,1),unloading(end,1),35);
        unloading_resampled(:,1) = resampling_time';
        unloading_resampled(:,2:6) = interp1(unloading(:,1),unloading(:,2:6),resampling_time');
        
        extracted_data_resampled = [loading_resampled;unloading_resampled];
        extracted_data_resampled = double(extracted_data_resampled);
    elseif(strcmp(data_type,'FL_static'))
        min_1 = min(extracted_data(1:round(end/2),5));
        min_2 = min(extracted_data(round(end/2):end,5));
        
        [locator_2] = find(extracted_data(1:end,5)<=max([min_1+min_2,0]));
        speed_locator = (locator_2(2:end)-locator_2(1:end-1));
        [max_jump,jump] = max(speed_locator);
        if(max_jump>1000)
            extracted_data = extracted_data(1+jump:1+(locator_2(jump+1)),:);
        else
            [locator_2] = find(extracted_data(round(end/2):end,5)<=max([min_2,0,0.1]));
            extracted_data = extracted_data(1:round(end/2)+(locator_2(1)),:);
        end
    
        resampling_time = linspace(extracted_data(1,1),extracted_data(end,1),70);
        extracted_data_resampled(:,1) = resampling_time';
        extracted_data_resampled(:,2:6) = interp1(extracted_data(:,1),extracted_data(:,2:6),resampling_time');
        extracted_data_resampled = double(extracted_data_resampled);
    else
        locator_2 = find(~isnan(extracted_data(:,4)));
        resampling_time = extracted_data(locator_2(1),1):0.001:extracted_data(locator_2(end),1);
%         resampling_time = linspace(extracted_data(locator(1),1),extracted_data(locator(end),1),2500);
        extracted_data_resampled(:,2:6) = interp1(extracted_data(locator_2,1),extracted_data(locator_2,2:6),resampling_time');
        extracted_data_resampled(:,1) = resampling_time;
        extracted_data_resampled = double(extracted_data_resampled);
    end
end

exported_data = array2table(extracted_data,'VariableNames',{'Time [s]','P1 [mmHg]','P2 [mmHg]','Diameter [um]','Force [g]','Stretch [-]'});
if(strcmp(data_type,'PD_static') || strcmp(data_type,'PD_dynamic'))
    filename = fullfile(path,PD{k});
else
    filename = fullfile(path,FL{k});
end
save(filename,'exported_data');

% location_D = ~isnan(extracted_data(:,3));
% 
% location_D_2 = find(location_D==1);
% 
% extracted_data = extracted_data(location_D_2,:)
% 
% extracted_data = extracted_data(1000:end,:);
% 
% [values,positions] = findpeaks(-sgolayfilt(double(extracted_data(:,1)),2,41),'MinPeakHeight',-125);
% extracted_data_loop = 0;
% for i = 1:length(positions)-1
%     extracted_data_loop = extracted_data_loop + extracted_data(positions(i):positions(i)+400,:);
% end
% extracted_data_loop = extracted_data_loop/(length(positions)-1);
% 
% 
% % extracted_data = interp1(1:size(extracted_data,1), extracted_data, linspace(1,size(extracted_data,1),200));
% 
% prompt = {'Enter new file name'};
% title = 'Saving results...';
% definput = {''};
% opts.Interpreter = 'tex';
% mouse_ID = inputdlg(prompt,title,[1 40],definput,opts);
% 
% save(sprintf('data_%s.mat',mouse_ID{1}), 'extracted_data')