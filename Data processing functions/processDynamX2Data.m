function [timeDAQV,pressure1V,pressure2V,interpolatedOuterDiameter,transducerForceV,interpolatedStretchV] = processDynamX2Data(path)
% This function is used to import and process the DynamX2 data (i.e., the
% video data, the DAQ data, and the axial motor data). This function is
% used only when the data had not been already processed previously either
% through this app or externally (importing data with GPU parallel
% computing is much faster). For more information, please contact Koen van
% der Laan.

    if(ismac)
        [A,~] = importdata([path,'/DAQ.txt'],'\t',5); % load DAQ data
    else
        [A,~] = importdata([path,'\DAQ.txt'],'\t',5); % load DAQ data
    end

    if isempty(strfind(A{5},'.'))
        sep=',';
    else
        sep='.';
    end
    
%% Load camera video
    if(ismac)
        mov = VideoReader([path '/CAM.avi']); % load video data
    else
        mov = VideoReader([path '\CAM.avi']); % load video data
    end
    frames = read(mov,[1,mov.NumFrames]);
    frames = squeeze(frames(:,:,1,:));
    clear mov 
    
    % find background intensity vallue to subtract from frames 
    int_count = histogram(frames(:),0:255);
    frames = frames-find(int_count==max(int_count),1)-1;
    
    % set D tracking parameters
    cutoff = 10;
    steep = 0.5;
    struc_s = ones(21,21,1);
    struc_l = ones(1,301,1);
    
%% Perform diameter tracking            
    %make psuedo binary grayscale image
    scaled = 1./(1+steep.*exp(cutoff-single(frames)));
    clear frames
    % erode small strucutre
    scaled_er = imopen(scaled,struc_s);
    clear scaled
    %close dark center of vessel
    scaled_dl = imclose(scaled_er,struc_l);
    clear scaled_er
    
    % sum horizontal psuedo-binary pixels to determine diameter 
    d = sum(scaled_dl,2);
    clear scaled_dl
    D = squeeze(mean(d,1));
    clear d
    
%% Load DAQ data
    if(ismac)
        DAQ = readmatrix([path '/DAQ.txt'],'NumHeaderLines',23,'DecimalSeparator',sep,'Delimiter','\t','ConsecutiveDelimitersRule', 'join');
    else
        DAQ = readmatrix([path '\DAQ.txt'],'NumHeaderLines',23,'DecimalSeparator',sep,'Delimiter','\t','ConsecutiveDelimitersRule', 'join');
    end
    timeDAQV = DAQ(:,1);
    transducerForceV = DAQ(:,2);
    pressure1V = DAQ(:,3);
    pressure2V = DAQ(:,4);
    flash = DAQ(:,9);
    
    % determine period/frequency of camera frames (DBA 2021 november)
    diff = (flash(2:end)-flash(1:end-1));
    % diff = (flash(51:end)-flash(50:end-1));
    edge_p = find(diff>2);
    periods = edge_p(2:end)-edge_p(1:end-1);
    cam_period = mode(periods)*(timeDAQV(2)-timeDAQV(1));
    
    % determine starting time of camera frames (DBA 2021 november)
    bin_flash = flash.*0;
    bin_flash(flash>2) = 1;
    a = mode(cumsum(bin_flash));
    time_cam_start = timeDAQV(find(cumsum(bin_flash)==a,1,'last')+1);
    
    % make time_cam from start and period
    time_cam = (0:length(D)-1).*cam_period+time_cam_start;
    
%% Load motors data
    if(ismac)
        motors = readmatrix([path '/motors.txt'],'NumHeaderLines',13,'DecimalSeparator',sep,'Delimiter','\t');
    else
        motors = readmatrix([path '\motors.txt'],'NumHeaderLines',13,'DecimalSeparator',sep,'Delimiter','\t');
    end
    timeMotorsV = motors(:,2);

    stretchV = motors(:,4);

%% Interpolation
    % camera diameter data to daq
    if(length(time_cam)<length(D))
       D = D(1:length(time_cam));
    else
       time_cam = time_cam(1:length(D));
    end
    interpolatedOuterDiameter = interp1(time_cam,D,timeDAQV);

    % motor data to daq data
    transducerForceV = interp1(timeMotorsV,transducerForceV,timeDAQV,'linear','extrap');
    interpolatedStretchV=interp1(timeMotorsV,stretchV,timeDAQV,'linear','extrap');
end
