function [dataMat] = importDataDynamX1(directoryName,test_name,passive,forceLengthCheck,dynamicCheck)

    Mtm = 1000;
    vesselVolume = pi*(((passive.OD/2/Mtm)^2-((passive.OD-passive.H)/2/Mtm)^2)*passive.L); % Volume of unloaded geometry [mm^3]
    sgolayK = 8;
    sgolayF = 51;
    additionalShiftSeconds = -0.002;
    syncMethod = 'electrical2';
    
    usedPressureChan = 10;
    usedForceChan = 14;
    usedElongationChan = 16;
    
    vidArtOrBart = 'vidBart';

 %% Extracting and processing TDMS data (pressure, transducer force and axial elongation)
    TDMSfileName = find_file('.tdms',0,test_name,directoryName);
    TDMSfileName = fullfile(directoryName, TDMSfileName);
    
    if(~dynamicCheck)
        TDMSPressureChannelTag = 'P1_SCALED'; %Note P1 is distal pressure
        TDMSAxialForceChannelTag = 'LOAD_SCALED_FILTERED';
        TDMSLengthstepChannelTag = 'ELONGATION_CORRECTED_MICRON'; 

        TDMS_data = TDMS_getStruct(TDMSfileName); 
        TDMSsyncIV = getSyncStarts(TDMS_data.Untitled.SYNC.data);
        
        pressureV = TDMS_data.Untitled.(TDMSPressureChannelTag).data(TDMSsyncIV)'; % pressure data
        transducerForceV = TDMS_data.Untitled.(TDMSAxialForceChannelTag).data(TDMSsyncIV)'; % axial force data
        axialLengthV = (TDMS_data.Untitled.(TDMSLengthstepChannelTag).data(TDMSsyncIV)/Mtm + (passive.L))'; % axial length data
        axialStretchV= axialLengthV/(passive.L); % axial stretch data
    end

%% Extracting and processing VEVO data (diameter)
    if(forceLengthCheck) % With the old DynamX, the diameter distension is not recorded
        % during the axial force sweeps. Hence, diameter vectors are
        % assigned NaN values
        outerDiameterV = NaN*ones(length(pressureV),1); % outer diameter data set to NaN
        innerDiameterV = NaN*ones(length(pressureV),1); % inner diameter data set to NaN
        timeV = nan(length(pressureV),1);
    elseif(~dynamicCheck)
        VEVOfileName = find_file('.raw.DI.mat',0,test_name,directoryName);
        VEVOfileName = fullfile(directoryName, VEVOfileName);
        i = strfind(lower(VEVOfileName), '.');
        i = i(end-1);
        VEVObase = VEVOfileName(1:i-1);
        VEVOPhysiolFileName = [VEVObase, '.physio'];
        VEVOBModeFileName =   [VEVObase, '.bmode'];

        inDataS = load(VEVOfileName);
    
        vidArtOrBart = 'vidBart';
        switch vidArtOrBart
            case 'vidArt' %Data from Arnold (vidArt) is in pixels
                rawInnerDiameterV = inDataS.vArt.pixscal * mean(inDataS.ddistr);
            case 'vidBart' %Data from Bart (vidBart) is in microns
                rawInnerDiameterV = 0.001*mean(inDataS.ddistr);
        end
    
        [ecgDataV, ~, ~, ~, sampleTimeStampsTicksV,~] = ...
        readVEVORawPhysio(VEVOPhysiolFileName);
        VEVOsyncIV = getSyncStarts(ecgDataV);
        VEVOsampleTimesTicks = sampleTimeStampsTicksV(VEVOsyncIV); 
        [~, frameTimeStampsTicksV, ~, ~, ~] = readVEVORawBMode(VEVOBModeFileName); 

        if (length(frameTimeStampsTicksV) ~= length(rawInnerDiameterV))
            disp('WARNING: Number of frames in diameter file smaller than number of frames in BMODE file.')
            disp('  Assuming that the diameter file corresponds to the FIRST frames of the raw data.')
            disp('  I.e., that in wall tracking, the END of the data was cut off.')
        end
    
        diamInterpolator = griddedInterpolant(double(frameTimeStampsTicksV(1:length(rawInnerDiameterV))),rawInnerDiameterV);
        innerDiameterV = diamInterpolator(double(VEVOsampleTimesTicks))*Mtm; % inner diameter data in um
            % Note: ultrasound tracks the vessel's inner diameter!!
        riV = innerDiameterV/2/Mtm; %inner radius in mm
        roV = radiusFromIincompressibility(riV,axialLengthV,vesselVolume,2);
        outerDiameterV = roV*2*Mtm; % outer diameter data in um
        timeV = nan(length(pressureV),1);
    else
        VEVOWallTrackFileName = find_file('.raw.DI.mat',1,test_name,VEVODir_save);
        VEVOWallTrackFileName = fullfile(VEVODir_save, VEVOWallTrackFileName);

        [timeV, innerDiameterV,tdmsM,~,~] = getSyncedPA(TDMSfileName,VEVOWallTrackFileName(1:end-7), ... %getSyncedPA doesn't want .di.mat extension...
                sgolayK, sgolayF, syncMethod, NaN, additionalShiftSeconds, false, vidArtOrBart);

        pressureV = tdmsM(:,usedPressureChan); % pressure data
        innerDiameterV = innerDiameterV*Mtm; % internal diameter data
            % Note: ultrasound tracks the vessel's inner diameter!!
        transducerForceV = tdmsM(:,usedForceChan); % axial force data
        elongV = tdmsM(:,usedElongationChan); % axial elongation
        axialStretchV = elongV/Mtm/L+1; % axial stretch
        axialLengthV = elongV/Mtm+L; % axial length
        riV = innerDiameterV/2/Mtm; %inner radius in mm
        roV = radiusFromIincompressibility(riV,axialLengthV,vesselVolume,2);
        outerDiameterV = roV*2*Mtm; % outer diameter data in um
    end

    dataMat = [timeV,pressureV,innerDiameterV,outerDiameterV,transducerForceV,axialStretchV,axialLengthV];
end