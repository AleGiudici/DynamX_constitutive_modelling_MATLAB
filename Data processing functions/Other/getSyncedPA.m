function [tV, diamV, tdmsM, beatStartIV, beatEndIV] = getSyncedPA(fnTDMS, fnVEVOBase, sgolayK, sgolayF, syncMethod, syncPressureChannel, additionalShiftSeconds, doPlot, vidArtOrBart)

  verbose = false;
  if ~exist('doPlot', 'var')
    doPlot = false;
  end
  
  if ~exist('additionalShiftSeconds','var')
    additionalShiftSeconds = 0;
  end
  
  if additionalShiftSeconds ~= 0
    if verbose
      disp(sprintf('Adding extra dt of %f seconds', additionalShiftSeconds))
    end
  end
  
  %syncMethod can be
  % 'electrical'
  % 'maxSecDer' %Maximum second derivative
  % 'electrical2' %E-sync looking for single, wide sync pulse.
  
  
  fsTDMS = 2500;
%   tdmsSyncChannel = 6;
  TDMSSyncChannelTag = 'SYNC';
  TDMSChannelTag=["TEMP1","TEMP2","P0","P1","AI","SYNC","P2","LOAD","P0_SCALED","P1_SCALED","P2_SCALED","AI_SCALED","LOAD_SCALED_UNFILTERED","LOAD_SCALED_FILTERED","STEPPER_LOCATION_MICRON","ELONGATION_CORRECTED_MICRON","DRIVING_PRESSURE_SETPOINT"];
  interpolationMethod = 'linear';
  expectedLongPulseDurationMilliseconds  = 50; %For 'electrical2' method
  expectedShortPulseDurationMilliseconds = 10; %For 'electrical2' method; in case no long pulses found
  
  
  fResampleFilter = 1000;
  fResampleOutput = 1000;
  
  derFilt = 10; %Derivative filtering cut-off
  
  
%   global dllPath hPath
  %% Read TDMS
  TDMS_data = TDMS_getStruct(fnTDMS);
  syncSignalTDMSV = TDMS_data.Untitled.(TDMSSyncChannelTag).data';
  
  for i=1:length(TDMSChannelTag)
      inDataM(1:length(syncSignalTDMSV),i)=TDMS_data.Untitled.(TDMSChannelTag{i}).data';
  end
  
%   [inDataM, channelNamesC] = readTDMS(fnTDMS, dllPath, hPath);
  tTDMSV = ( (1/fsTDMS)* ( 0:(size(TDMS_data.Untitled.(TDMSSyncChannelTag).data',1)-1) ) )';
  
%   syncSignalTDMSV = inDataM(:,tdmsSyncChannel);
  
  switch syncMethod
    case 'electrical2'
      longSyncPulseStartTimeTDMSMilliseconds  = findWideSyncPulse(syncSignalTDMSV, tTDMSV*1000, expectedLongPulseDurationMilliseconds, true);
      shortSyncPulseStartTimeTDMSMilliseconds = findWideSyncPulse(syncSignalTDMSV, tTDMSV*1000, expectedShortPulseDurationMilliseconds, true);
      if isscalar(longSyncPulseStartTimeTDMSMilliseconds) && isnan(longSyncPulseStartTimeTDMSMilliseconds)
        warning('No long sync pulse detected in TDMS, now syncing based on short pulses')
        useShortInsteadOfLongPulse = true;
      else
        useShortInsteadOfLongPulse = false;
      end
    otherwise
      syncStartIndicesTDMSV = getSyncStarts(syncSignalTDMSV);
      syncStartTimesTDMSV = tTDMSV(syncStartIndicesTDMSV);
  end
  
  TDMSRawDataFC = arrayfun(@(dataI) griddedInterpolant(tTDMSV, inDataM(:,dataI), interpolationMethod, 'none'), 1:size(inDataM,2), 'UniformOutput', false);
  
  TDMSFilterTV = ( tTDMSV(1): 1/fResampleFilter : tTDMSV(end) )';
  TDMSPreFilterDVC  = arrayfun(@(dataI) TDMSRawDataFC{dataI}(TDMSFilterTV), 1:size(inDataM,2), 'UniformOutput', false);
  TDMSPostFilterDVC = arrayfun(@(dataI) sgolayfiltValid(TDMSPreFilterDVC{dataI},sgolayK,sgolayF), 1:size(inDataM,2), 'UniformOutput', false);
  TDMSDataFC = arrayfun(@(dataI) griddedInterpolant(TDMSFilterTV, TDMSPostFilterDVC{dataI}, interpolationMethod, 'none'), 1:size(inDataM,2), 'UniformOutput', false);
  
  %% Process distension MAT and VEVO physio
  distDataS = load([fnVEVOBase '.DI.mat'], '-mat'); %Distension data for 1000 frames.
  distV = distDataS.vArt.pixscal*mean(distDataS.ddistr(2:end-2,:))'; %Compute mean distension, exclude outer lines.
  
  switch vidArtOrBart
    case 'vidArt' %Data from Arnold (vidArt) is in pixels
      distV = distDataS.vArt.pixscal*mean(distDataS.ddistr(2:end-2,:))'; %Compute mean distension, exclude outer lines.
    case 'vidBart' %Data from Bart (vidBart) is in microns
      distV = 0.001*mean(distDataS.ddistr(2:end-2,:))'; %Compute mean distension, exclude outer lines.
  end
  
  
  bnd = distDataS.vArt.frsel; %See whether all frames were analysed using vidDist, or only a subset (for whatever reason).
  %Generally (if all frames analysed): bnd = [1 1000];
  
  %Find time points that correspond to these 1000 frames
  [~, frameTimeStampsTicksV, frameTimeStampsMilliSecondsV, ~, sysParsVEVOS] = readVEVORawBMode([fnVEVOBase, '.bmode']);
  
  frameTimeStampsTicksV = frameTimeStampsTicksV(bnd(1):bnd(2));
  frameTimeStampsMilliSecondsV = frameTimeStampsMilliSecondsV(bnd(1):bnd(2));
  frameTimeStampsSecondsV = frameTimeStampsMilliSecondsV*0.001;
  
  %Create interpolation function for raw distension data.
  distDataRawF = griddedInterpolant(frameTimeStampsSecondsV,distV, interpolationMethod, 'none');
  
  %Filter distension data and create post-filtering interpolation function.
  distDataFilterTV = ( frameTimeStampsSecondsV(1):1/fResampleFilter:frameTimeStampsSecondsV(end) )';
  distDataPreFilterDV = distDataRawF(distDataFilterTV);
  distDataPostFilterDV = sgolayfiltValid(distDataPreFilterDV,sgolayK,sgolayF);
  distDataF = griddedInterpolant(distDataFilterTV,distDataPostFilterDV, interpolationMethod, 'none');
  
  %Get physio data
  fNameRawPhysio = [fnVEVOBase, '.physio'];
  [ecgDataV, ~, ~, ~, ...
    sampleTimeStampsTicksV, sampleTimeStampsMillisecondsV] = ...
    readVEVORawPhysio(fNameRawPhysio);
  
  switch syncMethod
    case 'electrical2'
      expectedLongPulseDurationMilliseconds = 50;
      
      longSyncPulseStartTimeVEVOMilliseconds = findWideSyncPulse(ecgDataV, sampleTimeStampsMillisecondsV, expectedLongPulseDurationMilliseconds, false);
      if (isscalar(longSyncPulseStartTimeVEVOMilliseconds) && isnan(longSyncPulseStartTimeVEVOMilliseconds))
        warning('No long sync pulse detected in VEVO signal, now syncing based on short pulses')
        useShortInsteadOfLongPulse = true;
      end
      if useShortInsteadOfLongPulse
        shortSyncPulseStartTimeVEVOMilliseconds = findWideSyncPulse(ecgDataV, sampleTimeStampsMillisecondsV, expectedShortPulseDurationMilliseconds, true);
        dtElectric = (shortSyncPulseStartTimeVEVOMilliseconds - shortSyncPulseStartTimeTDMSMilliseconds)/1000;
      else
        dtElectric = (longSyncPulseStartTimeVEVOMilliseconds - longSyncPulseStartTimeTDMSMilliseconds)/1000;
      end
      
      totalDt = dtElectric + additionalShiftSeconds;
      
      tV = (frameTimeStampsSecondsV(1):1/fResampleOutput:frameTimeStampsSecondsV(end))';
      diamV = distDataF(tV);
    otherwise
      syncStartIndicesVEVOV = getSyncStarts(ecgDataV);
      syncStartTimesVEVOV = 0.001*sampleTimeStampsMillisecondsV(syncStartIndicesVEVOV);
      
      %%
      % The VEVO acquisition window is the limiting factor (1000 frames
      % typically). Find the sync pulses that are in this window.
      t1 = frameTimeStampsSecondsV(1);
      t2 = frameTimeStampsSecondsV(end);
      tempV = (syncStartTimesVEVOV>=t1) & (syncStartTimesVEVOV<=t2);
      syncStartTimesVEVOTrimmedV   = syncStartTimesVEVOV(tempV);
      %syncStartIndicesVEVOTrimmedV = syncStartIndicesVEVOV(tempV);
      
      nSyncPulsesUsed = size(syncStartTimesVEVOTrimmedV,1);
      dtV = syncStartTimesVEVOTrimmedV - syncStartTimesTDMSV(end-nSyncPulsesUsed+1:end);
      dtElectric = mean(dtV); %This is the time (in seconds) that the VEVO time is LATER than the TDMS time.
      if max(abs(diff(dtV))) > 0.050
        warning('Large dt dispersion!')
      end
      
      if verbose
        fprintf(['Electrical dt: %f' char(13)], dtElectric)
      end
      
      tV = (t1:1/fResampleOutput:t2)';
      diamV = distDataF(tV);
      
      switch syncMethod
        case 'electrical'
          totalDt = dtElectric;
        case 'maxSecDer'
          [bHigh,aHigh] = butter(1,derFilt/(fResampleOutput/2),'high'); %filter characteristics
          diamDerV = -1*nanfiltfilt(bHigh,aHigh,diamV);                                  %second derivative
          presDerV = -1*nanfiltfilt(bHigh,aHigh,TDMSPostFilterDVC{syncPressureChannel}); %second derivative
          
          %Find the indices of the ELECTRICAL sync pulses, in both diameter and
          %pressure vectors
          [~, syncDiamIV] = arrayfun(@(i) min(abs(syncStartTimesVEVOTrimmedV(i) - tV)), (1:size(syncStartTimesVEVOTrimmedV,1))');
          [~, syncPresIV] = arrayfun(@(i) min(abs(syncStartTimesVEVOTrimmedV(i) - dtElectric - TDMSFilterTV)), (1:size(syncStartTimesVEVOTrimmedV,1))');
          
          nFullBeats = size(syncDiamIV,1)-1;
          
          [maxDiamDersV, relMaxIDiamV] = arrayfun(@(i) max(diamDerV(syncDiamIV(i):syncDiamIV(i+1)), [], 'includenan'),(1:nFullBeats)');
          relMaxIDiamV(isnan(maxDiamDersV)) = NaN;
          maxIDiamV = relMaxIDiamV + syncDiamIV(1:end-1) - 1;
          [maxPresDersV, relMaxIPresV] = arrayfun(@(i) max(presDerV(syncPresIV(i):syncPresIV(i+1)), [], 'includenan'),(1:nFullBeats)');
          relMaxIPresV(isnan(maxPresDersV)) = NaN;
          maxIPresV = relMaxIPresV + syncPresIV(1:end-1) - 1;
          
          if doPlot
            figure
            a1 = subplot(2,1,1);
            plot(TDMSFilterTV, presDerV)
            hold on
            plot(indnan(TDMSFilterTV,maxIPresV), indnan(presDerV,maxIPresV), 'or')
            ylabel('Pressure second derivative')
            xlabel('Time')
            
            a2 = subplot(2,1,2);
            plot(tV, diamDerV)
            hold on
            plot(indnan(tV,maxIDiamV), indnan(diamDerV,maxIDiamV), 'or')
            ylabel('Diameter second derivative')
            xlabel('Time')
            totalDt = nanmean(indnan(tV,maxIDiamV) - indnan(TDMSFilterTV,maxIPresV));
          end
        otherwise
          error('Invalid syncMethod specified!')
      end
      
  end
  
  if verbose
    fprintf(['Total (used) dt: %f' char(13)], totalDt)
    if strcmp('syncMethod', 'maxSecDer')
      fprintf(['Extra dt due to use of second derivative alignment: %f' char(13)], totalDt-dtElectric)
    end
  end
  
  tdmsM = arrayfun(@(a) TDMSDataFC{a}(tV-totalDt), 1:size(inDataM,2), 'UniformOutput', false);
  tdmsM = cell2mat(tdmsM);
  
  if doPlot
    figure('Name', 'Synchronised signals')
    a1 = subplot(2,1,1);
    plot(tV, tdmsM(:,syncPressureChannel))
    ylabel(['Pressure (' channelNamesC{syncPressureChannel} ')'])
    xlabel('Time')
    
    a2 = subplot(2,1,2);
    plot(tV, diamV)
    ylabel('Diameter')
    xlabel('Time')
    
    linkaxes([a1 a2], 'x')
  end
  
  if strcmpi(syncMethod, 'maxSecDer')
    [bHigh,aHigh] = butter(1,derFilt/(fResampleOutput/2),'high'); %filter characteristics
    diamDerV = -1*nanfiltfilt(bHigh,aHigh,diamV);          %second derivative
    presDerV = -1*nanfiltfilt(bHigh,aHigh,tdmsM(:, syncPressureChannel));          %second derivative
    
    if doPlot
      a3 = subplot(4,1,3);
      plot(tV, presDerV)
      ylabel('Pressure second derivative')
      xlabel('Time')
      
      a4 = subplot(4,1,4);
      plot(tV, diamDerV)
      ylabel('Diameter second derivative')
      xlabel('Time')
      linkaxes([a1 a2 a3 a4], 'x')
    end
    
    %Redo detection for plotting
    [~, syncDiamIV] = arrayfun(@(i) min(abs(syncStartTimesVEVOTrimmedV(i) - tV)), (1:size(syncStartTimesVEVOTrimmedV,1))');
    nFullBeats = size(syncDiamIV,1)-1;
    [~, relMaxIDiamV] = arrayfun(@(i) max(diamDerV(syncDiamIV(i):syncDiamIV(i+1))),(1:nFullBeats)');
    maxIDiamV = relMaxIDiamV + syncDiamIV(1:end-1) - 1;
    [~, relMaxIPresV] = arrayfun(@(i) max(presDerV(syncDiamIV(i):syncDiamIV(i+1))),(1:nFullBeats)');
    maxIPresV = relMaxIPresV + syncDiamIV(1:end-1) - 1;
    
    if doPlot
      axes(a3)
      hold on
      plot(tV(maxIPresV), presDerV(maxIPresV), 'or')
      
      axes(a4)
      hold on
      plot(tV(maxIDiamV), diamDerV(maxIDiamV), 'or')
    end
    
    beatStartIV = maxIPresV(1:end-1);
    beatEndIV   = maxIPresV(2:end  )-1;
  else
    beatStartIV = [];
    beatEndIV = [];
  end
end

