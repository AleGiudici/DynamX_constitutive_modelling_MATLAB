function [diamStepsV, pressureStepsV, elongationStepsV, forceStepsV] = getStaticPA(fnTDMS, fnVEVOBase, vidArtOrBart)
  VEVOWallTrackFileName = [fnVEVOBase, '.DI.mat'];
  VEVOPhysiolFileName =   [fnVEVOBase, '.physio'];
  VEVOBModeFileName =     [fnVEVOBase, '.bmode'];
  
  TDMSPressureChannelTag = 'P1_SCALED'; %Note, normally P1 is distal pressure
  %Can be 'P0_SCALED', 'P1_SCALED', 'P2_SCALED', 'DRIVING_PRESSURE_SETPOINT'
  
  TDMSFs = 2500;
  
  %% TDMS processing
  global dllPath
  global hPath
  [TDMSDataM,TDMSChannelNamesC]=readTDMS(fnTDMS, dllPath, hPath);
  TDMSPressureChannelNumber = find(strcmp(TDMSChannelNamesC, TDMSPressureChannelTag));
  TDMSTimeV = (0:(size(TDMSDataM,1)-1))/TDMSFs;
  
  TDMSsyncIV = getSyncStarts(TDMSDataM(:,6));
  if length(TDMSsyncIV) ~= 81
    warning('Number of sync pulses in TDMS unequal to 81!')
  end
  pressureStepsV = TDMSDataM(TDMSsyncIV,TDMSPressureChannelNumber);
  elongationStepsV = TDMSDataM(TDMSsyncIV, ...
    find(strcmp(TDMSChannelNamesC, 'ELONGATION_CORRECTED_MICRON')));
  forceStepsV = TDMSDataM(TDMSsyncIV, ...
    find(strcmp(TDMSChannelNamesC, 'LOAD_SCALED_FILTERED')));
  
  %% VEVO processing
  inDataS = load(VEVOWallTrackFileName);
  switch vidArtOrBart
    case 'vidArt' %Data from Arnold (vidArt) is in pixels
      diamV = inDataS.vArt.pixscal * mean(inDataS.ddistr(2:end-1,:));
    case 'vidBart' %Data from Bart (vidBart) is in microns
      diamV = mean(inDataS.ddistr(2:end-1,:))/1000;
  end
      
  
  [ecgDataV, ~, ~, ~, sampleTimeStampsTicksV, sampleTimeStampsMillisecondsV] = ...
    readVEVORawPhysio(VEVOPhysiolFileName);
  VEVOsyncIV = getSyncStarts(ecgDataV);
  VEVOsampleTimesTicks = sampleTimeStampsTicksV(VEVOsyncIV); %This variable
  %contains the timing (in ticks) of the sync pulses
  
  [~, frameTimeStampsTicksV, ~, ~, ~] = readVEVORawBMode(VEVOBModeFileName);
  %Timing (in ticks) of the FRAMES, i.e., of the diameters in diamV
  
  if (length(frameTimeStampsTicksV) ~= length(diamV))
    disp('WARNING: Number of frames in diameter file smaller than number of frames in BMODE file.')
    disp('  Assuming that the diameter file corresponds to the FIRST frames of the raw data.')
    disp('  I.e., that in wall tracking, the END of the data was cut off.')
  end
  
  diamInterpolator = griddedInterpolant(double(frameTimeStampsTicksV(1:length(diamV))),diamV);
  diamStepsV = diamInterpolator(double(VEVOsampleTimesTicks));
  
%   %% Plotting
%   plot(0.25*pi*diamStepsV.^2, pressureStepsV)
%   xlabel('Cross-sectional area [mm^2]')
%   ylabel('Pressure [mmHg]')
  
end