function [forceStepsV, elongationStepsV, pressureStepsV] = getStaticFL(fnTDMS)
  %, fnVEVOBase
  %VEVOWallTrackFileName = [fnVEVOBase, '.DI.mat'];
  %VEVOPhysiolFileName =   [fnVEVOBase, '.physio'];
  %VEVOBModeFileName =     [fnVEVOBase, '.bmode'];
  
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
  pressureStepsV   = TDMSDataM(TDMSsyncIV, ...
    find(strcmp(TDMSChannelNamesC, 'P1_SCALED')));
  forceStepsV      = TDMSDataM(TDMSsyncIV, ...
    find(strcmp(TDMSChannelNamesC, 'LOAD_SCALED_FILTERED')));
  elongationStepsV = TDMSDataM(TDMSsyncIV, ...
    find(strcmp(TDMSChannelNamesC, 'ELONGATION_CORRECTED_MICRON')));
  
 
end