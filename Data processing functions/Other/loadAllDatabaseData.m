function [dynamicDataT, staticPDDataT, staticFLDataT] = loadAllDatabaseData(databaseT, folderName, linesToLoad)
  if ~exist('linesToLoad', 'var')
    linesToLoad = 1:size(databaseT,1);
  end
  usedPressureChan = 10;
  dynamicProtocolPressuresV = [80;120;120;120;160];
  dynamicProtocolFrequenciesV = [5;2.5;5;10;5];
  sgolayK = 8;
  sgolayF = 51;
  additionalShiftSeconds = -0.004;
  syncMethod = 'electrical2';
  
  vidArtOrBart = 'vidBart';
  
  staticPDProtocolStretchesV = [0.95;1.00;1.05];
  staticFLProtocolPressuresV = [10; 60; 100; 140; 200];
  

  
  
  dynamicDataT = table();
  staticPDDataT = table();
  staticFLDataT = table();
  
  for iDataLine = linesToLoad
    fprintf('--==Loading mouse ID %s, sample ID %s, experiment ID %s==--\n', databaseT.MouseID{iDataLine}, databaseT.SampleID{iDataLine}, databaseT.ExperimentID{iDataLine})
    
    %--DYNAMIC DATA LOADING--
    %Determine folder/file name base part
    fullDynamicLVFileNameBase = fullfile(folderName, ...
      [databaseT.BaseFolderName{iDataLine} '_Labview'], ...
      databaseT.DynamicLVFileNameBase{iDataLine});
    fullDynamicVEVOFileNameBase = fullfile(folderName, ...
      [databaseT.BaseFolderName{iDataLine} '_VEVO'], ...
      databaseT.DynamicVEVOFileNameBase{iDataLine});
    %Load dynamic data for all (typically five) protocol steps
    dynamicDataSingleSampleT = loadSingleSampleDynamicData(fullDynamicLVFileNameBase, fullDynamicVEVOFileNameBase, dynamicProtocolPressuresV, dynamicProtocolFrequenciesV, sgolayK, sgolayF, syncMethod, additionalShiftSeconds, usedPressureChan, vidArtOrBart);
    dynamicDataSingleSampleT.MouseID      = repmat(databaseT.MouseID(iDataLine), length(dynamicProtocolPressuresV), 1);
    dynamicDataSingleSampleT.ExperimentID = repmat(databaseT.ExperimentID(iDataLine), length(dynamicProtocolPressuresV), 1);
    dynamicDataSingleSampleT.SampleID     = repmat(databaseT.SampleID(iDataLine), length(dynamicProtocolPressuresV), 1);
    dynamicDataT = [dynamicDataT;dynamicDataSingleSampleT];
    
    
    %--STATIC PD DATA LOADING--
    %Determine folder/file name base part
    fullStaticPDLVFileNameBase = fullfile(folderName, ...
      [databaseT.BaseFolderName{iDataLine} '_Labview'], ...
      databaseT.PressureSweepLVNameBase{iDataLine});
    fullStaticPDVEVOFileNameBase = fullfile(folderName, ...
      [databaseT.BaseFolderName{iDataLine} '_VEVO'], ...
      databaseT.PressureSweepVEVONameBase{iDataLine});
    
    
    staticPDDataSingleSampleT = loadSingleSampleStaticPDData(fullStaticPDLVFileNameBase, fullStaticPDVEVOFileNameBase, staticPDProtocolStretchesV, vidArtOrBart);
    staticPDDataSingleSampleT.MouseID      = repmat(databaseT.MouseID(iDataLine), length(staticPDProtocolStretchesV), 1);
    staticPDDataSingleSampleT.ExperimentID = repmat(databaseT.ExperimentID(iDataLine), length(staticPDProtocolStretchesV), 1);
    staticPDDataSingleSampleT.SampleID     = repmat(databaseT.SampleID(iDataLine), length(staticPDProtocolStretchesV), 1);
    staticPDDataSingleSampleT.UnloadedLength = repmat(databaseT.UnloadedLength_mm_(iDataLine), length(staticPDProtocolStretchesV), 1);
    staticPDDataSingleSampleT.InVivoLength   = repmat(databaseT.InVivoLength_mm_(iDataLine),   length(staticPDProtocolStretchesV), 1);
    l0 = databaseT.UnloadedLength_mm_(iDataLine);
    staticPDDataSingleSampleT.AxialStretchV  = cellfun(@(a) (a/1000+l0)./l0, staticPDDataSingleSampleT.ElongationV, 'UniformOutput', false);
    staticPDDataSingleSampleT.InVivoStretch  = staticPDDataSingleSampleT.InVivoLength./staticPDDataSingleSampleT.UnloadedLength;
    
    staticPDDataT = [staticPDDataT;staticPDDataSingleSampleT];
    
    %--STATIC FL DATA LOADING--
    fullStaticFLLVFileNameBase = fullfile(folderName, ...
      [databaseT.BaseFolderName{iDataLine} '_Labview'], ...
      databaseT.ForceSweepLVFileNameBase{iDataLine});
    
    staticFLDataSingleSampleT = loadSingleSampleStaticFLData(fullStaticFLLVFileNameBase, staticFLProtocolPressuresV);
    staticFLDataSingleSampleT.MouseID        = repmat(databaseT.MouseID(iDataLine),            length(staticFLProtocolPressuresV), 1);
    staticFLDataSingleSampleT.ExperimentID   = repmat(databaseT.ExperimentID(iDataLine),       length(staticFLProtocolPressuresV), 1);
    staticFLDataSingleSampleT.SampleID       = repmat(databaseT.SampleID(iDataLine),         length(staticFLProtocolPressuresV), 1);
    staticFLDataSingleSampleT.UnloadedLength = repmat(databaseT.UnloadedLength_mm_(iDataLine), length(staticFLProtocolPressuresV), 1);
    staticFLDataSingleSampleT.InVivoLength   = repmat(databaseT.InVivoLength_mm_(iDataLine),   length(staticFLProtocolPressuresV), 1);
    l0 = databaseT.UnloadedLength_mm_(iDataLine);
    staticFLDataSingleSampleT.AxialStretchV  = cellfun(@(a) (a/1000+l0)./l0, staticFLDataSingleSampleT.ElongationV, 'UniformOutput', false);
    staticFLDataSingleSampleT.InVivoStretch  = staticFLDataSingleSampleT.InVivoLength./staticFLDataSingleSampleT.UnloadedLength;
    
    staticFLDataT = [staticFLDataT;staticFLDataSingleSampleT];
  end
end

