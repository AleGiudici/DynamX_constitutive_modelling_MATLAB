function dynamicDataT = loadSingleSampleDynamicData(fullDynamicLVFileNameBase, fullDynamicVEVOFileNameBase, pressuresV, frequenciesV, sgolayK, sgolayF, syncMethod, additionalShiftSeconds, usedPressureChan, vidArtOrBart)
  if length(pressuresV) ~= length(frequenciesV)
    error('Number of pressures should equal number of stretches')
  end
  
  nSteps = length(frequenciesV);
  tVC          = cell(nSteps, 1);
  diamVC       = cell(nSteps, 1);
  presVC       = cell(nSteps, 1);
  beatStartIVC = cell(nSteps, 1);
  beatEndIVC   = cell(nSteps, 1);
  for iProtocol = 1:nSteps %Typically 1:5
    p = pressuresV(iProtocol);
    f = frequenciesV(iProtocol);
    fprintf('P=%.0f, f=%.1f\n',p,f)
    fileNameLV =   [fullDynamicLVFileNameBase, 'p', num2str(p), '_f', num2str(f), '.tdms'];
    
    %Find complete VEVO filename
    %VEVO filenames have dates appended... find the actual file name
    fileNameVEVO = [fullDynamicVEVOFileNameBase, 'p', num2str(p), '_f', num2str(f)];
    
    dirC = dir([fileNameVEVO '*.DI.mat']);
    proceed = true;
    if length(dirC) > 1
      warning('Multiple files starting with %s exist!', fileNameVEVO)
      proceed = false;
    elseif isempty(dirC)
      warning('No file found starting with %s!', fileNameVEVO)
      proceed = false;
    end
    [vevoPath,~,~] = fileparts(fileNameVEVO);
    fileNameVEVO = fullfile(vevoPath, dirC.name);
    

    if ~exist(fileNameLV, 'file')
      warning('LabVIEW File not found: %s, fileNameLV')
      proceed = false;
    end
    if ~exist(fileNameVEVO, 'file')
      warning('VEVO file not found: %s, fileNameVEVO')
      proceed = false;
    end
    
    if proceed
      [tVC{iProtocol}, diamVC{iProtocol}, tdmsM, beatStartIVC{iProtocol}, beatEndIVC{iProtocol}] = getSyncedPA(...
        fileNameLV, ...
        fileNameVEVO(1:end-7), ... %getSyncedPA doesn't want .di.mat extension...
        sgolayK, sgolayF, syncMethod, NaN, additionalShiftSeconds, false, vidArtOrBart);
      presVC{iProtocol} = tdmsM(:,usedPressureChan);
    else
      tVC{iProtocol} = NaN;
      diamVC{iProtocol} = NaN;
      beatStartIVC{iProtocol} = NaN;
      beatEndIVC{iProtocol} = NaN;
      presVC{iProtocol} = NaN;
    end
  end
  
  dynamicDataT = table(pressuresV, frequenciesV, tVC, diamVC, presVC, beatStartIVC, beatEndIVC, ...
    'VariableNames', {'ProtocolPressure', 'ProtocolFrequency', 'TimeV', 'DiameterV', 'PressureV', 'BeatStartIV', 'BeatEndIV'});
  

end