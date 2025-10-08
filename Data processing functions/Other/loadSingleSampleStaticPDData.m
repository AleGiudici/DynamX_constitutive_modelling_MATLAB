function staticPDDataT = loadSingleSampleStaticPDData(fullStaticPDLVFileNameBase, fullStaticPDVEVOFileNameBase, stretchesV, vidArtOrBart)
  
  nSteps = length(stretchesV);
  diamVC       = cell(nSteps, 1);
  presVC       = cell(nSteps, 1);
  elongationVC = cell(nSteps, 1);
  forceVC = cell(nSteps, 1);
  for iProtocol = 1:nSteps %Typically 1:3
    stretch = stretchesV(iProtocol);
    fprintf('lambda=%.2f\n', stretch)
    fileNameLV =   [fullStaticPDLVFileNameBase, sprintf('%.2f', stretch), '.tdms'];
    
    %Find complete VEVO filename
    %VEVO filenames have dates appended... find the actual file name
    fileNameVEVO = [fullStaticPDVEVOFileNameBase, sprintf('%.2f', stretch)];
    
    dirC = dir([fileNameVEVO '*.DI.mat']);
    proceed = true;    
    if length(dirC) > 1
      warning('Multiple files starting with %s exist! Skipping...', fileNameVEVO)
      proceed = false;
    elseif isempty(dirC)
      warning('No .di.mat VEVO file starting with %s found!', fileNameVEVO)
      proceed = false;
    end
    [vevoPath,~,~] = fileparts(fullStaticPDVEVOFileNameBase);
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
      [diamVC{iProtocol}, presVC{iProtocol}, elongationVC{iProtocol}, forceVC{iProtocol}] = getStaticPA(fileNameLV, fileNameVEVO(1:end-7), ... %getSyncedPA doesn't want .di.mat extension...
        vidArtOrBart);
    end
    
    if length(diamVC{iProtocol}) ~= 81
      warning('Static diameterV does not have 81 samples!')
      proceed = false;
    end
    if length(presVC{iProtocol}) ~= 81
      warning('Static pressureV does not have 81 samples!')
      proceed = false;
    end    
    
    if ~proceed
      diamVC{iProtocol} = NaN;
      presVC{iProtocol} = NaN;
      elongationVC{iProtocol} = NaN;
      forceVC{iProtocol} = NaN;
    end
  end
  
  staticPDDataT = table(stretchesV, diamVC, presVC, elongationVC, forceVC, ...
    'VariableNames', {'ProtocolAxialStretch', 'DiameterV', 'PressureV', 'ElongationV', 'ForceV'});
end