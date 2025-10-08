function staticFLDataT = loadSingleSampleStaticFLData(fullStaticFLLVFileNameBase, pressuresV, usedPressureChan)
  % fullStaticFLVEVOFileNameBase to be added when needed
  
  nSteps = length(pressuresV);
  forceVC  = cell(nSteps, 1);
  elongationVC = cell(nSteps, 1);
  presVC = cell(nSteps, 1);
  for iProtocol = 1:nSteps %Typically 1:5 (pressures 10 60 100 140 200)
    p = pressuresV(iProtocol);
    fprintf('P=%.0f\n', p)
    fileNameLV =   [fullStaticFLLVFileNameBase, 'p', num2str(p), '.tdms'];
    
%     %Find complete VEVO filename
%     %VEVO filenames have dates appended... find the actual file name
%     fileNameVEVO = [fullStaticFLVEVOFileNameBase, sprintf('%.2f', p)];
%     
%     dirC = dir([fileNameVEVO '*.DI.mat']);
%     proceed = true;    
%     if length(dirC) > 1
%       warning('Multiple files starting with %s exist! Skipping...', fileNameVEVO)
%       proceed = false;
%     elseif isempty(dirC)
%       warning('No .di.mat VEVO file starting with %s found!', fileNameVEVO)
%       proceed = false;
%     end
%     [vevoPath,~,~] = fileparts(fullStaticFLVEVOFileNameBase);
%     fileNameVEVO = fullfile(vevoPath, dirC.name);
    
    proceed = true;
    if ~exist(fileNameLV, 'file')
      warning('LabVIEW File not found: %s, fileNameLV')
      proceed = false;
    end
%     if ~exist(fileNameVEVO, 'file')
%       warning('VEVO file not found: %s, fileNameVEVO')
%       proceed = false;
%     end
    
    if proceed
      %[forceVC{iProtocol}, elongationVC{iProtocol}] = getStaticFL(fileNameLV);
      [forceVC{iProtocol}, elongationVC{iProtocol}, presVC{iProtocol}] = getStaticFL(fileNameLV);
    else
      forceVC{iProtocol} = NaN;
      elongationVC{iProtocol} = NaN;
      presVC{iProtocol} = NaN;
    end
  end
  
  staticFLDataT = table(pressuresV, forceVC, elongationVC, presVC, ...
    'VariableNames', {'ProtocolPressure', 'ForceV', 'ElongationV', 'PressureV'});
end