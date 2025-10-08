function createStaticFrameSet(sourceFolderName,destFolderName,fileNamePrefix,minPressure,maxPressure,stepPressure,goBack)
  %UNTITLED2 Summary of this function goes here
  %   Detailed explanation goes here
  d = dir([sourceFolderName '\*.raw.bmode']);
  fileNamesListC = {d.name}';
  pV = minPressure:stepPressure:maxPressure;
  pC = arrayfun(@(a) sprintf('%i',a),pV,'UniformOutput', false);
  
  if goBack
    p2V = (maxPressure-stepPressure):-stepPressure:minPressure;
    p2C = arrayfun(@(a) sprintf('b%i',a),p2V,'UniformOutput', false);
  end
  
  pFullC = [pC p2C];
  
  outMIsInitialised = false;
  for i = 1:length(pFullC)
    p = pFullC{i};
    fileNameSought = sprintf('%s%s-', fileNamePrefix, p);
    l = length(fileNameSought);
    
    correspondingIndices = cellfun(@(a) cmp(a, l, fileNameSought), fileNamesListC);
    if sum(correspondingIndices) == 0
      error(['File ' fileNameSought '* not found'])
    elseif sum(correspondingIndices) > 1
      disp(['Found multiple files for ' fileNameSought '*:'])
      error('Handling of that not implemented yet')
    else
      fileName = fileNamesListC{correspondingIndices};
    end
    
    %fullFileNameBMode = fullfile(sourceFolderName,fileName);
    fileNameNoExt = fileName(1:end-6);
    
    extC = {'bmode', 'xml', 'event', 'physio'};
    
    for j = 1:length(extC)
      fullFileNameSource = fullfile(sourceFolderName, [fileNameNoExt '.' extC{j}]);
      fullFileNameDest = fullfile(destFolderName, sprintf('%03i %s.%s',i, fileNameNoExt, extC{j}));
      copyfile(fullFileNameSource, fullFileNameDest)
    end
  end
end

function c = cmp(a, lMin, b)
  if length(a) >= lMin
    c = strcmp(a(1:lMin),b);
  else
    c = false;
  end
end