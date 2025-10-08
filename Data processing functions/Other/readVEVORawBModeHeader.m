function [nLines, nSamples, nFrames, pixelSizeV, acquisitionDate] = readVEVORawBModeHeader(fNameRawBMode)
  %READVEVORAWBMODEHEADER Function to read VEVO 8-bit RAW B-mode header info.
  % [nLines, nSamples, nFrames, pixelSizeV, acquisitionDate] = ...
  %   READVEVORAWBMODEHEADER(fNameRawBMode) Read and parse header from VEVO
  %   .raw.bmode file and the corresponding .raw.xml file. The full
  %   filename of the .raw.bmode file is specified in fNameRawBMode.
  % nLines is a double scalar specifying the number of US lines.
  % nSamples is a double scalar specifying the number of samples along an
  %   US line.
  % nFrames is the number of acquired frames.
  % pixelSizeV is a 2*1 double vector that denotes the pixel size in
  %   micrometer of the images in frameData3M.
  %   pixelSizeV(1) is the pixel size in depth direction,
  %   pixelSizeV(2) is the pixel size in horizontal direction.
  %   Note that the images in frameData3M are of the dimensions
  %   samples*lines, i.e., displaying them without scaling will cause the
  %   images to be stretched.
  % acquisitionDate is a character array specifying the acquisition date.
  %
  % V1.0 by Bart Spronck.
  
  
  %% Parse XML data and get number of lines and samples per line.
  [a,b] = fileparts(fNameRawBMode);
  fNameRawXML = fullfile(a,[b '.xml']);
  if ~exist(fNameRawXML, 'file')
    error('.raw.xml file was not found at %s', fNameRawXML)
  end
  
  sysParsS = readVEVORawXML(fNameRawXML);
  
  namesC = {sysParsS.name};   %Create cell array with parameter names.
  valuesC = {sysParsS.value}; %Create cell array with parameter values.
  nLines =   str2double(valuesC{strcmp(namesC,'B-Mode/Lines'        )});
  nSamples = str2double(valuesC{strcmp(namesC,'B-Mode/Samples'      )});
  maxDepth = str2double(valuesC{strcmp(namesC,'B-Mode/Depth'        )});
  minDepth = str2double(valuesC{strcmp(namesC,'B-Mode/Depth-Offset' )});
  width    = str2double(valuesC{strcmp(namesC,'B-Mode/Width'        )});
  acquisitionDate =     valuesC{strcmp(namesC,'Acquired-Date'        )};
  
  % Calculate pixel size
  pixelSizeV = zeros(2,1);
  pixelSizeV(1) = (maxDepth - minDepth) * 1000 / nSamples;
  pixelSizeV(2) = width * 1000 / nLines;
  
  %% Read binary data header
  fid = fopen(fNameRawBMode);
  
  %Read file header
  dwVersion   = fread(fid,1,'uint32'); %Version number of raw data file.
  %Current version is 3.
  dwNumFrames = fread(fid,1,'uint32'); %Number of frames in this raw data
  %file.
  dwInfo      = fread(fid,1,'uint32'); %Information bitfield used to identify
  %the type of data in the frame
  %Only dwInfo==8 is supported by this
  %script, corresponding to 8-bit RAW
  %data.
  temp        = fread(fid,7,'uint32'); %Unused bytes for future use.
  
  if dwVersion ~= 3
    error('Expected file version 3!')
  end
  if dwInfo ~= 8
    error('Expected 8-bit raw file! Found dwInfo==%i', dwInfo)
  end
  
  nFrames = double(dwNumFrames);
  fclose(fid);  
end