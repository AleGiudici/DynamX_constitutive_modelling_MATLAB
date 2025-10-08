function [frameData3M, frameTimeStampsTicksV, frameTimeStampsMilliSecondsV, pixelSizeV, sysParsS] = readVEVORawBMode(fNameRawBMode, nFrames)
  %READVEVORAWBMODE Function to read VEVO *.raw.bmode files.
  % [frameData3M, frameTimeStampsTicksV, ...
  %  frameTimeStampsMilliSecondsV, sysParsS] = ...
  %    READVEVORAWBMODE(fNameRawBMode, fNameRawXML) Read data from a VEVO
  %      .raw.bmode file whose full filename is specified in
  %      fNameRawBMode. 
  %      nFrames is an optional parameter specifying the number of frames
  %      to read. If it is omitted, all frames are read.
  % frameData3M is an nSamples*nLines*nFrames uint8 matrix of image frames.
  % frameTimeStampsTicksV is an nFrames*1 uint32 vector of frame times in
  %   ticks of a 400,000 Hz clock.
  % frameTimeStampsMilliSecondsV is an nFrames*1 double vector of frame
  %   times in milliseconds.
  % pixelSizeV is a 2*1 double vector that denotes the pixel size in
  %   micrometer of the images in frameData3M.
  %   pixelSizeV(1) is the pixel size in depth direction,
  %   pixelSizeV(2) is the pixel size in horizontal direction.
  %   Note that the images in frameData3M are of the dimensions
  %   samples*lines, i.e., displaying them without scaling will cause the
  %   images to be stretched.
  % sysParsS is a structure of system parameters obtained from the .raw.xml
  %   file (for specification see the readVevoRawXML function).
  %
  % V1.1 by Bart Spronck.
  
  
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
  
  % Calculate pixel size
  pixelSizeV = zeros(2,1);
  pixelSizeV(1) = (maxDepth - minDepth) * 1000 / nSamples;
  pixelSizeV(2) = width * 1000 / nLines;
  
  %% Read binary data
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
  
  if exist('nFrames', 'var')
    if nFrames > dwNumFrames
      fprintf(['The requested number of frames (%i) is larger than the available number of frames (%i). Only using available frames.' char(13)], nFrames, dwNumFrames)
      nFrames = dwNumFrames;
    end
  else
    nFrames = dwNumFrames;
    %fprintf(['Requested number of frames is not specified. Reading all frames (%i).' char(13)], nFrames)
  end  
  
  %Initialise arrays to save frame data
  frameData3M = zeros(nSamples, nLines, nFrames, 'uint8');
  frameTimeStampsTicksV = zeros(nFrames,1, 'uint32');
  frameTimeStampsMilliSecondsV = zeros(nFrames,1, 'double');
  

  
  for nFrame = 1:nFrames
    %Read frame header
    dwTimeStamp = fread(fid,1,'uint32');   %Hardware time stamp counted as ticks of a 400,000 Hz clock
    dbTimeStamp = fread(fid,1,'double');   %Hardware time stamp counted in ms. Calculated as dwTimeStamp * 1000 / 400000
    dwFrameNumber = fread(fid,1,'uint32'); %Frame number, starting from 0
    dwInfo = fread(fid,1,'uint32');        %Is 0 if frame contains valid data
    dwPacketSize = fread(fid,1,'uint32');  %Size in bytes of frame data (not including header)
    temp = fread(fid,8,'uint32');
    
    if dwInfo ~= 0
      error('Frame does not contain valid data')
    end
    
    if dwPacketSize ~= nLines*nSamples
      error('Corrupt frame found')
    end
    
    %Read frame image data into an uint8 1D array
    frameData = uint8(fread(fid,dwPacketSize,'uint8'));
    
    %Process data
    frameDataM = reshape(frameData, [nSamples nLines]); %Reshape to frame
    %dimensions.
    
    frameData3M(:,:,nFrame) = frameDataM;
    frameTimeStampsTicksV(nFrame)        = dwTimeStamp;
    frameTimeStampsMilliSecondsV(nFrame) = dbTimeStamp;
  end
  
  fclose(fid);
  
end