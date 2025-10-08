function [ecgDataV, respirationDataV, temperatureDataV, bloodPressureDataV, sampleTimeStampsTicksV, sampleTimeStampsMillisecondsV] = ...
    readVEVORawPhysio(fNameRawPhysio)
  %READVEVORAWPHYSIO Function to read VEVO *.raw.physio files.
  % [ecgDataV, respirationDataV, temperatureDataV, bloodPressureDataV, ...
  %  sampleTimeStampTicksV, sampleTimeStampMillisecondsV] = ...
  %    READVEVORAWPHYSIO(fNameRawPhysio) Read data from a VEVO
  %       .raw.physio file whose full filename is specified in
  %       fNameRawPhysio.
  % ecgDataV, respirationDataV, temperatureDataV, and bloodPressureDataV
  %   are nSamples*1 int16 arrays that contain (unscaled) ECG, respiration,
  %   temperature and blood pressure values, respectively.
  % sampleTimeStampsTicksV is an nSamples*1 uint32 vector of frame times in
  %   ticks of a 400,000 Hz clock.
  % sampleTimeStampsMillisecondsV is an nSamples*1 double vector of frame
  %   times in milliseconds.
  %
  % V1.0 by Bart Spronck.

  fid = fopen(fNameRawPhysio);
  
  %Read file header
  dwVersion   = fread(fid,1,'uint32'); %Version number of raw data file.
  %                                     Current version is 3.
  dwNumFrames = fread(fid,1,'uint32'); %Number of frames in this raw data
  %                                     file.
  dwInfo      = fread(fid,1,'uint32');  %Information bitfield used to identify
  %                                     the type of data in the frame
  temp        = fread(fid,7,'uint32'); %Unused bytes for future use.
  
  if dwVersion ~= 3
    error('Expected file version 3!')
  end
  %disp(['Number of "frames" in physio data file: ' num2str(dwNumFrames)])
  if dwInfo ~= 0
    error('Expected dwInfo==0')
  end
  
  %Initialise arrays to save data
  %In physio data, the number of samples per frame may vary. Therefore, we
  %do not know the final array size on forehand.
  %Since MATLAB appending (e.g. a = [a; newSample]) is very slow, the
  %functions fastAppendIni, fastAppendAdd and fastAppendFinish are used,
  %which pre-allocate an array.
  
  fastAppendBlockSize = 256*1000;
  [physioDataM, physioDataCountI]   = fastAppendIni(fastAppendBlockSize, 4, 'int16');
  [sampleTimeStampsTicksV, stCountI] = fastAppendIni(fastAppendBlockSize, 1, 'uint32');
  
  for nFrame = 1:double(dwNumFrames)
    %Read frame header
    dwTimeStamp = fread(fid,1,'uint32');   %Hardware time stamp counted as ticks of a 400,000 Hz clock
    dbTimeStamp = fread(fid,1,'double');   %Hardware time stamp counted in ms. Calculated as dwTimeStamp * 1000 / 400000
    dwFrameNumber = fread(fid,1,'uint32'); %Frame number, starting from 0
    dwInfo = fread(fid,1,'uint32');        %Is 0 if frame contains valid data
    dwPacketSize = fread(fid,1,'uint32');  %Size in bytes of frame data (not including header)
    temp = fread(fid,8,'uint32');
    
    nSamplesPerChannel = dwPacketSize/2/4;
    
    ecgFrameDataV = fread(fid,nSamplesPerChannel,'int16');
    respirationFrameDataV = fread(fid,nSamplesPerChannel,'int16');
    temperatureFrameDataV = fread(fid,nSamplesPerChannel,'int16');
    bloodPressureFrameDataV = fread(fid,nSamplesPerChannel,'int16');
    
    frameSampleTimeStampTicksV = uint32(  (0:(nSamplesPerChannel-1))*50 + dwTimeStamp  )';
    
    physioFrameDataM = [ecgFrameDataV, respirationFrameDataV, temperatureFrameDataV, bloodPressureFrameDataV];
    
    
    [physioDataM, physioDataCountI]   = fastAppendAdd(physioDataM,physioDataCountI,fastAppendBlockSize,physioFrameDataM);
    [sampleTimeStampsTicksV, stCountI] = fastAppendAdd(sampleTimeStampsTicksV, stCountI, fastAppendBlockSize, frameSampleTimeStampTicksV);
  end
  
  fclose(fid);
  physioDataM = fastAppendFinish(physioDataM,physioDataCountI);
  ecgDataV           = physioDataM(:,1);
  respirationDataV   = physioDataM(:,2);
  temperatureDataV   = physioDataM(:,3);
  bloodPressureDataV = physioDataM(:,4);
  sampleTimeStampsTicksV = fastAppendFinish(sampleTimeStampsTicksV,stCountI);
  sampleTimeStampsMillisecondsV = double(sampleTimeStampsTicksV) * 1000/400000;
end