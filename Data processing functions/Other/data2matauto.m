function data2matauto(bOnlyEEG,bOnlyBlood,bOnlyBP,dFILENAME,dPATHNAME,LengthOfDataSetInMinutes)
%Convert .DATA format from M-PAQ to .MAT file for enabling processing.
%March 3rd 2005
%Rik Vullings
%Last change:

%function DataToMat

%Variable names based on "IdeeQ recording file format" by Iwan de Jong.
%Most of algorithm by P. Aelen (change is the fact that the data file only
%consists of active channels, so at first is determined which channels are
%active and only for those channels the data is existing and converted).

%To prevent the computer from running out of memory
% CurrentDirectory=cd;
% DataDirectory='F:\AZM Schapenstudie\Reint';
% cd(DataDirectory);
%Temp
%DataDirectory='D:\Measurements\Data Judith';

OpenFileName=[dPATHNAME '\' dFILENAME];
SaveFileName=[dPATHNAME '\' dFILENAME(1:end-4),'MAT']; %minus 4 because of DATA
if bOnlyEEG==true
    SaveFileName=[dPATHNAME '\' dFILENAME(1:end-5),'_EEG.MAT'];
end
if bOnlyBlood==true
    SaveFileName=[dPATHNAME '\' dFILENAME(1:end-5),'_BLOOD.MAT'];
end
if bOnlyBP==true
    SaveFileName=[dPATHNAME '\' dFILENAME(1:end-5),'_BP.MAT'];
end
%Open file to convert

% if exist(SaveFileName)~=0
%     fprintf('File already exists. ');
%     overwrite=input('Overwrite? [y/n]: ','s');
%     if strcmp(overwrite,'n')==1
%         return
%     else
%         delete(SaveFileName)
%     end
% end

fid=fopen(OpenFileName,'r');
% cd(CurrentDirectory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Read header file for information      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sVersion=fread(fid,21,'char')'; sVersion=char(sVersion(2:end));
sIdentification=fread(fid,11,'char')'; sIdentification=char(sIdentification(2:end));
sDescription=fread(fid,41,'char')'; sDescription=char(sDescription(2:end));
sRemarks=fread(fid,101,'char')'; sRemarks=char(sRemarks(2:end));
sDate=fread(fid,11,'char')'; sDate=char(sDate(2:end));
sTime=fread(fid,9,'char')'; sTime=char(sTime(2:end));
sNumChannels=fread(fid,4,'int8')'; sNumChannels=str2double(char(sNumChannels(2:end)));
sHdrSize=fread(fid,7,'char')'; sHdrSize=str2double(char(sHdrSize(2:end)));
sNumBlocks=fread(fid,7,'char')'; sNumBlocks=str2double(char(sNumBlocks(2:end)));
sBlockSize=fread(fid,7,'char')'; sBlockSize=str2double(char(sBlockSize(2:end)));
bReserved0=fread(fid,37,'int8')';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Read record data header       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sName=zeros(sNumChannels,21);
sSampleRate=zeros(sNumChannels,1);
sNumSamplesInBlocks=zeros(sNumChannels,1);
sNumSamples=zeros(sNumChannels,1);
bReserved1=zeros(sNumChannels,202);
for h=1:sNumChannels
    sName(h,:)=fread(fid,21,'char')'; sName(h,:)=char(sName(h,:));
    temp=fread(fid,11,'char')'; sSampleRate(h)=str2double(char(temp(2:end)));
    temp=fread(fid,11,'char')'; sNumSamplesInBlocks(h)=str2double(char(temp(2:end)));
    temp=fread(fid,11,'char')'; sNumSamples(h)=str2double(char(temp(2:end)));     
    bReseverd1(h,:)=fread(fid,202,'char')';
end
sName=char(sName); sName(:,1)=[]; %clear first bit as this contains no information


MaximumAmountOfBlocksInOneDataSet=ceil(60*LengthOfDataSetInMinutes*max(sSampleRate)/max(sNumSamplesInBlocks));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Input options recording       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bOn=zeros(sNumChannels,1);
sChannelName=zeros(sNumChannels,20);
iChannelNameLength=zeros(sNumChannels,1);
iScanRate=zeros(sNumChannels,1);
iUnitLength=zeros(sNumChannels,1);
sUnit=zeros(sNumChannels,8);
fHardwareAmplifier=zeros(sNumChannels,1);
fGain=zeros(sNumChannels,1);
bAutoOffset=zeros(sNumChannels,1);
fOffset=zeros(sNumChannels,1);
fLPF=zeros(sNumChannels,1);
fHPF=zeros(sNumChannels,1);
%sInputOptionsRecording=zeros(sNumChannels,511);
for h=1:sNumChannels
    bOn(h)=fread(fid,1,'int8')';
    if bOnlyEEG==true
        %if (bOn(h)==1) && (isempty(strfind(sName(h,:),'Left'))==0 ||
        %isempty(strfind(sName(h,:),'left'))==0 || isempty(strfind(sName(h,:),'Right'))==0 || isempty(strfind(sName(h,:),'right'))==0)
        if (bOn(h)==1) && (isempty(strfind(sName(h,:),'EEG RA'))==0 || isempty(strfind(sName(h,:),'EEG RP'))==0 || isempty(strfind(sName(h,:),'EEG LA'))==0 || isempty(strfind(sName(h,:),'EEG LP'))==0 || isempty(strfind(sName(h,:),'EEG R ant'))==0 || isempty(strfind(sName(h,:),'EEG R post'))==0 || isempty(strfind(sName(h,:),'EEG L ant'))==0 || isempty(strfind(sName(h,:),'EEG L post'))==0)
            bOn(h)=1;
%         elseif (bOn(h)==1) && (isempty(strfind(sName(h,:),'EEG left'))==0 ||isempty(strfind(sName(h,:),'EEG right'))==0)
%             bOn(h)=1;
        else
            bOn(h)=0;
        end
    end
    if bOnlyBlood==true
        %if (bOn(h)==1) && (isempty(strfind(sName(h,:),'Left'))==0 || isempty(strfind(sName(h,:),'left'))==0 || isempty(strfind(sName(h,:),'Right'))==0 || isempty(strfind(sName(h,:),'right'))==0)
        if (bOn(h)==1) && (isempty(strfind(sName(h,:),'ECG'))==0 || isempty(strfind(sName(h,:),'corrected'))==0 || isempty(strfind(sName(h,:),'arterial'))==0 || isempty(strfind(sName(h,:),'amniotic'))==0)
        %if (bOn(h)==1) && (isempty(strfind(sName(h,:),'arterial'))==0 || isempty(strfind(sName(h,:),'amniotic'))==0)
            bOn(h)=1;
        else
            bOn(h)=0;
        end
    end
    if bOnlyBP==true     
        if (bOn(h)==1) && (isempty(strfind(sName(h,:),'corrected'))==0 || isempty(strfind(sName(h,:),'arterial'))==0 || isempty(strfind(sName(h,:),'amniotic'))==0)
            bOn(h)=1;
        else
            bOn(h)=0;
        end
    end 
        
    iChannelNameLength(h,1)=fread(fid,1,'int8'); sChannelNameTemp=[];
    sChannelNameTemp=fread(fid,20,'char')';sChannelName(h,1:iChannelNameLength(h,1))=sChannelNameTemp(1:iChannelNameLength(h,1));
    iIgnore=fread(fid,2,'int8'); %Ignore two bytes
    iScanRate(h,1)=fread(fid,1,'int32'); 
    iIgnore=fread(fid,1,'int8');
    iUnitLength(h,1)=fread(fid,1,'int8'); sUnitTemp=[];
    sUnitTemp=fread(fid,8,'char')'; sUnit(h,1:iUnitLength(h,1))=sUnitTemp(1:iUnitLength(h,1));
    iIgnore=fread(fid,2,'int8'); %Ignore two bytes
    fHardwareAmplifier(h,1)=fread(fid,1,'float64');
    fGain(h,1)=fread(fid,1,'float64');
    bAutoOffset(h,1)=fread(fid,1,'int8');
    iIgnore=fread(fid,7,'int8');
    fOffset(h,1)=fread(fid,1,'float64');
    fLPF(h,1)=fread(fid,1,'float64');
    fHPF(h,1)=fread(fid,1,'float64');
    iIgnore=fread(fid,424,'int8');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        AVI options recording       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sAviOptionsRecord=fread(fid,512,'char')';
bReserved2=fread(fid,7*512,'char')';
Temp=fread(fid,512,'char')'; %these are additional bits and not mentioned in the file specification



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        Read in data        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iActiveChannels=find(bOn(:)==1);
sChannelName=sChannelName(iActiveChannels(:),:);
iChannelNameLength=iChannelNameLength(iActiveChannels(:));
iScanRate=iScanRate(iActiveChannels(:));
iUnitLength=iUnitLength(iActiveChannels(:));
sUnit=sUnit(iActiveChannels(:),:);
fHardwareAmplifier=fHardwareAmplifier(iActiveChannels(:));
fGain=fGain(iActiveChannels(:));
bAutoOffset=bAutoOffset(iActiveChannels(:));
fOffset=fOffset(iActiveChannels(:));
fLPF=fLPF(iActiveChannels(:));
fHPF=fHPF(iActiveChannels(:));

fTimeBeginSubset=0;
if sNumBlocks<=MaximumAmountOfBlocksInOneDataSet
    %Determine number of active channels
    NumberOfChannelsActive=sum(bOn);
    %x=zeros(NumberOfChannelsActive,max(sNumSamples));
    x=zeros(max(sNumSamples),NumberOfChannelsActive);
    position=1;

    for h=1:sNumBlocks
        channel=1;
        for k=1:sNumChannels
            if bOn(k)==1;
                    data=fread(fid,sNumSamplesInBlocks(k),'float32')'; %The last data block might not be filled completely, but this will result in zeros added to the end.
                    data=rot90(data,3);
                    %x(channel,position:position+sNumSamplesInBlocks(k)-1)=data; %analog amplification of factor 10 for channels 1-8: not in software, so has to be rescaled                    
                    x(position:position+sNumSamplesInBlocks(k)-1,channel)=data;
                    channel=channel+1;                 
            elseif sNumSamples(k)~=0 %Channels present but empty would give problems in data processing.
                    data=fread(fid,sNumSamplesInBlocks(k),'float32')';  
            end
        end
        position=position+sNumSamplesInBlocks(k);
    end
    
    %Find correct channels names
    ChannelsUsed=zeros(NumberOfChannelsActive,size(sName,2));
    channel=1;
    for h=1:sNumChannels
        if bOn(h)==1
            ChannelsUsed(channel,:)=sName(h,:);
            channel=channel+1;
        end
    end
    ChannelsUsed=char(ChannelsUsed);
    
    %Remove possible zeros at end of x
    %x=x(:,1:max(sNumSamples));
    x=x(1:max(sNumSamples),:);
    Fs=sSampleRate(1); %Sample rate is the same for all channels
    

    
    
    %Save data
    SaveFileName=[SaveFileName(1:end-4) '_1.MAT'];
save(SaveFileName,'fTimeBeginSubset','sIdentification','sDate','sTime','ChannelsUsed','Fs','x','iScanRate','sUnit','fHardwareAmplifier','fGain','bAutoOffset','fOffset','fLPF','fHPF');
    
    fprintf('File converted successfully\n');
    %fclose(fid);

else
    OriginalSaveFileName=SaveFileName;
        
    %Determine number of data sets to divide the complete measurement into
    NumberOfDataSets=ceil(sNumBlocks/MaximumAmountOfBlocksInOneDataSet);

    for SetCounter=1:NumberOfDataSets
        position=1;
        
        %if SetCounter~=NumberOfDataSets
            %Determine number of active channels
            NumberOfChannelsActive=sum(bOn);
            %x=zeros(NumberOfChannelsActive,MaximumAmountOfBlocksInOneDataSet*max(sNumSamplesInBlocks));            
            x=zeros(MaximumAmountOfBlocksInOneDataSet*max(sNumSamplesInBlocks),NumberOfChannelsActive);            
        %else
        %    %For the last dataset x can be smaller (in the second
        %    %dimension). Has to be changed in code though
        %    NumberOfChannelsActive=sum(bOn);
        %    %x=zeros(NumberOfChannelsActive,MaximumAmountOfBlocksInOneDataSet*max(sNumSamplesInBlocks));            
        %    x=zeros(MaximumAmountOfBlocksInOneDataSet*max(sNumSamplesInBlocks),NumberOfChannelsActive);            
        %end           
        
        for h=(SetCounter-1)*MaximumAmountOfBlocksInOneDataSet+1:min(SetCounter*MaximumAmountOfBlocksInOneDataSet,sNumBlocks)
            channel=1;

            for k=1:sNumChannels
                if bOn(k)==1;                    
                    data=fread(fid,sNumSamplesInBlocks(k),'float32')'; %The last data block might not be filled completely, but this will result in zeros added to the end.
                    data=rot90(data,3);
                    %x(channel,position:position+sNumSamplesInBlocks(k)-1)=data; %analog amplification of factor 10 for channels 1-8: not in software, so has to be rescaled
                    x(position:position+sNumSamplesInBlocks(k)-1,channel)=data;
                    channel=channel+1;
                elseif sNumSamples(k)~=0
                    data=fread(fid,sNumSamplesInBlocks(k),'float32')';
                end
            end
            position=position+sNumSamplesInBlocks(k);
        end
        
        %Find correct channels names
        ChannelsUsed=zeros(NumberOfChannelsActive,size(sName,2));
        channel=1;
        for h=1:sNumChannels
            if bOn(h)==1
                ChannelsUsed(channel,:)=sName(h,:);
                channel=channel+1;
            end
        end
        ChannelsUsed=char(ChannelsUsed);
        
        %The last data set can be zero-padded. To save disk space, clear
        %this
        if SetCounter==NumberOfDataSets
            LengthOfLastDataSet=max(sNumSamples)-((NumberOfDataSets-1)*MaximumAmountOfBlocksInOneDataSet*max(sNumSamplesInBlocks)+1);
            %x=x(:,1:LengthOfLastDataSet);
            x=x(1:LengthOfLastDataSet,:);
        end
        
                
        Fs=sSampleRate(1); %Sample rate is the same for all channels
        
        %Save data
		
        SaveFileName=[OriginalSaveFileName(1:end-4) '_' num2str(SetCounter) '.MAT'];
        save(SaveFileName,'fTimeBeginSubset','sIdentification','sDate','sTime','ChannelsUsed','Fs','x','iScanRate','sUnit','fHardwareAmplifier','fGain','bAutoOffset','fOffset','fLPF','fHPF');
        
        fTimeBeginSubset=fTimeBeginSubset+(position-1)/Fs;
        
        fprintf('Dataset %2.0f converted successfully\n',SetCounter);           
    
    end
end
fclose(fid);