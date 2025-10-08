function [averagedData] = ensambleAveragingDynamicData(rawDataMat,HR)
% ensambleAveragingDynamicData(rawDataMat,HR) creates the ensable averaged
% version of the experimental (harmonic oscillatory test) data in the input
% matrix "rawDataMat". The input heart rate "HR" is used to ensure that
% irrespectively of the loading frequency, the ensable averaged matrix is
% made of a fixed number of datapoints.

    samplingTime = rawDataMat(2,1)-rawDataMat(1,1);

%% Segment heartbeats
    [~,sinePeaksLocator]=findpeaks(max(rawDataMat(:,3))-rawDataMat(:,3),'MinPeakHeight',...
                0.85*max(max(rawDataMat(:,3))-rawDataMat(:,3)),'MinPeakDistance',0.5/HR/samplingTime);

%% Ensable averaging
    duration=diff(sinePeaksLocator);
    dataMat=zeros(min(duration),7);

    nBeats=length(sinePeaksLocator)-1;
    for i=1:nBeats
        dataMat=dataMat+rawDataMat(sinePeaksLocator(i):sinePeaksLocator(i)+min(duration)-1,:)/nBeats;
    end
    dataMat(:,1)=rawDataMat(sinePeaksLocator(1):sinePeaksLocator(1)+min(duration)-1,1);

%% Resampling to the same number of datapoints
    resampling_time = (dataMat(1,1):samplingTime*10/HR*10:dataMat(end,1))';

    averagedData = [resampling_time interp1(dataMat(:,1),dataMat(:,2:end),resampling_time)];
    
    averagedData(:,1) = averagedData(:,1)-averagedData(1,1);
end