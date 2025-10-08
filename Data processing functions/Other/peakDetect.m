function [maxPosV, minPosV] = peakDetect(dataV)
  thr1 = nanmean(dataV);
  thr2 = min(dataV);
  thr = thr2 + (thr1-thr2)*0.5;
  
  lV = double(dataV > thr); %lV: largerThanVector
  difV = diff(lV);
  upslopePoints = find(difV == 1);
  
  %Check whether this are true upslopes.
  %Once, an upslope was detected on a noisy diastolic decay (HR 2 Hz).
  %Make sure we're really going up.
  
  deltaUpslopeVerification = 10; %Check of a point that is 10 points behind the detected upslope is indeed higher than the upslope point.
  %If that point is unavailable, check whether the point that is 10 points
  %before the detected upslope is lower.
  
  upslopeDeltasFirstPoints = upslopePoints;
  upslopeDeltasLastPoints = upslopePoints + deltaUpslopeVerification;
  
  if upslopeDeltasLastPoints(end) > length(dataV)
    upslopeDeltasFirstPoints(end) = upslopePoints(end) - deltaUpslopeVerification;
    upslopeDeltasLastPoints(end) = upslopePoints(end);
  end
  
  trueUpslopes = (dataV(upslopeDeltasLastPoints) - ...
    dataV(upslopeDeltasFirstPoints)) > 0;
  
  upslopePoints = upslopePoints(trueUpslopes);
  
  minPosV = arrayfun(@(a,b) a-1+minindex(dataV(a:b)),upslopePoints(1:end-1), upslopePoints(2:end));
  maxPosV = arrayfun(@(a,b) a-1+maxindex(dataV(a:b)),upslopePoints(1:end-1), upslopePoints(2:end));
  
  doPlot = false;
  if doPlot
    figure
    plot(dataV)
    hold on
    plot(upslopePoints, dataV(upslopePoints),'xk')
    plot(minPosV, dataV(minPosV),'ok')
    plot(maxPosV, dataV(maxPosV),'or')
    plot([1 length(dataV)], [1 1]*thr)
  end
end
function y = minindex(x)
  [~,y] = min(x);
end
function y = maxindex(x)
  [~,y] = max(x);
end