function longSyncPulseStartTimeMilliseconds = findWideSyncPulse(syncSignalV, sampleTimeStampsMillisecondsV, expectedPulseDurationMilliseconds, allowMultipleLongPulses)
  % If allowMultipleLongPulses==true, the LAST long sync pulse will be
  % returned.
  syncThreshold = max(syncSignalV)/2;
  
  aV = double(syncSignalV > syncThreshold);
  dV = diff(aV);
  syncStartIndicesV = find(dV==1);
  syncStartIndicesV = syncStartIndicesV + 1; %The sync index is the index
  
  syncEndIndicesV = find(dV==-1);
  syncEndIndicesV = syncEndIndicesV + 1; %The sync index is the index
  %of the first value that exceeds the threshold.
  
  if syncEndIndicesV(1) < syncStartIndicesV(1)
    syncEndIndicesV = syncEndIndicesV(2:end);
  end
  
  if length(syncStartIndicesV) > length(syncEndIndicesV)
    syncStartIndicesV = syncStartIndicesV(1:end-1);
  end
  
  if length(syncStartIndicesV) ~= length(syncEndIndicesV)
    error('This can''t be true...')
  end
  
  syncStartTimesV = sampleTimeStampsMillisecondsV(syncStartIndicesV);
  syncEndTimesV = sampleTimeStampsMillisecondsV(syncEndIndicesV);
  
  syncDurationsV = syncEndTimesV - syncStartTimesV;
  
  iV = find((syncDurationsV > 0.90*expectedPulseDurationMilliseconds) & (syncDurationsV < 1.10*expectedPulseDurationMilliseconds));
  if allowMultipleLongPulses
    if length(iV) > 1
      disp('WARNING: Multiple long sync pulses detected, using the LAST one!')
      iV = iV(end);
    end
  else
    if length(iV) > 1
      error('Multiple sync pulses detected!')
    end
  end
  if isempty(iV)
    warning('No long sync pulses detected!')
    longSyncPulseStartTimeMilliseconds = NaN;
  else
    longSyncPulseStartTimeMilliseconds = sampleTimeStampsMillisecondsV(syncStartIndicesV(iV));
  end
  
end

