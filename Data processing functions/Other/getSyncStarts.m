function syncStartIndicesV = getSyncStarts(syncSignalV, varargin)
  if nargin == 1
    syncThreshold = max(syncSignalV)/2;
  else
    syncThreshold = varargin{1};
  end
  
  aV = double(syncSignalV > syncThreshold);
  dV = diff(aV);
  syncStartIndicesV = find(dV==1);
  syncStartIndicesV = syncStartIndicesV + 1; %The sync index is the index
  %of the first value that exceeds the threshold.
end

