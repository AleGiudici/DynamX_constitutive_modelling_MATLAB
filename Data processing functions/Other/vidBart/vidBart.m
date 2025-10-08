function vidBart()
  %VIDBART Function to perform wall tracking in VEVO .raw.bmode files.
  % Code based on concepts from vidArt by Arnold Hoeks.
  %
  % V1.0 by Bart Spronck.
  close all
  
  parS.relativeDepthWindow = 0.2;  %Defines height of tracking window, below and above previously detected position, as a fraction of vessel diameter in frame 0 (total height is thus 2*parS.relativeDepthWindow).
  parS.wallThreshold = 0.8;        %Relative to maximum intensity
  parS.trackingColumnWidth = 300;  %Tracking column width [micrometers]
  parS.columnOverlap = 0.5;        %0.5 indicates half-overlapping columns
  parS.interpFactor = 100;         %Depth interpolation factor. Note: much higher than Arnolds, but yields much better tracking as well.
  parS.interpMethod = 'spline';    %Can be 'linear', 'spline', ... (all options supported by MATLAB's interp1 function).
  
  parS.alsoTrackOuterWall = false; %Also track outer wall. Note: the tracking window is still located around the INNER wall.
  parS.plotEveryFrame = false;     %Update plot during tracking. Takes ages, but may come handy for debugging.
  quickForTest = false; %Only for debugging
  
  [f, p] = uigetfile('*.raw.bmode', 'Select a .raw.bmode VEVO file');
  if ~ischar(f)
    error('No file selected')
  end
  VEVOBModeFileName = fullfile(p,f);

  %Override file name here if needed:
  %fName = 'F:\Bart\Work\Projects\TPLSM-HFUS\Data\VEVO\2018_10_26\2018-10-26-14-10-44_2018.10.26 - Dynamic, 1.0_test2-2018-10-26-14-02-22_1.raw.bmode';
  
  wallColoursM = lines(2); %Define colours to use for plotting near and far walls
  
  %% Read VEVO data
  [frameData3M, frameTimeStampsTicksV, frameTimeStampsMilliSecondsV, pixelSizeV, sysParsS] = readVEVORawBMode(VEVOBModeFileName);
  %pixelSizeV contains pixel size in [depth;horizontal] directions

  nFrames = size(frameData3M,3);
  
  %% Get rectangle for ROI selection and crop data
  figure
  figPosV = get(gcf, 'Position');
  imshow(frameData3M(:,:,1), [0 255])
  axH = gca;
  set(axH,'DataAspectRatio', [pixelSizeV;1]);
  set(axH,'Position',[0,0,1,1]);
  set(gcf, 'Name', 'Select ROI, then double click inside it')
  %set(gcf, 'Units', 'normalized')
  %set(gcf, 'Position', [0.05 0.05 0.90 0.80])
  set(gcf, 'Position', figPosV)
  
  h = imrect(axH, [20, 20, 50, 150]);
  if ~quickForTest
    rectPos = wait(h); %User has to double-click rect
  end
  rectPos = round(rectPos);
  frameDataCropped3M = frameData3M(rectPos(2) + (0:rectPos(4)), rectPos(1) + (0:rectPos(3)), :);
  
  ROIHeight = size(frameDataCropped3M,1); %in pixels
  ROIWidth  = size(frameDataCropped3M,2); %in pixels
  
  %Calculate the number of tracking columns corresponding to the ROI width
  nTrackingColumns = round((ROIWidth*pixelSizeV(2) - parS.trackingColumnWidth)/(parS.trackingColumnWidth*(1-parS.columnOverlap)) + 1);
  fprintf('Number of tracking columns: %.0f\n', nTrackingColumns)
  
  %% Initialise tracking (four points) (similar to Arnold's IniLP(ima))
  figH = figure;
  figPosV = get(gcf, 'Position');
  imH = imshow(frameDataCropped3M(:,:,1), [0 255]);
  set(gcf, 'Position', figPosV)
  hold on
  
  axH = gca;
  set(axH,'DataAspectRatio', [pixelSizeV;1]);
  set(axH,'Position',[0,0,1,1]);
  set(figH,'Pointer','crosshair')
  set(imH, 'HitTest', 'off');
  set(axH,'ButtonDownFcn',@getFourClickedPoints, 'PickableParts', 'all');
  
  setappdata(figH, 'ClickedPointCount', 0);
  setappdata(figH, 'ClickedPointLocationsV', zeros(4,2));
  uiwait();
  clickedPointLocationsV = getappdata(figH, 'ClickedPointLocationsV');
  %We now have the locations of our four clicked points.
  
  %Sort clickedPointLocationsV
  clickedPointLocationsV=sortrows(clickedPointLocationsV);
  for wxp=2:2:4           % organize as ant/post (left) and ant/post (right)
    if clickedPointLocationsV(wxp,2)<clickedPointLocationsV(wxp-1,2)
      tmp=clickedPointLocationsV(wxp,2);
      clickedPointLocationsV(wxp,2)=clickedPointLocationsV(wxp-1,2);
      clickedPointLocationsV(wxp-1,2)=tmp;
    end
  end
  wPosEachLineM=zeros(2,ROIWidth); %Wall position for EACH (ROI) IMAGE LINE (no columning yet!) for FIRST FRAME
  for wall=1:2        % fit straight lines through anterior/posterior points
    xp1=clickedPointLocationsV(wall,1);
    yp1=clickedPointLocationsV(wall,2);
    xp2=clickedPointLocationsV(wall+2,1);
    yp2=clickedPointLocationsV(wall+2,2);
    slope=(yp1-yp2)/(xp1-xp2);
    offs=(yp2*xp1-yp1*xp2)/(xp1-xp2);
    
    lineV = 1:ROIWidth;
    wPosEachLineM(wall,:,1)=slope*(lineV-1)+offs;
  end
  
  %Find approximate diameter from these points
  initialDiameterPixels = mean(diff(wPosEachLineM));
  initialDiameterMicrons = initialDiameterPixels*pixelSizeV(1);
  
  %Determine depthWindow using this the estimated diameter
  depthWindow = round(initialDiameterPixels*parS.relativeDepthWindow*parS.interpFactor);  %in INTERPOLATED (vertical) pixels
  bgWindow = round(initialDiameterPixels*parS.relativeDepthWindow*parS.interpFactor*0.3); %in INTERPOLATED (vertical) pixels
  
  %Plot interpolated 'guesses'
  for iWall = 1:2
    plot(lineV,wPosEachLineM(iWall,:), 'LineWidth', 3, 'Color', wallColoursM(iWall,:))
  end
  set(gcf, 'Name', [sprintf('Average diameter: %.1f',initialDiameterMicrons) ' \mum'])
  trackingColumns3M = zeros(ROIHeight, nTrackingColumns, nFrames);
  wPosColumnsIniM = zeros(2,nTrackingColumns); %2 for 2 walls
  
  %% Generate 'image' sequence with columns for tracking
  nFullColumnsThatFitInROI = (1-parS.columnOverlap)*(nTrackingColumns-1)+1;
  actualColumnWidth = ROIWidth / nFullColumnsThatFitInROI;
  actualColumnWidthMicrons = actualColumnWidth*pixelSizeV(2);
  fprintf('Actual tracking column width: %f µm; parameter setting: %f\n', actualColumnWidthMicrons, parS.trackingColumnWidth)
  
  for iCol = 1:nTrackingColumns
    iLineMin = round(1 + (1-parS.columnOverlap)*(iCol-1)*actualColumnWidth);
    iLineMax = round(    (1-parS.columnOverlap)*(iCol-1)*actualColumnWidth + actualColumnWidth);
    
    trackingColumns3M(:,iCol,:) = mean(frameDataCropped3M(:,iLineMin:iLineMax,:),2);
    wPosColumnsIniM(1,iCol) = mean(wPosEachLineM(1,iLineMin:iLineMax));
    wPosColumnsIniM(2,iCol) = mean(wPosEachLineM(2,iLineMin:iLineMax));
  end
  
  figure
  figPosV = get(gcf, 'Position');
  imshow(trackingColumns3M(:,:,1), [0 255])
  aspectV = [pixelSizeV(1);pixelSizeV(2)*actualColumnWidth*(1-parS.columnOverlap);1];
  axH = gca;
  set(axH,'DataAspectRatio', aspectV);
  set(axH,'Position',[0,0,1,1])
  set(gcf, 'Position', figPosV)
  hold on
  plot(1:nTrackingColumns,wPosColumnsIniM, 'LineWidth', 3)
  set(gcf, 'Name', 'Initial guesses for tracking')
  
  %% Wall positions are now initialized
  %This is where, in arnolds code, these starting positions can be adjusted.
  %Can be implemented in future, but is less critical for us since our
  %vessel is straight.
  
  %% Let's track some walls
  %Detection steps:
  %Find local maximum
  %Find lumen value ('baseline' or 'background')
  %Determine where relative intensity surpasses relative threshold
  %(check for double echo)
  
  innerWallPosColumnsM = zeros(2,nTrackingColumns,nFrames);
  if parS.alsoTrackOuterWall
    outerWallPosColumnsM = zeros(2,nTrackingColumns,nFrames);
  else
    outerWallPosColumnsM = [];
  end
  
  if parS.plotEveryFrame
    figH = figure;
      figPosV = get(gcf, 'Position');
  end
  
  disp('Tracking...')
  for iFrame = 1:nFrames
    %Interpolate data for tracking
    interpolatedFrameM = interp1(1:ROIHeight, trackingColumns3M(:,:,iFrame), 1:1/parS.interpFactor:ROIHeight, parS.interpMethod);
    wPosColumnsInterpIniM = (wPosColumnsIniM-1)*parS.interpFactor+1;
    
    if parS.plotEveryFrame %Plot frame image
      figure(figH)
      hold off
      imshow(interpolatedFrameM, [0 255])
      aspectV = [pixelSizeV(1);pixelSizeV(2)*actualColumnWidth*(1-parS.columnOverlap)*parS.interpFactor;1];
      axH = gca;
      set(axH,'DataAspectRatio', aspectV);
      set(axH,'Position',[0,0,1,1])
      set(gcf, 'Position', figPosV)  

      hold on
    end
    
    for iWall = 1:2
      if iFrame == 1
        searchWindowMinIV = round(wPosColumnsInterpIniM(iWall,:) - depthWindow);
        searchWindowMaxIV = round(wPosColumnsInterpIniM(iWall,:) + depthWindow);
      else
        searchWindowMinIV = round(innerWallPosColumnsM(iWall,:,iFrame-1) - depthWindow);
        searchWindowMaxIV = round(innerWallPosColumnsM(iWall,:,iFrame-1) + depthWindow);
      end
      
      %Check that search window does not extend outside ROI
      if (min(searchWindowMinIV) < 1) || (max(searchWindowMaxIV)>parS.interpFactor*(ROIHeight-1)+1)
        error('Search window extends outside ROI. Please choose a (vertically) larger ROI.')
      end
      
      if parS.plotEveryFrame %Plot detection window borders
        plot(1:nTrackingColumns,searchWindowMinIV, 'LineWidth', 1, 'Color', wallColoursM(iWall,:))
        plot(1:nTrackingColumns,searchWindowMaxIV, 'LineWidth', 1, 'Color', wallColoursM(iWall,:))
      end
      
      dataM = zeros(2*depthWindow+1,nTrackingColumns);
      for iCol = 1:nTrackingColumns
        dataM(:,iCol) = interpolatedFrameM(searchWindowMinIV(iCol):searchWindowMaxIV(iCol),iCol);
      end
      
      if iWall == 1
        bg = mean(mean(dataM(end-bgWindow+1:end,:)));
      else
        bg = mean(mean(dataM(1:bgWindow,:)));
      end
      maxV = max(dataM,[],1);
      
      A=(dataM-bg)>parS.wallThreshold*repmat(maxV-bg,2*depthWindow+1,1);
      for iCol = 1:nTrackingColumns
        if iWall == 1
          innerWallPosColumnsM(iWall,iCol,iFrame) = find(A(:,iCol),1,'last') + searchWindowMinIV(iCol);
          outerWallPosColumnsM(iWall,iCol,iFrame) = find(A(:,iCol),1,'first') + searchWindowMinIV(iCol);
        else
          innerWallPosColumnsM(iWall,iCol,iFrame) = find(A(:,iCol),1,'first') + searchWindowMinIV(iCol);
          outerWallPosColumnsM(iWall,iCol,iFrame) = find(A(:,iCol),1,'last') + searchWindowMinIV(iCol);
        end
      end
      if parS.plotEveryFrame %Plot detected edge
        plot(1:nTrackingColumns,innerWallPosColumnsM(iWall,:,iFrame), 'LineWidth', 3, 'Color', wallColoursM(iWall,:))
        if parS.alsoTrackOuterWall
          plot(1:nTrackingColumns,outerWallPosColumnsM(iWall,:,iFrame), '--', 'LineWidth', 3, 'Color', wallColoursM(iWall,:))
        end
      end
    end %for iWall
    drawnow
  end %For iFrame
  disp('...done!')
  
  %% Plot diameter waveforms for all lines
  figure
  set(gcf, 'Name', 'Diameter tracking results')
  if parS.alsoTrackOuterWall
    s1H = subplot(3,1,1);
  end
  
  iFrameV = 1:nFrames;
  innerDiameterM = permute(diff(innerWallPosColumnsM,1,1),[2,3,1]);
  innerDiameterM = innerDiameterM*pixelSizeV(1)/parS.interpFactor;
  plot(iFrameV,innerDiameterM)
  hold on
  plot(iFrameV,mean(innerDiameterM), 'LineWidth', 3)
  xlabel('Frame number')
  ylabel('Inner diameter [\mum]')
  
  if parS.alsoTrackOuterWall
    s2H = subplot(3,1,2);
    outerDiameterM = permute(diff(outerWallPosColumnsM,1,1),[2,3,1]);
    outerDiameterM = outerDiameterM*pixelSizeV(1)/parS.interpFactor;
    plot(iFrameV,outerDiameterM)
    hold on
    plot(iFrameV,mean(outerDiameterM), 'LineWidth', 3)
    xlabel('Frame number')
    ylabel('Outer diameter [\mum]')
    
    s3H = subplot(3,1,3);
    thicknessNearWallM = innerWallPosColumnsM(1,:,:) - outerWallPosColumnsM(1,:,:);
    thicknessNearWallM = permute(thicknessNearWallM,[2,3,1])*pixelSizeV(1)/parS.interpFactor;
    thicknessFarWallM = outerWallPosColumnsM(2,:,:) - innerWallPosColumnsM(2,:,:);
    thicknessFarWallM = permute(thicknessFarWallM,[2,3,1])*pixelSizeV(1)/parS.interpFactor;
    
    plot(iFrameV,mean(thicknessNearWallM), 'LineWidth', 3, 'Color', wallColoursM(1,:))
    hold on
    plot(iFrameV,mean(thicknessFarWallM), 'LineWidth', 3, 'Color', wallColoursM(2,:))
    
    xlabel('Frame number')
    ylabel('Wall thickness [\mum]')
    legend('Near wall', 'Far wall')
    
    linkaxes([s1H, s2H, s3H], 'x')
  end
  
  %% Make interactive window to inspect fitting
  fig2H = figure;
  axH = axes('Parent', fig2H, 'Position', [0.01 0.07 0.98 0.92]);
  
  b = uicontrol('Parent',fig2H,'Style','slider', 'Units', 'normalized', ...  
  'Position',[0.01 0.01 0.98 0.05], 'value',1, 'min',1, 'max', nFrames, ...
  'SliderStep', [0.001 0.01]);
  b.Callback = @(hObj,eventData) updateInteractivePlotWindow(...
    hObj, wPosColumnsInterpIniM, trackingColumns3M, innerWallPosColumnsM, ...
    outerWallPosColumnsM, axH, actualColumnWidth, pixelSizeV, nTrackingColumns, wallColoursM, depthWindow, parS);
  b.Callback(b,[]); %Run callback so plot shows.
  
  %% Save output
  i = strfind(lower(VEVOBModeFileName), '.raw.bmode');
  i = i(1);
  VEVObase = VEVOBModeFileName(1:i-1);
  VEVOWallTrackFileName = [VEVObase '.raw.DI.mat'];
  
  disp('Saving file...')
  saveOutput(VEVOWallTrackFileName, pixelSizeV, innerDiameterM, parS)
  disp('...done!')
end

function updateInteractivePlotWindow(hObj, wPosColumnsInterpIniM, trackingColumns3M, innerWallPosColumnsM, outerWallPosColumnsM, axH, columnWidth, pixelSizeV, nTrackingColumns, wallColoursM, depthWindow, parS)
  iFrame = round(hObj.Value); %Get frame number
  
  ROIHeight = size(trackingColumns3M,1);
  interpolatedFrameM = interp1(1:ROIHeight, trackingColumns3M(:,:,iFrame), 1:1/parS.interpFactor:ROIHeight, parS.interpMethod); %Generate interpolated data
  
  %Show image
  axes(axH)
  hold off
  imshow(interpolatedFrameM, [0 255])
  aspectV = [pixelSizeV(1);pixelSizeV(2)*columnWidth*(1-parS.columnOverlap)*parS.interpFactor;1];
  set(axH,'DataAspectRatio', aspectV);
  hold on
  
  for iWall = 1:2
    %Re-compute tracking window
    if iFrame == 1
      searchWindowMinIV = round(wPosColumnsInterpIniM(iWall,:) - depthWindow);
      searchWindowMaxIV = round(wPosColumnsInterpIniM(iWall,:) + depthWindow);
    else
      searchWindowMinIV = round(innerWallPosColumnsM(iWall,:,iFrame-1) - depthWindow);
      searchWindowMaxIV = round(innerWallPosColumnsM(iWall,:,iFrame-1) + depthWindow);
    end
    
    %Show detected wall
    plot(1:nTrackingColumns,innerWallPosColumnsM(iWall,:,iFrame), 'LineWidth', 3, 'Color', wallColoursM(iWall,:))
    if parS.alsoTrackOuterWall
      plot(1:nTrackingColumns,outerWallPosColumnsM(iWall,:,iFrame), '--', 'LineWidth', 3, 'Color', wallColoursM(iWall,:))
    end
    
    %Show detection window
    plot(1:nTrackingColumns,searchWindowMinIV, 'LineWidth', 1, 'Color', wallColoursM(iWall,:))
    plot(1:nTrackingColumns,searchWindowMaxIV, 'LineWidth', 1, 'Color', wallColoursM(iWall,:))
  end
  
  plot(1:nTrackingColumns,searchWindowMinIV, 'LineWidth', 1, 'Color', wallColoursM(iWall,:))
  plot(1:nTrackingColumns,searchWindowMaxIV, 'LineWidth', 1, 'Color', wallColoursM(iWall,:))
  
  figH = ancestor(axH, 'figure');
  set(figH, 'Name', sprintf('Frame number: %.0f', iFrame))
end

function saveOutput(VEVOWallTrackFileName, pixelSizeV, innerDiameterM, parS)
  %Save data in the same format as Arnolds data.
  vArt = struct();
  vArt.pixscal = pixelSizeV(1)/1000; %in mm!!
  vArt.frsel = [1 size(innerDiameterM,2)];
  ddistr = innerDiameterM; %nLines x nFrames diameter matrix (in pixels)
  save(VEVOWallTrackFileName, 'vArt', 'ddistr', 'parS', '-mat')
end

function getFourClickedPoints(hObj,~)
  figH = ancestor(gca, 'figure');
  pCnt = getappdata(figH, 'ClickedPointCount');
  pCnt = pCnt + 1;
  setappdata(figH, 'ClickedPointCount', pCnt) %Increment clicked point counter
  
  clickedPointLocationsV = getappdata(figH, 'ClickedPointLocationsV');
  cp = get(hObj,'CurrentPoint');
  xp=cp(1,1);
  yp=cp(1,2);
  clickedPointLocationsV(pCnt,:) = [xp yp];
  setappdata(figH, 'ClickedPointLocationsV', clickedPointLocationsV)
  
  plot(xp,yp,'o','MarkerSize',5,'color','r','MarkerFaceColor','w');
  
  if pCnt >=4
    uiresume();
  end
end