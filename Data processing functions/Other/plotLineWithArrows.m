function plotLineWithArrows(xV, yV, arrowInterval, arrowLengthScaling, arrowHeadStyle, arrowHeadSize)
  
  intervalMode = 'time'; %distance: equidistant on plot, or 'time': at equal time intervals
  
  %Delete NaN elements
  iOK = ~isnan(xV) & ~isnan(yV);
  xV = xV(iOK);
  yV = yV(iOK);
  
  pH = plot(xV, yV);
  hold on
  
  tV = (1:length(xV))';
  xInt = griddedInterpolant(tV, xV);
  yInt = griddedInterpolant(tV, yV);
  
  dt = 0.00001;
  xDiff = @(t) (xInt(t+dt) - xInt(t))/dt;
  yDiff = @(t) (yInt(t+dt) - yInt(t))/dt;
  
  %arrowScaling = 0.001;
  
  %arrowInterval = 20;
  tArrowV = 1:arrowInterval:tV(end);
  nArrows = length(tArrowV);
  
  dxdtArrV = zeros(nArrows,1);
  dydtArrV = zeros(nArrows,1);
  xArrV = zeros(nArrows,1);
  yArrV = zeros(nArrows,1);
  
  for iArrow = 1:nArrows
    tArrow = tArrowV(iArrow);
    dxdtArrV(iArrow) = xDiff(tArrow);
    dydtArrV(iArrow) = yDiff(tArrow);
    xArrV(iArrow) = xInt(tArrow);
    yArrV(iArrow) = yInt(tArrow);
  end
  
  %Plot two invisible dots to ensure that axes limits are large enough to
  %plot all arrows
  plot([min([xArrV; xArrV+dxdtArrV*arrowLengthScaling]) max([xArrV; xArrV+dxdtArrV*arrowLengthScaling])], ...
    [min([yArrV; yArrV+dydtArrV*arrowLengthScaling]) max([yArrV; yArrV+dydtArrV*arrowLengthScaling])], '.w') %Plot invisible dots
  
  %Actually plot the arrows:
  for iArrow = 1:nArrows
    xInDataUnits = [xArrV(iArrow) xArrV(iArrow)+dxdtArrV(iArrow)*arrowLengthScaling];
    yInDataUnits = [yArrV(iArrow) yArrV(iArrow)+dydtArrV(iArrow)*arrowLengthScaling];
    arrowH = drawArrow(xInDataUnits, yInDataUnits, gca);
    arrowH.Color = pH.Color;
    arrowH.HeadStyle = arrowHeadStyle;
    arrowH.HeadLength = arrowHeadSize;
    arrowH.HeadWidth = arrowHeadSize;
  end
end