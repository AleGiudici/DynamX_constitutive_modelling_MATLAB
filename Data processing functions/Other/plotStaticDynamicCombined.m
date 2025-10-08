function plotStaticDynamicCombined(staticPDDataT, dynamicDataT, mouseID, sampleID, experimentID, protocolPressuresV, frequenciesV, plotArrows, arrowInterval, arrowLengthScaling, arrowHeadStyle, arrowHeadSize)
  %Get data for this experiment only
  staticPDDataSingleExperimentT = staticPDDataT( ...
    strcmp(staticPDDataT.MouseID, mouseID) & ...
    strcmp(staticPDDataT.SampleID, sampleID) & ...
    strcmp(staticPDDataT.ExperimentID, experimentID) ...
    ,:);
  
  dynamicDataSingleExperimentT = dynamicDataT( ...
    strcmp(dynamicDataT.MouseID, mouseID) & ...
    strcmp(dynamicDataT.SampleID, sampleID) & ...
    strcmp(dynamicDataT.ExperimentID, experimentID) ...
    ,:);
  
  staticDataPointsToUseForDynamicComparison = 81:-1:41; %Use deflation points of static curve
  
  iV = staticPDDataSingleExperimentT.ProtocolAxialStretch == 1; %Which row corresponds to a stretch of 1.0x in vivo?
  %Typically iV = [0;1;0]
  
  pStaticV = staticPDDataSingleExperimentT(iV,:).PressureV{1}(staticDataPointsToUseForDynamicComparison);
  dStaticV = staticPDDataSingleExperimentT(iV,:).DiameterV{1}(staticDataPointsToUseForDynamicComparison);
  
  plot(dStaticV*1000, pStaticV, '-k', 'LineWidth', 1)
  hold on
  
  %Plot three pressures for f=5
  
  for iPlot = 1:length(protocolPressuresV)
    protocolPressure = protocolPressuresV(iPlot);
    frequency = frequenciesV(iPlot);
    iV = (dynamicDataSingleExperimentT.ProtocolPressure == protocolPressure) & (dynamicDataSingleExperimentT.ProtocolFrequency == frequency);
    pDynamicV = dynamicDataSingleExperimentT(iV,:).PressureV{1};
    dDynamicV = dynamicDataSingleExperimentT(iV,:).DiameterV{1};
    if plotArrows
      %arrowInterval = 30;
      %arrowLengthScaling = 0.00001;
      plotLineWithArrows(dDynamicV*1000, pDynamicV, arrowInterval, arrowLengthScaling, arrowHeadStyle, arrowHeadSize)
    else
      plot(dDynamicV*1000, pDynamicV)
    end
  end
  xlabel('Diameter [\mum]')
  ylabel('Pressure [mmHg]')
  
  legendC = arrayfun(@(p, f) sprintf('Dynamic, {\\itP}_{protocol}=%.0f mmHg, {\\itf}=%.1f Hz', p, f), protocolPressuresV, frequenciesV, 'UniformOutput', false);
  legendC = ['Quasi-static' legendC];
  legend(legendC, 'Location', 'NorthWest')
end

