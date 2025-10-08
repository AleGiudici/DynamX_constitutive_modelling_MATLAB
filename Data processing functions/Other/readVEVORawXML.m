function sysParsS = readVEVORawXML(fNameRawXML)
  % READVEVORAWXML Function to read VEVO *.raw.xml files that contain imaging settings
  %   sysParsS = READVEVORAWXML(fNameRawXML) Reads a .raw.xml file and
  %     parses its contents into the structure sysParsS.
  %     sysParsS is an nx1 structure with n the number of parameters
  %     specified in the XML file (typically 53).
  %     For each parameter, sysParsS contains three fields which are all
  %     character arrays:
  %       sysParsS(i).name:  parameter name (e.g. 'Time-Stamp-Clock')
  %       sysParsS(i).value: parameter value (e.g. '400000')
  %       sysParsS(i).units: parameter units (e.g. 'Hz')
  %     The units field may be empty for dimensionless parameters.
  %
  % V1.0 by Bart Spronck.
  
  xmlS = parseXML(fNameRawXML);
  
  sysParsS = struct(); %Empty struct to save parameters
  
  okIV = false(length(xmlS.Children),1);
  for i = 1:length(xmlS.Children)
    if ~isempty(xmlS.Children(i).Attributes)
      okIV(i) = true;
      for j = 1:length(xmlS.Children(i).Attributes)
        sysParsS(i,1).(xmlS.Children(i).Attributes(j).Name) = xmlS.Children(i).Attributes(j).Value;
      end
    end
  end
  sysParsS = sysParsS(okIV); %Remove empty fields.
end

