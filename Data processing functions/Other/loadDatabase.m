function databaseT = loadDatabase(folderName, databaseFileName)
  fileName = fullfile(folderName, databaseFileName);
  
  varNamesAndTypesC = { ... 
    'MouseID'                   'char'                 
    'SampleID'                  'char'               
    'ExperimentID'              'char'            
    'Date'                      'datetime'      
    'UnloadedLength_mm_'        'double'      
    'InVivoLength_mm_'          'double'      
    'BaseFolderName'            'char'
    'DynamicLVFileNameBase'     'char'
    'ForceSweepLVFileNameBase'  'char'
    'PressureSweepLVNameBase'   'char'
    'DynamicVEVOFileNameBase'   'char'
    'ForceVEVOFileNameBase'     'char'
    'PressureSweepVEVONameBase' 'char'
    };
  
  opts = detectImportOptions(fileName);
  %opts = spreadsheetImportOptions;
  opts.VariableNames = varNamesAndTypesC(:,1);
  opts.DataRange = 'A2'; %Discard first row
  opts = setvartype(opts, varNamesAndTypesC(:,2));
  databaseT = readtable(fileName, opts);
end

