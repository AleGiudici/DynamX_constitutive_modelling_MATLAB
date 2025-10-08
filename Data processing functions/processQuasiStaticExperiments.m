function [passive, directoryName] = processQuasiStaticExperiments(app)

    % initialisation
    passive.analysis_date = datetime("today"); % save the date of the analysis
    Mtm = 1000; % mm to micrometer
    
    %% Get unloaded configuration from GUI
    R0 = str2double(app.unloaded_Ro.Value); % unloaded outer radius [um]
    Ri = str2double(app.unloaded_Ri.Value); % unloaded inner radius [um]
    
    passive.OD = R0*2; % unloaded intact outer diameter [um]
    passive.L= str2double(app.unloaded_lz.Value); % unloaded intact length [mm]
    passive.H = R0-Ri; % unloaded intact wall thickness [um]
    
    app.axial_stretch_menu.Value = 'Raw data';
    
%% Protocol data
    [PD_stretch_passive,FL_pressure_passive,~,~,~,~,~,~,~,~,~,~,numberTests,correctingAxialStretch] = retriveProtocol(app);
    
%% Processing of automatically generated static P-A sweeps, with sync pulse alignment
    % if you worked on Dynamice 1, BEFORE RUNNING THIS, open vidDist and process the VEVO .raw.bmode file.
    % This yields a .rawdi.mat file in the same folder as the .raw.bmode file
    
%% Importing force sweep data
    dataType = 'FL_static';
    forceLengthCheck = 1;
    dynamicCheck = 0;
    
    k=1; % reset counter
    answer = 'Yes'; % to start the while cycle
    colour_scheme = {'r','b','g','m','k'}; % colour scheme for plots
    
    FL = strings(1,length(FL_pressure_passive));
    FL_legend = strings(1,length(FL_pressure_passive));
    FL_pressure = strings(1,length(FL_pressure_passive));
    
    while(strcmp(answer,'Yes'))
        % LOCATING FILE DIRECTORY
        if(~exist('directoryName','var')) % on the first iteration allows selecting the directory
            directoryName = uigetdir('Select the folder with the experimental data');
        end
    
        if(prod(directoryName == 0))
            break
        end
    
        FL_pressure{k} = FL_pressure_passive{k};
    
        FL{k} = ['FL_' FL_pressure{k}]; % experiment name
        FL_legend{k} = [FL_pressure{k} ' mmHg']; % entry for plot legend
        
        % IMPORTING DATA
        if(any(size(dir([directoryName '/*.raw.DI.mat' ]),1))) 
            dataMat = importDataDynamX1(directoryName,FL{k},passive,1);
        else
            dataMat = importDataDynamX2(app,directoryName,dataType,FL{k},passive,forceLengthCheck,dynamicCheck);
        end
        % dataMat is [timeV,pressureV,innerDiameterV,outerDiameterV,transducerForceV,axialStretchV,axialLengthV];
        %               1       2           3               4              5              6              7
        
        % Correcting for sample elongation in second leg of incubation
        % experiments
        dataMat(:,6) = dataMat(:,6)/correctingAxialStretch;

        %PLOTTING
        if(k~=1) % do not hold current plots at the first iteration
            hold (app.axes3_dat_proc_tab,'on')
            hold (app.axes4_dat_proc_tab,'on')
        end
        
        % outer diameter-axial stretch plot
        plot(app.axes3_dat_proc_tab,dataMat(:,6),dataMat(:,3)/Mtm,colour_scheme{k},'DisplayName',FL_legend{k}); 
        ylabel(app.axes3_dat_proc_tab,'Inner diameter [mm]','FontSize',16)
        xlabel(app.axes3_dat_proc_tab,'Axial stretch [-]','FontSize',16)
        hold (app.axes3_dat_proc_tab,'off')
    
        
        % axial force-axial stretch plot
        plot(app.axes4_dat_proc_tab,dataMat(:,6),dataMat(:,5),colour_scheme{k},'DisplayName',FL_legend{k}); 
        ylabel(app.axes4_dat_proc_tab,'Axial force [g]','FontSize',16)
        xlabel(app.axes4_dat_proc_tab,'Axial stretch [-]','FontSize',16)
        if(k==1)
            legend(app.axes4_dat_proc_tab,FL_legend{1},'Location','NorthWest','FontSize',13);
        end
        hold (app.axes4_dat_proc_tab,'off')
        
        % SAVING
        % splitting "loading" and "unloading" phases and saving in the
        % "passive" structure
        full = size(dataMat,1);
        div = round(full/2);
    
        if(round(full/2)*2 ~= full)
            additional_step = 0;
        else
            additional_step = 1;
        end
    
        headings = ["Time [s]","Pressure [mmHg]","Outer diameter [um]","Transducer axial force [g]","Axial stretch [-]","Axial length [mm]","Inner diameter [um]"];
        % Loading data
        test_data = [dataMat(1:div,1:2) dataMat(1:div,4:7) dataMat(1:div,3)];
        passive.static_data.loading.biaxial_Fl.(FL{k}) = array2table(test_data,'VariableNames',headings); % loading
        
        % Unloading data
        test_data = [dataMat(div+additional_step:full,1:2) dataMat(div+additional_step:full,4:7) dataMat(div+additional_step:full,3)];
        passive.static_data.unloading.biaxial_Fl.(FL{k}) = array2table(test_data,'VariableNames',headings); % unloading
    
        k = k+1; % moving to next experiment if answer = Yes
        app.data_proc_dial.Value = (k-1)/numberTests*100;
        if(k==size(FL_pressure_passive,2)+1)
            answer='No';
        end
    
        pause(0.01);
    end
    
    
%% Importing pressure sweep data
    dataType = 'PD_static';
    forceLengthCheck = 0;
    dynamicCheck = 0;
    
    k=1; % counter
    answer = 'Yes'; % to start the while cycle
    
    PD = strings(1,length(PD_stretch_passive));
    PD_legend = strings(1,length(PD_stretch_passive));
    PD_stretch = strings(1,length(PD_stretch_passive));
    
    while(strcmp(answer,'Yes'))
        PD_stretch{k} = PD_stretch_passive{k};
        PD{k} = ['PD_' PD_stretch{k}]; % Name of the pressure sweep experiment  
        PD_legend{k} = ['\lambda_z = ' PD_stretch{k} '% in vivo']; % plot legend entry for the current test
        
        % Identify the Dynamice set-up used (i.e. ultrasound vs video tracking)
        % on the basis of the files in the selected folder.
        if(any(size(dir([directoryName '/*.raw.DI.mat' ]),1))) 
            dataMat = importDataDynamX1(directoryName,PD{k},passive,1);
        else
            dataMat = importDataDynamX2(app,directoryName,dataType,PD{k},passive,forceLengthCheck,dynamicCheck);
        end

        % Correcting for sample elongation in second leg of incubation
        % experiments
        dataMat(:,6) = dataMat(:,6)/correctingAxialStretch;
    
        % PLOTTING
        if(strcmp(app.static_dyna_menu_data_proc.Value,'Static'))
            if(k~=1) % do not hold current plots at the first iteration
                hold (app.axes1_dat_proc_tab,'on')
                hold (app.axes2_dat_proc_tab,'on')
            end
            % dataMat is [timeV,pressureV,innerDiameterV,outerDiameterV,transducerForceV,axialStretchV,axialLengthV];
            %               1       2           3               4              5              6              7
            
            % pressure-diameter plot
            plot(app.axes1_dat_proc_tab,dataMat(:,3)/Mtm,dataMat(:,2),colour_scheme{k},'DisplayName',PD_legend{k})
            xlabel(app.axes1_dat_proc_tab,'Inner diameter [mm]','FontSize',16)
            ylabel(app.axes1_dat_proc_tab,'Pressure [mmHg]','FontSize',16)
            if(k==1)
                legend(app.axes1_dat_proc_tab,PD_legend{1},'Location','NorthWest','FontSize',13);
            end
            hold (app.axes1_dat_proc_tab,'off')
        
            % force-pressure plot
            plot(app.axes2_dat_proc_tab,dataMat(:,2),dataMat(:,5),colour_scheme{k},'DisplayName',PD_legend{k}); 
            xlabel(app.axes2_dat_proc_tab,'Pressure [mmHg]','FontSize',16)
            ylabel(app.axes2_dat_proc_tab,'Axial force [g]','FontSize',16)
            hold (app.axes2_dat_proc_tab,'off')
    
        elseif(strcmp(app.static_dyna_menu_data_proc.Value(1:7),'Dynamic'))
            if(str2double(PD_stretch{k})==100)
                % pressure-diameter plot
                plot(app.axes1_dat_proc_tab,dataMat(:,3)/Mtm, dataMat(:,2),'k--','DisplayName',PD_legend{k})
                xlabel(app.axes1_dat_proc_tab,'Inner diameter [mm]','FontSize',16)
                ylabel(app.axes1_dat_proc_tab,'Pressure [mmHg]','FontSize',16)
                legend(app.axes1_dat_proc_tab,PD_legend{k},'Location','NorthWest','FontSize',13);
    
                % force-pressure plot
                plot(app.axes2_dat_proc_tab,dataMat(:,2),dataMat(:,5),'k--','DisplayName',PD_legend{k}); 
                xlabel(app.axes2_dat_proc_tab,'Pressure [mmHg]','FontSize',16)
                ylabel(app.axes2_dat_proc_tab,'Axial force [g]','FontSize',16)
                hold (app.axes2_dat_proc_tab,'off')
            end
        end
        
        % SAVING
    
        % split relationship into "loading" and "unloading" phases and save in
        % the "passive" structure
        full = size(dataMat,1);
        div = round(full/2);
    
        if(round(full/2)*2 ~= full)
            additional_step = 0;
        else
            additional_step = 1;
        end
    
        headings = {'Time [s]','Pressure [mmHg]','Outer diameter [um]','Transducer axial force [g]','Axial stretch [-]','Axial length [mm]','Inner diameter [um]'};
        % Loading data
        test_data = [dataMat(1:div,1:2) dataMat(1:div,4:7) dataMat(1:div,3)];
        passive.static_data.loading.biaxial_Pd.(PD{k}) = array2table(test_data,'VariableNames',headings); % loading
        
        % Unloading data
        test_data = [dataMat(div+additional_step:full,1:2) dataMat(div+additional_step:full,4:7) dataMat(div+additional_step:full,3)];
        passive.static_data.unloading.biaxial_Pd.(PD{k}) = array2table(test_data,'VariableNames',headings); % unloading
       
        k=k+1; % move to next experiments if answer = 'Yes', exit otherwise.
        app.data_proc_dial.Value = (k-1+size(FL_pressure_passive,2))/numberTests*100;
        if(k==size(PD_stretch_passive,2)+1)
            answer='No';
        end
    
        pause(0.01);
    end
    
    if(prod(directoryName~=0))
        legend(app.axes4_dat_proc_tab,FL_legend{1:k-1},'Location','NorthWest','FontSize',13);
    
        passive.static_data.PD = PD; 
        passive.static_data.FL = FL;
    else
        passive = [];
    end
    pause(2);

%% Estimating in vivo axial stretch from intersection point between axial
%  force sweeps at 60, 100 and 140 mmHg.
    [passive,app] = findInVivoAxialStretch(app,passive,FL,colour_scheme);
    
%% Calculating axial structural stiffness at the in vivo axial stretch
    
    x_lb = app.axes4_dat_proc_tab.XLim(1);
    x_ub = app.axes4_dat_proc_tab.XLim(2);
    y_lb = app.axes4_dat_proc_tab.YLim(1);
    y_ub = app.axes4_dat_proc_tab.YLim(2);
    
    % [passive] = calculateAxialStructStiffness(FL,passive);

    axis(app.axes4_dat_proc_tab, [x_lb x_ub y_lb y_ub]);
end