function passive = processDynamicExperiments(passive, app, directoryName)
% This function allows processing and formatting the experimental data
% acquired with Dynamice 1 and 2 during dynamic pressure wavefroms
% experiments. The function receives as input the structure 'passive' that
% contains already the processed data of the quasi-static pressure sweeps
% and force sweeps experiments. Then, this structure is further populated
% with the processed data of the dynamic experiments and given as output.

%% Protocol data
    [PD_stretch_passive,FL_pressure_passive,PD_DBP_passive,PD_SBP_passive,PD_Hz_passive,PD_label_passive,...
        ~,~,~,~,~,~,numberTests,correctingAxialStretch] = retriveProtocol(app);

%% Importing reference quasi-static pressure sweep at in vivo axial length
    static_ref_data = (table2array(passive.static_data.loading.biaxial_Pd.PD_100)...
        +table2array(flip(passive.static_data.unloading.biaxial_Pd.PD_100)))/2;

    static_P_ref = static_ref_data(:,2);
    static_ri_ref = static_ref_data(:,7)/2;
    
    if(static_P_ref(1)==static_P_ref(2))
        static_P_ref = static_P_ref(2:end);
        static_ri_ref = static_ri_ref(2:end);
    end
    
    static_P_ref_res = 15:0.1:175;
    static_ri_ref_res = interp1(static_P_ref,static_ri_ref,static_P_ref_res);

%% Importing dynamic data
% if you worked on Dynamice 1, BEFORE RUNNING THIS, open vidDist and process the VEVO .raw.bmode file.
% This yields a .rawdi.mat file in thesame folder as the .raw.bmode file

    colour_scheme = {'r','b','g','m','k','y'}; % for plots
    iTest = 1;
    answer = 'Yes'; % default answer to start the while loop
    dataType = 'PD_dynamic';
    dynamicCheck = 1;
    forceLengthCheck = 0;
    PD = strings(1,length(PD_DBP_passive));
    legend_entry = strings(1,length(PD_DBP_passive));
    haemodynamic_var = zeros(8,length(PD_DBP_passive));
    
    while(strcmp(answer,'Yes'))
        PD{iTest} = ['PD_' PD_DBP_passive{iTest} '_' PD_SBP_passive{iTest} '_' PD_Hz_passive{iTest} '_' PD_label_passive{iTest}]; % experiment name
    
        if(any(size(dir([directoryName '/*.raw.DI.mat' ]),1))) 
            dataMat = importDataDynamX1(directoryName,PD,passive,forceLengthCheck,dynamicCheck);
        else
            dataMat = importDataDynamX2(app,directoryName,dataType,PD{iTest},passive,forceLengthCheck,dynamicCheck);
        end
        % dataMat is [timeV,pressureV,innerDiameterV,outerDiameterV,transducerForceV,axialStretchV,axialLengthV];
        %               1       2           3               4              5              6              7
        
        % Correcting for sample elongation in second leg of incubation
        % experiments
        dataMat(:,6) = dataMat(:,6)/correctingAxialStretch;
        
        if(strcmp(app.protocol_menu.Value,'New') || strcmp(app.protocol_menu.Value,'Custom'))
            legend_entry{iTest} = [PD_DBP_passive{iTest} '-' PD_SBP_passive{iTest} ' mmHg, ' num2str(round(str2double(PD_Hz_passive{iTest})/100)/10) ' Hz, ' PD_label_passive{iTest}]; % legend entry string
        else
            legend_entry{iTest} = [PD_DBP_passive{iTest} '-' PD_SBP_passive{iTest} ' mmHg, ' num2str(round(str2double(PD_Hz_passive{iTest})/100)/10) ' Hz']; % legend entry string
        end
        
        j = 1;
    
        % PLOTTING
        if(strcmp(app.protocol_menu.Value,'New') || strcmp(app.protocol_menu.Value,'Custom'))
            if(strcmp(app.static_dyna_menu_data_proc.Value,'Dynamic (physiological)') && strcmp(PD_label_passive{iTest},'physio'))
                % pessure-diameter plot
                hold (app.axes1_dat_proc_tab, 'on');
                plot(app.axes1_dat_proc_tab, dataMat(:,3)/10^3, dataMat(:,2), colour_scheme{iTest}, 'DisplayName',legend_entry{iTest})
                hold (app.axes1_dat_proc_tab, 'off');
                
                % force-pressure plot
                hold (app.axes2_dat_proc_tab, 'on');
                plot(app.axes2_dat_proc_tab,dataMat(:,2), dataMat(:,5), colour_scheme{iTest}); 
                hold (app.axes2_dat_proc_tab, 'off');
            elseif(strcmp(app.static_dyna_menu_data_proc.Value,'Dynamic (sine)') && strcmp(PD_label_passive{iTest},'sine'))
                % pessure-diameter plot
                hold (app.axes1_dat_proc_tab, 'on');
                plot(app.axes1_dat_proc_tab, dataMat(:,3)/10^3, dataMat(:,2), colour_scheme{j}, 'DisplayName',legend_entry{iTest})
                hold (app.axes1_dat_proc_tab, 'off');
                
                % force-pressure plot
                hold (app.axes2_dat_proc_tab, 'on');
                plot(app.axes2_dat_proc_tab, dataMat(:,2), dataMat(:,5), colour_scheme{j}); 
                hold (app.axes2_dat_proc_tab, 'off');
                
                if(iTest<length(PD_SBP_passive))
                    if(strcmp(PD_SBP_passive{iTest},PD_SBP_passive{iTest+1}))
                        j = j+1;
                    else
                        j = 1;
                    end
                end
            end
        else
            app.static_dyna_menu_data_proc.Value = 'Dynamic (physiological)';
            % pessure-diameter plot
            hold (app.axes1_dat_proc_tab, 'on');
            plot(app.axes1_dat_proc_tab, dataMat(:,3)/10^3, dataMat(:,2), colour_scheme{iTest}, 'DisplayName',legend_entry{iTest})
            hold (app.axes1_dat_proc_tab, 'off');
            
            % force-pressure plot
            hold (app.axes2_dat_proc_tab, 'on');
            plot(app.axes2_dat_proc_tab, dataMat(:,2), dataMat(:,5), colour_scheme{iTest}); 
            hold (app.axes2_dat_proc_tab, 'off');
        end
    
        % CALCULATING HAEMODYNAMIC METRICS
        haemodynamic_var(1:8,iTest) = calculateCircStructStiffness(dataMat,static_P_ref_res,static_ri_ref_res);
            
        % SAVING
        headings = ["Time [s]","Pressure [mmHg]","Outer diameter [um]","Transducer axial force [g]","Axial stretch [-]","Axial length [mm]","Inner diameter [um]"];
        test_data = [dataMat(:,1:2) dataMat(:,4:7) dataMat(:,3)];
        passive.dynamic_data.dynamic_PD.(PD{iTest}) = array2table(test_data,'VariableNames',headings);
        
        iTest = iTest+1; % move to next experiments if answer = 'Yes', exit otherwise.
        app.data_proc_dial.Value = (iTest-1+size(FL_pressure_passive,2)+size(PD_stretch_passive,2))/numberTests*100;
        if(iTest==size(PD_DBP_passive,2)+1)
            answer='No';
        end

        pause(0.01);
    end
    
    row_headings = ["SBP [mmHg]"; "DBP [mmHg]"; "PWV [m/s]"; "PWV (quasi-static) [m/s]"; "Distensibility [1/MPa]";...
        "Distensibility (quasi-static) [1/MPa]"; "Compliance [mm^4/N]"; "Compliance (quasi-static) [mm^4/Pa]"];
    
    passive.preliminary_analysis.haemodynamic_variables = array2table(haemodynamic_var,'RowNames',row_headings,'VariableNames',legend_entry);
    passive.dynamic_data.PD = PD;
end