function app = plotGraphsModellingTab(app,fit_data)

    cla(app.axes1_mod_tab);
    cla(app.axes2_mod_tab);

    if(strcmp(app.graph_choice_mod_tab_2.Value,'Dynamic loops (passive)'))
        set(app.graph_choice_mod_tab,'Enable','off');
    else
        set(app.graph_choice_mod_tab,'Enable','on');
    end

    if((~isfield(fit_data.fittedCurves,'dynamicFit')) && (strcmp(app.graph_choice_mod_tab_2.Value,'Dynamic loops (passive)')))
        app.graph_choice_mod_tab_2.Value='Pressure sweeps';

        app = appErrorDialogBox(app,'off','off','on','ERROR: missing dynamic fit data',...
                    'No dynamic fitting data available in the selected file.');
    end

    hold (app.axes1_mod_tab, 'on');
    hold (app.axes2_mod_tab, 'on');

    if(isfield(fit_data.fittedCurves,'quasiStaticFit'))
        fieldName = 'quasiStaticFit';
    else
        fieldName = 'staticFit';
    end
    
    if(strcmp(app.graph_choice_mod_tab_2.Value,'Pressure sweeps (passive)'))
        PD = fit_data.PD_static;
    
        lb_x = 10^7;
        ub_x = 0;
        ub_y = 0;
    
        color_sequence = {'r','b','g'};
        color_sybol_sequence = {'ro','bo','go'};
        legend_entry = strings(1,length(PD)*2);

        for i = 1:length(PD)            
            dM = fit_data.fittedCurves.(fieldName).(PD{i});            
            if(strcmp(app.graph_choice_mod_tab.Value,'Pressure and force'))            
                % hold (app.axes1_mod_tab, 'on')
                plot(app.axes1_mod_tab,dM.("Outer diameter [mm]"), dM.("Exp. pressure [mmHg]"),color_sybol_sequence{i}) % Plot experimental P-D relationship
                plot(app.axes1_mod_tab,dM.("Outer diameter [mm]"), dM.("Modelled pressure [mmHg]"),color_sequence{i}) % Plot fitted P-D relationship
                % hold (app.axes1_mod_tab, 'off')

                xlabel(app.axes1_mod_tab, 'Outer diameter [mm]','FontSize',16)
                ylabel(app.axes1_mod_tab, 'Pressure [mmHg]','FontSize',16)
    
                lb_x_prev = lb_x;
                lb_x = 0.95*min(dM.("Outer diameter [mm]"));
                if(lb_x >= lb_x_prev)
                    lb_x = lb_x_prev;
                end
    
                ub_x_prev = ub_x;
                ub_x = 1.05*max(dM.("Outer diameter [mm]"));
                if(ub_x <= ub_x_prev)
                    ub_x = ub_x_prev;
                end            
                axis(app.axes1_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 200]);
    
                % hold (app.axes2_mod_tab, 'on');
                plot(app.axes2_mod_tab,dM.("Exp. pressure [mmHg]"), dM.("Exp. transducer force [g]"),color_sybol_sequence{i}) % Plot experimental fT-P relationship
                plot(app.axes2_mod_tab,dM.("Modelled pressure [mmHg]"), dM.("Modelled transducer force [g]"),color_sequence{i}) % Plot fitted fT-P relationship
                % hold (app.axes2_mod_tab, 'off');

                xlabel(app.axes2_mod_tab, 'Pressure [mmHg]','FontSize',16)
                ylabel(app.axes2_mod_tab, 'Axial force [g]','FontSize',16)
               
                axis(app.axes2_mod_tab,[0 200 0 1.05*max(max(dM.("Exp. transducer force [g]")),max(dM.("Modelled transducer force [g]")))]);
    
            elseif(strcmp(app.graph_choice_mod_tab.Value,'Stresses'))
                % hold (app.axes1_mod_tab, 'on')
                plot(app.axes1_mod_tab,dM.("Circ. stretch [-]"), dM.("Exp. circ. stress [MPa]"),color_sybol_sequence{i}) % Plot experimental P-D relationship
                plot(app.axes1_mod_tab,dM.("Circ. stretch [-]"), dM.("Modelled circ. stress [MPa]"),color_sequence{i}) % Plot fitted P-D relationship
                % hold (app.axes1_mod_tab, 'off')

                xlabel(app.axes1_mod_tab, 'Circumferential stretch [-]','FontSize',16)
                ylabel(app.axes1_mod_tab, 'Circumferential stress [MPa]','FontSize',16)
    
                lb_x_prev = lb_x;
                lb_x = 0.95*min(dM.("Circ. stretch [-]"));
                if(lb_x >= lb_x_prev)
                    lb_x = lb_x_prev;
                end
    
                ub_x_prev = ub_x;
                ub_x = 1.05*max(dM.("Circ. stretch [-]"));
                if(ub_x <= ub_x_prev)
                    ub_x = ub_x_prev;
                end
    
                ub_y_prev = ub_y;
                ub_y = max(max(dM.("Exp. circ. stress [MPa]")),max(dM.("Exp. axial stress [MPa]")));
                if(ub_y <= ub_y_prev)
                    ub_y = ub_y_prev;
                end
    
                axis(app.axes1_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 round(ub_y*10)/10]);
    
                % hold (app.axes2_mod_tab, 'on');
                plot(app.axes2_mod_tab,dM.("Circ. stretch [-]"), dM.("Exp. axial stress [MPa]"),color_sybol_sequence{i}) % Plot experimental fT-P relationship
                plot(app.axes2_mod_tab,dM.("Circ. stretch [-]"), dM.("Modelled axial stress [MPa]"),color_sequence{i}) % Plot fitted fT-P relationship
                % hold (app.axes2_mod_tab, 'off');

                xlabel(app.axes2_mod_tab, 'Circumferential stretch [-]','FontSize',16)
                ylabel(app.axes2_mod_tab, 'Axial stress [MPa]','FontSize',16)
    
                lb_x_prev = lb_x;
                lb_x = 0.95*min(dM.("Circ. stretch [-]"));
                if(lb_x >= lb_x_prev)
                    lb_x = lb_x_prev;
                end
    
                ub_x_prev = ub_x;
                ub_x = 1.05*max(dM.("Circ. stretch [-]"));
                if(ub_x <= ub_x_prev)
                    ub_x = ub_x_prev;
                end
    
                axis(app.axes2_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 round(ub_y*10)/10]);            
            end
            legend_entry{2*i-1} = ['\lambda_z = ' PD{i}(4:end) '% in vivo, exp'];
            legend_entry{2*i} = ['\lambda_z = ' PD{i}(4:end) '% in vivo, mod'];
        end
        legend(app.axes1_mod_tab,legend_entry{1:i*2},'Location','NorthWest','FontSize',13);
        legend(app.axes2_mod_tab,'off');
               
    elseif(strcmp(app.graph_choice_mod_tab_2.Value,'Dynamic loops (passive)') && isfield(fit_data,'PD_dynamic'))
        colour_scheme = {'r','b','g','m','y','c'};
        stiffening_ratio_exp = table2array(fit_data.D_to_QS_stiffeningRatio.passive.Experimental);
        stiffening_ratio_mod = table2array(fit_data.D_to_QS_stiffeningRatio.passive.Modelled);
        
        freq = zeros(1,size(fit_data.PD_dynamic,2)/3);
        for iTest = 1:size(fit_data.PD_dynamic,2)/3
            loc_underscore = find(fit_data.PD_dynamic{iTest} == '_');
            freq(iTest) = str2double(fit_data.PD_dynamic{iTest}(loc_underscore(end-1)+1:loc_underscore(end)-1))/1000;
        end

        % hold(app.axes1_mod_tab,'on');
        % hold(app.axes2_mod_tab,'on');
        for i = 1:size(stiffening_ratio_exp,1)
            plot(app.axes1_mod_tab,freq,stiffening_ratio_exp(i,:),"Color",colour_scheme{i},"Marker","o","LineStyle","-.");
            plot(app.axes2_mod_tab,freq,stiffening_ratio_mod(i,:),"Color",colour_scheme{i},"Marker","o","LineStyle","-.");
        end
        ylabel(app.axes1_mod_tab,'Dynamic-to-QS stiffening ratio [-]','FontSize',16);
        xlabel(app.axes1_mod_tab,'Loading frequency [Hz]','FontSize',16);
        ylabel(app.axes2_mod_tab,'Dynamic-to-QS stiffening ratio [-]','FontSize',16);
        xlabel(app.axes2_mod_tab,'Loading frequency [Hz]','FontSize',16);
        legend(app.axes1_mod_tab,{'40-80 mmHg, exp.' '80-120 mmHg, exp.' '120-160 mmHg, exp.'},'FontSize',14);
        legend(app.axes2_mod_tab,{'40-80 mmHg, mod.' '80-120 mmHg, mod.' '120-160 mmHg, mod.'},'FontSize',14);
        axis(app.axes1_mod_tab,[0 max(freq)+2 0.8 round(max(max(stiffening_ratio_exp))*1.1*10)/10]);
        axis(app.axes2_mod_tab,[0 max(freq)+2 0.8 round(max(max(stiffening_ratio_exp))*1.1*10)/10]);

        % hold(app.axes1_mod_tab,'off');
        % hold(app.axes2_mod_tab,'off');
    elseif(strcmp(app.graph_choice_mod_tab_2.Value,'Pressure sweeps (active)') && fit_data.modelFormulation.activeModel.activator == 1)
        contractionStates = fieldnames(fit_data.fittedCurves.activeFit);
    
        lb_x = 10^7;
        ub_x = 0;
        ub_y = 0;

        PD = 'PD_100';
    
        color_sequence = {'r','b','g'};
        color_sybol_sequence = {'ro','bo','go'};
        legend_entry = strings(1,length(contractionStates)*2);

        for i = 1:length(contractionStates)    
            dM = fit_data.fittedCurves.activeFit.(contractionStates{i}).(PD);
    
            if(strcmp(app.graph_choice_mod_tab.Value,'Pressure and force'))
                % hold (app.axes1_mod_tab, 'on')
                plot(app.axes1_mod_tab,dM.("Outer diameter [mm]"), dM.("Exp. pressure [mmHg]"),color_sybol_sequence{i}) % Plot experimental P-D relationship
                plot(app.axes1_mod_tab,dM.("Outer diameter [mm]"), dM.("Modelled pressure [mmHg]"),color_sequence{i}) % Plot fitted P-D relationship
                % hold (app.axes1_mod_tab, 'off')

                xlabel(app.axes1_mod_tab, 'Outer diameter [mm]','FontSize',16)
                ylabel(app.axes1_mod_tab, 'Pressure [mmHg]','FontSize',16)
    
                lb_x_prev = lb_x;
                lb_x = 0.95*min(dM.("Outer diameter [mm]"));
                if(lb_x >= lb_x_prev)
                    lb_x = lb_x_prev;
                end
    
                ub_x_prev = ub_x;
                ub_x = 1.05*max(dM.("Outer diameter [mm]"));
                if(ub_x <= ub_x_prev)
                    ub_x = ub_x_prev;
                end
    
                axis(app.axes1_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 200]);

                % hold (app.axes2_mod_tab, 'on');
                plot(app.axes2_mod_tab,dM.("Exp. pressure [mmHg]"), dM.("Exp. transducer force [g]"),color_sybol_sequence{i}) % Plot experimental fT-P relationship
                plot(app.axes2_mod_tab,dM.("Modelled pressure [mmHg]"), dM.("Modelled transducer force [g]"),color_sequence{i}) % Plot fitted fT-P relationship
                % hold (app.axes2_mod_tab, 'off');

                xlabel(app.axes2_mod_tab, 'Pressure [mmHg]','FontSize',16)
                ylabel(app.axes2_mod_tab, 'Axial force [g]','FontSize',16)
               
                axis(app.axes2_mod_tab,[0 200 0 1.05*max(max(dM.("Exp. transducer force [g]")),max(dM.("Modelled transducer force [g]")))]);
    
            elseif(strcmp(app.graph_choice_mod_tab.Value,'Stresses'))    
                % hold (app.axes1_mod_tab, 'on')
                plot(app.axes1_mod_tab,dM.("Circ. stretch [-]"), dM.("Exp. circ. stress [MPa]"),color_sybol_sequence{i}) % Plot experimental P-D relationship
                plot(app.axes1_mod_tab,dM.("Circ. stretch [-]"), dM.("Modelled circ. stress [MPa]"),color_sequence{i}) % Plot fitted P-D relationship
                % hold (app.axes1_mod_tab, 'off')

                xlabel(app.axes1_mod_tab, 'Circumferential stretch [-]','FontSize',16)
                ylabel(app.axes1_mod_tab, 'Circumferential stress [MPa]','FontSize',16)
    
                lb_x_prev = lb_x;
                lb_x = 0.95*min(dM.("Circ. stretch [-]"));
                if(lb_x >= lb_x_prev)
                    lb_x = lb_x_prev;
                end
    
                ub_x_prev = ub_x;
                ub_x = 1.05*max(dM.("Circ. stretch [-]"));
                if(ub_x <= ub_x_prev)
                    ub_x = ub_x_prev;
                end
    
                ub_y_prev = ub_y;
                ub_y = max(max(dM.("Exp. circ. stress [MPa]")),max(dM.("Exp. axial stress [MPa]")));
                if(ub_y <= ub_y_prev)
                    ub_y = ub_y_prev;
                end
    
                axis(app.axes1_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 round(ub_y*10)/10]);

                % hold (app.axes2_mod_tab, 'on');
                plot(app.axes2_mod_tab,dM.("Circ. stretch [-]"), dM.("Exp. axial stress [MPa]"),color_sybol_sequence{i}) % Plot experimental fT-P relationship
                plot(app.axes2_mod_tab,dM.("Circ. stretch [-]"), dM.("Modelled axial stress [MPa]"),color_sequence{i}) % Plot fitted fT-P relationship
                % hold (app.axes2_mod_tab, 'off');

                xlabel(app.axes2_mod_tab, 'Circumferential stretch [-]','FontSize',16)
                ylabel(app.axes2_mod_tab, 'Axial stress [MPa]','FontSize',16)
    
                lb_x_prev = lb_x;
                lb_x = 0.95*min(dM.("Circ. stretch [-]"));
                if(lb_x >= lb_x_prev)
                    lb_x = lb_x_prev;
                end
    
                ub_x_prev = ub_x;
                ub_x = 1.05*max(dM.("Circ. stretch [-]"));
                if(ub_x <= ub_x_prev)
                    ub_x = ub_x_prev;
                end
    
                axis(app.axes2_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 round(ub_y*10)/10]);
    
            end
            legend_entry{2*i-1} = ['\lambda_z = ' PD(4:end) '% in vivo, exp'];
            legend_entry{2*i} = ['\lambda_z = ' PD(4:end) '% in vivo, mod'];
        end
    
        legend(app.axes1_mod_tab,legend_entry{1:i*2},'Location','NorthWest','FontSize',13);
        legend(app.axes2_mod_tab,'off');
    elseif(strcmp(app.graph_choice_mod_tab_2.Value,'Axial length sweeps'))
        if(isfield(fit_data,'FL'))
            FL = fit_data.FL;
                    
            lb_x_1 = 10^7;
            ub_x_1 = 0;
            ub_y_1 = 0;
            lb_x_2 = 10^7;
            ub_x_2 = 0;
            ub_y_2 = 0;
        
            color_sequence = {'r','b','g','m','y','c'};
            color_sybol_sequence = {'ro','bo','go','mo','yo','co'};
            legend_entry = strings(1,length(FL)*2);

            for i = 1:length(FL)
        
                dM = fit_data.fittedCurves.(fieldName).(FL{i});
        
                if(strcmp(app.graph_choice_mod_tab.Value,'Pressure and force'))
                    % hold (app.axes1_mod_tab, 'on')
                    plot(app.axes1_mod_tab,dM.("Outer diameter [mm]"), dM.("Exp. pressure [mmHg]"),color_sybol_sequence{i}) % Plot experimental P-D relationship
                    plot(app.axes1_mod_tab,dM.("Outer diameter [mm]"), dM.("Modelled pressure [mmHg]"),color_sequence{i}) % Plot fitted P-D relationship
                    % hold (app.axes1_mod_tab, 'off')

                    xlabel(app.axes1_mod_tab, 'Outer diameter [mm]','FontSize',16)
                    ylabel(app.axes1_mod_tab, 'Pressure [mmHg]','FontSize',16)
        
                    lb_x_prev_1 = lb_x_1;
                    lb_x_1 = 0.95*min(dM.("Outer diameter [mm]"));
                    if(lb_x_1 >= lb_x_prev_1)
                        lb_x_1 = lb_x_prev_1;
                    end
        
                    ub_x_prev_1 = ub_x_1;
                    ub_x_1 = 1.05*max(dM.("Outer diameter [mm]"));
                    if(ub_x_1 <= ub_x_prev_1)
                        ub_x_1 = ub_x_prev_1;
                    end
        
                    axis(app.axes1_mod_tab,[round(lb_x_1*10)/10 round(ub_x_1*10)/10 0 200]);
        
                    % hold (app.axes2_mod_tab, 'on');
                    plot(app.axes2_mod_tab,dM.("Exp. pressure [mmHg]"), dM.("Exp. transducer force [g]"),color_sybol_sequence{i}) % Plot experimental fT-P relationship
                    plot(app.axes2_mod_tab,dM.("Modelled pressure [mmHg]"),dM.("Modelled transducer force [g]"),color_sequence{i}) % Plot fitted fT-P relationship
                    % hold (app.axes2_mod_tab, 'off');

                    xlabel(app.axes2_mod_tab, 'Pressure [mmHg]','FontSize',16)
                    ylabel(app.axes2_mod_tab, 'Axial force [g]','FontSize',16)
        
                    lb_x_prev_2 = lb_x_2;
                    lb_x_2 = 0.95*min(dM.("Exp. pressure [mmHg]"));
                    if(lb_x_2 >= lb_x_prev_2)
                        lb_x_2 = lb_x_prev_2;
                    end
        
                    ub_x_prev_2 = ub_x_2;
                    ub_x_2 = 1.05*max(dM.("Exp. pressure [mmHg]"));
                    if(ub_x_2 <= ub_x_prev_2)
                        ub_x_2 = ub_x_prev_2;
                    end
        
                    axis(app.axes2_mod_tab,[round(lb_x_2*10)/10 round(ub_x_2*10)/10 0 1.05*max(max(dM.("Exp. transducer force [g]")),max(dM.("Modelled transducer force [g]")))]);
        
                elseif(strcmp(app.graph_choice_mod_tab.Value,'Stresses'))
                    % hold (app.axes1_mod_tab, 'on')
                    plot(app.axes1_mod_tab,dM.("Circ. stretch [-]"), dM.("Exp. circ. stress [MPa]"),color_sybol_sequence{i}) % Plot experimental P-D relationship
                    plot(app.axes1_mod_tab,dM.("Circ. stretch [-]"), dM.("Modelled circ. stress [MPa]"),color_sequence{i}) % Plot fitted P-D relationship
                    % hold (app.axes1_mod_tab, 'off')

                    xlabel(app.axes1_mod_tab, 'Circumferential stretch [-]','FontSize',16)
                    ylabel(app.axes1_mod_tab, 'Circumferential stress [MPa]','FontSize',16)
        
                    lb_x_prev_1 = lb_x_1;
                    lb_x_1 = 0.95*min(dM.("Circ. stretch [-]"));
                    if(lb_x_1 >= lb_x_prev_1)
                        lb_x_1 = lb_x_prev_1;
                    end
        
                    ub_x_prev_1 = ub_x_1;
                    ub_x_1 = 1.05*max(dM.("Circ. stretch [-]"));
                    if(ub_x_1 <= ub_x_prev_1)
                        ub_x_1 = ub_x_prev_1;
                    end
        
                    ub_y_prev_1 = ub_y_1;
                    ub_y_1 = max(max(dM.("Exp. circ. stress [MPa]")),max(dM.("Exp. axial stress [MPa]")));
                    if(ub_y_1 <= ub_y_prev_1)
                        ub_y_1 = ub_y_prev_1;
                    end
        
                    axis(app.axes1_mod_tab,[round(lb_x_1*10)/10 round(ub_x_1*10)/10 0 round(ub_y_1*10)/10]);

                    % hold (app.axes2_mod_tab, 'on');
                    plot(app.axes2_mod_tab,dM.("Axial stretch [-]"), dM.("Exp. axial stress [MPa]"),color_sybol_sequence{i}) % Plot experimental fT-P relationship
                    plot(app.axes2_mod_tab,dM.("Axial stretch [-]"), dM.("Modelled axial stress [MPa]"),color_sequence{i}) % Plot fitted fT-P relationship
                    % hold (app.axes2_mod_tab, 'off');

                    xlabel(app.axes2_mod_tab, 'Axial stretch [-]','FontSize',16)
                    ylabel(app.axes2_mod_tab, 'Axial stress [MPa]','FontSize',16)
        
                    lb_x_prev_2 = lb_x_2;
                    lb_x_2 = 0.95*min(dM.("Axial stretch [-]"));
                    if(lb_x_2 >= lb_x_prev_2)
                        lb_x_2 = lb_x_prev_2;
                    end
        
                    ub_x_prev_2 = ub_x_2;
                    ub_x_2 = 1.05*max(dM.("Axial stretch [-]"));
                    if(ub_x_2 <= ub_x_prev_2)
                        ub_x_2 = ub_x_prev_2;
                    end

                    ub_y_prev_2 = ub_y_2;
                    ub_y_2 = max(max(dM.("Exp. circ. stress [MPa]")),max(dM.("Exp. axial stress [MPa]")));
                    if(ub_y_2 <= ub_y_prev_2)
                        ub_y_2 = ub_y_prev_2;
                    end
        
                    axis(app.axes2_mod_tab,[round(lb_x_2*10)/10 round(ub_x_2*10)/10 0 round(ub_y_2*10)/10]);
        
                end
                legend_entry{2*i-1} = ['P = ' FL{i}(4:end) ' mmHg, exp'];
                legend_entry{2*i} = ['P = ' FL{i}(4:end) ' mmHg, mod'];
            end
        
            legend(app.axes1_mod_tab,legend_entry{1:i*2},'Location','NorthWest','FontSize',13);
        else
            app = appErrorDialogBox(app,'off','off','on','ERROR: Missing data',...
                    'Totally or partially missing data of the axial stretch sweeps/dynamic loops. Axial stretch sweeps data will not have diameter tracking if if data were acquired using DynamX-s1.');

            app.graph_choice_mod_tab_2.Value = 'Pressure sweeps';
        end
        hold (app.axes1_mod_tab, 'off');
        hold (app.axes2_mod_tab, 'off');
    end