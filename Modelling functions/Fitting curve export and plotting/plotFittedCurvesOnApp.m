function plotFittedCurvesOnApp(app,passive,PD,FL,viscousModel,fit_data)
% This function is used to plot the fitted curves on the app.
% Input parameters are the app data "app", the structure with passive
% experimental data "passive", the list of quasi-static pressure ("PD") or
% axial force ("FL") sweeps, the formulation of the viscous relaxation
% function ("viscousModel"), and the result structure ("fit_data").

    if(viscousModel.activator == 1)
        PD_dynamic = passive.dynamic_data.PD(4:end);
        test_list = [PD, FL, PD_dynamic];
    else
        test_list = [PD, FL];
    end

    if(isfield(fit_data.fittedCurves,'quasiStaticFit'))
        fieldName = 'quasiStaticFit';
    else
        fieldName = 'staticFit';
    end
    
    lb_x = 10^7; %for plots
    ub_x = 0; %for plots
    
    for i = 1:length(test_list)
        if(i<=length(PD))
            dataMat = table2array(fit_data.fittedCurves.(fieldName).(PD{i}));
        elseif(i>length(PD) && i<=length(PD)+length(FL))
            dataMat = table2array(fit_data.fittedCurves.(fieldName).(FL{i-length(PD)}));
        else
            dataMat = table2array(fit_data.fittedCurves.dynamicFit.(PD_dynamic{i-length(PD)-length(FL)}));
        end
        
        color_sequence = {'r','b','g','m','k','c'};
        color_symbol_sequence = {'ro','bo','go','mo','ko','co'};
        
        if((i<=length(PD)) && strcmp(app.graph_choice_mod_tab_2.Value,'Pressure sweeps'))
    
            plot_function_GUI(app, dataMat, color_symbol_sequence, color_sequence,i);
    
            legend_entry{2*i-1} = ['\lambda_z = ' PD{i}(4:end) '% in vivo, exp'];
            legend_entry{2*i} = ['\lambda_z = ' PD{i}(4:end) '% in vivo, mod'];
    
            legend(app.axes1_mod_tab,legend_entry,'Location','NorthWest','FontSize',13);
    
        elseif(strcmp(app.graph_choice_mod_tab_2.Value,'Axial length sweeps'))
            if(isempty('FL'))
                app.graph_choice_mod_tab_2.Value = 'Pressure sweeps';
                errordlg('Totally or partially missing data of the axial stretch sweeps. This is likely the case if data were acquired using DynamX-1','ERROR: Missing data.');
                plot_function_GUI(app,dataMat, color_symbol_sequence, color_sequence,i);
    
            elseif(i>length(PD) && i<=length(PD)+length(FL))
                plot_function_GUI(app, dataMat, color_symbol_sequence, color_sequence,i-length(PD));
                legend_entry{2*i-1-2*length(PD)} = ['P = ' FL{i-length(PD)}(4:end) ' mmHg, exp'];
                legend_entry{2*i-2*length(PD)} = ['P = ' FL{i-length(PD)}(4:end) ' mmHg, mod'];
                legend(app.axes1_mod_tab,legend_entry,'Location','NorthWest','FontSize',13);
            end
        
        elseif(strcmp(app.graph_choice_mod_tab_2.Value(1:13),'Dynamic loops'))
            if(isempty('PD_dynanamic'))
                app.graph_choice_mod_tab_2.Value = 'Pressure sweeps';
                errordlg('Totally or partially missing dynamic loops data. Check the choosen fitting options and the imported data.','ERROR: Missing data.');
                plot_function_GUI(app, dataMat, color_symbol_sequence, color_sequence,i-length(PD)-length(FL));
    
            elseif(i<=length(PD) && strcmp(PD{i},'PD_100'))
                plot_function_GUI(app, dataMat, {'--k'}, {'k'},1);
                lb_x_prev = lb_x;
                lb_x = 0.95*min(dataMat(:,4));
                if(lb_x >= lb_x_prev)
                    lb_x = lb_x_prev;
                end
                
                ub_x_prev = ub_x;
                ub_x = 1.05*max(dataMat(:,4));
                if(ub_x <= ub_x_prev)
                    ub_x = ub_x_prev;
                end
    
            elseif(i>length(PD)+length(FL))
                if(strcmp(app.graph_choice_mod_tab_2.Value,'Dynamic loops (physiological)') && strcmp(PD_dynamic{i-length(PD)-length(FL)}(end-5:end),'physio'))
                    plot_function_GUI(app, dataMat, color_symbol_sequence, color_sequence,i-length(PD)-length(FL));
                    axis(app.axes1_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 200]);
    
                elseif(strcmp(app.graph_choice_mod_tab_2.Value,'Dynamic loops (sine)') && strcmp(PD_dynamic{i-length(PD)-length(FL)}(end-3:end),'sine'))
                    if(i-length(PD)-length(FL) == 1)
                        colour_code = 1;
                    elseif(strcmp(PD_dynamic{i-length(PD)-length(FL)}(3:6), PD_dynamic{i-length(PD)-length(FL)-1}(3:6)))
                        colour_code = colour_code+1;
                    else
                        colour_code = 1;
                    end
                    plot_function_GUI(app, dataMat, color_symbol_sequence, color_sequence,colour_code,2);
                    axis(app.axes1_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 200]);
                end
            end
        end
    end
end