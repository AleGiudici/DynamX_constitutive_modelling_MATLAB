function plot_function_GUI(app,dM,color_symbol_sequence,color_sequence,i,counter)
% This function is used to plot the fitted curves on the app.

if(~exist('counter'))
    counter = i;
end

lb_x = 10^7; %for plots
ub_x = 0; %for plots
ub_y = 0; %for plots

if(strcmp(app.graph_choice_mod_tab.Value,'Pressure and force'))
        
    if(counter==1)
        plot(app.axes1_mod_tab,dM(:,4), dM(:,2),color_symbol_sequence{i}) % Plot experimental P-D relationship
        hold (app.axes1_mod_tab, 'on')
        plot(app.axes1_mod_tab,dM(:,4), dM(:,3),color_sequence{i}) % Plot fitted P-D relationship
        hold (app.axes1_mod_tab, 'off')
    else
        hold (app.axes1_mod_tab, 'on')
        plot(app.axes1_mod_tab,dM(:,4), dM(:,2),color_symbol_sequence{i}) % Plot experimental P-D relationship
        plot(app.axes1_mod_tab,dM(:,4), dM(:,3),color_sequence{i}) % Plot fitted P-D relationship
        hold (app.axes1_mod_tab, 'off')
    end
    xlabel(app.axes1_mod_tab, 'Outer diameter [mm]','FontSize',14)
    ylabel(app.axes1_mod_tab, 'Pressure [mmHg]','FontSize',14)
    
    lb_x_prev = lb_x;
    lb_x = 0.95*min(dM(:,4));
    if(lb_x >= lb_x_prev)
        lb_x = lb_x_prev;
    end
    
    ub_x_prev = ub_x;
    ub_x = 1.05*max(dM(:,4));
    if(ub_x <= ub_x_prev)
        ub_x = ub_x_prev;
    end
    
    axis(app.axes1_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 200]);

    if(counter==1)
        plot(app.axes2_mod_tab,dM(:,2), dM(:,9),color_symbol_sequence{i}) % Plot experimental fT-P relationship
        hold (app.axes2_mod_tab, 'on');
        plot(app.axes2_mod_tab,dM(:,3), dM(:,10),color_sequence{i}) % Plot fitted fT-P relationship
        hold (app.axes2_mod_tab, 'off');
    else
        hold (app.axes2_mod_tab, 'on');
        plot(app.axes2_mod_tab,dM(:,2), dM(:,9),color_symbol_sequence{i}) % Plot experimental fT-P relationship
        plot(app.axes2_mod_tab,dM(:,3),dM(:,10),color_sequence{i}) % Plot fitted fT-P relationship
        hold (app.axes2_mod_tab, 'off');
    end
    xlabel(app.axes2_mod_tab, 'Pressure [mmHg]','FontSize',14)
    ylabel(app.axes2_mod_tab, 'Axial force [g]','FontSize',14)
    
    lb_x_prev = lb_x;
    lb_x = 0.95*min(dM(:,4));
    if(lb_x >= lb_x_prev)
        lb_x = lb_x_prev;
    end
    
    ub_x_prev = ub_x;
    ub_x = 1.05*max(dM(:,4));
    if(ub_x <= ub_x_prev)
        ub_x = ub_x_prev;
    end
    
    axis(app.axes2_mod_tab,[0 200 0 max(max(dM(:,9)),max(dM(:,10)))]);

else

    if(counter==1)
        plot(app.axes1_mod_tab,dM(:,6), dM(:,7),color_symbol_sequence{i}) % Plot experimental P-D relationship
        hold (app.axes1_mod_tab, 'on')
        plot(app.axes1_mod_tab,dM(:,6), dM(:,8),color_sequence{i}) % Plot fitted P-D relationship
        hold (app.axes1_mod_tab, 'off')
    else
        hold (app.axes1_mod_tab, 'on')
        plot(app.axes1_mod_tab,dM(:,6), dM(:,7),color_symbol_sequence{i}) % Plot experimental P-D relationship
        plot(app.axes1_mod_tab,dM(:,6), dM(:,8),color_sequence{i}) % Plot fitted P-D relationship
        hold (app.axes1_mod_tab, 'off')
    end
    xlabel(app.axes1_mod_tab, 'Outer diameter [mm]','FontSize',14)
    ylabel(app.axes1_mod_tab, 'Pressure [mmHg]','FontSize',14)
    
    lb_x_prev = lb_x;
    lb_x = 0.95*min(dM(:,6));
    if(lb_x >= lb_x_prev)
        lb_x = lb_x_prev;
    end
    
    ub_x_prev = ub_x;
    ub_x = 1.05*max(dM(:,6));
    if(ub_x <= ub_x_prev)
        ub_x = ub_x_prev;
    end
    
    ub_y_prev = ub_y;
    ub_y = max(max(dM(:,7)),max(dM(:,12)));
    if(ub_y <= ub_y_prev)
        ub_y = ub_y_prev;
    end
    
    axis(app.axes1_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 round(ub_y*10)/10]);

    if(counter==1)
        plot(app.axes2_mod_tab,dM(:,6), dM(:,12),color_symbol_sequence{i}) % Plot experimental fT-P relationship
        hold (app.axes2_mod_tab, 'on');
        plot(app.axes2_mod_tab,dM(:,6), dM(:,13),color_sequence{i}) % Plot fitted fT-P relationship
        hold (app.axes2_mod_tab, 'off');
    else
        hold (app.axes2_mod_tab, 'on');
        plot(app.axes2_mod_tab,dM(:,6), dM(:,12),color_symbol_sequence{i}) % Plot experimental fT-P relationship
        plot(app.axes2_mod_tab,dM(:,6), dM(:,13),color_sequence{i}) % Plot fitted fT-P relationship
        hold (app.axes2_mod_tab, 'off');
    end
    xlabel(app.axes2_mod_tab, 'Pressure [mmHg]','FontSize',14)
    ylabel(app.axes2_mod_tab, 'Axial force [g]','FontSize',14)
    
    lb_x_prev = lb_x;
    lb_x = 0.95*min(dM(:,5));
    if(lb_x >= lb_x_prev)
        lb_x = lb_x_prev;
    end
    
    ub_x_prev = ub_x;
    ub_x = 1.05*max(dM(:,5));
    if(ub_x <= ub_x_prev)
        ub_x = ub_x_prev;
    end
    
    axis(app.axes2_mod_tab,[round(lb_x*10)/10 round(ub_x*10)/10 0 round(ub_y*10)/10]);
end
end