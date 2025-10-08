function [passive,app] = findInVivoAxialStretch(app,passive,FLexperiments,colour_scheme)
% This function is used to find the intersection point between the
% force sweep experiments at 60, 100 and 140 mmHg and estimate the in vivo
% axial stretch. The relevant information is then saved in the structure
% "passive" (field preliminary_analysis.in_vivo_axial_stretch").

    app.axial_stretch_menu.Value = 'In vivo axial stretch';

%% Importing and averaging loading/unloading of force sweep experiments    
    % Defining axial force and axial stretch tensors
    axialStretchV = struct;
    axialForceV = struct;
    
    legend_entry_FL = strings(3,1);
    
    for i = 2:4 % excluding 1st and 5th test - 10mmHg and 200 mmHg pressures
        dM = (table2array(passive.static_data.unloading.biaxial_Fl.(FLexperiments{i}))+...
            flip(table2array(passive.static_data.loading.biaxial_Fl.(FLexperiments{i}))))/2;
        axialStretchV(i).a = dM(:,5);
        axialForceV(i).b = dM(:,4);
    
        locator_force = find(axialStretchV(i).a(2:end)<axialStretchV(i).a(1:end-1));
        
        if(i>2)
            hold(app.axes4_dat_proc_tab, 'on');
        end
        axialStretchV(i).a = axialStretchV(i).a(locator_force);
        axialForceV(i).b = axialForceV(i).b(locator_force);
    
        plot(app.axes4_dat_proc_tab,axialStretchV(i).a,axialForceV(i).b,colour_scheme{i-1});
        legend_entry_FL{i-1} = [FLexperiments{i}(4:end) ' mmHg'];
    
        passive.preliminary_analysis.in_vivo_axial_stretch.curves.(FLexperiments{i}) = array2table([dM(:,5),dM(:,4)],'VariableName',["Axial stretch [-]","Force [g]"]);
    end
    
%% Finding intersection point
    combinationsV = nchoosek(2:4,2);
    intersection = zeros(3,2);
    
    axStretchRangeV = 0.4:0.0001:3;
    axStretchFoundV = zeros(size(combinationsV,1),1);
    for iComb = 1:size(combinationsV,1)
        i1 = combinationsV(iComb, 1);
        i2 = combinationsV(iComb, 2);
        axForce1V = axialForceV(i1).b;
        axForce2V = axialForceV(i2).b;
        axStretch1V = axialStretchV(i1).a;
        axStretch2V = axialStretchV(i2).a;
        
        axForceInt1V = interp1(axStretch1V, axForce1V, axStretchRangeV);
        axForceInt2V = interp1(axStretch2V, axForce2V, axStretchRangeV);
        
        axForceDiffV = axForceInt2V-axForceInt1V;
        [~,iMin] = min(abs(axForceDiffV));
        axStretchFound = axStretchRangeV(iMin);
        axStretchFoundV(iComb) = axStretchFound; 
        
        intersection(iComb,1:2) = [axStretchFound mean([axForceInt1V(iMin),axForceInt2V(iMin)])];
        
        xlabel(app.axes4_dat_proc_tab,'Axial stretch [-]','FontSize',16);
        ylabel(app.axes4_dat_proc_tab,'Axial force [g]','FontSize',16);
    end

    plot(app.axes4_dat_proc_tab,intersection(:,1), intersection(:,2),'kx','LineWidth',3)
        
    legend_entry_FL{length(legend_entry_FL)+1} = 'Intersection points';
    legend(app.axes4_dat_proc_tab, legend_entry_FL(1:i),'Location','NorthWest');
    
    hold(app.axes4_dat_proc_tab, 'off');

%% Saving to "passive" structure
    passive.preliminary_analysis.in_vivo_axial_stretch.intersection_points = array2table(intersection,'VariableName',["Axial stretch [-]","Force [g]"]);
    passive.preliminary_analysis.in_vivo_axial_stretch.value = mean(axStretchFoundV);
end