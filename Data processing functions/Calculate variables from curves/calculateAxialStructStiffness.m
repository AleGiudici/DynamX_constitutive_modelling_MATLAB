function [passive] = calculateAxialStructStiffness(FLexperiments,passive)
% This function calculate a structural stiffness coefficient from the
% experimental force sweep data as the slope of the transducer force-axial
% stretch relationship at the in vivo axial stretch ([g]).
%
% The axial stiffness coefficients are then saved in the "passive"
% structure (field "preliminary_analysis.axial_stiffness").

    axial_stiffness_coeff = zeros(1,length(FLexperiments));
    FL_legend = strings(1,length(FLexperiments));

    for iTest = 1:size(FLexperiments,2)
        dM = (table2array(passive.static_data.unloading.biaxial_Fl.(FLexperiments{iTest}))+...
            flip(table2array(passive.static_data.loading.biaxial_Fl.(FLexperiments{iTest}))))/2;
    
        [~,min_Locator] = min(abs(dM(:,5)-passive.preliminary_analysis.in_vivo_axial_stretch.value));
        
        if(min_Locator>2 && min_Locator<size(dM,1)-1)
            y_vect = dM(min_Locator-2:min_Locator+2,4);
            x_vect = dM(min_Locator-2:min_Locator+2,5);
        else
            y_vect = dM(min_Locator-1:min_Locator+1,4);
            x_vect = dM(min_Locator-1:min_Locator+1,5);
        end
        
        poly_coeff = polyfit(x_vect,y_vect,1);
        axial_stiffness_coeff(iTest) = poly_coeff(1);    

        FL_legend{iTest} = [FLexperiments{iTest}(4:end) ' mmHg'];
    end
    
    row_headings = {'Axial stiffness coefficient [g]'};
    passive.preliminary_analysis.axial_stiffness = array2table(axial_stiffness_coeff,'RowNames',row_headings,'VariableNames',FL_legend);
end