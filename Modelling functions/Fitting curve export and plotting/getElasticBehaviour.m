function [staticFitting] = getElasticBehaviour(inputMatElasticFitting,parameterV,elastic)
% This function returns strcuture of tables in which each table reports the
% comparison between measured and modelled elastic behaviours. 
% The input parameters are the matrix with the experimental data
% "inputMatElasticFitting", the vector of SEF model parameters "parameterV"
% and the chosen elastic model (i.e., SEF and related details) "elastic".

    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    gtN = 0.0098; %scaling factor to convert g to N

    [sigmatt_rr_compV,sigmazz_rr_compV,~,~,~]=elastic.fun(parameterV,inputMatElasticFitting(:,7),inputMatElasticFitting(:,8),inputMatElasticFitting(:,9));
    
    [pressureV_comp,transducerForceV_comp] = stress2pressure_force(sigmatt_rr_compV,sigmazz_rr_compV,...
        (inputMatElasticFitting(:,2)+inputMatElasticFitting(:,3))/2,inputMatElasticFitting(:,4));

    sigmatt_compV = sigmatt_rr_compV-pressureV_comp/2*mtN; % Lagrange multiplier
    sigmazz_compV = sigmazz_rr_compV-pressureV_comp/2*mtN; % Lagrange multiplier
    
    staticFittingMat = [zeros(length(sigmatt_compV),1), inputMatElasticFitting(:,1)/mtN, pressureV_comp, inputMatElasticFitting(:,2)*2, inputMatElasticFitting(:,3)*2,...
        inputMatElasticFitting(:,7), inputMatElasticFitting(:,10), sigmatt_compV, inputMatElasticFitting(:,5)/gtN,...
         transducerForceV_comp, inputMatElasticFitting(:,8), inputMatElasticFitting(:,11), sigmazz_compV];

    separator = find(inputMatElasticFitting(:,end));
    
    test_ID = ["PD_95", "PD_100", "PD_105", "FL_10", "FL_60", "FL_100", "FL_140", "FL_180"];
    headings = ["Time [s]", "Exp. pressure [mmHg]", "Modelled pressure [mmHg]", "Outer diameter [mm]", "Inner diameter [mm]"...
            "Circ. stretch [-]", "Exp. circ. stress [MPa]", "Modelled circ. stress [MPa]", "Exp. transducer force [g]", "Modelled transducer force [g]"...
            "Axial stretch [-]", "Exp. axial stress [MPa]", "Modelled axial stress [MPa]"];
    for i = 1:length(separator)
        if(i==length(separator))
            staticFitting.(test_ID{i}) = array2table(flip(staticFittingMat(separator(i):end,:)),"VariableNames",headings);
        else
            staticFitting.(test_ID{i}) = array2table(flip(staticFittingMat(separator(i):separator(i+1)-1,:)),"VariableNames",headings);
        end
    end
end