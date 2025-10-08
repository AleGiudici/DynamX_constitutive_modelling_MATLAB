% --------------------------------------------
% DYNAMX Constitutive Modelling Function
% Version: 2.2
%
% Alessandro Giudici & Shaiv Parikh
% Reesink and Spronck's lab
%
% Department Biomedical Engineering
% CARIM - School for Cardiovascular Diseases
% Maastricht University
% 
% The function follows this steps:
% 1) Model selection
% 2) First round of parameter estimation (here the wall is modelled as a
%    purely elastic material) --> This yields and estimate of collagen's
%    and elastin deposition stretches.
% 3) Second round of parameter estimation (here the wall is modelled as a
%    passive viscoelastic material) --> This step yields an estimate of the
%    viscous (relaxation function) parameters.
% 4) Third round of parameter estimation (here the wall is modelled as a
%    passive viscoelastic material) --> This step yields the final passive
%    strain energy function model parameters
% 5) Fourth round of parameter estimation (here the wall is modelled as an
%    active viscoelatic material) --> This step yields the active VSMC
%    model parameters
% 6) Plotting and saving results.
%
% NOTE: steps 3-5 are optional.
% -------------------------------------------

function [app,fit_data,elastic_save,visco_save,active_save,in_vivo_ref_save,axial_sweeps_check_save,PD_save] = ParameterEstimationFun(app,event,queuing,FileName,data,elastic_save,visco_save,active_save,in_vivo_ref_save,axial_sweeps_check_save,PD_save,i_sample)


%% Clearing GUI output

tic;

app.static_R2_value.Text = '....';

app.dynamic_R2_value.Text = '....';
app.memory_fun.Text = '....';
app.N_param_dynamic.Text = '....';


app.in_vivo_axial_stretch_mod.Text = '....';
app.in_vivo_axial_stretch_exp.Text = '....';

app.SEF.Text = '....';
app.N_param_static.Text = '....';
app.N_param_constrained.Text = '....';
app.ref_conf.Text = '....';

app.status_txt.Text = 'Initialising...';
pause(0.1);

%% Loading experimental data
load(FileName);
if(~isfield(data,'passive'))
    passive = experimental_data.passive;
end

%% Scaling factors
mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
mtm = 0.001; %scaling factor to convert mircometer to mm
gtN = 0.0098; %scaling factor to convert g to N
rtd = 180/pi; %scaling factor to convert radians to degrees

%% Geometrical parameters from unloaded intact configuration
Ro = ((passive.OD)/2)*mtm;%outer radius (mm)
H = (passive.H)*mtm; %wall thickness (mm)
Ri = Ro-H; %inner radius (mm)
L = passive.L; %unloaded intact length (mm)
Vol = pi*((Ro^2-Ri^2)*L); %volume of segment (mm^3)

unloaded_configuration = [Ri Ro H L];

%% Choice of the modelling approach
if(queuing == 0)
    [elasticModel,viscousModel,activeModel,in_vivo_ref,PD,axial_sweeps_check] = select_constitutive_model_GUI(passive,app);
else
    clear fit_data;
    if(i_sample == 1)
        [elasticModel,viscousModel,activeModel,in_vivo_ref,PD,axial_sweeps_check] = select_constitutive_model_GUI(passive,app);
        elastic_save = elasticModel;
        visco_save = viscousModel;
        active_save = activeModel;
        in_vivo_ref_save = in_vivo_ref;
        axial_sweeps_check_save = axial_sweeps_check;
        PD_save = PD;
    else
        elasticModel = elastic_save;
        viscousModel = visco_save;
        activeModel = active_save;
        in_vivo_ref = in_vivo_ref_save;
        axial_sweeps_check = axial_sweeps_check_save;
        PD = PD_save;
        elasticModel.imposedParam(11) = passive.preliminary_analysis.in_vivo_axial_stretch;%passive.preliminary_analysis.in_vivo_axial_stretch.value;
    end
end

if(activeModel.activator == 1)
    elastic_for_active = elasticModel;
end

fit_data.modelFormulation.elasticModel = elasticModel;
fit_data.modelFormulation.activeModel = activeModel;
fit_data.modelFormulation.viscousModel = viscousModel;

% fitting options
% elastic --> information on the quasi-static (elastic) modelling choices
% viscous --> information on the dynamic (viscous) modelling choices
%   in_vivo_ref --> sets the reference configuration to the unloaded
%   configuration or "in vivo" (100 mmHg, in vivo axial stretch)
%   configuration
% PD --> list of quasi-static pressure-diameter experiments

if(axial_sweeps_check)
    FL = passive.static_data.FL;
else
    FL = {};
end

%% Setting the reference configuration
[l_ref,ri_ref,h_ref,ro_ref,~,~,~] = setReferenceConfiguration(passive,in_vivo_ref);
reference_configuration = [ri_ref,ro_ref,h_ref,l_ref];

clear ri_ref ro_ref h_ref l_ref;

%% Importing data quasi-static experiments and creating input for elastic fitting

app.status_txt.Text = 'Loading data...';
pause(0.1);

[inputMatElasticFitting,weights_vector] = createInputMatrixElasticFitting(passive,PD,axial_sweeps_check,reference_configuration);

%% Optmisation steps

app.status_txt.Text = 'Estimating elastic model parameters...';
pause(0.1);

% set the upper bound of the collagen prestretch to the stretch between
% unloaded and reference configuration
% if(in_vivo_ref)
%     elasticModel.ub(end-1) = (1/inputMatElasticFitting(105,7));
%     elasticModel.lb(end-1) = (1/inputMatElasticFitting(105,7));
% end

check = 0;
while(check == 0)
    [parameterV_elastic_static,R2_static,RMSE_static] = runElasticFitting(inputMatElasticFitting,elasticModel);
        
    test_vector = parameterV_elastic_static-elasticModel.normV_baseline.*elasticModel.ub_baseline;
    locator_boundary_found = find(test_vector == 0);

    if(isempty(locator_boundary_found))
        check = 1;
    else
        if(locator_boundary_found == 9)
            check = 1;
        else
            elasticModel.normV_baseline(locator_boundary_found) = elasticModel.normV_baseline(locator_boundary_found)*1.2;
        end
    end
end

staticFitting = getElasticBehaviour(inputMatElasticFitting,parameterV_elastic_static,elasticModel);

if(viscousModel.activator == 0)
    fit_data.fittedCurves.staticFit = staticFitting;
end

elasticModel.imposedParam(end) = parameterV_elastic_static(end);
elasticModel.imposedParam(end-2) = parameterV_elastic_static(end-2);
elasticModel.normV = elasticModel.normV(1:end-2);
elasticModel.lb = elasticModel.lb(1:end-2);
elasticModel.ub = elasticModel.ub(1:end-2);
elasticModel.nParam = elasticModel.nParam-2;

clear inputMatElasticFitting test_vector locator_boundary_found;

if(viscousModel.activator==1)
    %% Importing data dynamic experiments
    [inputMatViscoelasticFitting1,inputMatViscoelasticFitting2] = createInputStructuresViscoelasticFitting(passive,PD,axial_sweeps_check,reference_configuration);    
              
    %% Optmisation steps - first round to determine the viscous parameters
    app.status_txt.Text = 'Estimating viscous model parameters...';
    pause(0.1);
    
    % Initialise the output data structure
    % fit_data = [];
    [parameterV_elastic_1stQLV,parameterV_viscous,RMSE_dynamic,R2_dynamic,fit_data] = runViscoelasticFitting(inputMatViscoelasticFitting2,elasticModel,viscousModel,ones(length(fieldnames(inputMatViscoelasticFitting2)),1),fit_data,1);
           
    %% Optmisation steps - second round to determine the elastic parameters
    app.status_txt.Text = 'Refining elastic model parameters...';
    pause(0.1);

    visco_2nd = viscousModel;

    viscous_normV = zeros(1,3*elasticModel.nConstituents);
    for i = 1:elasticModel.nConstituents
        viscous_normV(1+(i-1)*3:i*3) = viscousModel.normV;
    end
    visco_2nd.lb = parameterV_viscous./viscous_normV;
    visco_2nd.ub = parameterV_viscous./viscous_normV;
    
    [parameterV_elastic_2ndQLV,parameterV_viscous,RMSE_quasi_static,R2_quasi_static,fit_data] = runViscoelasticFitting(inputMatViscoelasticFitting1,elasticModel,visco_2nd,weights_vector,fit_data,2);
end

if(activeModel.activator==1)
    app.status_txt.Text = 'Estimating active model parameters...';
    pause(0.1);

    [inputMatActiveFitting,n_dataPoints,reference_configuration_for_active] = createInputStructureActiveFitting(experimental_data,reference_configuration);
    
    if(strcmp(activeModel.modelName,'Franchini 2022'))
        activeModel.lb_baseline(5) = (reference_configuration_for_active(1)+reference_configuration_for_active(2))/...
            (unloaded_configuration(1)+unloaded_configuration(2))/2;
        activeModel.ub_baseline(5) = activeModel.lb_baseline(5);
        activeModel.lb_baseline(6) = reference_configuration_for_active(4)/unloaded_configuration(4);
        activeModel.ub_baseline(6) = activeModel.lb_baseline(6);
    elseif(strcmp(activeModel.modelName,'Gaussian model'))
        activeModel.normV_baseline(5) = passive.preliminary_analysis.in_vivo_axial_stretch.value*passive.L/reference_configuration_for_active(4);
        activeModel.normV_baseline(4) = 1.1;
    end

    [parameterV_active,parameterV_viscous_active,R2_active,RMSE_active,fit_data] = runActiveFitting(inputMatActiveFitting,elasticModel,viscousModel,activeModel,parameterV_elastic_2ndQLV,...
        parameterV_viscous,reference_configuration_for_active,unloaded_configuration,n_dataPoints,ones(2,7),fit_data);
end

toc
%% Plots
%Plots from experimental data vs fitted curves
plotFittedCurvesOnApp(app,passive,PD,FL,viscousModel,fit_data);

%% Creating output structure
if(viscousModel.activator == 1)
    paramMat = [parameterV_elastic_2ndQLV',parameterV_elastic_1stQLV',parameterV_elastic_static',elasticModel.lb_baseline'.*elasticModel.normV_baseline',elasticModel.ub_baseline'.*elasticModel.normV_baseline',elasticModel.imposedParam'.*elasticModel.normV_baseline'];
    fit_data.modelParameters.staticSEFparameters = array2table(paramMat,'RowNames',elasticModel.paramName,'VariableNames',["Final QLV values","First QLV values","First QS values","Lower boundary","Upper boundary","Constrained values"]);
    if(activeModel.activator == 1)
        fit_data.fittingQuality = array2table([R2_dynamic, RMSE_dynamic; R2_quasi_static, RMSE_quasi_static; R2_static, RMSE_static; R2_active, RMSE_active],'VariableNames',{'R2', 'RMSE'},'RowNames',{'Dynamic', 'Quasi-static', 'Static', 'Active'});
    else
        fit_data.fittingQuality = array2table([R2_dynamic, RMSE_dynamic; R2_quasi_static, RMSE_quasi_static; R2_static, RMSE_static],'VariableNames',{'R2', 'RMSE'},'RowNames',{'Dynamic', 'Quasi-static', 'Static'});
    end
else
    paramMat = [parameterV_elastic_static',elasticModel.lb_baseline'.*elasticModel.normV_baseline',elasticModel.ub_baseline'.*elasticModel.normV_baseline',elasticModel.imposedParam'.*elasticModel.normV_baseline'];
    fit_data.modelParameters.staticSEFparameters = array2table(paramMat,'RowNames',elasticModel.paramName,'VariableNames',["First QS values","Lower boundary","Upper boundary","Constrained values"]);
    fit_data.fittingQuality = array2table([R2_static, RMSE_static],'VariableNames',{'R2', 'RMSE'},'RowNames',{'Static'});
end
% fit_data.staticR2 = R2_static;
if(isfield(passive.preliminary_analysis.in_vivo_axial_stretch,'value'))
    fit_data.in_vivo_axial_stretch = passive.preliminary_analysis.in_vivo_axial_stretch.value*passive.L/reference_configuration(4);
else
    fit_data.in_vivo_axial_stretch = passive.preliminary_analysis.in_vivo_axial_stretch*passive.L/reference_configuration(4);
end
headings_geometry = ["Inner radius [mm]","Outer radius [mm]","Thickness [mm]","Axial length [mm]"];
fit_data.vesselConfigurations.unloaded_configuration = array2table(unloaded_configuration,"VariableNames",headings_geometry);

fit_data.vesselConfigurations.reference_configuration = array2table(reference_configuration,"VariableNames",headings_geometry);

if(activeModel.activator == 1)
    fit_data.vesselConfigurations.reference_configuration_active = array2table(reference_configuration_for_active,"VariableNames",headings_geometry);
end

experimental_data.passive = passive;

fit_data.PD_static = PD;
fit_data.FL = FL;
if(viscousModel.activator == 1)
    fit_data.PD_dynamic = experimental_data.passive.dynamic_data.PD(4:end);
end

if(activeModel.activator == 1)      
    parameterV_active_converted(4) = parameterV_active(4)*(reference_configuration_for_active(1)+reference_configuration_for_active(2))/...
        (reference_configuration(1)+reference_configuration(2));   
    parameterV_active_converted(5) = parameterV_active(5)*reference_configuration_for_active(4)/reference_configuration(4);
    parameterV_active_converted(1:3) = parameterV_active(1:3);
    fit_data.modelParameters.activeSEFparameters = array2table([parameterV_active_converted',parameterV_active',(activeModel.lb_baseline.*activeModel.normV_baseline)',(activeModel.ub_baseline.*activeModel.normV_baseline)'],...
        'VariableNames',["Final values","Fitted values","Lower boundary","Upper boundary"],'RowNames',activeModel.paramName);
end

if(viscousModel.activator == 1 && activeModel.activator == 0)
    for i = 1:elasticModel.nConstituents
        lb_viscous(1+(i-1)*3:i*3) = viscousModel.lb;
        ub_viscous(1+(i-1)*3:i*3) = viscousModel.ub;
        imposed_viscous(1+(i-1)*3:i*3) = [NaN 0.001 NaN];
    end

    row_headings = {'Elastin viscous gain [-]','Elastin time constant 1 [s]','Elastin time constant 2 [s]'...
        'Collagen viscous gain [-]','Collagen time constant 1 [s]','Collagen time constant 2 [s]'};
    fit_data.modelParameters.viscousParameters=array2table([parameterV_viscous',lb_viscous',ub_viscous',imposed_viscous'],'VariableNames',...
        ["Values","Lower boundary","Upper boundary","Imposed parameters"],'RowNames',row_headings);
    fit_data.PD_dynamic = passive.dynamic_data.PD(4:end);
elseif(viscousModel.activator == 1 && activeModel.activator == 1)
    for i = 1:elasticModel.nConstituents+activeModel.nConstituents
        lb_viscous(1+(i-1)*3:i*3) = viscousModel.lb;
        ub_viscous(1+(i-1)*3:i*3) = viscousModel.ub;
        imposed_viscous(1+(i-1)*3:i*3) = [NaN 0.001 NaN];
    end
    parameterV_viscous_active(1:elasticModel.nConstituents*viscousModel.nParam) = parameterV_viscous;

    row_headings = {'Elastin viscous gain [-]','Elastin time constant 1 [s]','Elastin time constant 2 [s]'...
        'Collagen viscous gain [-]','Collagen time constant 1 [s]','Collagen time constant 2 [s]'...
        'VSMC viscous gain [-]','VSMC time constant 1 [s]','VSMC time constant 2 [s]'};

    fit_data.modelParameters.viscousParameters = array2table([parameterV_viscous_active',lb_viscous',ub_viscous',imposed_viscous'],'VariableNames',...
        ["Values","Lower boundary","Upper boundary","Imposed parameters"],'RowNames',row_headings);
end


%% P vs. D by back substitution
fit_data = simulateStaticBehaviour(fit_data);

app.status_txt.Text = 'Done!';
pause(0.1);

%% Printing output on GUI

if(exist('R2_quasi_static','var'))
    app.static_R2_value.Text = num2str(round(R2_quasi_static*100)/100);
else
    app.static_R2_value.Text = num2str(round(R2_static*100)/100);
end

if(viscousModel.activator == 1)
    app.dynamic_R2_value.Text = num2str(round(R2_dynamic*100)/100);
    app.memory_fun.Text = fit_data.modelFormulation.viscousModel.funName;
    app.N_param_dynamic.Text = num2str(fit_data.modelFormulation.viscousModel.nParam*2);
end

if(isfield(fit_data.fittedCurves,'quasiStaticFit'))
    fieldName = 'quasiStaticFit';
else
    fieldName = 'staticFit';
end

app.in_vivo_axial_stretch_mod.Text = num2str(round(fit_data.vesselConfigurations.reference_configuration.("Axial length [mm]")/...
    fit_data.vesselConfigurations.unloaded_configuration.("Axial length [mm]")*mean(fit_data.modBehaviour.PD_100.("Axial stretch [-]"))*100)/100);
app.in_vivo_axial_stretch_exp.Text = num2str(round(fit_data.vesselConfigurations.reference_configuration.("Axial length [mm]")/...
    fit_data.vesselConfigurations.unloaded_configuration.("Axial length [mm]")*mean(fit_data.fittedCurves.(fieldName).PD_100.("Axial stretch [-]"))*100)/100);

app.SEF.Text = fit_data.modelFormulation.elasticModel.modelName;
app.N_param_static.Text = num2str(length(table2array(fit_data.modelParameters.staticSEFparameters)));
app.N_param_constrained.Text = num2str(length(table2array(fit_data.modelParameters.staticSEFparameters))-length(fit_data.modelFormulation.elasticModel.lb));

if(fit_data.vesselConfigurations.reference_configuration.("Axial length [mm]")==fit_data.vesselConfigurations.unloaded_configuration.("Axial length [mm]"))
    app.ref_conf.Text = 'Unloaded vessel.';
else
    app.ref_conf.Text = 'in vivo (exp. axial stretch and 100 mmHg).';
end

%% Saving data
app.text_error_box.Value = '';
if(queuing == 0)
    while(strcmp(app.text_error_box.Value,''))
        app.error_box_tab.Visible = 'on';
        app.No_button_error.Visible = 'off';
        app.Yes_button_error.Visible = 'off';
        app.OK_button_error.Visible = 'on';
        
        app.No_button_error.Value = 0;
        app.Yes_button_error.Value = 0;
        app.OK_button_error.Value = 0;
        
        app.error_box_tab.Title = 'Please enter the desired file name:';
        
        while(app.OK_button_error.Value == 0)
            pause(0.1);
        end 
        app.error_box_tab.Visible = 'off';
    
        app.OK_button_error.Value = 0;
    end

    fit_data.file_name_exp_data = FileName;

    save(sprintf('fit_%s.mat',app.text_error_box.Value{1}), 'fit_data')

    if(viscousModel.activator == 1)
        app.status_txt.Text = 'Exporting pressure driven curves...';
        pause(0.1);
        
        fit_data = pressureDrivenFittedCurves(fit_data);
    end

    save(sprintf('fit_%s.mat',app.text_error_box.Value{1}), 'fit_data')
else
    fit_data.file_name_exp_data = FileName;

    save(sprintf('fit_%s.mat',data.file_list{i_sample}(1:end-4)), 'fit_data')

    if(viscousModel.activator == 1)
        fit_data = pressureDrivenFittedCurves(fit_data);
    end

    save(sprintf('fit_%s.mat',data.file_list{i_sample}(1:end-4)), 'fit_data')
end

end