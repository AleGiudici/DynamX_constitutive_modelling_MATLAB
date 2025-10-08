function [elasticModel,viscousModel,activeModel,in_vivo_ref,PD,axial_sweeps_check]=select_constitutive_model_GUI(passive,app)
% This function allows the user to customise the model. The function
% returns 4 output parameters:
%
% 1) elasticModel --> it's a structure containing the data/info necessessary to
%   formulate the quasi-static (elastic) passive strain energy function
%   (based on the user's choices). These include: the chosen strain energy 
%   function (.fun), the number of parameters of the model (.nParam_baseline), 
%   the number of modelled wall constituents (.nConstituents), the baseline 
%   normalisation vector to bound all parameters to the ~same range 
%   (.normV_baseline), a list of the parameter names (.paramName), the 
%   baseline upper and lower bound vector for the parameter search 
%   (.ub_baseline and .lb_baseline), the number of parameters included in 
%   the fitting routine (i.e., excluding fixed ones, .nParam), upper and 
%   lower bound vectors for the parameter search (i.e., excluding 
%   constrained parameters, .ub and .lb), and the final normalisation
%   vector (i.e., excluding constrained parameters, .normV).
%
%   The parameters of the elastic model can be constrained to user defined
%   values by typing these values in the parameter list.
%
% 2) viscousModel --> it's a structure containing the data/info necessessary to
%   formulate the dynamic (viscoelastic) model (based on the user's choices).
%   These include: an activator (if 0, the viscoelastic modelling is not
%   performed), the memory function to be used in Fung's model (.fun), the
%   number of viscous parameters (.nParam), and the normalisation vector to
%   bound all parameters to the same range (.normV).
%
% 3) activeModel --> it's a structure containing the data/info necessessary to
%   formulate the quasi-static (elastic) active strain energy function
%   (based on the user's choices). These include: the chosen active strain 
%   energy function (.fun), the number of parameters of the model 
%   (.nParam_baseline), the number of modelled active wall constituents 
%   (.nConstituents), the baseline normalisation vector to bound 
%   all parameters to the ~same range (.normV_baseline), a list of the
%   parameter names (.paramName), the baseline upper and lower bound vector 
%   for the parameter search (.ub_baseline and .lb_baseline).
% 
% 4) in_vivo_ref --> it's set to 1 if the user decides to fix the reference
%   configuration to the in vivo state (100 mmHg of pressure and experimentally 
%   estimated in vivo axial stretch). If set to 0, the reference is the
%   unloaded state.
%
% 5) PD --> list of quasi-static pressure sweep measurements to be used in
%   the parameter fitting.
%
% 6) axial_seeps_check --> determines whether axial sweep experiments are
%   included in the parameter estimation.
%
% Input parameters are the app data "app" and the structure with the
% experimental data of the passive experiments ("passive").

app.mod_settings_tab.Visible = 'on';

while(app.EnterButton.Value == 0)
    pause(0.1)
end

app.mod_settings_tab.Visible = 'off';
app.EnterButton.Value = 0;

elastic_model = app.SEF_menu.Value;
memory_function = app.memory_fun_menu.Value;
fitting_data_choice = app.input_data_menu.Value;
ref_conf = app.ref_conf_menu.Value;
active_model = app.active_model.Value;

% [elastic_model,memory_function,ref_conf,fitting_data_choice]=choose_constitutive_model;

list_models = ["4-fibre families", "2-fibre families", "2-fibre families with dispersion", "Zulliger's SEF","Zulliger's SEF (2-fibre)","4-fibre families (not neo-Hookean)","Bilayered model","UFD model (not neo-Hookean)"];
list_functions = ["None", "Fung", "Prony series"];
list_active_models = ["None", "Constant 2PK stress", "Constant 1PK stress", "Constant Cauchy stress", "Zulliger 2004","Rachev 1999","Franchini 2022","Gaussian model"];
list_reference = ["Unloaded vessel", "In vivo (100 mmHg and in vivo axial length)"];
list_fitting_data_options = ["All (pressure and axial force sweeps)", "Pressure sweeps", "Pressure sweep at in vivo length"];

if(strcmp(ref_conf,list_reference{1,1}))
    in_vivo_ref = false;
else
    in_vivo_ref = true;
end

if(strcmp(fitting_data_choice,list_fitting_data_options{1,1}))
    PD = passive.static_data.PD;
    axial_sweeps_check = 1;
elseif(strcmp(fitting_data_choice,list_fitting_data_options{1,2}))
    PD = passive.static_data.PD;
    axial_sweeps_check = 0;
else
    PD = {'PD_100'};
    axial_sweeps_check = 0;
end

if(strcmp(elastic_model,list_models{1,1}))
    elasticModel.modelName = '4-fiber families';
    elasticModel.fun = @FourFibreModel;
    elasticModel.nParam_baseline = 12;
    elasticModel.nConstituents = 2;
    elasticModel.normV_baseline = [0.2, 0.2, 100, 0.2, pi/2, 100, 0.2, 100, 50, 1, 1, 1];
    elasticModel.paramName = {'Elastin stiffness-like parameter [MPa]','Circumferential collagen stiffness-like parameter [MPa]',...
        'Circumferential collagen non-linearity parameter [-]','Diagonal collagen stiffness-like parameter [MPa]',...
        'Diagonal collagen fibre angle [deg]','Diagonal collagen non-linearity parameter [-]',...
        'Axial collagen stiffness-like parameter [MPa]','Axial collagen non-linearity parameter [-]',...
        'Collagen non-linearity parameter for compression [-]','Elastin circumferential deposition stretch [-]',...
        'Elastin axial deposition stretch [-]','Collagen deposition stretch [-]'};
    elasticModel.imposedParam = {'','','','','','','','','','','',''};
    if(in_vivo_ref)
        elasticModel.ub_baseline = [1,1,1,1,1,1,1,1,1,3,3,1.5];
        elasticModel.lb_baseline = [0,0,0,0,0,0,0,0,0,1,1,0.5];
    else
        elasticModel.ub_baseline = [1,1,1,1,1,1,1,1,1,2,2,1];
        elasticModel.lb_baseline = [0,0,0,0,0,0,0,0,0,1,1,0.5];
    end
elseif(strcmp(elastic_model,list_models{1,2}))
    elasticModel.modelName = '2-fibre families';
    elasticModel.fun = @TwoFibreModel;
    elasticModel.nParam_baseline = 8;
    elasticModel.normV_baseline = [1, 1, 20, pi/2, 20, 1, 1, 1];
    elasticModel.paramName = {'Elastin stiffness-like parameter [MPa]','Collagen stiffness-like parameter [MPa]',...
        'Collagen non-linearity parameter [-]','Collagen fibre angle [deg]'...
        'Collagen non-linearity parameter for compression [-]','Elastin circumferential deposition stretch [-]',...
        'Elastin axial deposition stretch [-]','Collagen deposition stretch [-]'};
    elasticModel.imposedParam = {'','','','','','','',''};
    elasticModel.nConstituents = 2;
    if(in_vivo_ref)
        elasticModel.ub_baseline = [1,1,1,1,1,3,2.5,2];
        elasticModel.lb_baseline = [0,0,0,0,0,1.5,1,1];
    else
        elasticModel.ub_baseline = [1,1,1,1,1,1,2,2,1];
        elasticModel.lb_baseline = [0,0,0,0,0,0,1,1,0.5];
    end
elseif(strcmp(elastic_model,list_models{1,3}))
    elasticModel.modelName = '2-fibre familis with fibre dispersion';
    elasticModel.fun = @TwoFibreModelDispersion;
    elasticModel.nParam_baseline = 9;
    elasticModel.normV_baseline = [1, 1, 150, pi/2, 1/3, 100, 1, 1, 1];
    elasticModel.paramName = {'Elastin stiffness-like parameter [MPa]','Collagen stiffness-like parameter [MPa]',...
        'Collagen non-linearity parameter [-]','Collagen fibre angle [deg]'...
        'Collagen fibre dispersion [-]','Collagen non-linearity parameter for compression [-]',...
        'Elastin circumferential deposition stretch [-]','Elastin axial deposition stretch [-]',...
        'Collagen deposition stretch [-]'};
    elasticModel.imposedParam = {'','','','','','','','',''};
    elasticModel.nConstituents = 2;
    if(in_vivo_ref)
        elasticModel.ub_baseline = [1,1,1,1,1,1,3,2.5,2];
        elasticModel.lb_baseline = [0,0,0,0,0,0,1.5,1,1];
    else
        elasticModel.ub_baseline = [1,1,1,1,1,1,2,2,1];
        elasticModel.lb_baseline = [0,0,0,0,0,0,1,1,0.5];
    end
elseif(strcmp(elastic_model,list_models{1,4}))
    elasticModel.modelName = 'Zulliger´s SEF';
    elasticModel.fun = @Zulliger;
    elasticModel.nParam_baseline = 8;
    elasticModel.nConstituents = 2;
    elasticModel.normV_baseline = [1, 200, 1, 2, 100, pi/2, 2, 2];
    elasticModel.paramName = {'Elastin stiffness-like parameter [MPa]','Collagen stiffness-like parameter [MPa]',...
        'Collagen siffness-like parameter for compression [MPa]','Collagen fibre engagment distribution scale parameter (b) [-]',...
        'Collagen fibre engagment distribution shape parameter (k) [-]','Diagonal collagen fibre angle [deg]',...
        'Elastin circumferential deposition stretch [-]','Elastin axial deposition stretch [-]'};
    elasticModel.imposedParam = {'','','','','','','',''};
    if(in_vivo_ref)
        errordlg('In-vivo reference configuration non implemented for this model','Error');
        [elasticModel,viscousModel,in_vivo_ref] = select_constitutive_model; % get the chosen constitutive function, number of parameters and the vector for their conversion to real values
    else
        elasticModel.ub_baseline = [1,0.4,1,1,1,1,1,1];
        elasticModel.lb_baseline = [0,0.01,0,0,0,0,0,0];
    end
elseif(strcmp(elastic_model,list_models{1,5}))
    elasticModel.modelName = '2-fibre familis Zulliger´s SEF';
    elasticModel.fun = @Zulliger2fibre;
    elasticModel.nParam_baseline = 12;
    elasticModel.nConstituents = 2;
    elasticModel.normV_baseline = [1, 200, 1, 2, 100, pi/2, 2, 2, 200, 1, 2, 100];
    elasticModel.paramName = {'Elastin stiffness-like parameter [MPa]','Collagen stiffness-like parameter [MPa]',...
        'Collagen siffness-like parameter for compression [MPa]','Collagen fibre engagment distribution scale parameter (b) [-]',...
        'Collagen fibre engagment distribution shape parameter (k) [-]','Diagonal collagen fibre angle [deg]',...
        'Elastin circumferential deposition stretch [-]','Elastin axial deposition stretch [-]',...
        'Collagen stiffness-like parameter [MPa]',...
        'Collagen siffness-like parameter for compression [MPa]','Collagen fibre engagment distribution scale parameter (b) [-]',...
        'Collagen fibre engagment distribution shape parameter (k) [-]'};
    elasticModel.imposedParam = {'','','','','','','','','','','',''};
    if(in_vivo_ref)
        errordlg('In-vivo reference configuration non implemented for this model','Error');
        [elasticModel,viscousModel,in_vivo_ref] = select_constitutive_model; % get the chosen constitutive function, number of parameters and the vector for their conversion to real values
    else
        elasticModel.ub_baseline = [1,0.4,1,1,1,1,1,1,0.4,1,1,1];
        elasticModel.lb_baseline = [0,0.01,0,0,0,0,0,0,0.01,0,0,0];
    end
elseif(strcmp(elastic_model,list_models{1,6}))
    elasticModel.modelName = '4-fiber families (not neo-Hookean)';
    elasticModel.fun = @FourFibreModelNNH;
    elasticModel.nParam_baseline = 12;
    elasticModel.nConstituents = 2;
    elasticModel.normV_baseline = [0.2, 0.2, 50, 0.2, pi/2, 50, 0.2, 50, 50, 1, 1, 1];
    elasticModel.paramName = {'Elastin stiffness-like parameter [MPa]','Circumferential collagen stiffness-like parameter [MPa]',...
        'Circumferential collagen non-linearity parameter [-]','Diagonal collagen stiffness-like parameter [MPa]',...
        'Diagonal collagen fibre angle [deg]','Diagonal collagen non-linearity parameter [-]',...
        'Axial collagen stiffness-like parameter [MPa]','Axial collagen non-linearity parameter [-]',...
        'Collagen non-linearity parameter for compression [-]','Elastin circumferential deposition stretch [-]',...
        'Elastin axial deposition stretch [-]','Collagen deposition stretch [-]'};
    elasticModel.imposedParam = {'','','','','','','','','','','',''};
    if(in_vivo_ref)
        elasticModel.ub_baseline = [1,1,1,1,1,1,1,1,0,3,3,1.8];
        elasticModel.lb_baseline = [0,0,0,0,0,0,0,0,0,1,1,0.5];
    else
        elasticModel.ub_baseline = [1,1,1,1,1,1,1,1,0,2,2,1];
        elasticModel.lb_baseline = [0,0,0,0,0,0,0,0,0,1,1,0.5];
    end
elseif(strcmp(elastic_model,list_models{1,7}))
    elasticModel.modelName = 'Bylayered model';
    elasticModel.fun = @biLayeredModel;
    elasticModel.nParam_baseline = 14;
    elasticModel.normV_baseline = [1, 1, 10, pi/2, 1/3, 1, 200, pi/4, 1/3, 100, 1, 1, 1, 1];
    elasticModel.paramName = {'Elastin stiffness-like parameter [MPa]','Medial collagen stiffness-like parameter [MPa]',...
        'Medial collagen non-linearity parameter [-]','Medial collagen fibre angle [deg]',...
        'Medial collagen fibre dispersion [-]','Adventitial collagen stiffness-like parameter [MPa]',...
        'Adventitial collagen non-linearity parameter multiplier [-]','Adventitial collagen fibre angle multiplier [-]',...
        'Adventitial collagen fibre dispersion [-]','Collagen non-linearity parameter for compression [-]',...
        'Elastin circumferential deposition stretch [-]','Elastin axial deposition stretch [-]',...
        'Medial collagen deposition stretch [-]', 'Adventitial collagen deposition stretch [-]'};
    elasticModel.imposedParam = {'','','','','','','','','','','','','',''};
    elasticModel.nConstituents = 2;
    if(in_vivo_ref)
        elasticModel.ub_baseline = [1,1,1,1,1,1,1,1.2,1,1,3,2.5,2,2];
        elasticModel.lb_baseline = [0,0,0,0.5,0,0,0,0.8,0,0,1.5,1,1,0.8];
    else
        elasticModel.ub_baseline = [1,1,1,1,1,1,1,1,1,1,2,2,1,1];
        elasticModel.lb_baseline = [0,0,0,0,0,0,0,0,0,0,1,1,0.5,0.5];
    end
elseif(strcmp(elastic_model,list_models{1,8}))
    elasticModel.modelName = 'UFD model (not Neo-Hookean)';
    elasticModel.fun = @UFDmodelNNH;
    elasticModel.nParam_baseline = 8;
    elasticModel.normV_baseline = [1, 1, 50, 1, 1, 1, 1, 1];
    elasticModel.paramName = {'Elastin stiffness-like parameter [MPa]','Collagen stiffness-like parameter [MPa]',...
        'Collagen non-linearity parameter [-]','Collagen fibre distribution [-]','Collagen non-linearity parameter for compression [-]',...
        'Elastin circumferential deposition stretch [-]','Elastin axial deposition stretch [-]',...
        'Collagen deposition stretch [-]'};
    elasticModel.imposedParam = {'','','','','','','',''};
    elasticModel.nConstituents = 2;
    if(in_vivo_ref)
        elasticModel.ub_baseline = [1,1,1,1,1,3,3,1.8];
        elasticModel.lb_baseline = [0,0,0,0,0,1,1,0.8];
    else
        elasticModel.ub_baseline = [1,1,1,1,1,2,2,1];
        elasticModel.lb_baseline = [0,0,0,0,0,1,0.5,0.5];
    end
end    

for i = 1:elasticModel.nParam_baseline
    parameterLabelName = ['Parameter' num2str(i) 'Label'];
    app.(parameterLabelName).Text = elasticModel.paramName{i};
    app.(parameterLabelName).Text = elasticModel.paramName{i};
    parameterName = ['param_' num2str(i)];
    app.(parameterName).Visible = "on";
%     eval(['app.Parameter' num2str(i) 'Label.Text = elasticModel.paramName{i};']);
%     eval(['app.Parameter' num2str(i) 'Label.Visible = "on";']);
%     eval(['app.param_' num2str(i) '.Visible = "on";']);
end
for i = elasticModel.nParam_baseline+1:14
    parameterLabelName = ['Parameter' num2str(i) 'Label'];
    app.(parameterLabelName).Visible = "off";
    parameterName = ['param_' num2str(i)];
    app.(parameterName).Visible = "off";
%     eval(['app.Parameter' num2str(i) 'Label.Visible = "off";']);
%     eval(['app.param_' num2str(i) '.Visible = "off";']);
end

app.model_param_tab.Visible = 'on';

while(app.OKButton.Value == 0)
    pause(0.1)
end

app.OKButton.Value = 0;
app.model_param_tab.Visible = 'off';

% elastic=imposeParam(elastic);
count=0;
for i=1:elasticModel.nParam_baseline
    parameterName = ['param_' num2str(i)];
    count = count+isempty(app.(parameterName).Value);
%     count = count+isempty(eval(['app.param_' num2str(i) '.Value']));
end

elasticModel.nParam=count;
elasticModel.lb=zeros(1,count);
elasticModel.ub=zeros(1,count);
elasticModel.normV=zeros(1,count);
elasticModel.imposedParam=zeros(1,count);

j=0;
for i=1:elasticModel.nParam_baseline
    parameterName = ['param_' num2str(i)];
    if(isempty(app.(parameterName).Value))
        j=j+1;
        elasticModel.lb(j)=elasticModel.lb_baseline(i);
        elasticModel.ub(j)=elasticModel.ub_baseline(i);
        elasticModel.normV(j)=elasticModel.normV_baseline(i);
        elasticModel.imposedParam(i) = NaN;
    else
        elasticModel.imposedParam(i) = str2double(app.(parameterName).Value);
    end
end

if(strcmp(active_model,list_active_models{1,1}))
    activeModel.activator = false;
elseif(strcmp(active_model,list_active_models{1,2}))
    activeModel.activator = true;
    activeModel.modelName = 'Constant 2PK stress';
    activeModel.fun = @Constant2PKmodel;
    activeModel.nParam_baseline = 2;
    activeModel.nConstituents = 1;
    activeModel.normV_baseline = [0.2, 1];
    activeModel.paramName = {'VSMC stiffness-like parameter [MPa]','VSMC circ.deposition stretch [-]'};

    if(in_vivo_ref)
        activeModel.ub_baseline = [1,1.8];
        activeModel.lb_baseline = [0,0.5];
    else
        activeModel.ub_baseline = [1,1];
        activeModel.lb_baseline = [0,0.5];
    end
elseif(strcmp(active_model,list_active_models{1,5}))
    activeModel.activator = true;
    activeModel.modelName = 'Zulliger 2004';
    activeModel.fun = @zulligerSMC;
    activeModel.nParam_baseline = 2;
    activeModel.nConstituents = 1;
    activeModel.normV_baseline = [0.2, 1];
    activeModel.paramName = {'VSMC stiffness-like parameter [MPa]','VSMC circ.deposition stretch [-]'};

    if(in_vivo_ref)
        activeModel.ub_baseline = [1,1.8];
        activeModel.lb_baseline = [0,0.5];
    else
        activeModel.ub_baseline = [1,1];
        activeModel.lb_baseline = [0,0.5];
    end
elseif(strcmp(active_model,list_active_models{1,6}))
    activeModel.activator = true;
    activeModel.modelName = 'Rachev 1999';
    activeModel.fun = @rachevSMC;
    activeModel.nParam_baseline = 3;
    activeModel.nConstituents = 1;
    activeModel.normV_baseline = [0.4,1,1];
    activeModel.paramName = {'VSMC stiffness-like parameter [MPa]','VSMC optimal circ. stretch [-]', 'VSMC circ. stretch width [-]'};
    
    if(in_vivo_ref)
        activeModel.ub_baseline = [1,2,0.8];
        activeModel.lb_baseline = [0,0.8,0];
    else
        activeModel.ub_baseline = [1,1.5,0.8];
        activeModel.lb_baseline = [0,0,0];
    end
elseif(strcmp(active_model,list_active_models{1,7}))
    activeModel.activator = true;
    activeModel.modelName = 'Franchini 2022';
    activeModel.fun = @franchiniSMC;
    activeModel.nParam_baseline = 6;
    activeModel.nConstituents = 1;
    activeModel.normV_baseline = [0.1,1/2,10,10,2,2];
    activeModel.paramName = {'VSMC stiffness-like parameter [MPa]','VSMC dispersion coefficient [-]', 'VSMC shape parameter 1 [-]', 'VSMC shape parameter 2 [-]', 'VSMC circ. deposition stretch [-]','VSMC axial deposition stretch [-]'};
    if(in_vivo_ref)
        activeModel.ub_baseline = [0.1,1,1,0,2,2];
        activeModel.lb_baseline = [0.1,0,0,-1,1,1];
    else
        activeModel.ub_baseline = [0.1,1,1,0,1,1];
        activeModel.lb_baseline = [0.1,0,0,-1,0.5,0.5];
    end
elseif(strcmp(active_model,list_active_models{1,8}))
    activeModel.activator = true;
    activeModel.modelName = 'Gaussian model';
    activeModel.fun = @AleSMC;
    activeModel.nParam_baseline = 5;
    activeModel.nConstituents = 1;
    activeModel.normV_baseline = [2, 1/2, 2, 1, 1];
    activeModel.paramName = {'VSMC peak stress [MPa]','VSMC dispersion coefficient [-]','Gaussian width [-]','Optimal circ. stretch [-]','Optimal axial stretch [-]'};
    activeModel.ub_baseline = [1,1,1,1,1];
    activeModel.lb_baseline = [0,0,0,0,1];
end

if(strcmp(memory_function,list_functions{1,1}))
    viscousModel.activator=false;
elseif(strcmp(memory_function,list_functions{1,2}))
    viscousModel.funName = 'Fung´s relaxation function';
    viscousModel.fun=@FungMemoryFun;
    viscousModel.nParam=3;
    viscousModel.normV=[0.2 0.001 100];
    viscousModel.lb = [0 0.001 10]./viscousModel.normV;
    viscousModel.ub = [0.2 0.001 100]./viscousModel.normV;
    viscousModel.activator=true;
    
elseif(strcmp(viscous_model,list_functions{1,2}))

end  



%% OLD ATTEMPT - 4-fibre model with power 4 stretch dependency (did't improve fitting)
% elseif(strcmp(elastic_model,list_models{1,5}))
%     elastic.modelName='FourFiberModel4';
%     elastic.fun=@FourFibreModel4;
%     elastic.nParam_baseline=13;
%     elastic.nConstituents=2;
%     elastic.normV_baseline=[1, 1, 50, 1, pi/2, 50, 1, 50, 50, 1, 1, 1, 1];
%     elastic.paramName={'Elastin stiffness-like parameter [MPa]','Circumferential collagen stiffness-like parameter [MPa]',...
%         'Circumferential collagen non-linearity parameter [-]','Diagonal collagen stiffness-like parameter [MPa]',...
%         'Diagonal collagen fibre angle [deg]','Diagonal collagen non-linearity parameter [-]',...
%         'Axial collagen stiffness-like parameter [MPa]','Axial collagen non-linearity parameter [-]',...
%         'Collagen non-linearity parameter for compression [-]','Elastin circumferential deposition stretch [-]',...
%         'Elastin axial deposition stretch [-]','Collagen deposition stretch [-]','Circumferential collagen deposition stretch [-]'};
%     elastic.imposedParam = {'','','','','','','','','','','','',''};
%     if(in_vivo_ref)
%         elastic.ub_baseline=[1,1,1,1,1,1,1,1,1,10,10,2,2];
%         elastic.lb_baseline=[0,0,0,0,0,0,0,0,0,1,1,0.8,0.8];
%     else
%         elastic.ub_baseline=[1,1,1,1,1,1,1,1,1,2,2,1,1];
%         elastic.lb_baseline=[0,0,0,0,0,0,0,0,0,1,1,0.5,0.5];
%     end
