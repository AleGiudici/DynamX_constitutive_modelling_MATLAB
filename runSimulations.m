clear all;
close all;
clc;

mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
gtN = 0.0098; %scaling factor to convert g to N

file_name = 'fit_Young_5_corrected_VSMC_model';

load(file_name);

parV_elastic = fit_data.modelParameters.staticSEFparameters.("Final QLV values");
parV_viscous = fit_data.modelParameters.viscousParameters.Values;
parV_active = fit_data.modelParameters.activeSEFparameters.("Final values");

if(fit_data.modelFormulation.viscousModel.activator == 1)
    if(strcmp(fit_data.modelFormulation.elasticModel.modelName,'4-fiber families 2') || strcmp(fit_data.modelFormulation.elasticModel.modelName,'4-fiber familiy model') || strcmp(fit_data.modelFormulation.elasticModel.modelName,'4-fiber families (not neo-Hookean)'))
        parV_elastic(1) = parV_elastic(1)/(1+parV_viscous(1,1)*log(parV_viscous(3,1)/parV_viscous(2,1)));
        parV_elastic(2) = parV_elastic(2)/(1+parV_viscous(4,1)*log(parV_viscous(6,1)/parV_viscous(5,1)));
        parV_elastic(4) = parV_elastic(4)/(1+parV_viscous(4,1)*log(parV_viscous(6,1)/parV_viscous(5,1)));
        parV_elastic(7) = parV_elastic(7)/(1+parV_viscous(4,1)*log(parV_viscous(6,1)/parV_viscous(5,1)));

    elseif(strcmp(fit_data.modelFormulation.elasticModel.modelName,'2-fibre family model') || strcmp(fit_data.modelFormulation.elasticModel.modelName,'2-fibre family model with dispersion'))
        parV_elastic(1) = parV_elastic(1)/(1+parV_viscous(1,1)*log(parV_viscous(3,1)/parV_viscous(2,1)));
        parV_elastic(2) = parV_elastic(2)/(1+parV_viscous(4,1)*log(parV_viscous(6,1)/parV_viscous(5,1)));

    elseif(strcmp(fit_data.modelFormulation.elasticModel.modelName,'Zulliger´s model'))
        parV_elastic(1) = parV_elastic(1)/(1+parV_viscous(1,1)*log(parV_viscous(3,1)/parV_viscous(2,1)));
        parV_elastic(2) = parV_elastic(2)/(1+parV_viscous(4,1)*log(parV_viscous(6,1)/parV_viscous(5,1)));
        parV_elastic(3) = parV_elastic(3)/(1+parV_viscous(4,1)*log(parV_viscous(6,1)/parV_viscous(5,1)));
    end
    parV_active(1) = parV_active(1)/(1+parV_viscous(7,1)*log(parV_viscous(9,1)/parV_viscous(8,1)));
end


elasticModel = fit_data.modelFormulation.elasticModel;
viscousModel = fit_data.modelFormulation.viscousModel;
activeModel = fit_data.modelFormulation.activeModel;

Ro = fit_data.vesselConfigurations.unloaded_configuration.("Outer radius [mm]");
Ri = fit_data.vesselConfigurations.unloaded_configuration.("Inner radius [mm]");
L = fit_data.vesselConfigurations.unloaded_configuration.("Axial length [mm]");

ro_ref = fit_data.vesselConfigurations.reference_configuration.("Outer radius [mm]");
ri_ref = fit_data.vesselConfigurations.reference_configuration.("Inner radius [mm]");
l_ref = fit_data.vesselConfigurations.reference_configuration.("Axial length [mm]");

ro_ref_active = fit_data.vesselConfigurations.reference_configuration_active.("Outer radius [mm]");
ri_ref_active = fit_data.vesselConfigurations.reference_configuration_active.("Inner radius [mm]");
l_ref_active = fit_data.vesselConfigurations.reference_configuration_active.("Axial length [mm]");

%% Quasi-static pressure sweep (uniaxial) - passive

target_pressure = 0:10:200;

lambdat_0 = (Ro+Ri)/(ro_ref+ri_ref);
lambdaz_0 = L/l_ref;

lambdazV = zeros(1,length(target_pressure));
lambdatV = zeros(1,length(target_pressure));

x0 = [1 1];
lb = [0.5 0.5];
ub = [1.5 1.5];

options = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 2000,'OptimalityTolerance',10^-6,'FunctionTolerance',10^-6); % setting max. no. of iterations to 8000

for iTime = 1:length(lambdatV)
    newLambdas = lsqnonlin(@(newLambdas) findUniaxialLambdazLambdatElastic(newLambdas,parV_elastic,parV_active,elasticModel,activeModel,...
        target_pressure(iTime),ro_ref-ri_ref,(ri_ref+ro_ref)/2,0),x0,lb,ub,options);
    lambdatV(iTime) = newLambdas(1);
    lambdazV(iTime) = newLambdas(2);
end

rm_ref = (ro_ref+ri_ref)/2;
h_ref = ro_ref-ri_ref;

rm = rm_ref*lambdatV;
ri = sqrt(rm.^2-(rm_ref^2-(rm_ref-h_ref/2)^2)./lambdazV);
ro = sqrt(rm.^2+((rm_ref+h_ref/2)^2-rm_ref^2)./lambdazV);
h = ro-ri;

Di = ri*2;
Do = ro*2;

[sigmatt,sigmazz] = elasticModel.fun(parV_elastic,lambdatV,lambdazV,1./lambdatV./lambdazV);
P = sigmatt.*h./rm/mtN;
fT = sigmazz.*(pi*h.*(2*rm))/gtN;
t = nan(1,length(h));

[sigmatt_elastin,sigmazz_elastin] = elasticModel.fun(parV_elastic,lambdatV,lambdazV,1./lambdatV./lambdazV,[1,0]);
circ_ratio = sigmatt_elastin./sigmatt*100;
axial_ratio = nan(1,length(circ_ratio));

dataMat = [t' P' Di' Do' fT' lambdatV' sigmatt' lambdazV' sigmazz' circ_ratio' axial_ratio'];

variable_names = ["Time [s]", "Pressure [mmHg]", "Inner diameter [mm]", "Outer diameter [mm]", "Transducer force [g]",...
                    "Circ. stretch [-]", "Circ. stress [MPa]", "Axial stretch [-]","Axial stress [MPa]",...
                    "Elastin circ. load-bearing [%]","Elastin axial load-bearing [%]"];

fit_data.user_simulation.UniaxialPsweep = array2table(dataMat,"VariableNames",variable_names);

lambdaz_uc = lambdazV'/lambdazV(6);


%% Sinusoidal stretch (uniaxial) - passive
frequency = 10;
t = 0:0.005:1;

variableNames = ["Time [s]", "Pressure [mmHg]", "Inner diameter [mm]", "Outer diameter [mm]", "Transducer force [g]",...
                        "Circ. stretch [-]", "Circ. stress [MPa]", "Axial stretch [-]","Axial stress [MPa]","Elastin circ. load bearing [%]",...
                        "Elastin axial load bearing [%]","Collagen circ. load bearing [%]","Collagen axial load bearing [%]",...
                        "VSMC circ. load bearing [%]", "VSMC axial load bearing [%]"];

for P_centre_loop = [70 90 110 130 150 170]
    target_pressure = P_centre_loop+20*sin(2*pi*frequency*t);
    
    tV = [0, 1000+t];
    lambdatV = zeros(length(tV),1);
    lambdazV = zeros(length(tV),1);
    lambdatV(1) = lambdat_0;
    lambdazV(1) = lambdaz_0;
    
    xV_elastic_ve = fit_data.modelParameters.staticSEFparameters.("Final QLV values");
    
    for iTime = 2:length(lambdatV)
        newLambdas = lsqnonlin(@(newLambdas) findLambdazLambdatUniaxialViscoelastic(newLambdas,lambdatV(1:iTime-1),lambdazV(1:iTime-1),tV(1:iTime),...
            xV_elastic_ve,parV_viscous,parV_active,elasticModel,viscousModel,activeModel,target_pressure(iTime-1),h_ref,rm_ref,0),x0,lb,ub,options);
        lambdatV(iTime) = newLambdas(1);
        lambdazV(iTime) = newLambdas(2);
    end
   
    rm = rm_ref*lambdatV;
    ri = sqrt(rm.^2-(rm_ref^2-(rm_ref-h_ref/2)^2)./lambdazV);
    ro = sqrt(rm.^2+((rm_ref+h_ref/2)^2-rm_ref^2)./lambdazV);
    h = ro-ri;
    
    [exportMat] = simulateUniaxialViscoelasticBehaviour(tV,lambdatV,lambdazV,fit_data.modelFormulation,fit_data.modelParameters,...
            fit_data.vesselConfigurations.reference_configuration,0);%(xV_elastic_ve,parV_viscous,tV,lambdatV,lambdazV,ri,h,elasticModel,viscousModel);   
    exportMat.("Transducer force [g]") = zeros(21,1);
    test_name = ['UniaxialSineWave_P' num2str(P_centre_loop-20) '_' num2str(P_centre_loop+20) '_f10000_SNP'];

    fit_data.user_simulation.(test_name) = exportMat;
end


%% Sinusoidal stretch (pressurisation at in vivo axial stretch) - passive

frequency = 10;
t = 0:0.005:1;
triggerV = zeros(length(t),1);
triggerV(end-19) = 1;

for P_centre_loop = [70 90 110 130 150 170]
    target_pressure = P_centre_loop+20*sin(2*pi*frequency*t);
   
    [exportMat] = runSimulation_script(fit_data,target_pressure,t,fit_data.in_vivo_axial_stretch,0,triggerV);
    
    test_name = ['SineWave_P' num2str(P_centre_loop-20) '_' num2str(P_centre_loop+20) '_f10000_SNP'];

    fit_data.user_simulation.(test_name) = array2table(exportMat(:,:),"VariableNames",variableNames);
end


%% Quasi-static pressure sweep (uniaxial) - active

target_pressure = 0:10:200;

lambdat_0 = (Ro+Ri)/(ro_ref+ri_ref);
lambdaz_0 = L/l_ref;

lambdazV = zeros(1,length(target_pressure));
lambdatV = zeros(1,length(target_pressure));

x0 = [1 1];
lb = [0.3 0.3];
ub = [1.5 1.5];

options = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 8000,'OptimalityTolerance',10^-6,'FunctionTolerance',10^-6); % setting max. no. of iterations to 8000

for iTime = 1:length(lambdatV)
    lambda_values = lsqnonlin(@(lambda_values) findUniaxialLambdazLambdatElastic(lambda_values,parV_elastic,parV_active,elasticModel,...
        activeModel,target_pressure(iTime),ro_ref-ri_ref,(ri_ref+ro_ref)/2,1),x0,lb,ub,options);
    lambdatV(iTime) = lambda_values(1);
    lambdazV(iTime) = lambda_values(2);
end

rm_ref = (ro_ref+ri_ref)/2;
h_ref = ro_ref-ri_ref;

rm = rm_ref*lambdatV;
ri = sqrt(rm.^2-(rm_ref^2-(rm_ref-h_ref/2)^2)./lambdazV);
ro = sqrt(rm.^2+((rm_ref+h_ref/2)^2-rm_ref^2)./lambdazV);
h = ro-ri;

Di = ri*2;
Do = ro*2;

[sigmatt_pas,sigmazz_pas] = elasticModel.fun(parV_elastic,lambdatV,lambdazV,1./lambdatV./lambdazV);
[sigmatt_act,sigmazz_act] = activeModel.fun(parV_active,lambdatV,lambdazV,1./lambdatV./lambdazV);

sigmatt = sigmatt_pas+sigmatt_act;
sigmazz = sigmazz_pas+sigmazz_act;

P = sigmatt.*h./rm/mtN;
fT = sigmazz.*(pi*h.*(2*rm))/gtN;
t = nan(1,length(h));

circ_ratio = sigmatt_act./sigmatt*100;
axial_ratio = nan(1,length(circ_ratio));

dataMat = [t' P' Di' Do' fT' lambdatV' sigmatt' lambdazV' sigmazz' circ_ratio' axial_ratio'];

variable_names = ["Time [s]", "Pressure [mmHg]", "Inner diameter [mm]", "Outer diameter [mm]", "Transducer force [g]",...
                    "Circ. stretch [-]", "Circ. stress [MPa]", "Axial stretch [-]","Axial stress [MPa]",...
                    "VSMC circ. load-bearing [%]","VSMC axial load-bearing [%]"];

fit_data.user_simulation.UniaxialPsweepLname = array2table(dataMat,"VariableNames",variable_names);


%% Quasi-static pressure sweep (inflation at in vivo axial stretch) - active

lambdatV = 0.4:0.0001:1.1;
lambdazV = fit_data.in_vivo_axial_stretch;

[sigmatt_rr_pas,sigmazz_rr_pas] = elasticModel.fun(parV_elastic,lambdatV,lambdazV,1./lambdatV./lambdazV);
[sigmatt_rr_act,sigmazz_rr_act] = activeModel.fun(parV_active,lambdatV,lambdazV,1./lambdatV./lambdazV);

sigmatt_rr = sigmatt_rr_pas+sigmatt_rr_act;
sigmazz_rr = sigmazz_rr_pas+sigmazz_rr_act;

rm_ref = (ro_ref+ri_ref)/2;
h_ref = ro_ref-ri_ref;

rm = rm_ref*lambdatV;
ri = sqrt(rm.^2-(rm_ref^2-(rm_ref-h_ref/2)^2)./lambdazV);
ro = sqrt(rm.^2+((rm_ref+h_ref/2)^2-rm_ref^2)./lambdazV);
h = ro-ri;

[PV,~] = stress2pressure_force(sigmatt_rr,sigmazz_rr,rm,h);

lambdatV_res = interp1(PV,lambdatV,0:10:200);

[sigmatt_rr_pas,sigmazz_rr_pas] = elasticModel.fun(parV_elastic,lambdatV_res,lambdazV,1./lambdatV_res./lambdazV);
[sigmatt_rr_act,sigmazz_rr_act] = activeModel.fun(parV_active,lambdatV_res,lambdazV,1./lambdatV_res./lambdazV);

sigmatt_rr = sigmatt_rr_pas+sigmatt_rr_act;
sigmazz_rr = sigmazz_rr_pas+sigmazz_rr_act;

rm_ref = (ro_ref+ri_ref)/2;
h_ref = ro_ref-ri_ref;

rm = rm_ref*lambdatV_res;
ri = sqrt(rm.^2-(rm_ref^2-(rm_ref-h_ref/2)^2)./lambdazV);
ro = sqrt(rm.^2+((rm_ref+h_ref/2)^2-rm_ref^2)./lambdazV);
h = ro-ri;

[PV,fTV] = stress2pressure_force(sigmatt_rr,sigmazz_rr,rm,h);

sigmatt = sigmatt_rr-PV/2*mtN;
sigmazz = sigmazz_rr-PV/2*mtN;

circ_ratio = sigmatt_rr_act./sigmatt_rr*100;
axial_ratio = sigmazz_rr_act./sigmazz_rr*100;

t = nan(1,length(h));

dataMat = [t' PV' 2*ri' 2*ro' fTV' lambdatV_res' sigmatt' lambdazV*ones(1,length(h))' sigmazz' circ_ratio' axial_ratio'];

variable_names = ["Time [s]", "Pressure [mmHg]", "Inner diameter [mm]", "Outer diameter [mm]", "Transducer force [g]",...
                    "Circ. stretch [-]", "Circ. stress [MPa]", "Axial stretch [-]","Axial stress [MPa]",...
                    "VSMC circ. load-bearing [%]","VSMC axial load-bearing [%]"];

fit_data.user_simulation.PsweepLname = array2table(dataMat,"VariableNames",variable_names);

%% Sinusoidal stretch (uniaxial) - active
frequency = 10;
t = 0:0.005:1;

variableNames = ["Time [s]", "Pressure [mmHg]", "Inner diameter [mm]", "Outer diameter [mm]", "Transducer force [g]",...
                        "Circ. stretch [-]", "Circ. stress [MPa]", "Axial stretch [-]","Axial stress [MPa]","Elastin circ. load bearing [%]",...
                        "Elastin axial load bearing [%]","Collagen circ. load bearing [%]","Collagen axial load bearing [%]",...
                        "VSMC circ. load bearing [%]", "VSMC axial load bearing [%]"];

for P_centre_loop = [70 90 110 130 150 170]
    target_pressure = P_centre_loop+20*sin(2*pi*frequency*t);
    
    tV = [0, 5000+t];
    lambdatV = zeros(length(tV),1);
    lambdazV = zeros(length(tV),1);
    lambdatV(1) = lambdat_0;
    lambdazV(1) = lambdaz_0;
    
    xV_elastic_ve = fit_data.modelParameters.staticSEFparameters.("Final QLV values");
    parV_active = fit_data.modelParameters.activeSEFparameters.("Final values");
    
    for iTime = 2:length(lambdatV)
        newLambdas = lsqnonlin(@(newLambdas) findLambdazLambdatUniaxialViscoelastic(newLambdas,lambdatV(1:iTime-1),lambdazV(1:iTime-1),tV(1:iTime),...
            xV_elastic_ve,parV_viscous,parV_active,elasticModel,viscousModel,activeModel,target_pressure(iTime-1),h_ref,rm_ref,1),x0,lb,ub,options);
        lambdatV(iTime) = newLambdas(1);
        lambdazV(iTime) = newLambdas(2);
    end
   
    rm = rm_ref*lambdatV;
    ri = sqrt(rm.^2-(rm_ref^2-(rm_ref-h_ref/2)^2)./lambdazV);
    ro = sqrt(rm.^2+((rm_ref+h_ref/2)^2-rm_ref^2)./lambdazV);
    h = ro-ri;
    
    [exportMat] = simulateUniaxialViscoelasticBehaviour(tV,lambdatV,lambdazV,fit_data.modelFormulation,fit_data.modelParameters,...
        fit_data.vesselConfigurations.reference_configuration,1);%(xV_elastic_ve,parV_viscous,tV,lambdatV,lambdazV,ri,h,elasticModel,viscousModel);
    exportMat.("Transducer force [g]") = zeros(21,1);
    test_name = ['UniaxialSineWave_P' num2str(P_centre_loop-20) '_' num2str(P_centre_loop+20) '_f10000_Lname'];

    fit_data.user_simulation.(test_name) = exportMat;
end


%% Sinusoidal stretch (pressurisation at in vivo axial stretch) - active

frequency = 10;
t = 0:0.005:1;
triggerV = zeros(length(t),1);
triggerV(end-19) = 1;

for P_centre_loop = [70 90 110 130 150 170]
    target_pressure = P_centre_loop+20*sin(2*pi*frequency*t);
   
    [exportMat] = runSimulation_script(fit_data,target_pressure,t,fit_data.in_vivo_axial_stretch,1,triggerV);
    
    test_name = ['SineWave_P' num2str(P_centre_loop-20) '_' num2str(P_centre_loop+20) '_f10000_Lname'];

    fit_data.user_simulation.(test_name) = array2table(exportMat(:,:),"VariableNames",variableNames);
end

save(file_name,'fit_data')