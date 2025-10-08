function [costV] = objectiveActive(parameterV_normalised,inputDataMatrixActive,elasticModel,viscousModel,activeModel,n_dataPoints,weight_vector)
% Objective function for the optimisation of the model parameters of a
% visceoelastic VSMC model. The error is calculated using the
% normalised differences between modelled and experimental pressure and
% modelled and experimental transducer axial force.
% 
% Input parameters are normalised parameter vector (parameterV_normalised), 
% the matrix of experimental data (inputDataMatrixActive), the passive 
% model formulation (elasticModel), the viscous model formulation
% elastic (viscousModel), the active VSMC contraction formulation (activeModel),
% the cumulative number of datapoints in all experiments (n_dataPoints),
% and a vecotor of wieghts that allows to weigh differently different
% experiments. Note that parameterV_norm contains all the model parameters
% (i.e., the concatenation of the passive, active and viscous model
% parameters).

% De-normalise elastic parameter vector and add imposed parameters
    parameterV_temp=parameterV_normalised(1:elasticModel.nParam).*elasticModel.normV; 
    
    j=0;
    parameterV_elastic=zeros(1,elasticModel.nParam_baseline);
    for i=1:elasticModel.nParam_baseline
        if(isnan(elasticModel.imposedParam(i)))
            j=j+1;
            parameterV_elastic(i)=parameterV_temp(j);
        elseif(isinf(elasticModel.imposedParam(i)))
            if(strcmp(elasticModel.modelName,'4-fiber families (not neo-Hookean)') || strcmp(elasticModel.modelName,'4-fiber families'))
                if(i==4 || i==7)
                    parameterV_elastic(i)=parameterV_temp(2);
                elseif(i==6 || i==8)
                    parameterV_elastic(i)=parameterV_temp(3);
                end
            end
        else
            parameterV_elastic(i)=elasticModel.imposedParam(i);
        end
    end

% De-normalise active parameter vector
    parameterV_active = parameterV_normalised(elasticModel.nParam+1:elasticModel.nParam+activeModel.nParam_baseline).*activeModel.normV_baseline;

% De-normalise viscous parameter vector
    visco_normV = zeros(1,(elasticModel.nConstituents+activeModel.nConstituents)*3);
    for i = 1:elasticModel.nConstituents+activeModel.nConstituents
        visco_normV(1+(i-1)*3:i*3) = viscousModel.normV;
    end
    parameterV_viscous = parameterV_normalised(elasticModel.nParam+activeModel.nParam_baseline+1:end).*visco_normV;

% initialise vectors
    sigmatt_ve_compV_final = zeros(n_dataPoints,1);
    sigmazz_ve_compV_final = zeros(n_dataPoints,1);
    sigmarr_ve_compV_final = zeros(n_dataPoints,1);

    sigmatt_rr_ve_compV_final = zeros(n_dataPoints,1);
    sigmazz_rr_ve_compV_final = zeros(n_dataPoints,1);

    P_compV_final = zeros(n_dataPoints,1);
    fT_compV_final = zeros(n_dataPoints,1);

    P_expV_final = zeros(n_dataPoints,1);
    fT_expV_final = zeros(n_dataPoints,1);
    lambdatV_final = zeros(n_dataPoints,1);
    lambdazV_final = zeros(n_dataPoints,1);
    weightttV_final = zeros(n_dataPoints,1);
    weightzzV_final = zeros(n_dataPoints,1);

    ratio_exp = zeros(7,1);
    ratio_comp = zeros(7,1);
    err_slope = zeros(7,1);

    % calculate modelled behaviour for each test
    counter = 0;
    counter_dynamic_loop = 0;
    
    contraction_states = fieldnames(inputDataMatrixActive);
    
    for iState = 2:length(contraction_states)
        test_list = fieldnames(inputDataMatrixActive.(contraction_states{iState}));
    
        for i = 1:length(test_list) % number of test
            num_underscore = length(find(test_list{i}=='_'));
            
            processed_data = inputDataMatrixActive.(contraction_states{iState}).(test_list{i}).processed_data;
            % Importing input data for current test
            trigger_test = processed_data(:,13);
        
            tV = processed_data(:,1)'; % time vector
            lambdatV = processed_data(:,8); % circumferential stretch
            lambdazV = processed_data(:,9); % axial stretch
            lambdarV = processed_data(:,10); % radial stretch
        
            riV = processed_data(:,4); % inner radius
            hV = processed_data(:,5); % wall thickness
        
            sigmatt_expV = processed_data(:,11); % experimental circumferential Cauchy stress
        
            PV = processed_data(:,2); % experimental pressure
            fV = processed_data(:,6); % experimental force
            
            % create a time interval matrix for relaxation function
            tV_inverted=ones(1,length(tV))*max(tV)-tV;
            tM_temp=ones(length(tV),length(tV)).*tV_inverted-tV_inverted';
            tM = (tM_temp(2:end,1:end-1)+tM_temp(1:end-1,1:end-1))/2;
            locator=tM<0;
            tM(locator) = inf;
            
            % initialise computational Cauchy extra stress vectors
            sigmatt_ve_compV=zeros(length(tV),1); % circumferential 
            sigmazz_ve_compV=zeros(length(tV),1); % axial
            sigmarr_ve_compV=zeros(length(tV),1); % radial
            
            % calculate and superimpose the viscoelastic response of all constituents
            for iConstituent=1:elasticModel.nConstituents+activeModel.nConstituents*(iState-1)
                if(iConstituent<=elasticModel.nConstituents)
                    activator=zeros(1,2);
                    activator(iConstituent)=1;
            
                    % Elastic 2ndPK extra stress in the three principal directions
                    [sigmatt_rr_eV,sigmazz_rr_eV,~,~,~]=elasticModel.fun(parameterV_elastic,lambdatV,lambdazV,lambdarV,activator); 
                else
                    [sigmatt_rr_eV,sigmazz_rr_eV,~,~,~]=activeModel.fun(parameterV_active,lambdatV,lambdazV,lambdarV);
                end
                
                P_eV = sigmatt_rr_eV.*hV./(riV+hV/2);
    
                sigmatt_eV = sigmatt_rr_eV-P_eV/2; % Lagrange multiplier to get actual circ. stress
                sigmazz_eV = sigmazz_rr_eV-P_eV/2; % Lagrange multiplier to get actual axial stress
                sigmarr_eV = -P_eV/2; % Lagrange multiplier to get actual radial stress
        
                Stt_eV = sigmatt_eV./lambdatV.^2;
                Srr_eV = sigmarr_eV./lambdarV.^2;
        
                delta_Stt_eV = Stt_eV(2:end)-Stt_eV(1:end-1); % Circ. elastic 2ndPK extra stress increments
                delta_Srr_eV = Srr_eV(2:end)-Srr_eV(1:end-1); % Radial elastic 2ndPK extra stress increments
                
                % Relaxation function calculated at each instant of the time interval matrix
                G_fun = (1+parameterV_viscous((iConstituent-1)*viscousModel.nParam+1)*(expint(tM/parameterV_viscous((iConstituent-1)*viscousModel.nParam+3))-...
                    expint(tM/parameterV_viscous((iConstituent-1)*viscousModel.nParam+2))))/(1+parameterV_viscous((iConstituent-1)*viscousModel.nParam+1)*...
                    log(parameterV_viscous((iConstituent-1)*viscousModel.nParam+3)/parameterV_viscous((iConstituent-1)*viscousModel.nParam+2)));
                locator = isnan(G_fun);
                G_fun(locator) = 0;
                
                % Viscoelastic 2ndPK circumferential extra stress
                Stt_veV = G_fun*delta_Stt_eV+Stt_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
                
                % Viscoelastic 2ndPK radial extra stress
                Srr_veV = G_fun*delta_Srr_eV+Srr_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
                
                % Viscoelastic Cauchy circumferential extra stress
                sigmatt_ve_compV(1) = sigmatt_ve_compV(1)+Stt_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdatV(1)^2;
                sigmarr_ve_compV(1) = sigmarr_ve_compV(1)+Srr_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdarV(1)^2;
                
                % Viscoelastic Cauchy radial extra stress
                sigmatt_ve_compV(2:end)=sigmatt_ve_compV(2:end)+Stt_veV.*lambdatV(2:end).^2; % Circumferential Cauchy viscous stress
                sigmarr_ve_compV(2:end)=sigmarr_ve_compV(2:end)+Srr_veV.*lambdarV(2:end).^2; % Circumferential Cauchy viscous stress
                
                if(num_underscore==1)
                    Szz_eV = sigmazz_eV./lambdazV.^2;
                    delta_Szz_eV = Szz_eV(2:end)-Szz_eV(1:end-1); % Axial elastic 2ndPK extra stress increments
                    
                    % Viscoelastic 2ndPK axial extra stress
                    Szz_veV = G_fun*delta_Szz_eV+Szz_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2))); % Axial 2nd Piola-Kirchhoff viscous stress
                    
                    % Viscoelastic Cauchy axial extra stress
                    sigmazz_ve_compV(1) = sigmazz_ve_compV(1)+Szz_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdazV(1)^2;
                    sigmazz_ve_compV(2:end)=sigmazz_ve_compV(2:end)+Szz_veV.*lambdazV(2:end).^2; % Axial Cauchy viscous stress
                end
            end
              
            % Isolate actual datapoints to include in the cost function
            locator_1 = find(trigger_test == 1);
            sigmatt_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) = sigmatt_ve_compV(locator_1:end,1);
            sigmarr_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) = sigmarr_ve_compV(locator_1:end,1);
    
            sigmatt_rr_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) =...
                sigmatt_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter)-...
                sigmarr_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter);
    
            P_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) = sigmatt_rr_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter).*...
                hV(locator_1:end)./(riV(locator_1:end)+hV(locator_1:end)/2);
            
            P_expV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) = PV(locator_1:end);
            lambdatV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) = lambdatV(locator_1:end);
            
            % Weighting system to give the same fitting weight to quasi-static
            % pressure sweeps and quasi-static axial force sweeps
            weightttV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) = weight_vector(iState,i);
        
            if(num_underscore==1)
                % Isolate actual datapoints to include in the cost function
                sigmazz_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) = sigmazz_ve_compV(locator_1:end,1);
                sigmazz_rr_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) =...
                    sigmazz_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter)-...
                    sigmarr_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter);
    
                fT_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) =...
                    (2*sigmazz_rr_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter)-...
                    sigmatt_rr_ve_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter))*pi.*hV(locator_1:end).*...
                    (riV(locator_1:end)+hV(locator_1:end)/2);
    
                fT_expV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) =fV(locator_1:end);
                lambdazV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) = lambdazV(locator_1:end);
                
                P_vector_mod = P_compV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter);
                
                area_mod = sum((P_vector_mod(2:end)+P_vector_mod(1:end-1))/2.*diff(lambdatV(locator_1:end)));
                area_exp = sum((PV(locator_1:end-1)+PV(locator_1+1:end))/2.*diff(lambdatV(locator_1:end)));
                err_area = (area_mod-area_exp)/area_exp/3;
                
                % Weighting system to give the same fitting weight to quasi-static
                % pressure sweeps and quasi-static axial force sweeps
                weightzzV_final(1+counter:length(sigmatt_ve_compV(locator_1:end,1))+counter) = weight_vector(iState,i);
        
                sigmatt_qs_exp = sigmatt_expV(locator_1:end);
                sigmatt_qs_comp = sigmatt_ve_compV(locator_1:end,1);
                lambdatt_qs = lambdatV(locator_1:end);
        
                [~,max_qs] = max(sigmatt_qs_exp);
        
                sigmatt_qs_exp_loading = sigmatt_qs_exp(1:max_qs);
                sigmatt_qs_exp_unloading = sigmatt_qs_exp(max_qs:end);
                sigmatt_qs_comp_loading = sigmatt_qs_comp(1:max_qs);
                sigmatt_qs_comp_unloading = sigmatt_qs_comp(max_qs:end);
                lambdatt_qs_loading = lambdatt_qs(1:max_qs);
                lambdatt_qs_unloading = lambdatt_qs(max_qs:end);
            else
                counter_dynamic_loop = counter_dynamic_loop+1;
        
                sigmazz_ve_compV_final = sigmazz_ve_compV_final(1:end-length(sigmatt_ve_compV(locator_1:end,1)));
                fT_compV_final = fT_compV_final(1:end-length(sigmatt_ve_compV(locator_1:end,1)));
    
                fT_expV_final = fT_expV_final(1:end-length(sigmatt_ve_compV(locator_1:end,1)));
                lambdazV_final = lambdazV_final(1:end-length(sigmatt_ve_compV(locator_1:end,1)));
                weightzzV_final = weightzzV_final(1:end-length(sigmatt_ve_compV(locator_1:end,1)));
        
                sigmatt_loop_exp = sigmatt_expV(locator_1:end);
                sigmatt_loop_comp = sigmatt_ve_compV(locator_1:end,1);
                lambdatt_loop = lambdatV(locator_1:end);
        
                [~,max_loop] = max(sigmatt_loop_exp);
        
                sigmatt_loop_exp_loading = sigmatt_loop_exp(1:max_loop);
                sigmatt_loop_exp_unloading = sigmatt_loop_exp(max_loop:end);
                sigmatt_loop_comp_loading = sigmatt_loop_comp(1:max_loop);
                sigmatt_loop_comp_unloading = sigmatt_loop_comp(max_loop:end);
                lambdatt_loop_loading = lambdatt_loop(1:max_loop);
                lambdatt_loop_unloading = lambdatt_loop(max_loop:end);
        
                sigmatt_qs_exp_res = interp1(lambdatt_qs_loading,sigmatt_qs_exp_loading,lambdatt_loop_loading,'spline');
                regression_coeff_1 = polyfit(lambdatt_loop_loading,sigmatt_qs_exp_res,1);
                regression_coeff_2 = polyfit(lambdatt_loop_loading,sigmatt_loop_exp_loading,1);
        
                sigmatt_qs_exp_res = interp1(lambdatt_qs_unloading,sigmatt_qs_exp_unloading,lambdatt_loop_unloading,'spline');
                regression_coeff_3 = polyfit(lambdatt_loop_unloading,sigmatt_qs_exp_res,1);
                regression_coeff_4 = polyfit(lambdatt_loop_unloading,sigmatt_loop_exp_unloading,1);
        
                ratio_exp(counter_dynamic_loop) = (regression_coeff_4(1)+regression_coeff_2(1))/(regression_coeff_3(1)+regression_coeff_1(1));
        
                sigmatt_qs_comp_res = interp1(lambdatt_qs_loading,sigmatt_qs_comp_loading,lambdatt_loop_loading,'spline');
                regression_coeff_1 = polyfit(lambdatt_loop_loading,sigmatt_qs_comp_res,1);
                regression_coeff_2 = polyfit(lambdatt_loop_loading,sigmatt_loop_comp_loading,1);
        
                sigmatt_qs_comp_res = interp1(lambdatt_qs_unloading,sigmatt_qs_comp_unloading,lambdatt_loop_unloading,'spline');
                regression_coeff_3 = polyfit(lambdatt_loop_unloading,sigmatt_qs_comp_res,1);
                regression_coeff_4 = polyfit(lambdatt_loop_unloading,sigmatt_loop_comp_unloading,1);
        
                ratio_comp(counter_dynamic_loop) = (regression_coeff_4(1)+regression_coeff_2(1))/(regression_coeff_3(1)+regression_coeff_1(1));
                
                err_slope(counter_dynamic_loop) = ratio_comp(counter_dynamic_loop)-ratio_exp(counter_dynamic_loop);
            end
        
            counter = counter+length(sigmatt_ve_compV(locator_1:end,1));
        end
    end
    
    cost_tt = (P_compV_final - P_expV_final)./mean(P_expV_final).*weightttV_final;
        
    cost_zz = (fT_compV_final - fT_expV_final)./mean(fT_expV_final).*weightzzV_final;
    
%     cost_tt = [(sigmatt_ve_compV_final-sigmatt_expV_final).*denom_ttV_final/mean(P_expV_final).*weightttV_final; sigmatt_rr_tf*H/Ri/mean(P_expV_final)];
%     cost_zz = [(sigmazz_ve_compV_final-sigmazz_expV_final).*denom_zzV_final/mean(fT_expV_final).*weightzzV_final; sigmazz_rr_tf*pi*H.*(2*Ri+H)/mean(fT_expV_final)];
    
    costV = cat(1,cost_tt,cost_zz,err_slope,err_area); %cost function for optimization
end