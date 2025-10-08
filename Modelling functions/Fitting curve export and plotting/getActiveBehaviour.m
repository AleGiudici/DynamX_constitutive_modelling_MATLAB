function [fit_data,P_final,fT_final] = getActiveBehaviour(parameterV_elastic,parameterV_active,parameterV_viscous,inputMatActiveFitting,...
    elasticModel,activeModel,viscousModel,n_dataPoints,fit_data)
% getActiveBehaviour exports the results of the fitting of the contraction
% (active) data into the result structure "fit_data". Additionally it
% creates two vectors with the concatenated experimental pressure and axial
% force of all experiments (these are used to calculate R2 and RMSE).
%
% The function inputs are the parameters and formulation of the passive strain
% energy function (parameterV_elastic and elasticModel), the parameters and 
% formulation of the viscous relaxation function (parameterV_viscous and
% viscousModel), the input matrix for the active fitting, the cumulative number
% of datapoints in all experiments and the result structure (fit_data).


    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    gtN = 0.0098; %scaling factor to convert g to N

    contraction_states = fieldnames(inputMatActiveFitting);

    P_final = zeros(n_dataPoints,1);
    fT_final = zeros(n_dataPoints,1);
    counter = 0;

    for iState = 1:length(contraction_states)
        test_list = fieldnames(inputMatActiveFitting.(contraction_states{iState}));

        for i = 1:length(test_list)
            locator = find(inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,end)==1);
            nPoints = size(inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(locator:end,2),1);
            
            P_final(counter+1:counter+nPoints) = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(locator:end,2);
            fT_final(counter+1:counter+nPoints) = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(locator:end,6);
    
            counter = counter+nPoints;
    
            trigger_test = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,13);
        
            tV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,1)';
            lambdatV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,8);
            lambdazV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,9);
            lambdarV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,10);
            
            riV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,4);
            hV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,5);
            
            sigmatt_expV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,11);
            sigmazz_expV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,12);
            
            PV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,2);
            fTV = inputMatActiveFitting.(contraction_states{iState}).(test_list{i}).processed_data(:,6);
                   
            tV_inverted=ones(1,length(tV))*max(tV)-tV;
            tM_temp=ones(length(tV),length(tV)).*tV_inverted-tV_inverted';
            tM = (tM_temp(2:end,1:end-1)+tM_temp(1:end-1,1:end-1))/2;
            locator=tM<0;
            tM(locator) = inf;
                    
            sigmatt_ve_compV=zeros(length(tV),1);
            sigmazz_ve_compV=zeros(length(tV),1);
            sigmarr_ve_compV=zeros(length(tV),1);
                
            for iConstituent=1:elasticModel.nConstituents+activeModel.nConstituents*(iState-1)
                if(iConstituent<=elasticModel.nConstituents)
                    activator=zeros(1,2);
                    activator(iConstituent)=1;
            
                    % Elastic 2ndPK extra stress in the three principal directions
                    [sigmatt_rr_eV,sigmazz_rr_eV,~,~,~]=elasticModel.fun(parameterV_elastic,lambdatV,lambdazV,lambdarV,activator); 
                else
                    [sigmatt_rr_eV,sigmazz_rr_eV,~,~,~]=activeModel.fun(parameterV_active,lambdatV,lambdazV,lambdarV);
                end 
        
                [PeV_compV,~] = stress2pressure_force(sigmatt_rr_eV,sigmazz_rr_eV,(riV+hV/2),hV);
                
                sigmatt_eV = sigmatt_rr_eV-PeV_compV/2*mtN; % Lagrange multiplier to get actual circ. stress
                sigmazz_eV = sigmazz_rr_eV-PeV_compV/2*mtN; % Lagrange multiplier to get actual axial stress
                sigmarr_eV = -PeV_compV/2*mtN; % Lagrange multiplier to get actual radial stress
                
                % Convert Cauchy stresses to 2ndPK stresses
                Stt_eV = sigmatt_eV./lambdatV.^2; 
                Srr_eV = sigmarr_eV./lambdarV.^2; 
                Szz_eV = sigmazz_eV./lambdazV.^2; 
                
                delta_Stt_eV = Stt_eV(2:end)-Stt_eV(1:end-1); % Circ. elastic 2ndPK extra stress increments
                delta_Srr_eV = Srr_eV(2:end)-Srr_eV(1:end-1); % Radial elastic 2ndPK extra stress increments
                delta_Szz_eV = Szz_eV(2:end)-Szz_eV(1:end-1); % Axial elastic 2ndPK extra stress increments
                
                % Relaxation function calculated at each instant of the time interval matrix
                G_fun = viscousModel.fun(parameterV_viscous((iConstituent-1)*viscousModel.nParam+1:(iConstituent-1)*viscousModel.nParam+3),tM,0);
                locator = isnan(G_fun);
                G_fun(locator) = 0;
                
                % Viscoelastic 2ndPK circumferential extra stress
                Stt_veV = G_fun*delta_Stt_eV+Stt_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
                
                % Viscoelastic 2ndPK radial extra stress
                Srr_veV = G_fun*delta_Srr_eV+Srr_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
            
                % Viscoelastic 2ndPK axial extra stress
                Szz_veV = G_fun*delta_Szz_eV+Szz_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)));
                
                % Viscoelastic Cauchy circumferential extra stress
                sigmatt_ve_compV(1) = sigmatt_ve_compV(1)+Stt_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdatV(1)^2;
                sigmarr_ve_compV(1) = sigmarr_ve_compV(1)+Srr_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdarV(1)^2;
                sigmazz_ve_compV(1) = sigmazz_ve_compV(1)+Szz_eV(1)/(1+parameterV_viscous((iConstituent-1)*3+1)*log(parameterV_viscous((iConstituent-1)*3+3)/parameterV_viscous((iConstituent-1)*3+2)))*lambdazV(1)^2;
            
                sigmatt_ve_compV(2:end)=sigmatt_ve_compV(2:end)+Stt_veV.*lambdatV(2:end).^2; % Circumferential
                sigmarr_ve_compV(2:end)=sigmarr_ve_compV(2:end)+Srr_veV.*lambdarV(2:end).^2; % Radial
                sigmazz_ve_compV(2:end)=sigmazz_ve_compV(2:end)+Szz_veV.*lambdazV(2:end).^2; % Axial
            end
                
            % Difference between viscoelastic Cauchy stress in the circ. and radial direction
            sigmatt_rr_ve_compV = sigmatt_ve_compV-sigmarr_ve_compV;
            
            % Difference between viscoelastic Cauchy stress in the axial and radial direction
            sigmazz_rr_ve_compV = sigmazz_ve_compV-sigmarr_ve_compV;
            
            % Computational pressure and total axial force
            [PV_compV,fT_compV] = stress2pressure_force(sigmatt_rr_ve_compV,sigmazz_rr_ve_compV,(riV+hV/2),hV);
            
            locator_1 = find(trigger_test == 1);
            
            experimentFit = [tV(locator_1:end)' PV(locator_1:end)/mtN PV_compV(locator_1:end) (riV(locator_1:end)+hV(locator_1:end))*2 ...
                riV(locator_1:end)*2 lambdatV(locator_1:end) sigmatt_expV(locator_1:end) sigmatt_ve_compV(locator_1:end)...
                fTV(locator_1:end)/gtN fT_compV(locator_1:end) lambdazV(locator_1:end) sigmazz_expV(locator_1:end)...
                sigmazz_ve_compV(locator_1:end)];
            
            headings = ["Time [s]", "Exp. pressure [mmHg]", "Modelled pressure [mmHg]", "Outer diameter [mm]", "Inner diameter [mm]"...
                "Circ. stretch [-]", "Exp. circ. stress [MPa]", "Modelled circ. stress [MPa]", "Exp. transducer force [g]", "Modelled transducer force [g]"...
                "Axial stretch [-]", "Exp. axial stress [MPa]", "Modelled axial stress [MPa]"];
            
            fit_data.fittedCurves.activeFit.(contraction_states{iState}).(test_list{i}) = array2table(experimentFit,"VariableNames",headings);
        end   
    end
end