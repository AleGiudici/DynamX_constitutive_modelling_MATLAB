function [fit_data,P_final,fT_final] = getViscoelasticBehaviour(xV_elastic,xV_visco,inputMatViscoelasticFitting,elastic,visco,n_dataPoints,fit_data,optimisation_step)
% This function returns strcuture of tables in which each table reports the
% comparison between measured and modelled viscoelastic behaviours. 
% The input parameters are the structure with the experimental data
% "inputMatElasticFitting", the vector of SEF model parameters 
% "parameterV_elastic", the vector of viscous model parameters 
% "parameterV_viscous" and the chosen elastic and viscous models "elastic" 
% and "viscous", the output matrix "fit_data" and the optimisation step 
% number "optimisation_step" (either 1 or 2).

    mtN = 0.000133; %scaling factor to convert mmHg to N/mm2
    gtN = 0.0098; %scaling factor to convert g to N

    test_list = fieldnames(inputMatViscoelasticFitting);

    P_final = zeros(n_dataPoints,1);
    fT_final = zeros(n_dataPoints,1);
    counter = 0;

    for i = 1:length(test_list)
        locator = find(inputMatViscoelasticFitting.(test_list{i}).processed_data(:,end)==1);
        nPoints = size(inputMatViscoelasticFitting.(test_list{i}).processed_data(locator:end,2),1);
        
        P_final(counter+1:counter+nPoints) = inputMatViscoelasticFitting.(test_list{i}).processed_data(locator:end,2);
        fT_final(counter+1:counter+nPoints) = inputMatViscoelasticFitting.(test_list{i}).processed_data(locator:end,6);

        counter = counter+nPoints;

        trigger_test = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,13);
    
        tV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,1)';
        lambdattV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,8);
        lambdazzV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,9);
        lambdarrV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,10);
        
        riV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,4);
        hV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,5);
        
        sigmatt_expV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,11);
        sigmazz_expV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,12);
        
        PV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,2);
        fTV = inputMatViscoelasticFitting.(test_list{i}).processed_data(:,6);
               
        tV_inverted=ones(1,length(tV))*max(tV)-tV;
        tM_temp=ones(length(tV),length(tV)).*tV_inverted-tV_inverted';
        tM = (tM_temp(2:end,1:end-1)+tM_temp(1:end-1,1:end-1))/2;
        locator=tM<0;
        tM(locator) = inf;
                
        sigmatt_ve_compV=zeros(length(tV),1);
        sigmazz_ve_compV=zeros(length(tV),1);
        sigmarr_ve_compV=zeros(length(tV),1);
            
        for iConstituent=1:elastic.nConstituents
            activator=zeros(1,2);
            activator(iConstituent)=1;
        
            % Elastic 2ndPK extra stress in the three principal directions
            [sigmatt_rr_eV,sigmazz_rr_eV,~,~,~]=elastic.fun(xV_elastic,lambdattV,lambdazzV,lambdarrV,activator); 
    
            [PeV_compV,~] = stress2pressure_force(sigmatt_rr_eV,sigmazz_rr_eV,(riV+hV/2),hV);
            
            sigmatt_eV = sigmatt_rr_eV-PeV_compV/2*mtN; % Lagrange multiplier to get actual circ. stress
            sigmazz_eV = sigmazz_rr_eV-PeV_compV/2*mtN; % Lagrange multiplier to get actual axial stress
            sigmarr_eV = -PeV_compV/2*mtN; % Lagrange multiplier to get actual radial stress
            
            % Convert Cauchy stresses to 2ndPK stresses
            Stt_eV = sigmatt_eV./lambdattV.^2; 
            Srr_eV = sigmarr_eV./lambdarrV.^2; 
            Szz_eV = sigmazz_eV./lambdazzV.^2; 
            
            delta_Stt_eV = Stt_eV(2:end)-Stt_eV(1:end-1); % Circ. elastic 2ndPK extra stress increments
            delta_Srr_eV = Srr_eV(2:end)-Srr_eV(1:end-1); % Radial elastic 2ndPK extra stress increments
            delta_Szz_eV = Szz_eV(2:end)-Szz_eV(1:end-1); % Axial elastic 2ndPK extra stress increments
            
            % Relaxation function calculated at each instant of the time interval matrix
            G_fun = visco.fun(xV_visco((iConstituent-1)*visco.nParam+1:(iConstituent-1)*visco.nParam+3),tM,0);
            locator = isnan(G_fun);
            G_fun(locator) = 0;
            
            % Viscoelastic 2ndPK circumferential extra stress
            Stt_veV = G_fun*delta_Stt_eV+Stt_eV(1)/(1+xV_visco((iConstituent-1)*3+1)*log(xV_visco((iConstituent-1)*3+3)/xV_visco((iConstituent-1)*3+2)));
            
            % Viscoelastic 2ndPK radial extra stress
            Srr_veV = G_fun*delta_Srr_eV+Srr_eV(1)/(1+xV_visco((iConstituent-1)*3+1)*log(xV_visco((iConstituent-1)*3+3)/xV_visco((iConstituent-1)*3+2)));
        
            % Viscoelastic 2ndPK axial extra stress
            Szz_veV = G_fun*delta_Szz_eV+Szz_eV(1)/(1+xV_visco((iConstituent-1)*3+1)*log(xV_visco((iConstituent-1)*3+3)/xV_visco((iConstituent-1)*3+2)));
            
            % Viscoelastic Cauchy circumferential extra stress
            sigmatt_ve_compV(1) = sigmatt_ve_compV(1)+Stt_eV(1)/(1+xV_visco((iConstituent-1)*3+1)*log(xV_visco((iConstituent-1)*3+3)/xV_visco((iConstituent-1)*3+2)))*lambdattV(1)^2;
            sigmarr_ve_compV(1) = sigmarr_ve_compV(1)+Srr_eV(1)/(1+xV_visco((iConstituent-1)*3+1)*log(xV_visco((iConstituent-1)*3+3)/xV_visco((iConstituent-1)*3+2)))*lambdarrV(1)^2;
            sigmazz_ve_compV(1) = sigmazz_ve_compV(1)+Szz_eV(1)/(1+xV_visco((iConstituent-1)*3+1)*log(xV_visco((iConstituent-1)*3+3)/xV_visco((iConstituent-1)*3+2)))*lambdazzV(1)^2;
        
            sigmatt_ve_compV(2:end)=sigmatt_ve_compV(2:end)+Stt_veV.*lambdattV(2:end).^2; % Circumferential
            sigmarr_ve_compV(2:end)=sigmarr_ve_compV(2:end)+Srr_veV.*lambdarrV(2:end).^2; % Radial
            sigmazz_ve_compV(2:end)=sigmazz_ve_compV(2:end)+Szz_veV.*lambdazzV(2:end).^2; % Axial
        end
            
        % Difference between viscoelastic Cauchy stress in the circ. and radial direction
        sigmatt_rr_ve_compV = sigmatt_ve_compV-sigmarr_ve_compV;
        
        % Difference between viscoelastic Cauchy stress in the axial and radial direction
        sigmazz_rr_ve_compV = sigmazz_ve_compV-sigmarr_ve_compV;
        
        % Computational pressure and total axial force
        [PV_compV,fT_compV] = stress2pressure_force(sigmatt_rr_ve_compV,sigmazz_rr_ve_compV,(riV+hV/2),hV);
        
        locator_1 = find(trigger_test == 1);
        
        experimentFit = [tV(locator_1:end)' PV(locator_1:end)/mtN PV_compV(locator_1:end) (riV(locator_1:end)+hV(locator_1:end))*2 ...
            riV(locator_1:end)*2 lambdattV(locator_1:end) sigmatt_expV(locator_1:end) sigmatt_ve_compV(locator_1:end)...
            fTV(locator_1:end)/gtN fT_compV(locator_1:end) lambdazzV(locator_1:end) sigmazz_expV(locator_1:end)...
            sigmazz_ve_compV(locator_1:end)];
        
        headings = ["Time [s]", "Exp. pressure [mmHg]", "Modelled pressure [mmHg]", "Outer diameter [mm]", "Inner diameter [mm]"...
            "Circ. stretch [-]", "Exp. circ. stress [MPa]", "Modelled circ. stress [MPa]", "Exp. transducer force [g]", "Modelled transducer force [g]"...
            "Axial stretch [-]", "Exp. axial stress [MPa]", "Modelled axial stress [MPa]"];
        
        if(optimisation_step == 1)
            fit_data.fittedCurves.dynamicFit.(test_list{i}) = array2table(experimentFit,"VariableNames",headings);
        else
            fit_data.fittedCurves.quasiStaticFit.(test_list{i}) = array2table(experimentFit,"VariableNames",headings);
        end
    end   
end