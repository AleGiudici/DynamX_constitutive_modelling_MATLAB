function experimental_data = processContractionData(passive, app, event, VEVODir_save)
    
    V2mN = 36.58; % Volt to mN conversion for rat aorta
    p2um = 1.638; % pixel to um conversion for rat aorta
    mN2g = 1/9.81; % mN to g comversion

    L = passive.L;

    if(ismac)
        folder_location = [VEVODir_save '/Contr/'];
    else
        folder_location = [VEVODir_save '\Contr\'];
    end

    protocol_type = app.protocol_menu.Value;

%% Protocol data
    [PD_stretch_passive,FL_pressure_passive,PD_DBP_passive,~,~,~,vasoactive_agents,PD_stretch_active,...
        PD_DBP_active,PD_SBP_active,PD_Hz_active,PD_label_active,numberTests] = retriveProtocol(app);

    vesselVolume = ((passive.OD/1000)^2/4-((passive.OD-2*passive.H)/1000)^2/4)*pi*passive.L;
    
    for iAgent = 1:size(vasoactive_agents,2)
        % if(iAgent == 1)
        %% Quasi-static pressure sweep (Fully relaxed)
        exp_name = [vasoactive_agents{iAgent} '_QS_P50_190'];

        if(ismac)
            path = [folder_location '/' exp_name];
            load([path '/Data.mat']);
            dataMat = matfile([path '/Data.mat']);
        else
            path = [folder_location '\' exp_name];
            load([path '\Data.mat']);
            dataMat = matfile([path '\Data.mat']);
        end
        
        name_fields = fieldnames(dataMat);
        time_daq = eval([name_fields{2,1} '.time_pressure']);
        p1 = eval([name_fields{2,1} '.p1']);
        p2 = eval([name_fields{2,1} '.p2']);
        force = eval([name_fields{2,1} '.load']);
        stretch = eval([name_fields{2,1} '.stretch']);
        D = eval([name_fields{2,1} '.D']);
        time_motors = eval([name_fields{2,1} '.time_stretch']);
        
        force = interp1(time_motors,force,time_daq,'linear','extrap');
        
        int_stretch=interp1(time_motors,stretch,time_daq,'linear','extrap');
        
        if(isnan(int_stretch(1)))
            int_stretch = ones(length(int_stretch),1)*1.6138;
        end
        
        target_P = 55:10:185;
        target_P_2 = 51:10:181;
        target_P_3 = 59:10:189;        
        
        for i = 1:14
            [locator] = find(p2<target_P(i));
            [~,jump] = max(diff(locator));
            
            p2_first_step = p2(locator(1:jump));
            time_daq_first_step = time_daq(locator(1:jump));

            [locator_2] = p2_first_step>target_P_2(i);
            [~,jump_2] = max(diff(locator_2));
            
            reg_coeff = polyfit(p2_first_step(jump_2+1:end),time_daq_first_step(jump_2+1:end),1);
        
            time_step = reg_coeff(1)*(target_P(i)-5)+reg_coeff(2);
            
            [~,locator2] = min(abs(time_step-time_daq));
        
            tV(i) = time_daq(locator2);
            DV(i) = D(locator2);
            p1V(i) = p1(locator2);
            p2V(i) = p2(locator2);
            lambdazV(i) = int_stretch(locator2);
            forceV(i) = force(locator2);
        end

                
        for i = 1:14
            [locator] = find(p2>target_P(15-i));
        
            p2_first_step = p2(locator);
            time_daq_first_step = time_daq(locator);

            [locator_2] = find(p2_first_step<target_P_3(end-i+1));
            [~,jump_2] = max(diff(locator_2));
        
            reg_coeff = polyfit(p2_first_step(locator_2(jump_2+1):end),time_daq_first_step(locator_2(jump_2+1):end),1);
        
            time_step = reg_coeff(1)*(target_P(15-i)+5)+reg_coeff(2);
        
            [~,locator2] = min(abs(time_step-time_daq));
        
            tV(14+i) = time_daq(locator2);
            DV(14+i) = D(locator2);
            p1V(14+i) = p1(locator2);
            p2V(14+i) = p2(locator2);
            lambdazV(14+i) = int_stretch(locator2);
            forceV(14+i) = force(locator2);
        end
        
        [locator] = find(p2<55);
        if(~isempty(locator))
            [~,jump] = max(diff(locator));
        end
        
        if(jump ~= length(locator)-1)
            p2_first_step = p2(locator(jump+1:end));
            time_daq_first_step = time_daq(locator(jump+1:end));

            [locator_2] = p2_first_step>51;
            [~,jump_2] = min(diff(locator_2));
            
            reg_coeff = polyfit(p2_first_step(1:jump_2),time_daq_first_step(1:jump_2),1);
            time_step = reg_coeff(1)*(50)+reg_coeff(2);
            
            [~,locator2] = min(abs(time_step-time_daq));
            
            tV(29) = time_daq(locator2);
            DV(29) = D(locator2);
            p1V(29) = p1(locator2);
            p2V(29) = p2(locator2);
            lambdazV(29) = int_stretch(locator2);
            forceV(29) = force(locator2);
        end
        
        Di = 2*radiusFromIincompressibility(DV/2/1000*p2um,lambdazV*L,vesselVolume,1)*1000;
    
        exp_name = ['PD_' PD_stretch_active{1}];
        
        column_heading = {'Time [s]' 'Pressure [mmHg]' 'Outer diameter [um]', 'Axial force [g]', 'Axial stretch [-]', 'Axial length [mm]', 'Inner diameter [um]'};
        if(iAgent == 1)
            active.passive_reference.(exp_name) = array2table([tV', p2V', DV'*p2um, forceV'*V2mN*mN2g, lambdazV', lambdazV'*L, Di'],"VariableNames",column_heading);
        else
            active.(vasoactive_agents{iAgent}).(exp_name) = array2table([tV', p2V', DV'*p2um, forceV'*V2mN*mN2g, lambdazV', lambdazV'*L, Di'],"VariableNames",column_heading);
        end
        % end
        
        %% Dynamic loops
        for i = 1:size(PD_label_active,2)
            list_exp = [vasoactive_agents{iAgent} '_Dyn_f' num2str(str2num(PD_Hz_active{i})/1000) '_P' PD_DBP_active{i} '_' PD_SBP_active{i}];
            exp_name = ['PD_' PD_SBP_active{i} '_' PD_Hz_active{i}];
    
            path = [folder_location list_exp];
            
            if(ismac)
                load([path '/Data.mat']);
                dataMat = matfile([path '/Data.mat']);
            else
                load([path '\Data.mat']);
                dataMat = matfile([path '\Data.mat']);
            end
            
            name_fields = fieldnames(dataMat);
            time_daq = eval([name_fields{2,1} '.time_pressure']);
            p1 = eval([name_fields{2,1} '.p1']);
            p2 = eval([name_fields{2,1} '.p2']);
            force = eval([name_fields{2,1} '.load']);
            stretch = eval([name_fields{2,1} '.stretch']);
            D = eval([name_fields{2,1} '.D']);
            time_motors = eval([name_fields{2,1} '.time_stretch']);
            
            force = interp1(time_motors,force,time_daq,'linear','extrap');
            
            int_stretch=interp1(time_motors,stretch,time_daq,'linear','extrap');
            
            if(isnan(int_stretch(1)))
                int_stretch = ones(length(int_stretch),1)*1.6138;
            end
        
            Di = 2*radiusFromIincompressibility(D/2/1000*p2um,int_stretch*L,vesselVolume,1)*1000;
        
            column_heading = {'Time [s]' 'Pressure [mmHg]' 'Outer diameter [um]', 'Axial force [g]', 'Axial stretch [-]', 'Axial length [mm]', 'Inner diameter [um]'};

            if(iAgent == 1)
                active.passive_reference.(exp_name) = array2table([time_daq, p2, D*p2um, force*V2mN*mN2g, int_stretch, int_stretch*L, Di],"VariableNames",column_heading);
            else
                active.(vasoactive_agents{iAgent}).(exp_name) = array2table([time_daq, p2, D*p2um, force*V2mN*mN2g, int_stretch, int_stretch*L, Di],"VariableNames",column_heading);
            end
        end
    end

    experimental_data.passive = passive;
    experimental_data.active = active;
end