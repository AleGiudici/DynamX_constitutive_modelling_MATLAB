function [PD_stretch_passive,FL_pressure_passive,... % Quasi-static experiments
    PD_DBP_passive,PD_SBP_passive,PD_Hz_passive,PD_label_passive,... % Passive dynamic experiments
    vasoactive_agents,PD_stretch_active,PD_DBP_active,PD_SBP_active,PD_Hz_active,PD_label_active,... % Active experiments
    numberTests,correctingAxialStretch] = retriveProtocol(app)
    
    correctingAxialStretch = 1;

    if(strcmp(app.protocol_menu.Value,'Old'))
        PD_stretch_passive ={'95','100','105'}; % List of guessed axial stretches for pressure sweep experiments
        FL_pressure_passive ={'10','60','100','140','200'}; % List of guessed pressures for force sweep experiments
        PD_SBP_passive = {'80','120','120','120','120','160'}; %Number of Dynamic PD Tests
        PD_Hz_passive = {'5000', '2500', '5000', '7500', '10000', '5000'};
        PD_DBP_passive = {'40','80','80','80','80','120'};
        PD_label_passive = {'physio','physio','physio','physio','physio','physio'};
        
        vasoactive_agents = {};
        PD_SBP_active = {};
        PD_DBP_active = {};
        PD_Hz_active = {};
        PD_label_active = {};
        PD_stretch_active = {};

        numberTests = length(PD_stretch_passive)+length(FL_pressure_passive)+length(PD_DBP_passive);
    
    elseif(strcmp(app.protocol_menu.Value,'New'))
        PD_stretch_passive ={'95','100','105'}; % List of guessed axial stretches for pressure sweep experiments
        FL_pressure_passive ={'10','60','100','140','180'}; % List of guessed pressures for force sweep experiments
        PD_SBP_passive = {'80','120','160','80','80','80','80','80','80','120','120','120','120','120','120','160','160','160','160','160', '160'}; %Number of Dynamic PD Tests
        PD_Hz_passive = {'10000', '10000', '10000', '625', '1250', '2500', '5000', '10000','20000','625','1250','2500', '5000', '10000','20000','625','1250','2500', '5000', '10000','20000'};
        PD_DBP_passive = {'40','80','120','40','40','40','40','40','40','80','80','80','80','80','80','120','120','120','120','120','120'};
        PD_label_passive = {'physio','physio','physio','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine'};
        
        vasoactive_agents = {};
        PD_SBP_active = {};
        PD_DBP_active = {};
        PD_Hz_active = {};
        PD_label_active = {};
        PD_stretch_active = {};
        
        numberTests = length(PD_stretch_passive)+length(FL_pressure_passive);
        if(strcmp(app.dynamic_menu.Value,'On'))
            numberTests = numberTests+length(PD_DBP_passive);
        end
    
    elseif(strcmp(app.protocol_menu.Value,'Custom'))
        load(app.load_data_dat_proc.UserData.file_protocol);
        PD_stretch_passive = protocolMat.passive.PD_stretch_guess;
        FL_pressure_passive = protocolMat.passive.FL_pressure_guess;
                
        numberTests = length(PD_stretch_passive)+length(FL_pressure_passive);
        if(isfield(protocolMat.passive,'PD_DBP_guess'))
            PD_SBP_passive = protocolMat.passive.PD_SBP_guess;
            PD_DBP_passive = protocolMat.passive.PD_DBP_guess;
            PD_Hz_passive = protocolMat.passive.PD_Hz_guess;
            PD_label_passive = protocolMat.passive.PD_label_guess;

            if(strcmp(app.dynamic_menu.Value,'On'))
                numberTests = numberTests+length(PD_DBP_passive);
            end
        end
        if(isfield(protocolMat.active, 'vasoactive_agents'))
            vasoactive_agents = protocolMat.active.vasoactive_agents;
            PD_SBP_active = protocolMat.active.PD_SBP_guess;
            PD_DBP_active = protocolMat.active.PD_DBP_guess;
            PD_Hz_active = protocolMat.active.PD_Hz_guess;
            PD_label_active= protocolMat.active.PD_label_guess;
            PD_stretch_active = protocolMat.active.PD_stretch_guess;
            
            if(strcmp(app.contraction_menu.Value,'On'))
                numberTests = numberTests+length(vasoactive_agents)*(length(PD_stretch_active)+...
                length(PD_DBP_active));
            end
        end
        if(isfield(protocolMat, 'correctingAxialStretch'))
            correctingAxialStretch = protocolMat.correctingAxialStretch;
        end
    end
end