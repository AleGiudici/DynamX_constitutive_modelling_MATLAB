function path = findFilePathDynamX2(app,directoryName,forceLengthCheck,dynamicCheck,testName)
% According to the experiment that is being processed, this function
% returns the path of the folder in which the relevant data is stored.
% The input parameters are the app data ("app"), the name of the main
% directory where the all the experimental data are stored in subfolders
% ("directoryName"), a dichotomous variable which is 1 if the experiment
% is a force-length sweep ("forceLengthCheck"), a dichotomous variable
% which is 1 when the experiment is dynamic ("dynamicCheck"), and the name
% of the esperiment ("testName").
    
    locatorUnderscore = find(testName == '_');

    if(forceLengthCheck)
        FL_pressure = testName(locatorUnderscore+1:end);
        if(ismac)
            if(strcmp(app.protocol_menu.Value,'New') || strcmp(app.protocol_menu.Value,'Custom'))
                path = [directoryName '/FL/P' FL_pressure '_2'];
            else
                path = [directoryName '/Force sweep/P' FL_pressure];
            end
        else
            if(strcmp(app.protocol_menu.Value,'New') || strcmp(app.protocol_menu.Value,'Custom'))
                path = [directoryName '\FL\P' FL_pressure '_2'];
            else
                path = [directoryName '\Force sweep\P' FL_pressure];
            end
        end
    elseif(~dynamicCheck)
        PD_stretch = testName(locatorUnderscore+1:end);
        if(strcmp(app.protocol_menu.Value,'New') || strcmp(app.protocol_menu.Value,'Custom'))
            if(length(PD_stretch)>2)
                path = [directoryName '/PD/L' PD_stretch '_2'];
            else
                path = [directoryName '/PD/L0' PD_stretch '_2'];
            end
        else
            if(length(PD_stretch)>2)
                path = [directoryName '/Pressure sweep/L' PD_stretch(1) '.' PD_stretch(2:3)];
            else
                path = [directoryName '/Pressure sweep/L0.' PD_stretch];
            end
        end
    else
        PD_label = testName(locatorUnderscore(end)+1:end);
        PD_Hz = testName(locatorUnderscore(end-1)+1:locatorUnderscore(end)-1);
        PD_SBP = testName(locatorUnderscore(end-2)+1:locatorUnderscore(end-1)-1);
        PD_DBP = testName(locatorUnderscore(end-3)+1:locatorUnderscore(end-2)-1);

        if(strcmp(app.protocol_menu.Value,'New') || strcmp(app.protocol_menu.Value,'Custom'))
            if(strcmp(PD_label,'physio'))
                if((round(str2double(PD_Hz)/1000)-floor(str2double(PD_Hz)/1000)) == 1)
                    path = [directoryName '/Dyn/PP_P' PD_DBP '_' PD_SBP '_f' num2str(floor(str2double(PD_Hz)/1000)) '_5'];
                else
                    path = [directoryName '/Dyn/PP_P' PD_DBP '_' PD_SBP '_f' num2str(floor(str2double(PD_Hz)/1000))];
                end
            else
                if((round(str2double(PD_Hz)/1000)-floor(str2double(PD_Hz)/1000)) == 1 && round(str2double(PD_Hz)/1000)>=2)
                    path = [directoryName '/Dyn/SP_P' PD_DBP '_' PD_SBP '_f' num2str(floor(str2double(PD_Hz)/1000)) '_5'];
                elseif(round(str2double(PD_Hz)/1000) == 1 && str2double(PD_Hz)/1000 > 1)
                    path = [directoryName '/Dyn/SP_P' PD_DBP '_' PD_SBP '_f1_25'];
                elseif(round(str2double(PD_Hz)/1000) == 1 && str2double(PD_Hz)/1000 < 1)
                    path = [directoryName '/Dyn/SP_P' PD_DBP '_' PD_SBP '_f0_625'];
                else
                    path = [directoryName '/Dyn/SP_P' PD_DBP '_' PD_SBP '_f' num2str(floor(str2double(PD_Hz)/1000))];
                end
            end
        else
            path = [VEVODir_save '/Dynamic/p' PD_SBP '-' PD_DBP ' f' num2str(str2double(PD_Hz)/1000)];
        end
    end
end