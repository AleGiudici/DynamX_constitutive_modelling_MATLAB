function vidArt
% main program to process automatically video sequences( DICOM or *.avi) 
% to get
% - R-top
% - lumen diameter (distribution)
% - distension (distribution)
% - IMT (distribution)
% It is assumed that the video sequence is made available via a DICOM 
% file (to be preferred) or a movie with standard *.avi format. 
% Each recording device has its own image size, affecting the distribution 
% of objects across the image (size and position B-mode  
% field, size and position sonogram, ECG), requiring readjustment 
% of the region of interest (ROI) for each type of video sequence.
% 
% The program is composed of the following functions (alphabetical order):
% BeatStat     - Extracts statistics per beat
% CheckStruct  - readjusts stored settings structure to current structure
% ClearXLsheet - Clears/removes sheet from Excel-workbook
% Dist         - to read/process a videoframe and store result
% GenFilName   - checks/modifies filenames if wildcard without extension
% GenSubPlot   - generates/scales  figure with subplots
% Get_avi_Frame- select movie segment
% GetFileList  - Note for files with .mat-extension prefix does not work
% GetIma       - reads video-image from file (FIRST FRAME DISCARDED!)
% getIMT       - initiate detection of intima position get IMT
% getInfo()    - provides list of abbreviations
% getMan_ECG   - Add/adjust R-top events manually
% getMWave     - get/process mean distension wave
% getPixScal   - get depth scaling (standard unit: mm/pixel)
% getRtop      - identifies start of cardiac cycle
% getSubROI    - Splits ROI in subparts
% getWaveStat  - Extract distension waveform statistics
% getXT        - get timepoint of level crossing
% imt_ln       - detect lumen-intima transition for current line
% IniLP        - manual anterior/posterior wall idenfication
% limmess      - message with limited presentation time
% linefit      - generates (curved) regression line through datapoints
% mod_MarkPos  - Manual modify/edit intima/adventitia position
% ProcWave     - actual distension waveform analysis
% Save2txt     - save in ASCII format
% Save2XL      - save detailed results to Excel file
% set_Broi	   - sets region of interest (ROI) for B-mode 
% wdet         - search for anterior/posterior walls
% 
% Calling convention (in order of appearance):
% line 261:  GetFileList(fext)
% line 316:  CheckStruct(vArt,defvArt)
% line 408:  Dist(filename)
% line 650:  GetIma()
% line 802:  getSubROI()
% line 829:  ProcWave(fn)
% line 975:  Save2XL(ddiam,DistHead)
% line 1096: Save2txt(vlog)
% line 1216: getMWave
% line 1148: linefit(data,acclevel,ord)
% line 1189: GenFilName(fn,fext,mode)
% line 1319: ECGtr=getRtop(diam)
% line 1363: ddiam=getMan_ECG(ddiam,mdiam,hsub)
% line 1447: res=getWaveStat(rtop,sdistr)
% line 1517: xt=getXT(segm,fps,xl)
% line 1538: set_Broi(ima) 
% line 1641: wpos=IniLP(ima)
% line 1745: wpos=wdet(ima,cblk,wall,disp)
% line 1826: [fh,hsub]=GenSubPlot(npl,fn,FOS)
% line 1866: frsel=Get_avi_Frame(ima,fn)
% line 1947: rprt=ClearXLsheet(xlsfn,xltab)
% line 1989: [wpos,imt]=getIMT(ima,wpos,imt,cblk) 
% line 2084: ipos=imt_ln(echoln,wp,cstr)
% line 2174: [wpos,IMTpos]=mod_MarkPos(wpos,IMTpos,cblk,cstr)
% line 2328: getPixScal()
% line 2382: limmess(str,dur)
% line 2391: [infotxt,abbr]=getInfo()
% 
% Most functions are embedded, so for compilation use:
%   mcc -m vidArt.m startup.m
% Arnold Hoeks, Biofysica/Biomed Eng
% Jeire Steinbuch, Biomed Eng
% Bart Spronck, Biomed Eng
% Maastricht, Jan 2015

% change v1511 compared to v1510
% In wall tracking results, clicking on the diameter lines (right pane) not
% only selects the corresponding *segment* in the left pane, but also the
% corresponding time point.
% Added vArt.excLimit option.

% change v1510 compared to 1509
% Added vArt fields:
%   vArt.flist
%   vArt.mergeFiles
%   vArt.mwexplFactor
%   vArt.wexplDiamScaleFactor

% change v1509 compared to v1501
% Access to VEVO video files 

% change v1501 compared to v1407
% Double wall threshold evaluation to exclude large intima reflections

%change v1407 compared to v1406
%Transducer info is stored in results by adding lines 530-532, 808,952
% re-evaluate/reproces: previous imt values are shown instead off new imt

%change v1406 compared to v1405
%ECG trigger from movie contains zero at first value, so zero will be
%removed
%Selecting frames also selects the ECG triggers, so frsel is used for
%changing ECG triggers
%wexpl is set in mm
%new variable vArt.lwwr and set to value 3 in stead of 4.

%change v1405 compared to v1404
%dicomread can not read certain compressed images in matlab 2010, therefore
%2013 is used.

% This is the start of the main function vidArt
clear all;
global vArt
close all;
warning off;
% default settings
vArt.fn='';                                                    % filename
vArt.FOS=12;                                              % text fontsize
vArt.version='Video Dist. V15.11';
vArt.fnr=0;                                         % current file number
vArt.nrfiles=0;          % total number of files to process in this batch
vArt.frbwind=[];               % image window diameter [X,Y,width,height]
vArt.MinRoiWidth=150;                  % minimal width Region of interest
vArt.ROIwidth=30;                      % default ROI width in terms of mm
vArt.dip=4;                                  % depth interpolation factor
vArt.mwexpl=0.8;                % wall exploration range peak value in mm
vArt.wexpl=[];               % wall exploration range peak value in pixel
vArt.wthr=0.65;              % wall threshold relative to peak adventitia
vArt.wthr2=(1+vArt.wthr)/2;                   % threshold for second echo
vArt.ithr=0.6;       % threshold to search for relevant intime derivative
vArt.dres=0.05;           % estimated depth US resolution (mm) within ROI
vArt.lres=4*vArt.dres;  % estimated lateral US resolution (mm) within ROI
vArt.lwwr=3;        % lateral window width in terms of lateral resolution
vArt.wspan=7;                        % wall smoothing span in terms of mm
vArt.blksize=40;% preferential width (pixels) imagesegment wall detection
vArt.wwidth=40;       % actual width (pixels) imagesegment wall detection
vArt.nblk=5;                % number of image segments for wall detection
vArt.delseg=[];             % number of ignored/deleted segments at edges
vArt.bshift=1/2;                       % 1/2 shift segments (50% overlap)
vArt.sgwind=0.2; % smoothing window (s) sgolay-filter distension waveform
vArt.rtop=[]; % diastolic time points, assumed to be coincident wit R-top
vArt.necg=[];    % number of R-tops (beats), one more than cardiac cycles
vArt.MinRelDis=2;                        % minimal relative distension [%]
vArt.imt_fb=0;            % graphical feedback on intima detection (if 1)
vArt.acqdate=[];                          % acquisition date (DICOM)-file
vArt.fps=25;                                  % default frames per second
vArt.nfr=100;                        % default (minimum) number of frames
vArt.tmax=6;       % long video clips will be truncated to 'tmax' seconds
vArt.rWidth=[];                         % number of echo lines within ROI
vArt.rDepth=[];                % number of videolines (pixels) within ROI
vArt.rtthr=[10 90];  % lower and upper distension thresholds for risetime
vArt.pixscal=0;              % pixel scaling dimension; pixelsize in MM!!
vArt.proc=0;     % (RE-)process only distension (0) or distension AND IMT
vArt.rproc=0;                   % 0 is normal analysis, else reprocessing
vArt.probe='';                                         % probe identifier
vArt.mfext='DI.mat';   % extension for mat file with intermediate results
vArt.source='*.avi';                % video source, either *.avi or DICOM
vArt.log='Distlog.xls';
vArt.flist=cell(0); %               list of static VEVO BMODE file names.
%                                        Used when vArt.mergeFiles==true.
vArt.mergeFiles = false; %  true when combining separate VEVO BMODE files
%                                     (one frame each) into one sequence.
vArt.mwexplFactor = 5;       % conversion factor between dres and mwexpl.
%                                                      Arnolds default: 3
vArt.wexplDiamScaleFactor = 2;         % scaling factor for diameter when
%                                  adjusting diameter. Arnolds default: 4
vArt.excLimit = 4;                 %'Excursion limit'. Arnolds default: 2
defvArt=vArt;                               % copy first default settings
cfold=cd;
try 
    load('vArtSett.mat','vArt');                    % load previous 'vArt'
    vArt=CheckStruct(vArt,defvArt);                     % update structure
    disp('vidArt settings loaded from file.')
catch
  disp('No vidArt settings loaded from file. Using defaults.')
end 
vArt.version=defvArt.version; 
vArt.rproc=0;
set(0, 'DefaultUICOntrolFontSize',vArt.FOS);          % enlarge font size
srdy=ModSett();
save([cfold,'\','vArtSett.mat'],'vArt');
vArt.mwexpl=vArt.mwexplFactor*vArt.dres;        % wall exploration range peak value in mm
while 2>1                                       % loop until program exit
    if vArt.proc==0
        pstring='REPROCESS (distension waveforms only)';
    else
        pstring='REPROCESS (both distension waveforms and IMT)';
    end
    butt=listdlg('Liststring',...
        {'Exit';...                                            % option 1
        'Process *.avi recording (mm scaling)';...             % option 2
        'Process *.dcm recording (DICOM)';...                  % option 3
        'Process DICOM recording (No extension) or VEVO';...   % option 4
        pstring;...                                            % option 5
        ['General information (Version ',vArt.version,')'];... % option 6
        ['Modify settings (font size, clip limit, filter, ',...
        'resolution)']},...                                     % option 7        
        'Selectionmode','Single',...
        'ListSize',[400,200],...
        'Promptstring','Distension waveforms from video clips',...
        'Name',['B-mode Analysis ',vArt.version]);
    if butt==1                                             % request exit
        break;
    elseif butt==3                                 % read video recording
        vArt.source='*.dcm';                        % process DICOM movie
    elseif butt==4
        vArt.source='*.*';                          % process DICOM movie
        reply=questdlg({['Though currently files with and without ',...
            'extensions may be listed only files without extension ',...
            'will eventually be considered.'];...
            'Select some or all files and let the program decide!!';...
            'While previewing a movie you may still skip a file.';...
            '';'Press OK to continue'},'','DICOM','VEVO','Cancel','DICOM');
        if strcmp(reply,'Cancel')
            butt=1;
        elseif strcmp(reply,'VEVO')
            vArt.source='*.bmode';
        end
    elseif butt==2                       % then it must be an *.avi movie
        vArt.source='*.avi';
        if vArt.pixscal==0
            vArt.pixscal=-1;            % request mm-calibration procedure
        end
    end
    if butt>1&&butt<5
        [flist,vArt.nrfiles,ferr]=GetFileList(vArt.source);
        if ferr==0                   % at least one valid file identified
          if (vArt.nrfiles > 1) && strcmp(vArt.source, '*.bmode')
            answer = questdlg('Multiple .bmode files were selected. Do these files represent single frames and should the sequence be processed as one?','','No (default)', 'Yes', 'No (default)');
            if strcmp(answer, 'Yes')
              mergeFiles = true;
            else
              mergeFiles = false;
            end
          else
            mergeFiles = false;
          end
          
          if mergeFiles
              vArt.fnr=1;
              vArt.fn=char(flist(vArt.fnr)); %Use first file name here.
              vArt.flist = flist;
              vArt.ROIwidth=30;
              vArt.delseg=[];
              vArt.necg=[];
              vArt.rtop=[];
              vArt.acqdate=[];           % acquisition date (DICOM)-file
              vArt.rWidth=[];          % number of echo lines within ROI
              vArt.rDepth=[]; % number of videolines (pixels) within ROI
              vArt.pixscal=0;
              vArt.frsel=[];
              vArt.iDepth=[];
              vArt.iWidth=[];
              vArt.probe=[];
              vArt.mergeFiles=true;
              Dist();              
          else
            for fnr=1:vArt.nrfiles               % start processing files
              vArt.fnr=fnr;
              vArt.fn=char(flist(vArt.fnr));
              vArt.ROIwidth=30;
              vArt.delseg=[];
              vArt.necg=[];
              vArt.rtop=[];
              vArt.acqdate=[];           % acquisition date (DICOM)-file
              vArt.rWidth=[];          % number of echo lines within ROI
              vArt.rDepth=[]; % number of videolines (pixels) within ROI
              vArt.pixscal=0;
              vArt.frsel=[];
              vArt.iDepth=[];
              vArt.iWidth=[];
              vArt.probe=[];
              vArt.mergeFiles=false;
              Dist();                      % go to main analysis program
            end
          end
        end
        srdy=0;
    elseif butt==5                 % recall and re-evaluate processed data
        vArt.rproc=1;
        if srdy==0
            srdy=ModSett();
        end
        if vArt.rproc==1
            [flist,nrfiles,ferr]=GetFileList(['*',vArt.mfext]);
        else
            ferr=-1;
        end;
        if ferr==0
            defvArt=vArt;                         % retain previous 'vArt'
            for fnr=1:nrfiles
                matfn=char(flist{fnr});
                if strfind(matfn,vArt.mfext)    % this is a proper matfile
                    load(matfn,'vArt','cblk','rima','wposT','wposR',...
                        'ddistr','sdistr','adiam','imt');
                    vArt=CheckStruct(vArt,defvArt);     % update structure
                    vArt.proc=defvArt.proc;              % use new setting
                    vArt.rproc=1;
                    vArt.sgwind=defvArt.sgwind;      
                    vArt.rtthr=defvArt.rtthr;
                    save(matfn,'vArt','cblk','rima','wposT','wposR',...
                        'ddistr','sdistr','adiam','imt');
                    ProcWave(ddistr,matfn);   % reprocess waveform (& IMT)
                    if vArt.proc==1            % analysis IMT distribution
                        [wposR,imt]=getIMT(rima,sdistr,wposR,imt,cblk);
                    end
                    save(matfn,'vArt','cblk','rima','wposT','wposR','ddistr',...
                        'sdistr','adiam','imt');
                    Save2XL(sdistr,imt);
                end
            end
        end
        vArt.rproc=0;
    end
    if butt==6             % provide general information about the program
        mstr={'The "vidArt" program analyses B-mode ultrasound (US) clips to extract:';...
            '(1): anterior and posterior media-adventita position';...
            '(2): artery diameter and distension';...
            '(3): intima position and, hence, intima-media thickness (IMT)';...
            '(4): onset of cardiac cycles, allowing synchronized parameter extraction';...
            'All data are gathered and presented as a spatial distribution (2D), allowing evaluation at end-diastole of the mean value as well as the spatial homogeneity (spatial standard deviation) of diameter, distension and wall thickness';...
            'B-mode clips should reveal artery (wall) details, i.e:';... 
            '-the best depth resolution -> high US bandwidth -> high frequency';...
            '-limited (fixed) depth range of 3-4 cm and NO echo overexposure.';...
            'Video-clips should have either DICOM or avi-format. DICOM is preferred because it supports a high frame-rate consistent with B-mode acquisition and it carries image details.';...
            'If video-clips were recorded on DVD, they should be converted to avi-format (suggestion: "Free Video Convertor" or "Any DVD Convertor"), preferrable with the frame-rate and image size of the original recording. Use as codec (image compressor) DivX or MPEG-4. DVD to avi-conversion might be complicated because of limited selectability of image clips (simple convertors accept only single clips), codecs and image specifics. This program will reject incompatible conversions.';...
            '';'Results for the processed video files are accumulated per video file folder in an Excel logfile "Distlog.xls". The user is referred to "VidArtInfo.txt" for a detailed description of the algorithms applied to assess adventitia and intima position.'};
        waitfor(questdlg(mstr,['General INFO vidArt ',vArt.version],'OK','OK'));
    elseif butt==7
        srdy=ModSett();
        save([cfold,'\','vArtSett.mat'],'vArt');
    end
end
cd(cfold);
save('vArtSett.mat','vArt');
set(0, 'DefaultUICOntrolFontSize', 8);        % restore default font size
close all
end

function [flist,nrfiles,ferr]=GetFileList(fext)
% Function searches for files with the specified extension 'fext' and
% returns the list 'flist' as a cell-structure (independent of the number
% of file-names 'nrfiles' returned). The error flag 'ferr' will be zero if
% at least one file-name was returned (otherwise 'flist' will be empty and
% 'nrfiles'=0). If 'ferr'=0, the folder with the file is made the current
% folder (directory).
% Because Matlab ignores the prefix 'ab.mat' an additional check is
% performed on strings leading the extension.
% Arnold Hoeks, Biomed. Eng, MU, Aug 2015
[fnlist,PathName,FileIndex]=uigetfile(fext,['Select single/multiple/all'...
    ' files.'],'MultiSelect','on');
ferr=-1;
nrfiles=0;
flist=[];
if FileIndex~=0                       % at least one valid file identified
    cd(PathName);                     % switch folder to identified folder
    if iscell(fnlist)==0             % single name is returned as a string
        flist=cellstr(fnlist);
    else
        flist=sortrows(fnlist');
    end 
    nrfiles=length(flist);
    ferr=0;
end
end

function sett=CheckStruct(sett,defsett)
% function checks whether the default structure 'defsett' deviates from a
% recalled structure 'sett' and modifies 'sett' accordingly by removing and
% deleting fields. In this way the fieldnames in the stored/recalled 
% structure remain automatically consistent with the default.
% Arnold Hoeks, Biomedical Engineering, Maastricht University, April 2015
fldname1=fieldnames(defsett);                % update 'sett'-structure
fldname2=fieldnames(sett); 
for ln=1:length(fldname1)
    if isempty(find(strcmp(fldname1(ln),fldname2)>0, 1))     % add field
        sett=setfield(sett,fldname1{ln},getfield(defsett,fldname1{ln}));
    end
end
fldname2=fieldnames(sett); 
for ln=1:length(fldname2)                     % remove superfluous fields
    if isempty(find(strcmp(fldname2(ln),fldname1)>0, 1))      % add field
        sett=rmfield(sett,fldname2{ln});
    end
end
end

function srdy=ModSett()
global vArt
butt=2;
valstr={'0.05';'0.1';'0.15';'0.2';'0.25';'0.30'};
value=[0.05, 0.1, 0.15, 0.2, 0.25, 0.30];
while butt>1
    dres=sprintf('%d',round(1000*vArt.dres));
    if vArt.proc==0
        pstring='Process distension waveforms only';
    else
        pstring='Process both distension waveforms and IMT';
    end
    butt=listdlg('Liststring',...
        {'Return to Main Menu';...                              % option 1
        pstring;...                                             % option 2
        ['Record limit (now ',num2str(vArt.tmax),' s)'];...     % option 3
        ['Distension rise time threshold (now ',...
            num2str(vArt.rtthr(1)),' %)'];...                   % option 4
        ['Distension waveform lowpass filter (now ',...         % option 5
            num2str(vArt.sgwind),' s)'];...
        ['Ultrasound depth resolution (now ',...
            dres,' um)'];...
        ['Font Size (currently ',num2str(vArt.FOS),')']},...    % option 7
        'Selectionmode','Single',...
        'ListSize',[375,160],...
        'Promptstring','Modify settings',...
        'Name',['B-mode Analysis ',vArt.version]);
    if butt==2
        vArt.proc=mod(vArt.proc+1,2);
    elseif butt==3
        rep=questdlg(['Select Clip Limit (currently truncated at ',...
            num2str(vArt.tmax),' s)'],'','6','9','12',num2str(vArt.tmax));
        if ~isempty(rep)
            vArt.tmax=sscanf(rep,'%d');
        end
    elseif butt==4                 % adjust setting for risetime threshold
        rep=questdlg(['Select distension risetime threshold (currently ',...
            num2str(vArt.rtthr(1)),' %)'],'','10','15','20',...
            num2str(vArt.rtthr(1)));
        if ~isempty(rep)
            vArt.rtthr(1)=sscanf(rep,'%d');% lower threshold  for risetime
            vArt.rtthr(2)=100-vArt.rtthr(1);% upper threshold for risetime
        end
    elseif butt==5               % set spanwidth distension lowpass filter
        rep=listdlg('liststring',valstr,'SelectionMode','single',...
            'InitialValue',find(value==vArt.sgwind),'PromptString',...
            'Spanwidth lowpass (s)','ListSize',[170,150],'Name','Select');
        if ~isempty(rep)
            vArt.sgwind=sscanf(valstr{rep},'%f');
        end
    elseif butt==6               % set ultrasound depth resolution
        rep=char(inputdlg_check(['Depth resolution in um (lateral res.',...
            ' 4 times larger)'],'US Resolution',1,...
            {sprintf('%d',1000*vArt.dres)}));
        if ~isempty(rep)
            vArt.dres=sscanf(rep,'%d')/1000;             % expressed in mm
            vArt.lres=4*vArt.dres;    % lateral resolution default 5 times
        end
    elseif  butt==7
        rep=questdlg(['Select Font Size (currently ',...
            num2str(vArt.FOS),')'],'','8','10','12',num2str(vArt.FOS));
        if ~isempty(rep)
            vArt.FOS=sscanf(rep,'%d');
        end
        set(0, 'DefaultUICOntrolFontSize',vArt.FOS);   % enlarge font size
    end
end
srdy=1;
end

function Dist()           %                 read/processes video sequences
global vArt hsub
vArt.fps=25;                                   % default frames per second
vArt.nfr=100;                         % default (minimum) number of frames
[ima,err]=GetIma();                                  % read video sequence
vArt.fps=round(10*vArt.fps)/10;                      % round to one deimal
if err<0
    return
end
rpcont=-1;                                 % skip/repeat/continue selector
while rpcont<0
    frsel=Get_avi_Frame(ima,vArt.fn);      % Get begin and end-frame movie
    if frsel(1)==0
        return                                   % current file is skipped
    end
    dell_first=find(vArt.rtop<frsel(1));
    dell_end=find(vArt.rtop>frsel(2));
    if ~isempty(dell_end)||~isempty(dell_first)
        vArt.rtop_org=vArt.rtop;
        if isempty(dell_end);
            dell_end=length(vArt.rtop)+1;
        end
        if ~isempty(dell_first);
            if dell_end(1)==dell_first(end)+1
                vArt.rtop=[];
            else
                vArt.rtop=vArt.rtop(dell_first(end)+1:dell_end(1)-1);
                vArt.rtop=vArt.rtop-frsel(1)+1;
            end
        else
            vArt.rtop=vArt.rtop(1:dell_end(1)-1);
        end
    end
    vArt.nfr=frsel(2)-frsel(1)+1;     % actual number of frames to process
    vArt.frsel=frsel;
    [wpos1,err,cblk]=set_Broi(squeeze(ima(:,:,frsel(1))),1);     % Set ROI
    wpos=wpos1;
    if err<0
        break;
    end
    frbwind=vArt.frbwind;
    lcp=(wpos(1,:)+wpos(2,:))/2;                   % lumen center position
    xh=1:vArt.rWidth;
    [fh, hsub]=GenSubPlot(2,vArt.fn,vArt.FOS);
    axp=get(hsub(2),'Position');
    axp=[axp(1),0.015,axp(3),0.035];
    ch2=uicontrol('style','text','units','normalized','position',axp,...
        'string','Left button down within figure to pause/continue');
    set(fh,'WindowButtonDownFcn',@DoPause);
    menustr={'Echo/diameter';'anterior wall';'posterior wall'};
    fsh=uicontrol('style','popupmenu','string',menustr,...
        'units','normalized','Position',[0.03 0.015 0.1 0.035],...
        'callback', @doChange);
    taxis=(0:vArt.nfr-1)/vArt.fps;
    tstr='movie time(s)';
    subplot(hsub(2));
    if vArt.pixscal>0
        plot(taxis,mean(lcp)*vArt.pixscal);  
        ylabel('diameter (mm)');
    else
        plot(taxis,mean(lcp));                           % use pixel-scale
        ylabel('diameter (pixel)');
    end
    xlabel(tstr);
    wdisp=0;             % show intermediate results wall-lumen detection
    lbp=2*vArt.wexpl; %        minimum/safe range to lumen-wall interface
    wrange=4*vArt.wexpl;  % data window range around lumen-wall interface
    tmp=zeros(wrange+1,vArt.nblk);
    wref=zeros(2,vArt.nblk);      % reference for estimated wall positions
    wposT=zeros(2,vArt.nblk,vArt.nfr);
    ddistr=zeros(vArt.nblk,vArt.nfr);
    ddistr(1:vArt.nblk,1:vArt.nfr)=mean(wpos(2,:)-wpos(1,:));
    ddistr(:,1)=(wpos(2,:)-wpos(1,:))';
    wspan_pix=round(vArt.wspan/vArt.pixscal);            % wspan in pixels
    shift=vArt.bshift*vArt.wwidth;
    wspan=round((wspan_pix+shift-vArt.wwidth)/shift);   % smoothing window
    wspan=max(3,wspan); % minimum smoothing window =3 for 2th order filter
    iwspan=wspan_pix;                         % smoothing window for image
    for fr=1:vArt.nfr
        if fr>vArt.nfr-2      % ensure for the last frame diameter display
            wdisp=0;
        end
        im1=imcrop(squeeze(ima(:,:,fr+frsel(1)-1)),frbwind);
        if fr<5     % set wall reference to start-up wall-lumen transition
            wref=(wref*(fr-1)+wpos(:,:))/fr;                    % average
            wref=min(frbwind(4)-(3*vArt.wexpl),...
                max(3*vArt.wexpl,wref));            % enforce some margin
        end
        for wall=1:2
            for ln=1:vArt.nblk        % compute local average skewed area!
                bp=round(wpos(wall,ln)-lbp);
                ep=min(frbwind(4)-2,bp+wrange);
                tmp(:,ln)=mean(im1(bp:ep,cblk(1,ln):cblk(3,ln)),2);
            end
            wpos(wall,:)=wpos(wall,:)+wdet(tmp,wall,wdisp);  % update wall
            wpos(wall,:)=smooth(wpos(wall,:)',wspan,'sgolay')';
            % restrict excursions with respect to 'wref' 
            wpos(wall,:)=min(wref(wall,:)+vArt.excLimit*vArt.wexpl,...
                max(wref(wall,:)-vArt.excLimit*vArt.wexpl,wpos(wall,:)));
        end
        ddistr(:,fr)=wpos(2,:)-wpos(1,:);
        lcp=lcp-(lcp-(wpos(1,:)+wpos(2,:))/2)/2;%LPF lumen center position
        wposT(:,:,fr)=wpos;
        % Now display results
        if wdisp==0         % diameter display requested
             figure(fh)
             plotImage()%(im1, cblk, wpos, lcp, fr)
            subplot(hsub(2));   % 2D-diameter distribution (time-position)
            hold off
            if vArt.pixscal>0
                rawlines = plot(taxis,ddistr*vArt.pixscal);
                ylabel('diameter (mm)');
            else
                rawlines = plot(taxis,ddistr);
                ylabel('diameter (pixel)');
            end
            axis tight;
            xlabel(tstr);
        end
        pause(0.02)
    end  
    set(fh,'WindowButtonDownFcn',[]);      % remove pause/continue option
    delete(fsh);    
    rpcont=1;
    uicontrol('style','popupmenu','units','normalized',...
        'Position',[0.03 0.015 0.15 0.035],'string',{'Select action:',...
        'Continue (wave processing)','Save figure',...
        'Repeat processing','Skip file'},'callback', @doCont); 
    set(ch2,'string',...
     'Select a line to show corresponding region in left figure'); 
    widthpixel=cell(2,vArt.nblk);
    for line=1:vArt.nblk       % link lines with the corresponding regions
        widthpixel{line} = [cblk(1,line),cblk(3,line)];
        setappdata(rawlines(line), 'messages', widthpixel{line});
    end 
    dcm_obj = datacursormode(fh);               % enable data cursor mode
    set(dcm_obj,'UpdateFcn',@myupdatefcn, ...
        'DisplayStyle','datatip','Enable','on');
    uiwait(fh);       % wait until "continue",skip or "repeat" is selected
    close all 
    if rpcont>0     % if rpcont<0 -> repeat; elseif rpcont==0 -> skip file
        if (rpcont>0)                   % 1: continue, 0: skip; -1: repeat
            [~, mname]=fileparts(vArt.fn);
            matfn=[mname,vArt.mfext];     % name MAT-file to store results
            [sdistr,adiam,rpcont]=ProcWave(ddistr,matfn);  % wave process.
            wposR=wposT(:,:,vArt.rtop);     % wall positions at R-top only
            frbwind=round(vArt.frbwind); 
            rima=zeros(frbwind(4)+1,frbwind(3)+1,vArt.necg,'uint8');
            for fr=1:vArt.necg
                rima(:,:,fr)=imcrop(squeeze(ima(:,:,vArt.rtop(fr)+frsel(1)-1)),frbwind);
            end
            imt=[];
            if vArt.proc==1                            % analysis IMT distribution
                [wposR,imt]=getIMT(rima,sdistr,wposR,[],cblk);
            end
            save(matfn,'vArt','cblk','rima','wposT','wposR','ddistr',...
                'sdistr','adiam','imt');
            Save2XL(sdistr,imt);
       elseif rpcont==0 % skip further processing
            break
        end
    end
end
% nested functions for callback routines
  function plotImage %Plot image in left subplot, together with tracked walls and centre line.
    subplot(hsub(1));
    imagesc(im1);                           % present B-mode image
    hold on;
    colormap gray;
    swp=zeros(2,vArt.rWidth);
    for wall=1:2
      swp(wall,:)=interp1(cblk(2,:),wpos(wall,:),...
        1:vArt.rWidth,'linear','extrap');
      swp(wall,:)=smooth(swp(wall,:)',iwspan,'sgolay')';
    end
    scp=interp1(cblk(2,:),lcp,1:vArt.rWidth,'linear','extrap');
    plot(xh,swp,'b',  'linewidth',2);                % smoothed wall
    plot(xh,scp,'--w','linewidth',2);       % smoothed center line
    title(sprintf(['Current frame %d (of %d frames), ',...
      'segmented in %d blocks'],fr,vArt.nfr,vArt.nblk));
    xlabel('width (pixel)');
    ylabel('depth (pixel)');
    hold off;
  end
  function txt=myupdatefcn(~,event_obj)
    %Update left (image) plot to clicked time position
    pos = get(event_obj,'Position');
    t = pos(1); %Time clicked in right plot.
    [~, fr] = min(abs(t-taxis)); %Find frame number from time
    im1=imcrop(squeeze(ima(:,:,fr+frsel(1)-1)),frbwind);
    wpos = wposT(:,:,fr);
    lcp=(wpos(1,:)+wpos(2,:))/2;
    plotImage()    
    
    %Do highlighting (Jeires work)
    set(rawlines,'LineWidth',1)
    width = getappdata(event_obj.Target, 'messages');
    txt   = sprintf('width (pixel): %d to %d',width(1),width(2));
    set(event_obj.Target,'LineWidth',2)
    subplot(hsub(1));
    hold on;
    try            % remove overlay line segment with increased thickness
      linehandle=findobj(gca, 'type', 'line', 'linewidth', 3);
      delete(linehandle);
    end
    
    %Find swp
    swp=zeros(2,vArt.rWidth);
    for wall=1:2
      swp(wall,:)=interp1(cblk(2,:),wpos(wall,:),...
        1:vArt.rWidth,'linear','extrap');
      swp(wall,:)=smooth(swp(wall,:)',iwspan,'sgolay')';
    end
    
    plot(width(1):width(2),swp(1:2,width(1):width(2)),...
      'linewidth',3,'Color','r');         % create overlay line segment
    hold off;
  end
function doChange(~,~)% select other presentation (interface <-> waveform)
    wdisp=get(fsh,'Value')-1;
end
function doCont(ch1,~)
    val=get(ch1,'Value');
        set(ch1,'Visible','on');
    if val==3                               % save image
        print(fh,'-djpeg',GenFilName(vArt.fn,'*.jpg',0));
        set(ch1,'Value',1);
    elseif val==4                                     % repeat processing
        rpcont=-1;
        %vArt.rtop=vArt.rtop_org;
    elseif val==5                                             % skip file
        rpcont=0;
    end
    if val==2||val==4||val==5         % Continue (continue, repeat, skip)
        delete(fh)
        uiresume
    end
end
function DoPause(fh,~)
    set(fh,'WindowButtonDownFcn',@DoResume);
    uiwait
end
function DoResume(fh,~)
    set(fh,'WindowButtonDownFcn',@DoPause);
    uiresume
end
end


function [ima,err]=GetIma()
% Function reads entire videofragment or video-images up to a maximum of
% 6.5 seconds from either DICOM or avi-source and returns result in 3D
% array 'ima' [y,x,fr]. Dicom movies may contain a list with 'rtop'. 'Err'
%  will be set to -1 if the function fails to read the specified segment.
global vArt;
rtv=[];
vArt.rtop=[];
ima=[];
err=-1;
datstr='';
vArt.probe='';
% First evaluate video file info
if strcmp(vArt.source,'*.*')||strcmp(vArt.source,'*.dcm')  % DICOM file
    try
        info=dicominfo(vArt.fn);
    catch
        wh=warndlg({'This is NOT a DICOM movie!';...
            'Try again, but select a proper file!';...
            '';'Press OK to continue'});
        uiwait(wh); 
        return                                     % abort/skip this file
    end
    if isfield(info,'StudyDate')
        datstr=info.StudyDate;
    end
    if isfield(info.SequenceOfUltrasoundRegions.Item_1,...
            'PhysicalDeltaY')  % 1 pixel=??mm (thus about 0.09 mm)
        vArt.pixscal=...
            10*info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY;
    end
    if ~isfield(info,'NumberOfFrames')||(info.NumberOfFrames<20) 
        return               % NOT a movie, probably single DICOM picture
    end
    vArt.nfr=info.NumberOfFrames;
    if isfield(info,'RWaveTimeVector');% Philips iu22 may store ECG R-top
        rtv=info.RWaveTimeVector;                           % R-top in ms
        if rtv(1)==0&&length(rtv)>1
            rtv=rtv(2:end);
        end
    end
    if isfield(info,'FrameTime')% Philips iu22 may include recording time
        vArt.fps=1000./info.FrameTime;
        if ~isempty(rtv)                   % express rtop as frame number
            vArt.rtop=round(double(rtv/info.FrameTime));
        end
    end
    if isfield(info,'TransducerData')
        vArt.probe=info.TransducerData(1:end-2);
    end
elseif strcmp(vArt.source,'*.avi')          % then it must be an avi movie
    try
        iobj=VideoReader(vArt.fn);
        vArt.nfr=get(iobj, 'NumberOfFrames');
        frh=get(iobj, 'Height');
        frw=get(iobj, 'Width');
        if isempty(vArt.nfr)||(vArt.nfr<20)  
            return                  % NOT a movie, probably single picture
        end
        vArt.fps=get(iobj,'FrameRate');
    catch
        wh=warndlg({'Cannot open/read this avi-file!';...
            'Will skip this file!';'';'Press OK to continue'});
        uiwait(wh);
        return
    end
elseif strcmp(vArt.source,'*.bmode')       % it is a VEVO B-mode data file
  try
    if vArt.mergeFiles %Merge multiple .bmode frames into one movie
      
      for fnr=1:vArt.nrfiles
        fn = vArt.flist{fnr};
        str=sprintf(['Reading VEVO file "%s".\n\nThis is file %d of ',...
          '%d files.\nB-mode image is also converted to square pixels.\n'],...
          fn, fnr, vArt.nrfiles);
        if fnr == 1
          wm = waitbar(0, str);
        else
          waitbar(fnr/vArt.nrfiles, wm, str);
        end
        [imaTemp,pixSize,ferr]=readVEVOBModeRAW(fn,-1);
        if (ferr<0)||(isempty(vArt.nfr))%||(vArt.nfr<20)
          return
        end
        if size(imaTemp,3) > 1
          error(['File ' fn ' contains >1 frames'])
        end
        if fnr == 1
          ima1 = zeros([size(imaTemp) vArt.nrfiles], class(imaTemp));
          ima1(:,:,1) = imaTemp;
        else
          ima1(:,:,fnr) = imaTemp;
        end
      end
      close(wm);
      vArt.fps=10;            % frames per second. Just an arbitrary number here, but may affect temporal filtering.
      %If vArt.fps is not overwritten, Inf will be there.
    else
      str=sprintf(['Reading VEVO file "%s".\n\nThis is file %d of ',...
        '%d files.\nB-mode image is also converted to square pixels.\n'],...
        vArt.fn,vArt.fnr,vArt.nrfiles);
      wm=msgbox({str,'Have some patience!'});
      [ima1,pixSize,ferr]=readVEVOBModeRAW(vArt.fn,-1);   % read all frames
      if (ferr<0)||(isempty(vArt.nfr)||(vArt.nfr<20))
        return
      end
      
      delete(wm);
    end
    vArt.nfr=size(ima1,3);
    
    [npy,npx,nfr]=size(ima1);
    pxratio=pixSize(1)/pixSize(2);          % ratio of depth and width
    npx=ceil(npx/pxratio);
    ima=zeros(npy,npx,nfr);
    for fr=1:nfr                             % convert to square pixel
      ima(:,:,fr)=imresize(ima1(:,:,fr),[npy npx]);
    end
    
  catch
    wh=warndlg({'Cannot open/read this VEVO Bmode-file!';...
      'Will skip this file!';'';'Press OK to continue'});
    uiwait(wh);
    return
  end
  err=0;
  return
end
% Now set (number of) frames to read
if vArt.nfr/vArt.fps>vArt.tmax       % recording longer than max. seconds
    frindex=1:floor(vArt.tmax*vArt.fps);
else
    frindex=[];
end
if length(datstr)<1           % no info on acquisition date; use file date
    prop=dir(vArt.fn);                               % get file properties
    datstr=datestr(prop.date, 'yyyymmdd');
end  
vArt.acqdate=sprintf('%s/%s/%s',datstr(1:4),datstr(5:6),datstr(7:8)); 
if strcmp(vArt.source,'*.*')||strcmp(vArt.source,'*.dcm')  % DICOM file
    try 
        str=sprintf('Reading DICOM file "%s".\n\nThis is file %d of %d files.\n',...
            vArt.fn,vArt.fnr,vArt.nrfiles);
        wm=msgbox({str,'Have some patience!'});
        if isempty(frindex)                             % read all frames
            ima1=dicomread(vArt.fn,'frames','all');
        else                                  % read only selected frames
            ima1=dicomread(vArt.fn,'frames',frindex);
        end
        delete(wm);
    catch                                  % catch possible error message
        wh=warndlg({'Cannot accomodate such a large file!';...
            'Try to reduce movie-size; Will skip this file!';...
            '';'Press OK to continue'});
        uiwait(wh);
        return
    end
elseif strcmp(vArt.source,'*.avi')          % then it must be an avi movie
    try
        str=sprintf('Reading *.avi file "%s".\n\nThis is file %d of %d files.\n',...
            vArt.fn,vArt.fnr,vArt.nrfiles);
        wm=msgbox({str,'Have some patience!'});
        if isempty(frindex)                             % read all frames
            ima1=read(iobj);
        else                                 % read only selected frames
            ima1=zeros(frh,frw,3,length(frindex),'uint8');
            for fr=1:length(frindex)
                ima1(:,:,:,fr)=read(iobj,frindex(fr));
            end
        end
        delete(wm);
    catch
        wh=warndlg({'Cannot accomodate such a large file!';...
            'Try to reduce movie-size; Will skip this file!';...
            '';'Press OK to continue'});
        uiwait(wh);
        return
    end
end
[npy,vArt.rWidth,~,vArt.nfr]=size(ima1);
ima=zeros(npy,vArt.rWidth,vArt.nfr-1,'uint8');  %  memory for entire video
for fr=2:vArt.nfr                  % conversion requires large memory size
    ima(:,:,fr-1)=rgb2gray(squeeze(ima1(:,:,:,fr)));
end
vArt.nfr =vArt.nfr-1;% First frame discarded because it might be corrupted
err=0;
end

function cblk=getSubROI()
% Function splits ROI in halfoverlapping subsegments (blocks), allocates
% memory for outputs: DistHead (cell), ddiam (diameter waveforms)
global vArt
blksize=vArt.lwwr*vArt.lres/vArt.pixscal; 
shift=vArt.bshift*blksize;
vArt.nblk=round((vArt.rWidth+shift-blksize)/shift);   % number of segments
vArt.nblk=2*floor(vArt.nblk/2)+1;                    % odd number segments
vArt.wwidth=vArt.rWidth/(vArt.bshift*(vArt.nblk-1)+1);
cblk=zeros(3,vArt.nblk);
if vArt.nblk>6                      % discard at both sides edge waveforms
    vArt.delseg=1;
    if vArt.nblk>8
        vArt.delseg=2;            % many waveforms: discard 2 at each side
    end
else
    vArt.delseg=0;
end
for blk=1:vArt.nblk                         % calculate limits of segments
    cblk(1,blk)=(blk-1)*vArt.rWidth/(vArt.nblk+1/vArt.bshift-1)+1;%left 
    cblk(2,blk)=cblk(1,blk)+vArt.wwidth/2-1;             % center position
    cblk(3,blk)=cblk(1,blk)+vArt.wwidth-1;                 % right extreme
end
cblk=round(cblk);
cblk(3,vArt.nblk)=vArt.rWidth;
end

function [sdistr,adiam,cont]=ProcWave(ddistr,matfn)      % wave processing
% Smooths the curves,
% deletes the curves at both ends, and extracts waveform characteristics.
% If waveforms are inconsistent, the user can call for repeated movie
% analysis with an adjusted window size (function returns "cont"=-1).
% 'cont' will be 0 for an insufficient number of cardiac cycles, otherwise
% 'cont' will be 1.
global vArt
tbase=((1:vArt.nfr)-1)/vArt.fps;                              % time base
[fh,hsub]=GenSubPlot(3,vArt.fn,vArt.FOS);        % figure with 3 subplots
set(fh,'NumberTitle','off','Name',['File: ',strrep(vArt.fn,'_','-')]);
if vArt.pixscal>0
    ylstr='Diameter (mm)';
    pixscal=vArt.pixscal;
else
    ylstr='Diameter (pixel)';
    pixscal=1;
end
subplot(hsub(1));
plot(tbase,ddistr*pixscal);
xlabel('Video time (s)');
ylabel(ylstr);
title('Unfiltered (raw) distension distribution');
axis tight
if vArt.rproc==0
    pstring='Edit manually R-top positions';
else
    pstring='Edit R-top positions blocked';
end
menustr={'Select',...
    'Repeat processing',...                                    % option 2
    pstring,...                                                % option 3
    'Continue (IMT without feedback)',...                      % option 4
    'Continue (IMT WITH feedback)',...                         % option 5
    'General Info regarding options'};                         % option 6
valstr={'Set lp-filter(s)','0.05','0.1','0.15','0.2','0.25','0.30'};
if vArt.proc==0
    menustr(4:5)={'Continue (next file)','Continue (next file)'};
end
cont=0;
ford=3;                                                     % filter order
rectime=vArt.nfr/vArt.fps;                        % total recording time
fspan=2*round(2.5*vArt.fps/2)+1;   % 2 seconds span-width highpass filter
sdistr=zeros(vArt.nblk,vArt.nfr);
ah=uicontrol('style','popupmenu','string',valstr,'units','normalized',...
    'Position',[0.14 0.015 0.1 0.035],'callback',@doAdjust);%setting filter
wh=uicontrol('style','popupmenu','string',menustr,'units','normalized',...
    'Position',[0.03 0.015 0.1 0.035],'callback', @doWait);         % menu
ch1=uicontrol('style','text','units','normalized','position',...
    [0.3 0.015 0.4 0.035]);
ch2=uicontrol('style','text','units','normalized','position',...
    [0.3 0.015 0.4 0.035],'Visible','off','BackgroundColor','r');
while cont==0                % loop until (temporal) smoothing is accepted
    sgwind=max(ford+1,2*round((vArt.sgwind*vArt.fps)/2)+1);  %smooth wind.
    for blk=1:vArt.nblk
        tmp=smooth(ddistr(blk,:),sgwind,'sgolay',ford);  % temporal smooth
        if rectime>fspan/vArt.fps
            tmp=tmp-smooth(tmp,fspan,'sgolay',2)+mean(tmp);  % hp-filter!!
        end
        sdistr(blk,:)=tmp;
    end
    subplot(hsub(2));
    hold off
    plot(tbase,sdistr(1+vArt.delseg:end-vArt.delseg,:)*pixscal);%center ROI
    hold on
    plot(tbase,sdistr(1:vArt.delseg,:)*pixscal,'r--');% left/right red
    plot(tbase,sdistr(end-vArt.delseg+1:end,:)*pixscal,'r--');
    mstr=sprintf(['%d-order Golay smoothing filter %4.2fs',...
        ' @%4.1f fps (N=%d)\nDashed lines are at left/right edges'],...
        ford,vArt.sgwind,round(10*vArt.fps)/10,vArt.nblk);
    title(mstr);
    xlabel('Video time (s)');
    ylabel(ylstr);
    axis tight
    getMWave(hsub(3),sdistr(vArt.delseg+1:end-vArt.delseg,:),pixscal);
    if length(vArt.rtop)<3
        beep
        set(ch1,'Visible','off');
        set(ch2,'Visible','on','string',...
            'Number of cardiac cycles is low: Edit R-top');
    else
        set(ch1,'Visible','on','string','Check onset cardiac cycle');
        set(ch2,'Visible','off');
    end
    uiwait
%     if vArt.rproc==1
%         cont=1;
%         pause(1);
%     else
%         uiwait;                          % loop is terminated if 'cont'~=0
%     end
end
set(fh,'PaperPositionMode','auto');           % save/replace current image
print(fh,GenFilName(vArt.fn,'*.png',1),'-dpng','-r150');
vArt.necg=length(vArt.rtop);
adiam=getMWave(hsub(3),sdistr(vArt.delseg+1:end-vArt.delseg,:),pixscal);
if vArt.necg<2 
    questdlg({'Insufficient number (<2) of cardiac cycles!';...
        'Impossible to perform statistical analysis.';...
        'Will skip further processing and save procedure.';...
        'If possible, return to manual adjustment R-top position';...
        'via repeated analysis of the video movie.';...
        ' ';'Press OK to continue'},'Statistical problem','OK','OK');
    cont=0;                                      % skip further processing
end 
close(fh);
function doAdjust(~,~)% adjust interactively span width temporal smoothing
    psel=get(ah,'Value');
    if psel>1
        vArt.sgwind=sscanf(char(valstr{psel}),'%f');
    end
    set(ah,'Value',1);
    uiresume
end
function doWait(~,~)                            % evaluate menu selection
    psel=get(wh,'Value');
    if psel==2
        cont=-1;            % request/initiate repeat processing of movie
        uiresume
    elseif psel==3&&(vArt.rproc==0) % modify/set manually positions R-top
        set([ah,wh],'Visible','off');
        set(ch1,'visible','on','string',...
            'Click adds, but within 0.3 sec removes, marker.');
        getMan_ECG(mean(sdistr(vArt.delseg+1:end-vArt.delseg,:)),...
            hsub(3),pixscal);
        set([ah,wh],'Visible','on');
    elseif psel==4                                             % continue
        cont=1;
        vArt.imt_fb=0;
        uiresume
    elseif psel==5                                             % continue
        cont=1;
        vArt.imt_fb=1;       % enable graphical feedback intima detection
        uiresume
    elseif psel==6      % provide general infoabout the available options
        mstr={'The current figure is saved automatically when leaving this window.';
            '';'Waveform inconsistency may be related to inadequate selection of the region of interest (ROI) or initial adventitia identification. It is the best to repeat movie processing and pay careful attention to the detection process and movie quality (motion artifacts?)';
            '';'Visual inspection of the average distension waveform quickly reveals missing beats. The option "Edit R-top" allows addition, removal and repositioning of any or all end-diastolic markers (should ALWAYS precede the end-diastolic minimum).';
            '';'The process of adventitia and intima detection is quite similar. However, the first one heavily depends on amplitude level while the second one considers the maximum of the first derivative. The problem is which is the true maximum (the first one inwards from the adventitia or a high one outwards from the lumen if the intima/media exhibits additional peaks?). The detection process can be followed (at a reduced processing rate) by activating the feedback option.';
            '';'In the subsequent processing step, manual editing of the anterior and/or posterior adventitia positioning is followed by intima processing, involving intima detection (if processing started with the evaluation of a video) and manual intima adjustment. However, intima detection is bypassed if the file was processed before (intermediate results were recalled from a "*.DI.mat"-file), unless "intima feedback" was activated. Editing of a previously stored intima without feedback avoids repetitive editing.'};
        questdlg(mstr,'Processing info','OK','OK');
    end
    set(wh,'Value',1);
end
end

function Save2XL(sdistr,imt)
global vArt
% And save the results to Excel logfile
ncycl=min(6,vArt.necg-1);
wah=warndlg({'Saving processing results to Excel-file';...
    'Wait with your next action until this warning disappears'});
vlog=cell(12,40);
vlog{1,1}=vArt.version;
vlog(2,2:3)={'TBD','TBD'};
vlog{3,3}=vArt.probe;
vlog(2,9:40)={'RR(ms)',...                                        % 1 par
    'diast.mean','diast.med','diast.hom','diast.A',...            % 4 par
    'syst.mean','syst.med','syst.hom','syst.A'...                 % 4 par
    'dist.mean','dist.med','dist.hom','dist.A',...                % 4 par
    'rdist.mean','rdist.median','rdist.hom','rdist.A(%)',...      % 4 par
    'RT.mean(ms)','RT.median(ms)','RT.hom(ms)','RT.A',...         % 4 par
    'IMT.mean','IMT.med','IMT.hom','IMT.A',...                    % 4 par
    'rIMT.mean(%)','rIMT.median','rIMT.hom','rIMT.A'...           % 4 par
    'IMTmax','TBD','TBD'};                                        % 4 par
vlog(1,9)=num2cell(vArt.rtop(1));
vlog(1,[10,14,18,22,26,30,34])={'diastole','systole','distension',...
    'relative distension(%)','risetime (ms)','IMT','relative IMT(%)'};
fn=[];
if (length(vArt.fn)>4)&&(strcmp(vArt.fn(2:4),'000'))
    str=['Replace current file name "',vArt.fn,'" in folder "',...
        cd,['" by a unique name (for current folder) ',...
        'containing patient-ID']];
    while isempty(fn)
        fn=char(inputdlg_check(str, 'Name modification'));
        str='You must provide a discriminating file name, dummy!';
    end
else
    fn=vArt.fn;
end
vlog(2:5,1)={'file';fn;fn;fn};
for row=6:6+ncycl-1
    vlog{row,8}=sprintf('beat %d',row-5);
end
vlog(3:5,8)={'mean';'median';'std'};
vlog(1:11,4)={'source';'fps(Hz)';'nfr';'rec.time';'Depth';'Width';...
    'PixSize';'blksize';'Nrblks';'blkshift';'RiseThreshold'};
vlog{1,5}=vArt.source;
if vArt.pixscal>0
    pixscal=1000*vArt.pixscal;% express distance scale as micrometer/pixel
else                
    pixscal=1;
end
vlog(2:11,5)=num2cell([round(10*vArt.fps)/10;vArt.nfr;...
    round(10*vArt.nfr/vArt.fps)/10;vArt.rDepth;vArt.rWidth;round(pixscal);...
    round(10*vArt.wwidth)/10;vArt.nblk;round(10*vArt.bshift)/10;vArt.rtthr(1)]);
vlog(1:8,6)={'dres';'lres';'interpol';'Wall Thr.';'Intima Thr.';...
    'sm.window(s)';'ECG-cycles';'ECG-bpm'};
vlog(1:8,7)=num2cell([round(10*[vArt.dres;vArt.lres])/10;vArt.dip;...
    round(100*[vArt.wthr;vArt.ithr;vArt.sgwind])/100;ncycl;...
    round(60*vArt.fps*ncycl/(vArt.rtop(end)-vArt.rtop(1)))]);
vlog(9,6:7)={'Acq.date',vArt.acqdate};
clcdate=clock;
dat=sprintf('%d/%02d/%0d', clcdate(1), clcdate(2),clcdate(3));
vlog(10,6:7)={'proc.date',dat};
if vArt.proc>0
    imt=squeeze(imt(:,:,3));            % select only column with distance
    wdat=getWaveStat(sdistr,imt);
else
    wdat=getWaveStat(sdistr,[]);
end
if pixscal==1
    vlog{1,11}='Distances in pixels';
    p2m=[4,8,12,16,24];
    wdat(:,p2m)=round(100*wdat(:,p2m))/100;                     % 2 digits
    p2m=[1:3,5:7,9:11,13:15,17,22:23,25:29,30];
    wdat(:,p2m)=round(10*wdat(:,p2m))/10;                   % single digit
    p2m=(18:21);
    wdat(:,p2m)=round(wdat(:,p2m));                             % no digit
else
    vlog{1,11}='Distances in micrometer';
    p2m=[2:13,22:25,30];                           % convert to micrometer
    wdat(:,p2m)=round(wdat(:,p2m)*pixscal);
    wdat=round(100*wdat)/100;                % 2 digits for [4,8,12,16,24]
    p2m=[1,14:15,17,26:27,29];
    wdat(:,p2m)=round(10*wdat(:,p2m))/10;                   % single digit
    p2m=[2:3,5:7,9:11,13,18:25,30];
    wdat(:,p2m)=round(wdat(:,p2m));                             % no digit
end
vlog(3:5+ncycl,9:38)=num2cell(wdat);   
if exist(vArt.log,'file')                             % existing log-file
    bpxls=check_bpxls(vArt.log,'logsheet',vArt.fn,12);
    newfil=0;
else
    bpxls=1;
    newfil=1;
end
xlswrite(vArt.log,vlog,'logsheet',sprintf('A%d',bpxls));
if newfil==1                % this is the first time for current directory
    ClearXLsheet(vArt.log);                         % remove empty sheets
    [~,abbr]=getInfo();
    xlswrite(vArt.log,abbr,'logsheet','AP1');  % write list abbreviations
end
Save2txt(vlog);                                      % save to ASCII file
if ishandle(wah)
    delete(wah)
end
function bpxls=check_bpxls(xlfn,xlsheet,fn,xlrow)  % check begin position
% Finds begin position current ExceL sheet to append data
bpxls=1; % default begin position
if exist(xlfn,'file')                                 % Excel-file exists
    [~,desc]=xlsfinfo(xlfn);
    if ~isempty(find(strcmp(xlsheet,desc))>0)   % contains sheet 'xlsheet'
        [~,txt]=xlsread(xlfn,xlsheet);
        ffnd=find(strcmp(txt(:,1),fn)==1);
        if ~isempty(ffnd)            % file 'fn' has been processed before
            bpxls=xlrow*floor(ffnd(end)/xlrow)+1;           % replace data
        else
            bpxls=length(txt(:,4))-find(cellfun('isempty',...
                flipud(txt(:,4)))==0,1);
            bpxls=xlrow*ceil(bpxls/xlrow)+1;                 % append data
        end
    end
end   
end
end

function Save2txt(vlog)
% Saves alle relevant data to an ASCII-file for conversion to a data-base
% Arnold Hoeks, Biomedical Engineering, May 2013
global vArt
ncell=126;
txtcell=cell(2,ncell);
txtcell(1:2,1)={'vidArt';vArt.version};
txtcell(1:2,2)={'Filename';vArt.fn};
splitp=regexp(vArt.fn,'_','split');
if (length(vArt.fn)>20 && length(splitp)==6)      % this is a ParisK file
    tmp=splitp{1}(end-2:end);
    txtcell(1:2,4)={'patientnr';tmp};       
    txtcell(1:2,5)=['ID';splitp(1)]; 
    txtcell(1:2,6)=['USmode';splitp(3)];
    txtcell(1:2,7)=['artery';splitp(4)];  
    txtcell(1:2,8)=['side';splitp(2)]; 
    txtcell(1:2,9)=['angle';splitp(5)]; 
    tmp=splitp{6}(1);
    txtcell(1:2,10)={'duplo';tmp}; 
else
    txtcell(1:2,4)={'patientnr';'NA'};       
    txtcell(1:2,5)={'ID';'NA'}; 
    txtcell(1:2,6)={'USmode';'NA'};
    txtcell(1:2,7)={'artery';'NA'};  
    txtcell(1:2,8)={'side';'NA'}; 
    txtcell(1:2,9)={'angle';'NA'}; 
    txtcell(1:2,10)={'duplo';'NA'}; 
end
txtcell(1:2,11)={'Rise Threshold';num2str(vArt.rtthr(1))};
txtcell(1:2,12)={'File/processing properties';''};
txtcell(1:2,13)={'Source';vArt.source};
txtcell(1:2,14)={'fps(Hz)';num2str(vlog{2,5})};
txtcell(1:2,15)={'nfr';num2str(vlog{3,5})};
txtcell(1:2,16)={'rec.time(s)';num2str(vlog{4,5})};
txtcell(1:2,17)={'ROIdepth';num2str(vlog{5,5})};
txtcell(1:2,18)={'ROIwidth';num2str(vlog{6,5})};
txtcell(1:2,19)={'PixSize';num2str(vlog{7,5})};
txtcell(1:2,20)={'blksize(l)';num2str(vlog{8,5})};
txtcell(1:2,21)={'Nrblks';num2str(vlog{9,5})};
txtcell(1:2,22)={'blkshift';num2str(vlog{10,5})};
txtcell(1:2,23)={'dres';num2str(vlog{1,7})};
txtcell(1:2,24)={'res';num2str(vlog{2,7})};
txtcell(1:2,25)={'dinterpol';num2str(vlog{3,7})};
txtcell(1:2,26)={'WallThr';num2str(vlog{4,7})};
txtcell(1:2,27)={'IntThr';num2str(vlog{5,7})};
txtcell(1:2,28)={'sm.wind(s)';num2str(vlog{6,7})};
txtcell(1:2,29)={'ECG-cycle';num2str(vlog{7,7})};
txtcell(1:2,30)={'ECG-bpm';num2str(vlog{8,7})};
txtcell(1:2,31)={'acq.date';num2str(vlog{9,7})};
txtcell(1:2,32)={'proc.date';num2str(vlog{10,7})};
txtcell(1:2,33)={'Probe';vlog{3,3}};
txtcell(1:2,36)={'Artery Characteristics';''};
txtcell(1:2,37)={'RR.tmn';num2str(vlog{3,9})};
txtcell(1:2,38)={'RR.tmed';num2str(vlog{4,9})};
txtcell(1:2,39)={'RR.tsd';num2str(vlog{5,9})};
% diameter spatial mean; temporal mean/median/sd
txtcell(1:2,40:42)=BeatStat(squeeze(vlog(3:5,10)),'Ddiast.smn');% diastole
txtcell(1:2,43:45)=BeatStat(squeeze(vlog(3:5,14)),'Dsyst.smn');  % systole
txtcell(1:2,46:48)=BeatStat(squeeze(vlog(3:5,18)),'Ddist.smn');%distension
% diameter spatial median; temporal mean/median/sd
txtcell(1:2,49:51)=BeatStat(squeeze(vlog(3:5,11)),'Ddiast.smed');%diastole
txtcell(1:2,52:54)=BeatStat(squeeze(vlog(3:5,15)),'Dsyst.smed'); % systole
txtcell(1:2,55:57)=BeatStat(squeeze(vlog(3:5,19)),'Ddist.smn');%distension
% diameter spatial homogeneity; temporal mean/median/sd
txtcell(1:2,58:60)=BeatStat(squeeze(vlog(3:5,12)),'Ddiast.shom');%diastole
txtcell(1:2,61:63)=BeatStat(squeeze(vlog(3:5,16)),'Dsyst.shom'); % systole
txtcell(1:2,64:66)=BeatStat(squeeze(vlog(3:5,20)),'Ddist.shom');%distension
% spatial average diameter waveform; temporal mean/median/sd
txtcell(1:2,67:69)=BeatStat(squeeze(vlog(3:5,13)),'Ddiast.A');  % diastole
txtcell(1:2,70:72)=BeatStat(squeeze(vlog(3:5,17)),'Dsyst.A');   % dsystole
txtcell(1:2,73:75)=BeatStat(squeeze(vlog(3:5,21)),'Ddist.A'); % distension
% relative distension spatial mean/median/sd; temporal mean/median/sd
txtcell(1:2,76:78)=BeatStat(squeeze(vlog(3:5,22)),'rDist.smn');
txtcell(1:2,79:81)=BeatStat(squeeze(vlog(3:5,23)),'rDist.smed');
txtcell(1:2,82:84)=BeatStat(squeeze(vlog(3:5,24)),'rDist.shom');
txtcell(1:2,85:87)=BeatStat(squeeze(vlog(3:5,25)),'rDist.A');
% risetime spatial mean/median/std; temporal mean/median/sd
txtcell(1:2,88:90)=BeatStat(squeeze(vlog(3:5,26)),'RT.smn');
txtcell(1:2,91:93)=BeatStat(squeeze(vlog(3:5,27)),'RT.smed');
txtcell(1:2,94:96)=BeatStat(squeeze(vlog(3:5,28)),'RT.shom');
txtcell(1:2,97:99)=BeatStat(squeeze(vlog(3:5,29)),'RT.A');
% IMT spatial mean/median/std; temporal mean/median/sd
txtcell(1:2,100:102)=BeatStat(squeeze(vlog(3:5,30)),'IMT.smn');
txtcell(1:2,103:105)=BeatStat(squeeze(vlog(3:5,31)),'IMT.smed');
txtcell(1:2,106:108)=BeatStat(squeeze(vlog(3:5,32)),'IMT.shom');
txtcell(1:2,109:111)=BeatStat(squeeze(vlog(3:5,33)),'IMT.A');
% relative IMT spatial mean/median/std; temporal mean/median/sd
txtcell(1:2,112:114)=BeatStat(squeeze(vlog(3:5,34)),'rIMT.smn');
txtcell(1:2,115:117)=BeatStat(squeeze(vlog(3:5,35)),'rIMT.smed');
txtcell(1:2,118:120)=BeatStat(squeeze(vlog(3:5,36)),'rIMT.shom');
txtcell(1:2,121:123)=BeatStat(squeeze(vlog(3:5,37)),'rIMT.A');
% and finally maximum IMT; temporal mean/median/sd
txtcell(1:2,124:126)=BeatStat(squeeze(vlog(3:5,38)),'IMT.max');

%  Now create long string with the contents of a row concatenated
txtstr=txtcell{1,1};
ncell=length(txtcell(1,:));
for entry=2:ncell
    txtstr=[txtstr,', ',txtcell{1,entry}];
end
txtstr=[txtstr,'\n'];
[~,fn,~]=fileparts(vArt.fn);
fn=[fn,'.txt']; % replace file extension
fp=fopen(fn,'wt');
fprintf(fp,txtstr);
fclose(fp);
txtstr=txtcell{2,1};
for entry=2:ncell
    txtstr=[txtstr,', ',txtcell{2,entry}];
end
fp=fopen(fn,'a');
fprintf(fp,txtstr);
fclose(fp);
end

function dattxt=BeatStat(dat,str)
dattxt=cell(2,3);
dattxt(1:2,1)={[str,'.tmn'];num2str(dat{1})}; 
dattxt(1:2,2)={[str,'.tmed'];num2str(dat{2})};
dattxt(1:2,3)={[str,'.tsd'];num2str(dat{3})};
end

function adiam=getMWave(fh,sdistr,pixscal)
global vArt
% Gets, evaluates mean distension wave and initiates detection cardiac
% cycle ("R-top"). R-top position is added to 'ddiam', 2nd column
tbase=((1:vArt.nfr)-1)/vArt.fps;                              % time base
adiam=mean(sdistr);
subplot(fh);
hold off;
plot(tbase,adiam*pixscal,'k');
title('Average distension (edge waveforms excluded)');
if isempty(vArt.rtop)   % no ECG-triggers included, evaluate mean diameter
    tmp=linefit(adiam,1,1);    % first order regression line; no exclusion
    ECGtr=getRtop(adiam-tmp);
    vArt.rtop=find(ECGtr>0);                            % R-top positions
end
hold on
plot(tbase(vArt.rtop),adiam(vArt.rtop)*pixscal,'o','MarkerSize',5,...
    'color','r','MarkerFaceColor','w');
axis tight
if vArt.pixscal>0
    ylabel('Mean Diameter (mm)');
else
    ylabel('Mean Diameter (pixel)');
end
xlabel('Video time (s)');
end

function [ln,ind,rms]=linefit(data,acclevel,ord)
% function evaluates regression line of order 'ord' through 'data' with 
% outlier exclusion in a recursive way (maximum 3 loops) where eventually 
% the fraction 'acclevel' of the  datapoints will be retained
% The result is returned in 'ln' with 'ind' containg the index of data
% points that are retained and 'rms' the error of the accepted data-points
% and the regression line through them
% Arnold Hoeks, Biomedical Engineering, Maastricht University, april 2011
npt=length(data);
if acclevel*npt<5  % insufficient number of data-points
    ln=[];
    ind=[];
    rms=0;
    return
end
% intermediate results are stored in 'regr(:,1:4)' with
% regr(:,1) index
% regr(:,2) data
% regr(:,3) regression line
% regr(:,4) deviation between data and regression line
regr=zeros(npt,4);
regr(:,1)=1:npt;
regr(:,2)=data;
discard=max(1,round((1-acclevel)*npt/3));
while 2>0
    [pol, S]=polyfit(regr(:,1),regr(:,2),ord); 
    regr(:,3)=polyval(pol,regr(:,1));           % generate regression line
    regr(:,4)=regr(:,2)-regr(:,3);                             % deviation
    rms=std(regr(:,4));                   % compute root-mean-square error
    regr=sortrows(regr(:,:),4);                        % sort on deviation
    if length(regr(:,1))<=acclevel*npt
        break
    end
    regr=regr(1+discard:length(regr(:,1))-discard,:);   % discard extremes
    regr=sortrows(regr(:,:),1);                        % sort on deviation
end
ind=sort(regr(:,1));        % remaining accepted values in ascending order
ln=polyval(pol,1:npt);                          % generate regression line
end


function gfn=GenFilName(fn,fext,mode)
% Generates filename including index number for files with wildcard 
% extension 'fext' (e.g. '*.jpg') to avoid that file is
% overwritten; allows therefore repeated save calls. If mode==1, the
% extension is attached to the filename, otherwise it is omitted.
flist=ls(fext);
[~,lfn]=fileparts(fn);
if isempty(flist)   % no conflict
    gfn=lfn;
else
    nmatch=0;
    [nfil,~]=size(flist);
    for fil=1:nfil
        if ~isempty(strfind(flist(fil,:),lfn))
            nmatch=nmatch+1;
        end
    end
    if nmatch==0
        gfn=lfn;
    else  
        gfn=[lfn,sprintf('_%d',nmatch+1)];
    end
end
if mode==1
    pp=strfind(fext,'.'); % find begin extension
    gfn=[gfn fext(pp:end)];
end
end


function ECGtr=getRtop(pdiam)
npt=length(pdiam);
ECGtr=zeros(1,npt);
tmp=abs(fft(pdiam));
[~,mind]=max(tmp(1:round(end/2)));
RR=round(npt/(mind-1));       % mean RR-interval in sample points (frames)
bp=1;
wind=round(0.7*RR);
wind2=round(0.8*wind);
[~,mind]=min(pdiam(bp:bp+wind));
if mind<3||(mind==1+wind)  % search for first local minimum to start from
    bp=bp+wind-1;
else
    ECGtr(mind-2)=1;
    bp=mind+wind;
end
while bp<npt-wind
    [~,mind]=min(pdiam(bp:bp+wind));
    if mind>wind-1                                 % shorten search window
        [~,mind]=min(pdiam(bp:bp+wind2));                 % wind2<0.8*wind
        if diff(pdiam(bp+mind:bp+mind+1)>0)% apparently an earlier minimum
            ECGtr(bp-3+mind)=1;
            bp=bp-3+mind+wind2;
        elseif bp<npt-wind-wind2-1        % continue with a forward search
            bp=bp+wind;
            [~,mind]=min(pdiam(bp:bp+wind2));
            ECGtr(bp-3+mind)=1;
            bp=bp-3+mind+wind2;
        else
            bp=bp-3+mind+wind2;
        end
    else
        ECGtr(bp-3+mind)=1;
        bp=bp-3+mind+wind;
    end
end
if bp<npt                          % search for the trailing local minimum
    [~,mind]=min(pdiam(bp:end));
    if mind+bp-1>bp&&mind+bp-1<npt
        ECGtr(mind+bp-3)=1;
    end
end
end

function getMan_ECG(mdiam,hsub,pixscal)
global vArt
subplot(hsub);
axlim=get(gca,'Position');                                  % normalized!
fglim=get(gcf,'Position');                                   % in pixels!
x1=axlim(1)*fglim(3)+fglim(1); 
x2=(axlim(1)+axlim(3))*fglim(3)+fglim(1);
y1=axlim(2)*fglim(4)+fglim(2);
y2=(axlim(2)+axlim(4))*fglim(4)+fglim(2);
xlabel('Identify/change onset cardiac cycles.'); 
tbase=((1:vArt.nfr)-1)/vArt.fps;                              % time base
hold off
mdiam=mdiam*pixscal;
plot(tbase,mdiam)
hold on
axis tight
str=sprintf('Add/remove R-top events\nClick "Continue" if satisfied');
title(str);
xlabel('Video time [s]');
xl=get(hsub,'XLim');
yl=get(hsub,'YLim');
set(gcf,'WindowButtonMotionFcn',@changepointer,...
    'WindowButtonDownFcn',@getPPos,'WindowButtonUpFcn',[]);
ch=uicontrol('style','pushbutton','Units','normalized','position',...
    [0.87,0.01,0.1,0.035],'String','Accept/Continue','callback',@doFinish);
fin=0;
while fin==0 
    if ~isempty(vArt.rtop)
        ind=length(vArt.rtop);
        eh=zeros(1,ind);
        for tr=1:ind
            eh(tr)=plot(tbase(vArt.rtop(tr)),mdiam(vArt.rtop(tr)),'o',...
                'MarkerSize',5,'color','r','MarkerFaceColor','w');
        end
    else
        ind=0;
    end
    neh=ind;
    uiwait()
    if (gca==hsub)&&(fin==0)
        prind=find(abs(vArt.rtop-xp)<0.3*vArt.fps);% distance 2 previous <0.3s?
        if (ind==0)||(isempty(prind))
            ind=ind+1;                                       % extra R-top
            vArt.rtop(ind)=xp; % units pixels!!
        elseif ~isempty(prind)          % remove previous marker frim list
            if length(vArt.rtop)>1
                    vArt.rtop=[vArt.rtop(1:prind-1),...
                        vArt.rtop(prind+1:end)];
            else
                vArt.rtop=[];
            end
        end
        if neh>0
            delete(eh);
        end
    end
    vArt.rtop=sort(vArt.rtop);
end
set(gcf,'WindowButtonMotionFcn',[],'WindowButtonDownFcn',[],...
    'WindowButtonUpFcn',[]);              % deactivate callback functions
title('Average distension (edge waveforms excluded)'); 
delete(ch);
function changepointer(~,~)
    pntr=get(0,'PointerLocation');         % returns coordinates in pixels
    if pntr(1)>x1&&pntr(1)<x2&&pntr(2)>y1&&pntr(2)<y2 % within subwindow??
        set(gcf,'Pointer','fullcrosshair')
    else
        set(gcf,'Pointer','arrow')
    end
end
function getPPos(~,~)
     pntr=get(hsub,'CurrentPoint');   % returns coordinates in plot units
     if pntr(1,1)>xl(1)&&pntr(1,1)<xl(2)&&pntr(1,2)>yl(1)&&pntr(1,2)<yl(2)
        xp=round(pntr(1,1)*vArt.fps+1);% only x-coordinate is of interest
        uiresume
     end
end
function doFinish(~,~)
    fin=1;
    uiresume;
end
end


function res=getWaveStat(dwave,imt)
% Extracts statistics from the set of distension waves contained in
% 'dwave', organized as [blk,time], with edge signals removed.
% 'res' contains the following info for each complete cardiac cycle:
% RR (ms), diast.(mean,median,sd), syst.diam(mean,median,sd), 
% pp.dist(mean,median,sd), risetime(mean,median,sd) in ms. The same
% morphological parameters are extracted for the average waveform.
% Moreover, the relative distension (%), the IMT (mean, media, std) 
% and the relative wall thickness are attached: total 29 columns
% Each column  is headed by the mean, median and standard deviation 
% of the rows below (beats).
global vArt
dwave=dwave(vArt.delseg+1:vArt.nblk-vArt.delseg,:);         % remove edges
if ~isempty(imt)
    imt=imt(vArt.delseg+1:vArt.nblk-vArt.delseg,:);         % remove edges
end
[nsig,~]=size(dwave);
ncycle=min(6,vArt.necg-1);          % number of complete cycles, maximal 6
value=zeros(nsig,ncycle,7);               % to store temp. wave parameters
ncol=30;
res=zeros(3+ncycle,ncol);
RRsp=mean(diff(vArt.rtop));            % mean RR interval in sample points
syswind=ceil(RRsp/2);          % crudely estimated duration systolic phase
for beat=1:ncycle                                % process segment signals
    res(3+beat,1)=1000*(vArt.rtop(beat+1)-vArt.rtop(beat))/vArt.fps;  % ms
    bp=vArt.rtop(beat);                              % location begin beat
    sigs=[];
    for sig=1:nsig 
        bsegm=dwave(sig,bp:vArt.rtop(beat+1));              % beat segment
        [~,indx]=max(bsegm(1:syswind));                          % systole
        ssegm=bsegm(1:indx);                      % early systolic segment
        value(sig,beat,2)=max(bsegm(1:syswind));                 % systole
        value(sig,beat,1)=min(bsegm(1:indx));                   % diastole
        value(sig,beat,3)=value(sig,beat,2)-value(sig,beat,1);      % dist
        value(sig,beat,4)=100*value(sig,beat,3)/value(sig,beat,1);%reldist
        value(sig,beat,5)=getXT(ssegm,vArt.rtthr(2))-getXT(ssegm,vArt.rtthr(1));
        if ~isempty(imt)
            value(sig,beat,6)=imt(sig,beat);                         % imt
            value(sig,beat,7)=100*value(sig,beat,6)/value(sig,beat,1);
        end
        value(sig,beat,8)=(3*value(sig,beat,1)^2)/(8*value(sig,beat,6)*value(sig,beat,3));%E modulus
        if value(sig,beat,4)>=vArt.MinRelDis
            sigs=[sigs,sig];
        end
    end
    for par=1:7
        col=4*(par-1)+1;
        res(beat+3,col+1)=mean(value(sigs,beat,par));
        res(beat+3,col+2)=median(value(sigs,beat,par));
        res(beat+3,col+3)=std(value(sigs,beat,par));
    end
end
for beat=1:ncycle                       % process spatial average waveform
    bp=vArt.rtop(beat);                              % location begin beat
    bsegm=mean(dwave(:,bp:vArt.rtop(beat+1)),1);            % average wave
    [~,indx]=max(bsegm(1:syswind));                              % systole
    ssegm=bsegm(1:indx);                          % early systolic segment
    res(beat+3,9)=max(bsegm(1:syswind));                         % systole
    res(beat+3,5)=min(bsegm(1:indx));                           % diastole
    res(beat+3,13)=res(beat+3,9)-res(beat+3,5);              % distension
    res(beat+3,17)=100*res(beat+3,13)/res(beat+3,5);  % rel. distension(%)
    res(beat+3,21)=getXT(ssegm,vArt.rtthr(2))-getXT(ssegm,vArt.rtthr(1));
    if ~isempty(imt)
        res(beat+3,25)=mean(imt(:,beat));                           % aIMT
        res(beat+3,29)=100*res(beat+3,25)/res(beat+3,5);   % rel. aIMT (%)
        res(beat+3,30)=max(imt(:,beat));                            % aIMT
    end
end
for col=1:ncol
    res(1,col)=mean(res(4:end,col));               % temporal (beats) mean
    res(2,col)=median(res(4:end,col));           % temporal (beats) median
    res(3,col)=std(res(4:end,col));       % temporal (beats) standard dev.
end
end


function xt=getXT(segm,xl)
% Function determines time at which "segm" crosses "xl" percent level 
% relative to peak-peak value using linear interpolation. "xt" is true 
% time, expressed in msec,  relative to the start of the waveform, 
% irrespective of the actual position of the diastolic minimum.
% The time difference between the upper and lower level crossings 
% (in procent) determines the rise-time (requires 2 calls).
global vArt
[syst,tsyst]=max(segm);  % systole
[dias,tdias]=min(segm(1:tsyst));  % diastole
thr=dias+xl*(syst-dias)/100; % threshold level
xt=find(segm(tdias:tsyst)>=thr,1)+tdias-1; % first relevant level crossing
if isempty(xt)||(xt<=1)||(segm(xt)-segm(xt-1)==0)
    xt=round(1000*(tdias+xl*(tsyst-tdias)/100)/vArt.fps);  
else
    xt=xt+(thr-segm(xt))/(segm(xt)-segm(xt-1));   % linear interpolation
    xt=1000*xt/vArt.fps;                               % convert to msec
end
end


function [wpos,err,cblk]=set_Broi(ima,ini)
global vArt modstop        % 'modstop' is set in 'mod_MarkPos' (if 2 redo)
[vArt.iDepth,vArt.iWidth]=size(ima);                      % raw image size
cblk=[];
[fh,~]=GenSubPlot(1,vArt.fn,vArt.FOS);
imagesc(ima);
colormap gray; 
if strcmp(vArt.source,'*.avi')
    getPixScal();                           % pixel size expressed in CM!!
end
vArt.wexpl=round(vArt.mwexpl/vArt.pixscal);% exploration range mm to pixel
str=sprintf('Depth scaling =%d um/pixel',round(1000*vArt.pixscal));
if vArt.pixscal>0
    str=sprintf('%s (%d pixel/mm)',str,round(1/vArt.pixscal));
end
title([str,'; Frame Rate=',num2str(vArt.fps),' Hz']);
xl=get(gca,'Xlim');
yl=get(gca,'Ylim');
if isempty(vArt.frbwind)||(sum(vArt.frbwind([1,3]))>xl(2)||sum(vArt.frbwind([2,4]))>yl(2))
    vArt.frbwind=round([xl(1)+(xl(2)-xl(1))/4,yl(1)+(yl(2)-yl(1))/5,...
        0.5*(xl(2)-xl(1)),0.5*(yl(2)-yl(1))]);
end
vArt.dblgap=round(vArt.dip*vArt.dres*.75/vArt.pixscal);    % gap echo
while true                % cycle until window size is sufficiently large
    th=uicontrol('style','text','units','normalized','string',...
        'Drag/resize rectangle to cover artery and click Accept',...
        'position',[0.2 0.015 0.55 0.035]);
    bh=uicontrol('style','pushb','units','normalized','string',...
        'Accept/Continue','Position',[0.77 0.01 0.2 0.05],...
        'callback',@doAccept);
    flag=1;
    hrect=imrect(gca,vArt.frbwind);
    % api.setResizable(false); % fixed size, comment for readjustment
    api=iptgetapi(hrect);
    while true
        pause(0.2);
        vArt.frbwind=round(api.getPosition());
        if flag<0
           break
        end
    end    
    delete(bh);delete(th);
    vArt.frbwind=round(api.getPosition());
    api.delete();                                       % remove rectangle
    if ini==0;                            % only readjustment position ROI
        close(fh);
        return
    end 
    [wpos,err]=IniLP(imcrop(ima,vArt.frbwind)); % rect=[xp yp xn yn]
    if err<0
        limmess(['Proper positioning ROI impossible. ',...
            'File processing aborted.'],2);
        close(fh)
        return
    elseif err==0                             % everything seems to be OK
        break;
    end
end
close(fh);
[fh,~]=GenSubPlot(1,vArt.fn,vArt.FOS);      % check initilization 
imagesc(imcrop(ima,vArt.frbwind));
colormap gray;
cblk=getSubROI();
wpos=wpos(:,cblk(2,:));
im1=imcrop(ima(:,:),vArt.frbwind);
tmp=zeros(4*vArt.wexpl+1,vArt.nblk);
wspan_pix=vArt.wspan/vArt.pixscal;
shift=vArt.bshift*vArt.wwidth;
wspan=round((wspan_pix+shift-vArt.wwidth)/shift);
wspan=max(3,wspan);       % minimal smoothing window for 2th order filter
modstop=2;
while modstop==2
    for wall=1:2
        for ln=1:vArt.nblk                                      % skewed area!
            bp=round(wpos(wall,ln,1)-2*vArt.wexpl);
            tmp(:,ln)=mean(im1(bp:bp+4*vArt.wexpl,cblk(1,ln):cblk(3,ln)),2);
        end
        wpos(wall,:,1)=wpos(wall,:,1)+wdet(tmp,wall,0);  
        wpos(wall,:,1)=smooth(wpos(wall,:,1)',wspan/vArt.nblk,'rloess')';
                           % update wall
    end
    [wpos,~]=mod_MarkPos(wpos,[],cblk(2,:),...
        'New initialization wall position'); % function modifies 'modstop'
end
close(fh);
function doAccept(~,~)  
    if vArt.frbwind(3)>vArt.MinRoiWidth
        vArt.frbwind(3)=2*round(vArt.frbwind(3)/2)+1;  % ensure odd width
        flag=-1;
    else
        beep;
        title(['Select wider window (>',num2str(vArt.MinRoiWidth),...
            ' pixels)']);
        pause(0.3);
        beep;
    end
end    
end
 
function [wpos,err]=IniLP(ima)
% Function coordinates ant/post wall identification and returns the 
% results in "wpos(wall,ln)". 
% 'err'<0 -> impossible to set properly ROI size/position
% 'err'=0 -> ROI size/position accepted
% 'err'>0 -> resize ROI
% Arnold Hoeks, Biomedical Engineering, Maastricht University, May 2013
global vArt
[vArt.rDepth,vArt.rWidth]=size(ima);
[fh,~]=GenSubPlot(1,vArt.fn,vArt.FOS);                  % create new image
hima=imagesc(ima);
colormap gray;
xlabel(sprintf(['Estimated resolution %4.1f (depth) by %4.1f ',...
    '(lateral) pixels'],round([vArt.dres,vArt.lres]/vArt.pixscal)));
posax=get(gca,'position');
set(gca,'position',[posax(1) posax(2)+0.05 posax(3) posax(4)-0.05]);
set(gcf,'WindowButtonMotionFcn',@changepointer); 
set(gca,'ButtonDownFcn',@getpoints);
hold on                                         % to add points to figure
set(hima,'Hittest','off');           % imagesc resets figure properties!!
pcnt=0;
inip=zeros(2,4);               % coordinates of identified wall positions
th=uicontrol('style','text','units','normalized','position',...
    [0.2 0.015, 0.6 0.035],'FontSize',vArt.FOS,'string',...
    'Select twice (left/right) anterior/posterior wall.');
while pcnt<4
    pause(0.2);
end
delete(th);
inip=sortrows(inip')';
for wxp=2:2:4           % organize as ant/post (left) and ant/post (right)
    if inip(2,wxp)<inip(2,wxp-1) 
        tmp=inip(2,wxp);
        inip(2,wxp)=inip(2,wxp-1);
        inip(2,wxp-1)=tmp;
    end
end
wpos=zeros(2,vArt.rWidth);
for wall=1:2        % fit straight lines through anterior/posterior points
    xp1=inip(1,wall);
    yp1=inip(2,wall);
    xp2=inip(1,wall+2);
    yp2=inip(2,wall+2);
    slope=(yp1-yp2)/(xp1-xp2); 
    offs=(yp2*xp1-yp1*xp2)/(xp1-xp2);
    for ln=1:vArt.rWidth
        wpos(wall,ln,1)=slope*(ln-1)+offs;
    end 
end
cdev=round((vArt.frbwind(4)-(max(wpos(2,:))+min(wpos(1,:))))/2);%deviation
wpos=wpos+cdev;                                             % center lumen
vArt.frbwind(2)=vArt.frbwind(2)-cdev;               % adjust ROI to center
vArt.wexpl=min(vArt.wexpl, round(mean(wpos(2,:)-wpos(1,:))/vArt.wexplDiamScaleFactor));% scale to diameter
thr=4*round(vArt.wexpl);
vshift=min(0,min(wpos(1,:))-thr); % check margin and expand ROI (vshift<0)
vArt.frbwind(2)=vArt.frbwind(2)+vshift;           % shift image vertically
vArt.frbwind(4)=vArt.frbwind(4)-2*vshift;                   % expand image
wpos=wpos-vshift;
err=0;
vArt.frbwind=round(vArt.frbwind);
if vArt.frbwind(2)<1
    err=-1;
end
delete(fh);
function getpoints(hObj,~)
    pcnt=pcnt+1;
    cp = get(hObj,'CurrentPoint');
    xp=cp(1,1);
    yp=cp(1,2);
    inip(1,pcnt)=xp;
    inip(2,pcnt)=yp;
    plot(xp,yp,'o','MarkerSize',5,'color','r','MarkerFaceColor','w');
end
function changepointer(~,~)  % controls crosshair/arrow appearance pointer
    [ppOK,~,~]=CheckPP();
    if ppOK                                                 % within axes
         set(gcf,'Pointer','crosshair')
    else
         set(gcf,'Pointer','arrow')
    end
end
end

function wdev=wdet(ima,wall,disp)
% Function to detect lumen wall transition in image 'ima' for either 
% anterior or posterior wall, starting from lumen center
% (anterior section is flipped). Realize that for a skewed.curved center
%  line, the selected part of the image is reshaped accordingly such that 
% 'ima' is always centered symmetrically around the fit through
% the wall-media position. 
% The detection is performed on average echo level across the image 
% of segments of the image as specified in "cblk" within 2 times the
% units for the exploration range (with wall-lumen at the center).
% Detected displacements are limited to one standard deviation from mean.
% The result of edge detection, based on crossing relative to 
% the local maximum is returned in 'wdev' for all vertical image lines 
% within 'ima', containing only the DEVIATION from the centerline; 
global vArt hsub                           % 'hsub' is handle of subimage
if wall==1                              % for anterior wall flip the image
    ima=flipud(ima);
end
[rDepth,~]=size(ima);                   % 'rDepth' is 4 times 'vArt.wexpl'
base=mean(mean(ima(1:3,:)));                      % baseline (lumen) value 
% atop=mean(max(ima(vArt.wexpl:3*vArt.wexpl,:)));         % mean top value
bp=vArt.dip*vArt.wexpl;                 % begin-position exploration range
ep=vArt.dip*(rDepth-vArt.wexpl);             % end position analysis range
cp=round(length(1:1/vArt.dip:rDepth)/2);% center of points after interpol.
wpos=zeros(1,vArt.nblk);               % detected wall positions per block
imai=interp1(1:rDepth,ima,1:1/vArt.dip:rDepth); % depth interpolated image 
gtop=0;                                      % global top adventitia value
shift=round(0.3*cp);
for blk=1:vArt.nblk
    Asig=imai(:,blk);
    [val,~]=max(Asig(bp:ep));         % local maximum value, peak position
    gtop=max(gtop,val);                    % absolute top value adventitia
    twp=find(Asig(bp:ep)-base>vArt.wthr*(val-base),1);% threshold wall pos
    if isempty(twp)
        wpos(blk)=cp+shift;                                % move outwards
    else
        twp2=find(Asig(bp:ep)-base>vArt.wthr2*(val-base),1);% double echo?
        if twp2-twp>vArt.dblgap        % minimum gap between double echoes
            twp=twp2;                                  % skip leading echo
        end
        wpos(blk)=twp+bp;                    % position threshold crossing
    end
end
% Limit detected displacements to one standard deviation from mean
ave=mean(wpos);                                                  % average
sd=std(wpos);                                         % standard deviation
wpos=min(ave+sd,max(ave-sd,wpos))/vArt.dip;             % limit excursions
if disp==wall  % presentation of either anterior or posterior wall results
    subplot(hsub(1));
    hold off
    imagesc(ima);
    colormap gray
    hold on
    plot(wpos,'Linewidth',3);
    title('Echo with detected wall-lumen transition');
    xlabel('Echo line');
    ylabel('Depth (pixels)');
    subplot(hsub(2));
    hold off
    plot(ima(:,:),'LineWidth',2,'LineStyle','-');
    hold on
    for blk=1:vArt.nblk
        plot(wpos(blk),imai(round(wpos(blk)*vArt.dip),blk),'o',...
            'MarkerSize',7,'color','k','MarkerFaceColor','k');
    end
    plot([vArt.wexpl,3*vArt.wexpl],[0,0],'LineWidth',5,'color','r');
    axis tight
    ylim(gca,[0 1.1*gtop]);
    title('Average echo level(depth) with detected position',...
        'FontSize',vArt.FOS);
    xlabel('Depth (pixels); red line is exploration range');
    ylabel('Echo level');
    pause(0.1);
end
if wall==1    % remove offset (centerline position), retain only deviation
    wdev=2*vArt.wexpl-wpos;                % flip result for anterior wall
else
    wdev=wpos-2*vArt.wexpl;
end
end

function [fh,hsub]=GenSubPlot(npl,fn,FOS)
% Function creates figure with 'npl' subplots on a single (if npl=1:3) row
% or on a double row (f npl>3) with a modest margin. The function returns 
% the handle of  the figure ('fh') as well as the handles ('hsub') of the 
% subplots. The fontsize 'FOS' affects the axis labels and the title
if nargin<3
    FOS=10;                                           % default FontSize
end
fh=figure('Toolbar','Figure','units','normalized','NumberTitle','off',...
    'Name',['File: ',strrep(fn,'_','-')]);
fpos=[0.2 0.2 0.6 0.7];
oldscr=get(0,'units');
set(0,'units','pixels');
scr=get(0,'screensize');
set(0,'units',oldscr);
xm=80;
ym=50;
xpmarg=xm/(scr(3)*fpos(3));
ypmarg=ym/(scr(4)*fpos(4));
if npl==1
    set(gcf,'Position',fpos);
    set(gca,'position',[xpmarg,2*ypmarg,1-2*xpmarg,1-3*ypmarg],'FontSize',FOS)
    hsub=[];
    return
elseif npl==2    
    fpos(1)=0.15;
    fpos(3)=0.7;    
elseif npl==3
    fpos(1)=0.1;
    fpos(3)=0.8;
elseif npl==4  % distribute over 2 rows
    fpos=[0.15 0.04 0.7 0.9];
else
    fpos=[0.05 0.04 0.9 0.9];
end
set(gcf,'Position',fpos);
nrow=1;
ncol=npl;
if npl>3
    nrow=2;
    ncol=ceil(npl/2);
end 
xpmarg=xm/(scr(3)*fpos(3));
ypmarg=ym/(scr(4)*fpos(4));
ww=(1-(ncol+1)*xpmarg)/ncol;
wh=(1-(1.5*(nrow+1)*ypmarg))/nrow;% reserve space for uicontrol buttons
cnt=0;
axp=[1,1-wh-ypmarg,ww,wh];
for row=1:nrow
    for col=1:ncol
        cnt=cnt+1;
        axp(1)=xpmarg*col+(col-1)*ww;
        hsub(cnt)=subplot(nrow,ncol,cnt);
        set(gca,'Position',axp,'FontSize',FOS);
    end
    axp(2)=axp(2)-wh-1.5*ypmarg; 
end
end

function frsel=Get_avi_Frame(ima,fn)
% Function to select interactively a frame sequence from a movie sequence. 
% The selected  range of frames is returned in 'frsel'[begin,end] 
% Arnold Hoeks, March 2012
% Biomedical Engineering, Maastricht University
global vArt
[~,~,nfr]=size(ima);
frsel=[1 nfr];
[fh,~]=GenSubPlot(1,fn,vArt.FOS);
posax=get(fh,'Position');
imagesc(squeeze(ima(:,:,1)));
colormap gray
hcnt=uicontrol('style','text','units','normalized','Position',...
    [0.11 0.015 0.1 0.035]);       % set user interface control buttons
hbar=uicontrol('style','slider','units','normalized','Position',...
    [0.23 0.015 0.3 0.035],'callback',@DoFrame);
set(hbar,'SliderStep',[1/nfr 10/nfr]);
playbut=uicontrol('style','togglebutton','string','Pause','units',...
    'normalized','Position',[0.63 0.015 0.1 0.035],...
    'callback',@DoPlay);
pbval=get(playbut,'Value');
menustr={'Select Action','First Frame','Last Frame','Continue','Skip this file'};
finbut=uicontrol('style','popupmenu','string',menustr,'units',...
    'normalized','Position',[0.75 0.015 0.2 0.035],'callback',@DoStop);
fr=1;                % first frame
avistop=0;
while avistop==0                      % start cycling until frame accepted
    menustr{2}=['First Frame: ',num2str(frsel(1))];
    menustr{3}=['Last Frame: ',num2str(frsel(2))];
    set(finbut,'String',menustr);
    imagesc(squeeze(ima(:,:,fr)));
    colormap gray
    pause(0.02);
    set(hcnt,'string',[num2str(fr),'/',num2str(nfr)]);
    set(hbar,'Value',fr/nfr);
    while (pbval>0)&&(avistop==0)       % loop until playbutton is paused
        pause(0.2)
    end
    fr=fr+1;
    if fr>nfr
        fr=1;                                     % restart at first frame
    end
end
delete(gcf);
function DoFrame(~,~)                                     % slider control
    frpos=get(hbar,'Value');               % retrieve current frame number
    fr=min(nfr,max(1,round(frpos*nfr)));
    set(hcnt,'string',[num2str(fr),'/',num2str(nfr)]);
    set(hbar,'Value',fr/nfr);
    imagesc(squeeze(ima(:,:,fr)));
    colormap gray;
end
function DoPlay(~,~)                            % controls play/pause mode
    pbval=get(playbut,'Value');
    if pbval==1
        set(playbut,'string','Play');
    else
        set(playbut,'string','Pause');
    end
end
function DoStop(~,~)                           % frame identified! Finish
    choice=get(finbut,'Value');
    if choice==2
        frsel(1)=fr;  % begin frame video segment
        frsel(1)=min(frsel(1),frsel(2));
    elseif choice==3
        frsel(2)=fr;  % last frame video segment
        frsel(2)=max(frsel(1),frsel(2));
    elseif choice==4
        avistop=1;
    elseif choice==5
        frsel(:)=0;
        avistop=1;
    end
    menustr{2}=['First Frame: ',num2str(frsel(1))];
    menustr{3}=['Last Frame: ',num2str(frsel(2))];
    set(finbut,'String',menustr);
    set(finbut,'Value',1);
end
end


function rprt=ClearXLsheet(xlsfn,xltab)
% Function to clear sheet "xltab" within Excel workbook "xlsfn"
% In the process, all empty sheets will be removed as well!
% If only "xlsfn' is provided, then the function deletes empty sheets!
% The number od cleared/deleted sheets is reported in "rprt".
% Arnold Hoeks, Maastricht, BME, March 2012
if isempty(strfind(xlsfn,'\'))    % full path for actxserver('Excel.Appl') 
    xlsfn = [cd '\' xlsfn];                         % construct full path
end 
rprt.deltd=0;
rprt.clrd=0;
XL=actxserver('Excel.Application'); 
set(XL,'Visible', 0);                   % Make the application invisible
set(XL,'DisplayAlerts',0);                % Make excel not display alerts
Workbook=XL.Workbooks.Open(xlsfn); % Open an Excel Workbook and activate it
WorkSheets=XL.Sheets;                    % Get the sheets in the Workbook
XL.EnableSound=false;
nshwb=WorkSheets.count;                     % number of sheets in workbook
if nargin==1
    xltab=[];
end
for sh=nshwb:-1:1
    current_sheet=get(WorkSheets,'Item',sh);
    if strcmp(WorkSheets.Item(sh).Name,xltab)
        WorkSheets.Item(sh).UsedRange.Delete;              % clear sheet
        rprt.clrd=rprt.clrd+1;
    % worksheets.Item(sheetIdx).UsedRange.Count is number of used cells.
    % This will be 1 for an empty sheet. It may be 1 for certain other
    % cases but then function will beep and not actually delete sheet.
    elseif (WorkSheets.count>1)&&(WorkSheets.Item(sh).UsedRange.Count==1)
        invoke(current_sheet,'Delete');              % delete empty sheet
        rprt.deltd=rprt.deltd+1;
    end
end
XL.EnableSound=true;
Workbook.Save;                                     % Now save the workbook
XL.Workbooks.Close;                                % Close the workbook
invoke(XL,'Quit');       
delete(XL);                      % Delete the handle to the ActiveX Object
end


function [wposR,imt]=getIMT(ima,sdistr,wposR,imt,cblk)
% input format 'ima' and 'wposR(1:2,1:nblk,1:necg)', at R-top
% all input, processing and output are in pixels
% IMT(1:iLine, 1:necg,1) line index for intima 
% IMT(1:iLine, 1:necg,2) intima position for valid intima (sp)
% IMT(1:iLine, 1:necg,3) intima media thickness (IMT,sp)
% Arnold Hoeks, University Maastricht, April 2013.
global vArt modstop                    % 'modstop is used in 'mod_MarkPos'
modstop=0;
[ny,nx,~]=size(ima);
if (isempty(imt))||(vArt.imt_fb~=0)
    imt=zeros(vArt.nblk,vArt.necg,3);
    for beat=1:vArt.necg
        imt(:,beat,1)=1:vArt.nblk;
    end
    imtflag=0;
else
    imtflag=1;
    imt(:,:,2)=(imt(:,:,2)-1)*vArt.dip+1;
end
wposR(:,:,:)=(wposR(:,:,:)-1)*vArt.dip+1;
[fh,~]=GenSubPlot(1,vArt.fn,vArt.FOS);
sh=uicontrol('style','pushbutton','string','Skip editing','units',...
    'normalized','Position',[0.05 0.015 0.2 0.035],'callback',@DoSkip);
sk=0;
sh=uicontrol('style','slider','units','normalized','position',...
    [0.93 0.2 0.04 0.6],'Value',0.5,'callback',@Bslider);
imagesc(ima(:,:,1));
axh=gca;
for beat=1:vArt.necg                                   % explore all beats
    cstr=sprintf(['Interpolated image (factor %d) for beat ',...
        '%d of %d beats'],vArt.dip,beat,vArt.necg);
    hold off
    bima=interp2(double(ima(:,:,beat)),1:nx,(1:1/vArt.dip:ny)','spline');
    figure(fh)
    if sk==0                                  % edit/modify wall positions
        imagesc(bima);                                % interpolated image 
        colormap gray
        pause(0.05);
        hold on                                  
        Bslider(sh,[]);    
        [wposR(:,:,beat),~]=mod_MarkPos(wposR(:,:,beat),[],cblk(2,:),cstr);
    end
    if imtflag==0                  % explore all lines for intima position
        for ln=1:vArt.nblk                          % select relevant data
            cstr=sprintf(['Current line %d of %d image segments '....
                '(beat %d of %d beats)'],ln,vArt.nblk,beat,vArt.necg);
            ep=wposR(2,ln,beat)+5*vArt.dip;
            npt=floor((wposR(2,ln,beat)-wposR(1,ln,beat))/2);% expl. range
            bp=ep-npt;
            mline=squeeze(mean(bima(bp:ep,cblk(1,ln):cblk(3,ln)),2));         
            imt(ln,beat,2)=imt_ln(mline,wposR(2,ln,beat)-bp,cstr)+bp;
        end
    end                                                  
    cstr=sprintf(['Interpolated image (factor %d) for beat ',...
        '%d of %d beats'],vArt.dip,beat,vArt.necg);
    if sk==0                                             % edit/modify IMT
        [~,imt(:,beat,2)]=mod_MarkPos(wposR(:,:,beat),imt(:,beat,2)',...
            cblk(2,:),cstr);
    end
    imt(:,beat,3)=wposR(2,:,beat)'-imt(:,beat,2);
end
wposR(:,:,:)=(wposR(:,:,:)-1)/vArt.dip+1; % remove interpolation in depth
for beat=1:vArt.necg-1 % correct smoothed diameter distribution
    for ln=1:vArt.nblk
        bp=vArt.rtop(beat);
        ep=vArt.rtop(beat+1)-1;
        dev=sdistr(ln,bp)-(wposR(2,ln,beat)-wposR(1,ln,beat));
        sdistr(ln,bp:ep)=sdistr(ln,bp:ep)-dev;
    end
end
imt(:,:,2:3)=(imt(:,:,2:3)-1)/vArt.dip+1;
Save2XL(sdistr,imt);
close(fh);
function DoSkip(~,~)
    delete(sh);
    sk=1;                            % skip editing wall position and IMT
    modstop=1;                                 % interrupt editing WP/IMT
    uiresume
end
function Bslider(sh,~)
    slpos=get(sh,'Value');
    set(axh,'cLimMode','auto');
    [blim]=get(axh,'cLim');
    % if position of slider is below center increase lower limit
    if slpos <= 0.5
        blim_new=[-1+blim(1)+2*(0.5-slpos)*(blim(2)-blim(1)),blim(2)];
    % if position of slider is above center decrease upper limit
    else
        blim_new=[blim(1),1+blim(2)-2*(slpos-0.5)*(blim(2)-blim(1))];
    end
    if blim_new(1)<blim_new(2)+2
        set(axh,'cLim',blim_new);
    end
end
end


function ipos=imt_ln(echoln,wp,cstr)
% To detect intima-media thickness (IMT) of posterior wall in end-diastole.
% Starting from the adventitia an interative search is made for a minimum
% (media?) in the first derivative followed by a backward search for the
% lumen-intima transition (a distinct peak in the first derivative)
% All input, processing and output are in pixels
global vArt
wp=round(wp);
wind=ceil(2.5*vArt.dip); % assumed echo resolution 2.5 videoline (0.25 mm)
xwin=3;       % one-sided expansion of evaluation window to ensure overlap
dln=diff(echoln);                               % compute first derivative
dln=[0;0;dln(2:end)];             % elminate leading artifacts and realign
[adv,~]=max(echoln(wp-wind:end));                  % maximum of adventitia
lumval=mean(echoln(1:30));
endvalue=0.3*(adv-lumval);
bpos=find(echoln-lumval>endvalue,1);
if max(echoln(1:wp-wind))<0.7*adv
    oldpos=length(echoln)-find(fliplr(echoln')<adv/2,1);
else
    oldpos=wp;
end
stint=max(1,oldpos-3*wind);
[~,minpos]=min(dln(stint:oldpos));
minpos=minpos+stint-1;                                    % position media
stint=minpos-wind;                       % reposition assumed start intima
[oldval,oldpos]=max(dln(stint:minpos));     % first guess derivative level
oldpos=oldpos+stint-1;
newpos=oldpos;
ipos=[];
mecho_old=mean(echoln(stint:minpos));
acount=0;
while stint>4*vArt.dip    % start from minimum, search inwards for maximum
    stint=stint-wind;
    [ival,newpos]=max(dln(stint:stint+2*wind)); 
    mecho_new=mean(echoln(stint:stint+2*wind));
    newpos=newpos+stint-1;
    if oldpos-newpos<xwin% local maximum, xwin avoids ambiguity on plateau
        if ~isempty(ipos)
            %imax=max(dln(max(1,newpos-3*wind):newpos));
            imax=max(dln(max(1,min(newpos-3*wind,bpos)):newpos));
            cthr=vArt.ithr*oldval;
            if ival<cthr && imax<cthr % new maximum substantially smaller
                break;       % maximum found at 'newpos', possibly intima
            end
        end
        oldval=max(ival,oldval);
        ipos=newpos;           % remember this location (possibly intima)
    end
    oldpos=newpos;
    if (mecho_new<0.4*mecho_old)&& acount>2
        break                    % safety valve, no distinct maximum found
    end
    mecho_old=mecho_new;
    acount=acount+1;
end 
if isempty(ipos)                    % safety valve: accept any large peak
    bp=max(1,newpos-wind);
    [~,ipos]=max(dln(bp:minpos));
    ipos=ipos+bp-1;
end 
if vArt.imt_fb~=0                 % graphical feedback on intima detection
    fh=figure;
    set(fh,'WindowButtonDownFcn',@DoPause);
    plot(echoln) 
    hold on
    plot(5*max(0,dln),'g');   % plot first derivative (scaled by factor 5)
    plot(ipos,echoln(ipos),'o','MarkerSize',6,'color','r',...
        'MarkerFaceColor','g');
    plot(wp,echoln(wp),'o','MarkerSize',6,'color','r',...
        'MarkerFaceColor','w');
    axis tight
    title(cstr);
    xstring=sprintf(['Averaged A-line (black) and derivative (green);',...
        '\nLeft button down to pause.']);
    xlabel(xstring);
    pause(1);
    delete(fh);
end
function DoPause(fh,~)
    set(fh,'WindowButtonDownFcn',@DoResume);
    xstring=sprintf(['Averaged A-line (black) and derivative (green);',...
        '\nLeft button down to continue.']);
    xlabel(xstring);
    uiwait
end
function DoResume(fh,~)
    set(fh,'WindowButtonDownFcn',@DoPause);
    uiresume
end
end
       
function [wpos,IMTpos]=mod_MarkPos(wpos,IMTpos,cblk,cstr)
% Allows interactive adjustment of wall and intima mark position. 
% Inputs: detected intima (IMTpos) and anterior/posterior adventitia 
% (wpos) positions at subsegments. Valid positions select manual 
% adjustment of either wall (IMTpos empty) or intima (both 'wpos' and 
% 'IMTpos' available) positions. It is assumed that display of B-mode 
% image for current beat is within current figure.
% Arnold Hoeks, MU, BME, April 2013
global vArt modstop          % 'modstop' is also set in function 'getIMT'
hold on
dscal=1000*vArt.pixscal/vArt.dip;%pixel scaling after depth interpolation
if ~isempty(strfind(cstr,'New'))
    cthr={'wall threshold';'0.5';'0.55';'0.6';'0.65';'0.7';'0.75';'0.8'};
    cthr{1}=sprintf('wall threshold=%4.2f',vArt.wthr);
    eh=uicontrol('style','popupmenu','units','normalized','string',...
    cthr,'Position',[0.03 0.01 0.2 0.05],'FontSize',vArt.FOS,'callback',@doThr);
    dscal=1000*vArt.pixscal;                        % uninterpolated image
end
finbut=uicontrol('style','pushbutton','string','Accept/Continue','units',...
    'normalized','Position',[0.75 0.01 0.2 0.05],'FontSize',vArt.FOS,...
    'callback',@DoStop);
modstop=0;
hold on
if isempty(IMTpos)                   % edit wall-lumen interface positions
    title(sprintf('Point with left button down at desired wall positions\n%s',...
        cstr));
    mode=0;
    swp=zeros(2,vArt.rWidth);
    hWP=zeros(2,vArt.nblk);
    for ln=1:vArt.nblk                      % plot wall adventitia markers
        cln=cblk(ln);
        hWP(1,ln)=plot(cln,wpos(1,ln),'o','MarkerSize',5,...
            'color','r','MarkerFaceColor','w');
        hWP(2,ln)=plot(cln,wpos(2,ln),'o','MarkerSize',5,...
            'color','r','MarkerFaceColor','w');
    end
elseif isempty(wpos)
    return
else
    title(sprintf('Point with left button down at desired intima positions\n%s',...
        cstr));
    mode=1;
    hIMT=zeros(1,vArt.nblk);
    for ln=1:vArt.nblk                                % plot intima makers
        cln=cblk(ln);
        hWP(1,ln)=plot(cln,wpos(1,ln),'o','MarkerSize',5,...
            'color','r','MarkerFaceColor','w');
        hWP(2,ln)=plot(cln,wpos(2,ln),'o','MarkerSize',5,...
            'color','r','MarkerFaceColor','w');
        hIMT(ln)=plot(cln,IMTpos(1,ln),'o','MarkerSize',5,...
            'color','r','MarkerFaceColor','g');
    end
end
iwspan=round(vArt.wspan/vArt.pixscal);
while modstop==0                                           % Start cycling
    set(gcf,'WindowButtonMotionFcn',@ChangePointer)
    if exist('slh')                                 % remove smoothed line
        delete(slh);
    end
    diam=wpos(2,:)-wpos(1,:);              % current diameter distribution
    for wall=1:2
        swp(wall,:)=interp1(cblk,wpos(wall,:),1:vArt.rWidth,...
            'linear','extrap');
        swp(wall,:)=smooth(swp(wall,:),iwspan/vArt.rWidth,'rloess');
    end
    slh=plot((1:vArt.rWidth),swp,'b');    
    if mode==0                         % present (smoothed) wall positions
        xlabel(sprintf('Diameter mean=%d um;  sd=%d um; (N=%d)',...
            round([mean(diam),std(diam)]*dscal),length(diam)));
    else                             % present (smoothed) intima positions
        invalid=find(wpos(2,:)-IMTpos<2.5);       % unacceptable (<200 um)
        if ~isempty(invalid)
            IMTpos(invalid)=wpos(2,invalid);
        end
        imt=wpos(2,:)-IMTpos;
        lumi=interp1(cblk,IMTpos,1:vArt.rWidth,'linear','extrap');
        lumi=smooth(lumi,iwspan/vArt.rWidth,'rloess');
        slh=plot(lumi,'g');
        xlabel(sprintf('IMT mean=%d um;  sd=%d um; (N=%d)',...
            round([mean(imt),std(imt)]*dscal),length(imt)));
    end
    uiwait
end
title(' ');
set(gcf,'pointer','arrow','WindowButtonMotionFcn',[],...
    'WindowButtonDownFcn',[],'WindowButtonUpFcn',[]);
delete(finbut);
function doThr(~,~)
    sel=get(eh,'value');
    if sel>1
        vArt.wthr=sscanf(char(cthr{sel}),'%f');
        cthr{1}=sprintf('wall threshold=%4.2f',vArt.wthr);
        set(eh,'value',1,'string',cthr);
        vArt.wthr2=(1+vArt.wthr)/2;        % set threshold for second echo
    end
    delete(eh);delete(hWP);delete(slh);
    modstop=2; 
    uiresume 
end
function ChangePointer(~,~)
    [ppOK,~,~]=CheckPP();
    if ppOK                                                 % within axes
        set(gcf,'Pointer','crosshair','WindowButtonDownFcn',@StartDis);
    else
        set(gcf,'Pointer','arrow','WindowButtonDownFcn',[]);
    end
end
function StartDis(~,~)
    [ppOK,~,~]=CheckPP();
    if ppOK                                                 % within axes
        set(gcf,'WindowButtonMotionFcn',@new_wmf,...
            'WindowButtonUpFcn',@StopDis);
    else
        set(gcf,'WindowButtonMotionFcn',@ChangePointer);
    end
end
function new_wmf(~,~)
[ppOK,xp,yp]=CheckPP();
if ppOK                                                 % within axes
    [~,xln]=min(abs(cblk-xp));   % look for closest x-position segment
    if mode==1                                   % modify IMT position
        IMTpos(xln)=yp; 
        set(hIMT(xln),'ydata',yp);       % move handle to new position
    else             % adjust either anterior or posterior wall marker
        [~,yln]=min(abs(wpos(:,xln)-yp));            % 'yln' is 1 or 2
        wpos(yln,xln)=yp;
        set(hWP(yln,xln),'ydata',yp);
    end
end
end
function StopDis(~,~)
    uiresume
end
function DoStop(~,~)
    modstop=1;
    uiresume
end
end

function [ppOK,xp,yp]=CheckPP()
start=get(gca,'Currentpoint');         % quoted in scaled axis units!!
xp=start(1,1);
yp=start(1,2);
xl=get(gca,'Xlim');
yl=get(gca,'Ylim');
if xp>xl(1)&&xp<xl(2)&&yp>yl(1)&&yp<yl(2)                % within axes
    ppOK=1;
else ppOK=0;
end
end
    
function getPixScal()
global vArt
ah=gca;
if vArt.pixscal>0
    str=sprintf('Current distance scaling is %6.2f pixel/mm',vArt.pixscal);
    str=sprintf('%s (%d pixel/mm).\n\nIs this factor correct?',...
        str,round(1/vArt.pixscal));
    reply=questdlg(str,'Calibration','No, recalibrate','Yes, accept',...
        'Yes, accept');
    if strcmp(reply,'Yes, accept')
        return;
    end
end
title('Distance calibration procedure');
th=uicontrol('style','text','units','normalized','position',...
    [0.2 0.015, 0.6 0.035],'FontSize',vArt.FOS,'string',...
    'Click pointer at 2 depthmarks exactly 2 cm apart.');
set(gcf,'WindowButtonDownFcn',@getDMark);
hold on;
pcnt=0;
inip=zeros(2,2);
uiwait
vArt.pixscal=20/abs(inip(2,1)-inip(2,2));                     % mm/pixel
set(gcf,'WindowButtonDownFcn',[]);
delete(th);
function getDMark(hobj,~)
    pcnt=pcnt+1;
    cp = get(ah,'CurrentPoint');
    xp=cp(1,1);
    yp=cp(1,2);
    inip(1,pcnt)=xp;
    inip(2,pcnt)=yp;
    dm(pcnt)=plot(xp,yp,'o','MarkerSize',5,'color','r',...
        'MarkerFaceColor','w');
    if pcnt==2
        uiresume
    end
end
end

function limmess(str,dur)
% function displays a message for 'dur' seconds 
% Unlike the Matlab message dialog, this function has a limited duration.
% Arnold Hoeks,MU, october 2010
if iscell(str)
    mstr=[str;' ';'Kill this message box if you are impatient!'];
else
    mstr={str;' ';'Kill this message box if you are impatient!'};
end
wh=waitbar(0,mstr);
nstep=100;
for wloop=0:dur/nstep:dur
    if (~ishandle(wh))
        break
    else        
        waitbar(wloop/dur)
    end
    pause(dur/nstep);
end
if ishandle(wh)
    delete(wh);
end
end

function [infotxt,abbr]=getInfo()
infotxt={'The vidDist-program processes B-mode ultrasound video-sequences, recorded either in avi- or in DICOM-mode, to extract the (change in) artery lumen diameter over a few cardiac cycles (typically 5-10).';...
    'It is assumed that transducer position is at the top of the image, hence superficial arteries will be positioned in the horizontal direction.';...
    'The vertically oriented A-lines are realigned such that within the converted subimage the wall-lumen transition approaches a straight (horizontal) line.';...
    'Wall detection is referenced to the short-time (25 frames) lowpass filtered wall-lumen position, smoothed across the image segments with a 2nd order Savitsky-Golay filter (span equals the number of segments).';...
    'The realignment allows calculation of the spatial average of a few adjacent A-mode lines within half overlapping image segments (~3 mm).';...
    'Detection of the lumen-wall transition at subsequently the anterior and posterior walls is based at a crossing of a relative threshold within a few mm of the wall reference line.';...
    'The incidental spatial diameter distribution is the instantaneous difference between the detected posterior and anterior wall positions across the image.';...
    'The subimage (within the main B-mode image) considered for wall detection is readjusted in the vertical direction to track the observed artery position';...
    'The diameter waveforms expressed in pixels are converted to mm (only for DICOM recordings)';... 
    'The observed diameter distribution as function of time is saved in a MatLab file together with details about video characteristics and processing conditions.';...
    'In the post-processing phase, the 2D diameter distribution is passed through a 2nd order Savitsky-Golay smoothing filter (span 0.25 s), followed by a second order Savitsky-Golay highpass filter (span 2.5 s).';...
    'The unfiltered and filtered distributions are displayed along side each other, allowing interactive selection of the lowpass span width.';...
    'If accepted, at both sides 2 waveforms are discarded since they may carry extremes caused by the spatial smoothing filter.';...
    'The smoothed 2D diameter waveforms (and the spatial average and standard deviation) are gathered in sheets of an Excel logfile.';...
    'In the last processing phase, the statistical characteristics of the diameter waveform distributions and their spatial average are extracted.';...
    'First the end-diastolic time points in the spatial average distension waveform are identified.';...
    'The signal point 3 sample points prior to the minimum is considered and marked as R-top.';...
    'Using these points as time reference, for each complete cardiac cycle the diastolic, systolic, peak-peak, rise-time and mean distension are extracted.';...
    'The rise time is the time elapsed between the 10% and 90% level crossing of the systolic upstroke.';...
    'Next, those characteristics are statistically (mean, standard deviation) evaluated over the cardiac cycles.';...
    'Note that the std of the 2D distension distribution deviates from std of its spatial average!';' ';...
    'Arnold Hoeks. Dept Biomedical Engineering Maastricht University. Version April2012.'}; 
abbr={'Abbreviations','video distension program';...
    'file','file name';...
    'Source','*.avi or *.dcm (DICOM) or no extension';...
    'fps(Hz)','frames per second';...
    'nfr','total number of frames within considered time fragment';...
    'rec.time','recording time in seconds';...
    'Depth','depth of subframe (pixels) with selected artery segment';...
    'Width','width of subframe (pixels) with selected artery segment';...
    'PixSize','pixel size (horizontal/vertical) in mm';...
    'blksize','blocksize of image segment wall identification(pixels)';...
    'Nrblks','Number of blocks in selected image segment';...
    'blkshift','block shift (0.5 means half-overlapping)';...
    'dres','assumed resolution in depth [mm]';....
    'lres','assumed lateral resolution [mm]';....
    'interpol','interpolation factor for wall detection';....
    'Wall thr','threshold wall detection relative to adventitia peak';....
    'intima thr','treshold relative to intima peak';....
    'sm.window(s)','length smoothing window Savitsky-Golay filter';...
    'ECG-cycles','number of complete cardiac cycles';...
    'ECG-bpm','Heart rate (beats per minute)';...
    'Acq.date','Acquisition date (only for Dicom images)';....
    'proc.date','data processing date';....
    'RR(ms)','Duration cardiac cycle [ms]';....
    'Diast.mean','mean (beats) of MEDIAN of diastolic diameters (ROI)';...
    'Diast.med','median (beats) of MEDIAN of diastolic diameters (ROI)';...
    'Diast.hom','mean homogeneity (beats) diastolic diameters (ROI)';...
    'Diast.A','mean (beats) of MEAN of diastolic diameters (ROI)';...
    'Syst.mean','mean (beats) of MEDIAN of systolic diameters (ROI)';...
    'Syst.med','median (beats) of MEDIAN of systolic diameters (ROI)';...
    'Syst.hom','mean homogeneity (beats) systolic diameters (ROI)';...
    'Syst.A','mean (beats) of MEAN of systolic diameters (ROI)';...
    'Dist.mean','mean (beats) of MEDIAN of distension (ROI)';...
    'Dist.med','median (beats) of MEDIAN of distension (ROI)';...
    'Dist.hom','mean homogeneity (beats) distension (ROI)';...
    'Dist.A','mean (beats) of MEAN of distension (ROI)';...
    'rdist.mean','mean (beats) of MEDIAN of relative distension (ROI)';...
    'rdist.med','median (beats) of MEDIAN of relative distension (ROI)';...
    'rdist.hom','mean homogeneity (beats) relative distension (ROI)';...
    'rdist.A(%)','mean (beats) of MEAN of relative distension (ROI)';...
    'RT.mean','mean (beats) of MEDIAN of rise time [ms] (ROI)';...
    'RT.med','median (beats) of MEDIAN of rise time (ROI)';...
    'RT.hom','mean homogeneity (beats) rise time (ROI)';...
    'RT.A','mean (beats) of MEAN of rise time (ROI)';...
    'IMT.mean','mean (beats) of MEDIAN of IMT [um] (ROI)';...
    'IMT.med','median (beats) of MEDIAN of IMT (ROI)';...
    'IMT.hom','mean homogeneity (beats) IMT (ROI)';...
    'IMT.A','mean (beats) of MEAN of IMT (ROI)';...
    'rIMT.mean','mean (beats) of MEDIAN of IMT/diam [%] (ROI)';...
    'rIMT.med','median (beats) of MEDIAN of IMT/diam [%] (ROI)';...
    'rIMT.hom','mean homogeneity (beats) IMT/diam [%] (ROI)';...
    'rIMT.A','mean (beats) of MEAN of IMT/diam [%] (ROI)';...
    'IMTmax','maximum intima-media wall thickness'};
end

function Answer=inputdlg_check(varargin)
  %Check whether xlate exists; otherwise inputdlg_mod does not function (as
  %is the case in MATLAB R2015b
  if exist('xlate')
    Answer=inputdlg_mod(varargin{:});
  else
    Answer=inputdlg(varargin{:});
  end
end

function Answer=inputdlg_mod(Prompt, Title, NumLines, DefAns, Resize)
%INPUTDLG Input dialog box.
%  ANSWER = INPUTDLG(PROMPT) creates a modal dialog box that returns user
%  input for multiple prompts in the cell array ANSWER. PROMPT is a cell
%  array containing the PROMPT strings.
%
%  INPUTDLG uses UIWAIT to suspend execution until the user responds.
%
%  ANSWER = INPUTDLG(PROMPT,NAME) specifies the title for the dialog.
%
%  ANSWER = INPUTDLG(PROMPT,NAME,NUMLINES) specifies the number of lines for
%  each answer in NUMLINES. NUMLINES may be a constant value or a column
%  vector having one element per PROMPT that specifies how many lines per
%  input field. NUMLINES may also be a matrix where the first column
%  specifies how many rows for the input field and the second column
%  specifies how many columns wide the input field should be.
%
%  ANSWER = INPUTDLG(PROMPT,NAME,NUMLINES,DEFAULTANSWER) specifies the
%  default answer to display for each PROMPT. DEFAULTANSWER must contain
%  the same number of elements as PROMPT and must be a cell array of
%  strings.
%
%  ANSWER = INPUTDLG(PROMPT,NAME,NUMLINES,DEFAULTANSWER,OPTIONS) specifies
%  additional options. If OPTIONS is the string 'on', the dialog is made
%  resizable. If OPTIONS is a structure, the fields Resize, WindowStyle, and
%  Interpreter are recognized. Resize can be either 'on' or
%  'off'. WindowStyle can be either 'normal' or 'modal'. Interpreter can be
%  either 'none' or 'tex'. If Interpreter is 'tex', the prompt strings are
%  rendered using LaTeX.
%
%  Examples:
%
%  prompt={'Enter the matrix size for x^2:','Enter the colormap name:'};
%  name='Input for Peaks function';
%  numlines=1;
%  defaultanswer={'20','hsv'};
%
%  answer=inputdlg(prompt,name,numlines,defaultanswer);
%
%  options.Resize='on';
%  options.WindowStyle='normal';
%  options.Interpreter='tex';
%
%  answer=inputdlg(prompt,name,numlines,defaultanswer,options);
%
%  See also DIALOG, ERRORDLG, HELPDLG, LISTDLG, MSGBOX,
%    QUESTDLG, TEXTWRAP, UIWAIT, WARNDLG .

%  Copyright 1994-2007 The MathWorks, Inc.
%  $Revision: 1.58.4.17 $

%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
error(nargchk(0,5,nargin));
error(nargoutchk(0,1,nargout));

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Handle Input Args %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
  Prompt='Input:';
end
if ~iscell(Prompt)
  Prompt={Prompt};
end
NumQuest=numel(Prompt);

if nargin<2,
  Title=' ';
end

if nargin<3
  NumLines=1;
end

if nargin<4
  DefAns=cell(NumQuest,1);
  for lp=1:NumQuest
    DefAns{lp}='';
  end
end

if nargin<5
  Resize = 'off';
end
WindowStyle='modal';
Interpreter='none';

Options = struct([]); %#ok
if nargin==5 && isstruct(Resize)
  Options = Resize;
  Resize  = 'off';
  if isfield(Options,'Resize'),      Resize=Options.Resize;           end
  if isfield(Options,'WindowStyle'), WindowStyle=Options.WindowStyle; end
  if isfield(Options,'Interpreter'), Interpreter=Options.Interpreter; end
end

[rw,cl]=size(NumLines);
OneVect = ones(NumQuest,1);
if (rw == 1 & cl == 2) %#ok Handle []
  NumLines=NumLines(OneVect,:);
elseif (rw == 1 & cl == 1) %#ok
  NumLines=NumLines(OneVect);
elseif (rw == 1 & cl == NumQuest) %#ok
  NumLines = NumLines';
elseif (rw ~= NumQuest | cl > 2) %#ok
  error('MATLAB:inputdlg:IncorrectSize', 'NumLines size is incorrect.')
end

if ~iscell(DefAns),
  error('MATLAB:inputdlg:InvalidDefaultAnswer', 'Default Answer must be a cell array of strings.');
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Create InputFig %%%
%%%%%%%%%%%%%%%%%%%%%%%
FigWidth=175;
FigHeight=100;
FigPos(3:4)=[FigWidth FigHeight];  %#ok
FigColor=get(0,'DefaultUicontrolBackgroundcolor');

InputFig=dialog(                     ...
  'Visible'          ,'off'      , ...
  'KeyPressFcn'      ,@doFigureKeyPress, ...
  'Name'             ,Title      , ...
  'Pointer'          ,'arrow'    , ...
  'Units'            ,'pixels'   , ...
  'UserData'         ,'Cancel'   , ...
  'Tag'              ,Title      , ...
  'HandleVisibility' ,'callback' , ...
  'Color'            ,FigColor   , ...
  'NextPlot'         ,'add'      , ...
  'WindowStyle'      ,WindowStyle, ...
  'Resize'           ,Resize       ...
  );
%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
DefOffset    = 5;
DefBtnWidth  = 53;
DefBtnHeight = 23;

TextInfo.Units              = 'pixels'   ;
TextInfo.FontSize           = get(0,'DefaultUICOntrolFontSize');
TextInfo.FontWeight         = get(InputFig,'DefaultTextFontWeight');
TextInfo.HorizontalAlignment= 'left'     ;
TextInfo.HandleVisibility   = 'callback' ;

StInfo=TextInfo;
StInfo.Style              = 'text'  ;
StInfo.BackgroundColor    = FigColor;

EdInfo=StInfo;
EdInfo.FontWeight      = get(InputFig,'DefaultUicontrolFontWeight');
EdInfo.Style           = 'edit';
EdInfo.BackgroundColor = 'white';

BtnInfo=StInfo;
BtnInfo.FontWeight          = get(InputFig,'DefaultUicontrolFontWeight');
BtnInfo.Style               = 'pushbutton';
BtnInfo.HorizontalAlignment = 'center';

% Add VerticalAlignment here as it is not applicable to the above.
TextInfo.VerticalAlignment  = 'bottom';
TextInfo.Color              = get(0,'FactoryUIControlForegroundColor');
% adjust button height and width
btnMargin=1.4;
ExtControl=uicontrol(InputFig   ,BtnInfo     , ...
  'String'   ,xlate('Cancel', '-s')        , ...
  'Visible'  ,'off'         ...
  );

% BtnYOffset  = DefOffset;
BtnExtent = get(ExtControl,'Extent');
BtnWidth  = max(DefBtnWidth,BtnExtent(3)+8);
BtnHeight = max(DefBtnHeight,BtnExtent(4)*btnMargin);
delete(ExtControl);

% Determine # of lines for all Prompts
TxtWidth=FigWidth-2*DefOffset;
ExtControl=uicontrol(InputFig   ,StInfo     , ...
  'String'   ,''         , ...
  'Position' ,[ DefOffset DefOffset 0.96*TxtWidth BtnHeight ] , ...
  'Visible'  ,'off'        ...
  );

WrapQuest=cell(NumQuest,1);
QuestPos=zeros(NumQuest,4);

for ExtLp=1:NumQuest
  if size(NumLines,2)==2
    [WrapQuest{ExtLp},QuestPos(ExtLp,1:4)]= ...
      textwrap(ExtControl,Prompt(ExtLp),NumLines(ExtLp,2));
  else
    [WrapQuest{ExtLp},QuestPos(ExtLp,1:4)]= ...
      textwrap(ExtControl,Prompt(ExtLp),80);
  end
end % for ExtLp

delete(ExtControl);
QuestWidth =QuestPos(:,3);
QuestHeight=QuestPos(:,4);

TxtHeight=QuestHeight(1)/size(WrapQuest{1,1},1);
EditHeight=TxtHeight*NumLines(:,1);
EditHeight(NumLines(:,1)==1)=EditHeight(NumLines(:,1)==1)+4;

FigHeight=(NumQuest+2)*DefOffset    + ...
  BtnHeight+sum(EditHeight) + ...
  sum(QuestHeight);

TxtXOffset=DefOffset;

QuestYOffset=zeros(NumQuest,1);
EditYOffset=zeros(NumQuest,1);
QuestYOffset(1)=FigHeight-DefOffset-QuestHeight(1);
EditYOffset(1)=QuestYOffset(1)-EditHeight(1);

for YOffLp=2:NumQuest,
  QuestYOffset(YOffLp)=EditYOffset(YOffLp-1)-QuestHeight(YOffLp)-DefOffset;
  EditYOffset(YOffLp)=QuestYOffset(YOffLp)-EditHeight(YOffLp);
end % for YOffLp

QuestHandle=[]; %#ok
EditHandle=[];
AxesHandle=axes('Parent',InputFig,'Position',[0 0 1 1],'Visible','off');
inputWidthSpecified = false;

for lp=1:NumQuest,
  if ~ischar(DefAns{lp}),
    delete(InputFig);
    error('MATLAB:inputdlg:InvalidInput', 'Default Answer must be a cell array of strings.');
  end

  EditHandle(lp)=uicontrol(InputFig    , ...
    EdInfo      , ...
    'Max'        ,NumLines(lp,1)       , ...
    'Position'   ,[ TxtXOffset EditYOffset(lp) TxtWidth EditHeight(lp) ], ...
    'String'     ,DefAns{lp}           , ...
    'Tag'        ,'Edit'                 ...
    );

  QuestHandle(lp)=text('Parent'     ,AxesHandle, ...
    TextInfo     , ...
    'Position'   ,[ TxtXOffset QuestYOffset(lp)], ...
    'String'     ,WrapQuest{lp}                 , ...
    'Interpreter',Interpreter                   , ...
    'Tag'        ,'Quest'                         ...
    );

  MinWidth = max(QuestWidth(:));
  if (size(NumLines,2) == 2)
    % input field width has been specified.
    inputWidthSpecified = true;
    EditWidth = setcolumnwidth(EditHandle(lp), NumLines(lp,1), NumLines(lp,2));
    MinWidth = max(MinWidth, EditWidth);
  end
  FigWidth=max(FigWidth, MinWidth+2*DefOffset);

end % for lp

% fig width may have changed, update the edit fields if they dont have user specified widths.
if ~inputWidthSpecified
  TxtWidth=FigWidth-2*DefOffset;
  for lp=1:NumQuest
    set(EditHandle(lp), 'Position', [TxtXOffset EditYOffset(lp) TxtWidth EditHeight(lp)]);
  end
end

FigPos=get(InputFig,'Position');
units_old=get(InputFig,'units');
FigWidth=max(FigWidth,2*(BtnWidth+DefOffset)+DefOffset);
FigPos(3)=FigWidth;
FigPos(4)=FigHeight;
set(InputFig,'Position',FigPos,'units','normalized');
FigPos=get(InputFig,'Position');
FigPos(1)=(1-FigPos(3))/2;
FigPos(2)=(1-FigPos(4))/2;
set(InputFig,'Position',FigPos,'units',units_old);

uicontrol(InputFig,...
  BtnInfo      , ...
  'Position'   ,[FigWidth-BtnWidth-2*DefOffset DefOffset BtnWidth BtnHeight ] , ...
  'KeyPressFcn',@doControlKeyPress , ...
  'String'     ,'OK'        ,...
  'Callback'   ,@doCallback ,...
  'BackGroundColor','green' ,...
  'Tag'        ,'OK'        ,...
  'UserData'   ,'OK'         ...
  );

handles = guihandles(InputFig);
handles.MinFigWidth = FigWidth;
handles.FigHeight   = FigHeight;
handles.TextMargin  = 2*DefOffset;
guidata(InputFig,handles);
set(InputFig,'ResizeFcn', {@doResize, inputWidthSpecified});

% make sure we are on screen
movegui(InputFig)

% if there is a figure out there and it's modal, we need to be modal too
if ~isempty(gcbf) && strcmp(get(gcbf,'WindowStyle'),'modal')
  set(InputFig,'WindowStyle','modal');
end

set(InputFig,'Visible','on');
drawnow;

if ~isempty(EditHandle)
  uicontrol(EditHandle(1));
end

if ishghandle(InputFig)
  % Go into uiwait if the figure handle is still valid.
  % This is mostly the case during regular use.
  uiwait(InputFig);
end

% Check handle validity again since we may be out of uiwait because the
% figure was deleted.
if ishghandle(InputFig)
  Answer={};
  if strcmp(get(InputFig,'UserData'),'OK'),
    Answer=cell(NumQuest,1);
    for lp=1:NumQuest,
      Answer(lp)=get(EditHandle(lp),{'String'});
    end
  end
  delete(InputFig);
else
  Answer={};
end
end

function doFigureKeyPress(obj, evd) %#ok
switch(evd.Key)
  case {'return','space'}
    set(gcbf,'UserData','OK');
    uiresume(gcbf);
  case {'escape'}
    delete(gcbf);
end
end

function doControlKeyPress(obj, evd) %#ok
switch(evd.Key)
  case {'return'}
    if ~strcmp(get(obj,'UserData'),'Cancel')
      set(gcbf,'UserData','OK');
      uiresume(gcbf);
    else
      delete(gcbf)
    end
  case 'escape'
    delete(gcbf)
end
end

function doCallback(obj, evd) %#ok
if ~strcmp(get(obj,'UserData'),'Cancel')
  set(gcbf,'UserData','OK');
  uiresume(gcbf);
else
  delete(gcbf)
end
end

function doResize(FigHandle, evd, multicolumn) %#ok
% TBD: Check difference in behavior w/ R13. May need to implement
% additional resize behavior/clean up.

Data=guidata(FigHandle);

resetPos = false;

FigPos = get(FigHandle,'Position');
FigWidth = FigPos(3);
FigHeight = FigPos(4);

if FigWidth < Data.MinFigWidth
  FigWidth  = Data.MinFigWidth;
  FigPos(3) = Data.MinFigWidth;
  resetPos = true;
end

% make sure edit fields use all available space if
% number of columns is not specified in dialog creation.
if ~multicolumn
  for lp = 1:length(Data.Edit)
    EditPos = get(Data.Edit(lp),'Position');
    EditPos(3) = FigWidth - Data.TextMargin;
    set(Data.Edit(lp),'Position',EditPos);
  end
end

if FigHeight ~= Data.FigHeight
  FigPos(4) = Data.FigHeight;
  resetPos = true;
end

if resetPos
  set(FigHandle,'Position',FigPos);
end
end

% set pixel width given the number of columns
function EditWidth = setcolumnwidth(object, rows, cols)
% Save current Units and String.
old_units = get(object, 'Units');
old_string = get(object, 'String');
old_position = get(object, 'Position');

set(object, 'Units', 'pixels')
set(object, 'String', char(ones(1,cols)*'x'));

new_extent = get(object,'Extent');
if (rows > 1)
  % For multiple rows, allow space for the scrollbar
  new_extent = new_extent + 19; % Width of the scrollbar
end
new_position = old_position;
new_position(3) = new_extent(3) + 1;
set(object, 'Position', new_position);

% reset string and units
set(object, 'String', old_string, 'Units', old_units);
EditWidth = new_extent(3);
end

function [frameData3M, pixelSizeV,ferr]=readVEVOBModeRAW(fNameRawBMode,nFrames)
% READVEVORAWBMODE Function to read VEVO *.raw.bmode files.
% nFrames specifies the number of frames to read or 
%   if nFrames==0 only specifics are read (no data), or
%   If nFrames==-1, all frames will be read.
% frameData3M is an nSamples*nLines*nFrames uint8 matrix of image frames.
% pixelSizeV is a 2*1 double vector that denotes the pixel size in
%   micrometer of the images in frameData3M.
%   pixelSizeV(1) is the pixel size in depth direction,
%   pixelSizeV(2) is the pixel size in horizontal direction.
%   Note that the images in frameData3M are of the dimensions
%   samples*lines, i.e., displaying them without scaling will cause the
%   images to be stretched.
% V1.1 by Bart Spronck.
%   Modified by Arnold Hoeks, august 2015
  
global vArt
ferr=0;
frameData3M=[];
pixelSizeV=[];
[a,b] = fileparts(fNameRawBMode);
fNameRawXML = fullfile(a,[b '.xml']);
if ~exist(fNameRawXML, 'file')
  ferr=-1;
  return
%     error('.raw.xml file was not found at %s', fNameRawXML)
end

sysParsS = readVEVORawXML(fNameRawXML);
namesC = {sysParsS.name};         % Create cell array with parameter names
valuesC = {sysParsS.value};      % Create cell array with parameter values
nLines =   str2double(valuesC{strcmp(namesC,'B-Mode/Lines'        )});
nSamples = str2double(valuesC{strcmp(namesC,'B-Mode/Samples'      )});
maxDepth = str2double(valuesC{strcmp(namesC,'B-Mode/Depth'        )});
minDepth = str2double(valuesC{strcmp(namesC,'B-Mode/Depth-Offset' )});
width    = str2double(valuesC{strcmp(namesC,'B-Mode/Width'        )});
vArt.acqdate = valuesC{strcmp(namesC,'Acquired-Date'        )};
datstr=datestr(vArt.acqdate, 'yyyymmdd');       % convert to yyyy/mm/dd
vArt.acqdate=sprintf('%s/%s/%s',datstr(1:4),datstr(5:6),datstr(7:8)); 

% Calculate pixel size
pixelSizeV = zeros(2,1);
pixelSizeV(1) = (maxDepth - minDepth) * 1000 / nSamples;
pixelSizeV(2) = width * 1000 / nLines;
vArt.pixscal=(maxDepth - minDepth)/nSamples;% pixel size (um per sample)

%% Read binary data
fid = fopen(fNameRawBMode);

%Read file header
dwVersion=fread(fid,1,'uint32');     % Version number (3) of raw data file
dwNumFrames=fread(fid,1,'uint32'); % Number of frames in this raw datafile
dwInfo=fread(fid,1,'uint32');%       Information bitfield used to identify
% the type of frame data. Only dwInfo==8 is supported, corresponding to
% 8-bit RAW data.

if (dwVersion~=3)||(dwInfo ~= 8)      % wrong version or wrong data format
    ferr=-1;
    fclose(fid);
    return
end

vArt.nfr=dwNumFrames;                         % number of available frames
if nFrames==0                                        % read only specifics
    fclose(fid);
    return
elseif nFrames>0
    nFrames=min(nFrames,dwNumFrames);    % read a limited number of frames
else
    nFrames=dwNumFrames;                    % else all frames will be read
end                                        
vArt.nfr=dwNumFrames; 
temp = fread(fid,7,'uint32');                % Unused bytes for future use

%Initialise arrays to save frame data
frameData3M = zeros(nSamples, nLines, nFrames, 'uint8');
frameTimeStampsTicksV = zeros(nFrames,1, 'uint32');
frameTimeMilliSecondsV = zeros(nFrames,1, 'double');

for fr = 1:nFrames                                     % Read frame header
    dwTimeStamp=fread(fid,1,'uint32');% Hardware time stamp counted as ticks of a 400 kHz clock
    dbTimeStamp=fread(fid,1,'double');% Hardware time stamp counted in ms. Calculated as dwTimeStamp * 1000 / 400000
    dwFrameNumber=fread(fid,1,'uint32'); % Frame number, starting from 0
    dwInfo=fread(fid,1,'uint32');    % Is 0 if frame contains valid data
    dwPacketSize=fread(fid,1,'uint32');  %Size in bytes of frame data (not including header)
    temp=fread(fid,8,'uint32');

    if dwInfo ~= 0
        ferr=-1;
        error('Frame does not contain valid data')
    end

    if dwPacketSize ~= nLines*nSamples
        ferr=-1;
        error('Corrupt frame found')
    end

    % Read frame image data into an uint8 1D array
    frameData=uint8(fread(fid,dwPacketSize,'uint8'));

    %Process data
    frameDataM=reshape(frameData, [nSamples nLines]);   % Reshape to frame
    %dimensions.

    frameData3M(:,:,fr)=frameDataM;
    frameTimeStampsTicksV(fr)       =dwTimeStamp;
    frameTimeMilliSecondsV(fr)=dbTimeStamp;
end

fclose(fid);
rectime=frameTimeMilliSecondsV(nFrames)-frameTimeMilliSecondsV(1);
vArt.fps=round((10000*nFrames/rectime))/10;            % frames per second
end

function sysParsS=readVEVORawXML(fNameRawXML)
  % READVEVORAWXML Function to read VEVO *.raw.xml files that contain imaging settings
  %   sysParsS=READVEVORAWXML(fNameRawXML) Reads a .raw.xml file and
  %     parses its contents into the structure sysParsS.
  %     sysParsS is an nx1 structure with n the number of parameters
  %     specified in the XML file (typically 53).
  %     For each parameter, sysParsS contains three fields which are all
  %     character arrays:
  %       sysParsS(i).name:  parameter name (e.g. 'Time-Stamp-Clock')
  %       sysParsS(i).value: parameter value (e.g. '400000')
  %       sysParsS(i).units: parameter units (e.g. 'Hz')
  %     The units field may be empty for dimensionless parameters.
  %
  % V1.0 by Bart Spronck.
  
  xmlS = parseXML(fNameRawXML);
  
  sysParsS = struct(); %Empty struct to save parameters
  
  okIV = false(length(xmlS.Children),1);
  for i = 1:length(xmlS.Children)
    if ~isempty(xmlS.Children(i).Attributes)
      okIV(i) = true;
      for j = 1:length(xmlS.Children(i).Attributes)
        sysParsS(i,1).(xmlS.Children(i).Attributes(j).Name) = xmlS.Children(i).Attributes(j).Value;
      end
    end
  end
  sysParsS = sysParsS(okIV); %Remove empty fields.
end

function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
% Copied from MATLAB R2015a documentation.
try
   tree = xmlread(filename);
catch
   error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems 
% with very deeply nested trees.
try
   theStruct = parseChildNodes(tree);
catch
   error('Unable to parse XML file %s.',filename);
end
end

% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end
end

% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.
attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end
end