clear all;
clc;
close all;

%% QUASI-STATIC PRESSURE SWEEPS
% Define the percentages of the in vivo axial stretch at which the quasi-
% static pressure sweeps were performed. E.g. write '95' for a test
% performed at 95% of the in vivo axial stretch.
protocolMat.passive.PD_stretch_guess = {'95','100','105'};

%% QUASI-STATIC AXIAL FORCE SWEEPS
% Define the pressures at which the quasi-static axial force sweeps
% were performed. E.g. write '10' for a axial froce sweep performed at the
% constant introluminar pressure of 19 mmHg.
protocolMat.passive.FL_pressure_guess ={'10','60','100','140','180'};

%% DYNAMIC EXPERIMENTS
% Define the systolic and diastolci pressures, the loading frequency (in
% mHz) and the waveform type ('physio' = physiological, 'sine' = sine wave)
% of each dynamic test. E.g. write '120' in protocolMat.PD_SBP_guess, '80'
% in protocolMat.PD_DBP_guess, '10000' in protocolMat.PD_Hz_guess and
% 'sine' in protocolMat.PD_label_guess for a dynamic test with sinusoidal
% pressure waveform between 120 and 80 mmHg at a frequency of 10 Hz.
protocolMat.passive.PD_SBP_guess = {'80','80','80','80','80','80','120','120','120','120','120','120','160','160','160','160','160','160'}; %Number of Dynamic PD Tests
protocolMat.passive.PD_DBP_guess = {'40','40','40','40','40','40','80','80','80','80','80','80','120','120','120','120','120','120'};
protocolMat.passive.PD_Hz_guess = {'625','1250','2500','5000','10000','20000','625','1250','2500','5000','10000','20000','625','1250','2500','5000','10000','20000'};
protocolMat.passive.PD_label_guess = {'sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine','sine'};

% CONTRACTION EXPERIMENTS
protocolMat.active.vasoactive_agents = {'SNP', 'Lname'};

protocolMat.active.PD_stretch_guess = {'100'};

protocolMat.active.PD_SBP_guess = {'90','110','130','150','170','190'}; %Number of Dynamic PD Tests
protocolMat.active.PD_DBP_guess = {'50','70','90','110','130','150'}; %Number of Dynamic PD Tests
protocolMat.active.PD_Hz_guess = {'10000','10000','10000','10000','10000','10000'}; %Number of Dynamic PD Tests
protocolMat.active.PD_label_guess = {'sine','sine','sine','sine','sine','sine'};

%% CALIBRATION COEFFICIENTS
% Define the calibration coefficients for the load cell and the camera.
protocolMat.p2um = 1.638; % pixel to um conversion for the camera
protocolMat.V2mN = 33.545; % Volt to mN conversion for load cell

%% AXIAL STRETCH CORRECTION COEFFICIENT
protocolMat.correctingAxialStretch = 1;
% This quantifies the axial elongation of the sample during incubation
% experiments. Leave = 1, unless you are working with a study where you
% characterise mechanics pre- and post-incubation.

%% SAVING THE PROTOCOL STRUCTURE
% Define the name of the protocol file
file_name = 'Margarita_incubation';

save(file_name,'protocolMat')
