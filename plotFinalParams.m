clear
close all 
% Plot predicted parameter distributions
folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select DIRECTORY/IES containing GA results') ; 
objective = 'single' ;      % GA objective type options: single, multi
logfactor = 2 ;
logplot = true ;
lim = 3 ;
data_model = 'experiment' ;     % options: experiment, Kernik19, Paci18, Tomek19
fit_model = 'Kernik19' ;    % options: Kernik19, Paci18

% function: f_plotDistribution(folders, groundTruthNorm, objective, logfactor, fp_out)
f_plotParameterEstimates(folders, objective, logfactor, logplot, lim, data_model, fit_model) ;
