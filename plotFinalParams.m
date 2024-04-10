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

% Plot Kernik baseline parameters / ground truth cell
% load([folders{1}, '/Details.mat'], 'cell_number') ;
% ground_truth = readmatrix('Pseudodataset/saved_data/ground_truth_conductances.xlsx') ;
% ground_truth = ground_truth(cell_number, :) ;  
% baseline = 1 ./ ground_truth ; 
% baseline_log = log2(baseline) ;
% plot(baseline_log, 'r^', "LineWidth", 2)
% print('-dpng', [fp_out, '_withBaseline']) ;
