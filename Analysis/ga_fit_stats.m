function [data] = ga_fit_stats(realdata, varargin)

%{
Function to calculate summary statistics for GA runs

Varargin: strings indicating which stats to calculate 
    'params': (logfactor.^minoptimparams).*baseline = predicted params
    'norm': params ./ ground truth
    'log2norm'
    'pcerror': 100*(predicted - ground truth) / ground truth
    'relerror': 100*|predicted - ground truth| / ground truth
    'abslog2norm'
    'sqerror': (predicted - ground truth)^2
    'normsqerror': (norm - 1)^2
Current: ga_fit_stats(false, 'params', 'norm', 'log2norm', 'abslog2norm') ;
    or ga_fit_stats(true, 'params', 'norm', 'log2norm', 'abslog2norm') ;
Output (optional):
    data: matrix of all output stats
%}

load GA/x_names.mat names

% baseline = [9.720613409241, 0.0435,  0.308027691379, 0.1178333333333, 0.0077,  0.218025, 0.133785777797606, 0.2625, 4.35e-04, 3.6704e-04, 1.105e-04, 12.5, 1100, 2.4761, 1.6e-06, 0.185] ;
% units = {'nS/pF', 'nS/pF', 'nS/pF', 'nS/pF', 'nS/pF', 'nS/pF', 'g_K1unit', 'A/F', 'nS/pF', 'nS/pF', 'mM/s', '/ms', 'A/F', 'A/F', '/ms', 'nS/pF'} ;

%% Inter-run and inter-cell optimized param variability from single GA runs
% Select directories with SGA run results
folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select DIRECTORIES containing GA results') ;
if ispc % Windows
    foldersplit = strsplit(folders{1}, '\') ;
else
    foldersplit = strsplit(folders{1}, '/') ;    
end
base = strjoin(foldersplit(1:end-1), '/') ; % experiment folder with all runs
fp_out = [base, '/stats_fitted_params.xlsx'] ;

% Run info: cell num, protocol sequence
load([folders{1}, '/Details.mat'], 'cell_number', 'protocol_number', 'isNormalized')

% sequence = regexprep(num2str(protocol_number), '\s*' , '-') ;
results = {} ;
logfactor = 2 ;

% Load ground truth parameters for each cell (or baseline if real data)
if realdata
    ground_truth = ones(1, length(names)) ;
else
    ground_truth = readmatrix('Pseudodataset/saved_data/ground_truth_conductances.xlsx') ;
    ground_truth = ground_truth(cell_number, :) ;
end
% ground_truth = ground_truth .* baseline ;

%% Calculate stats

params = zeros(length(results), length(names)) ;
norm_params = zeros(length(results), length(names)) ;
log_norm_params = zeros(length(results), length(names)) ;
percent_errors = zeros(length(results), length(names)) ;
rel_errors = zeros(length(results), length(names)) ;
log_abs_errors = zeros(length(results), length(names)) ;
sq_errors = zeros(length(results), length(names)) ;
norm_sq_errors = zeros(length(results), length(names)) ;

for i=1:length(folders)
    load([folders{i}, '/minoptimparams.mat'], 'minoptimparams')
    result = (logfactor.^minoptimparams) ;
    params(i, :) = result ;
    norm_params(i, :) = result ./ ground_truth ;
    log_norm_params(i, :) = log2(result ./ ground_truth) ;
    percent_errors(i, :) = 100*(result - ground_truth) ./ ground_truth ;
    rel_errors(i, :) = 100*(abs(result - ground_truth) ./ ground_truth) ;
    log_abs_errors(i, :) = abs(log_norm_params(i,:)) ;
    sq_errors(i, :) = (result - ground_truth).^2 ;
    norm_sq_errors(i, :) = (norm_params(i, :) - 1).^2 ;
end
save([base, '/stats_data.mat'], 'params', 'norm_params', 'log_norm_params', ...
    'percent_errors', 'rel_errors', 'log_abs_errors', 'sq_errors', 'norm_sq_errors')

%% Save summary stats, based on input
stats = {'mean','std','var','sem','min','max','meanci', 'range'} ; % stats to store
sz = [length(names), length(stats)+1] ; % mean CI is 2 values: lb and ub
data = cell(1, length(varargin)) ;
data_idx = 1 ;

% Parameter stats
if any(strcmp(varargin, 'params'))
    statarray = zeros(sz) ;
    for i=1:size(params, 2) % calculate stats for each parameter
        param = params(:, i) ;
        [mean,std,var,sem,min,max,meanci,range] = grpstats(param,[],stats) ;
        statarray(i, :) = [mean,std,var,sem,min,max,meanci,range] ;
    end
    colnames = {'mean','std','var','SEM','min','max','meanCI_lower', 'meanCI_upper', 'range'} ;
    tbl_params = array2table(statarray, 'VariableNames', colnames, 'RowNames', names) ;
    writetable(tbl_params, fp_out, 'WriteRowNames', true, 'Sheet','ParamStats')
    data{data_idx} = tbl_params ;
    data_idx = data_idx + 1 ;
end
    
% Ground truth-normalized parameters
if any(strcmp(varargin, 'norm'))
    statarray = zeros(sz) ;
    for i=1:size(norm_params, 2) % calculate stats for each parameter
        param = norm_params(:, i) ;
        [mean,std,var,sem,min,max,meanci,range] = grpstats(param,[],stats) ;
        statarray(i, :) = [mean,std,var,sem,min,max,meanci,range] ;
    end
    colnames = {'mean','std','var','SEM','min','max','meanCI_lower', 'meanCI_upper', 'range'} ;
    tbl_normparams = array2table(statarray, 'VariableNames', colnames, 'RowNames', names) ;
    writetable(tbl_normparams, fp_out, 'WriteRowNames', true, 'Sheet','NormParamStats')
    data{data_idx} = tbl_normparams ;
    data_idx = data_idx + 1 ;
end

% Log ground truth-normalized parameters
if any(strcmp(varargin, 'log2norm'))
    statarray = zeros(sz) ;
    for i=1:size(log_norm_params, 2) % calculate stats for each parameter
        param = log_norm_params(:, i) ;
        [mean,std,var,sem,min,max,meanci,range] = grpstats(param,[],stats) ;
        statarray(i, :) = [mean,std,var,sem,min,max,meanci,range] ;
    end
    colnames = {'mean','std','var','SEM','min','max','meanCI_lower', 'meanCI_upper', 'range'} ;
    tbl_lognormparams = array2table(statarray, 'VariableNames', colnames, 'RowNames', names) ;
    writetable(tbl_lognormparams, fp_out, 'WriteRowNames', true, 'Sheet','Log2NormParamStats')
    data{data_idx} = tbl_lognormparams ;
    data_idx = data_idx + 1 ;
end

% Percent change
if any(strcmp(varargin, 'pcerror'))
    statarray = zeros(sz) ;
    for i=1:size(percent_errors, 2) % calculate stats for each parameter
        pc_error = percent_errors(:, i) ;
        [mean,std,var,sem,min,max,meanci,range] = grpstats(pc_error,[],stats) ;
        statarray(i, :) = [mean,std,var,sem,min,max,meanci,range] ;
    end
    tbl_percent_error = array2table(statarray, 'VariableNames', colnames, 'RowNames', names) ;
    writetable(tbl_percent_error, fp_out, 'WriteRowNames', true, 'Sheet','PercentError')
    data{data_idx} = tbl_percent_error ;
    data_idx = data_idx + 1 ;    
end

% Relative error
if any(strcmp(varargin, 'relerror'))
    statarray = zeros(sz) ;
    for i=1:size(rel_errors, 2) % calculate stats for each parameter
        rel_error = rel_errors(:, i) ;
        [mean,std,var,sem,min,max,meanci,range] = grpstats(rel_error,[],stats) ;
        statarray(i, :) = [mean,std,var,sem,min,max,meanci,range] ;
    end
    tbl_rel_error = array2table(statarray, 'VariableNames', colnames, 'RowNames', names) ;
    writetable(tbl_rel_error, fp_out, 'WriteRowNames', true, 'Sheet','RelError')
    data{data_idx} = tbl_rel_error ;
    data_idx = data_idx + 1 ;  
end

% Log norm params error (absolute value)
if any(strcmp(varargin, 'abslog2norm'))
    statarray = zeros(sz) ;
    for i=1:size(log_abs_errors, 2) % calculate stats for each parameter
        rel_error = log_abs_errors(:, i) ;
        [mean,std,var,sem,min,max,meanci,range] = grpstats(rel_error,[],stats) ;
        statarray(i, :) = [mean,std,var,sem,min,max,meanci,range] ;
    end
    tbl_log_abs_error = array2table(statarray, 'VariableNames', colnames, 'RowNames', names) ;
    writetable(tbl_log_abs_error, fp_out, 'WriteRowNames', true, 'Sheet','LogAbsError')
    data{data_idx} = tbl_log_abs_error ;
    data_idx = data_idx + 1 ;  
end

% Squared errors
if any(strcmp(varargin, 'sqerror'))
    statarray = zeros(sz) ;
    for i=1:size(sq_errors, 2) % calculate stats for each parameter
        rel_error = sq_errors(:, i) ;
        [mean,std,var,sem,min,max,meanci,range] = grpstats(rel_error,[],stats) ;
        statarray(i, :) = [mean,std,var,sem,min,max,meanci,range] ;
    end
    tbl_sq_error = array2table(statarray, 'VariableNames', colnames, 'RowNames', names) ;
    writetable(tbl_sq_error, fp_out, 'WriteRowNames', true, 'Sheet','SquaredError')
    data{data_idx} = tbl_sq_error ;
    data_idx = data_idx + 1 ;  
end

% Predicted - ground truth squared errors
if any(strcmp(varargin, 'normsqerror'))
    statarray = zeros(sz) ;
    for i=1:size(norm_sq_errors, 2) % calculate stats for each parameter
        rel_error = norm_sq_errors(:, i) ;
        [mean,std,var,sem,min,max,meanci,range] = grpstats(rel_error,[],stats) ;
        statarray(i, :) = [mean,std,var,sem,min,max,meanci,range] ;
    end
    tbl_norm_sq_error = array2table(statarray, 'VariableNames', colnames, 'RowNames', names) ;
    writetable(tbl_norm_sq_error, fp_out, 'WriteRowNames', true, 'Sheet','NormSquaredError')
    data{data_idx} = tbl_norm_sq_error ;
    data_idx = data_idx + 1 ;  
end

% Ground truth
tbl_groundtruth = array2table(ground_truth', 'RowNames', names);
writetable(tbl_groundtruth, fp_out, 'WriteRowNames', true, 'Sheet','GroundTruth')
data{data_idx} = tbl_groundtruth ;

end
   
