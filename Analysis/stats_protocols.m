clear

%% Two-sample t-tests for comparing metrics between 2 protocols, 1 cell
addpath('../Other_helper_functions')

folders = uigetfile_n_dir([pwd, 'GA/Results'], 'Select DIRECTORY/IES containing GA results') ;
n_runs = 10 ; 
actual_IKr = 0.43438 ;

% Two sample t-test
if length(folders) ~= 2
    error("This test compares exactly 2 protocols.")
else
    param_MAE = zeros(n_runs, length(folders)) ;
    param_SD = zeros(n_runs, length(folders)) ;
    predicted_IKr = zeros(n_runs, length(folders)) ;
    for i=1:length(folders)
        load([folders{i}, '/thresholds_IKr.mat'], 'predictedThreshold')
        load([folders{i}, '/stats_data.mat'], 'log_abs_errors', 'log_norm_params')

        param_MAE(:, i) = mean(log_abs_errors,2) ;
        param_SD(:, i) = std(log_norm_params, 0, 2) ;
        predicted_IKr(:, i) = predictedThreshold' ;
    end
    IKr_errors = (predicted_IKr - actual_IKr).^2 ;
    
    % Test log-normalized parameter absolute errors
    disp("Comparing parameter errors for 2 protocols:")
    [h,p,ci,stats] = ttest2(param_MAE(:,1), param_MAE(:,2))
    
    % Test log-normalized parameter SD
    disp("Comparing parameter SDs for 2 protocols:")
    [h,p,ci,stats] = ttest2(param_SD(:,1), param_SD(:,2))    
    
    % Test predicted vs actual IKr threshold
    disp("Comparing IKr threshold errors for 2 protocols:")
    [h,p,ci,stats] = ttest2(IKr_errors(:,1), IKr_errors(:,2))      
    
end

%% Two-sample t-test: 2 protocols, multiple cells
clear

name = '1cond_APvCaT' ;
files = uigetfile_n_dir([pwd, '/Analysis/AggregateStats'], 'Select FILES containing compiled metrics') ;
n_cells = 4 ;
n_runs = 10 ; 
% actual_IKr = 0.43438 ;

% Two sample t-test
if length(files) ~= 2
    error("This test compares exactly 2 protocols.")
else
    param_errors = zeros(n_cells*n_runs, length(files)) ;
    param_MAE = zeros(n_cells*n_runs, length(files)) ;
    param_SD = zeros(n_cells*n_runs, length(files)) ;
    valid_errors = zeros(n_cells*n_runs, length(files)) ;
    valid_SD = zeros(n_cells*n_runs, length(files)) ;
    valid_acc = zeros(n_cells*n_runs, length(files)) ;

    for i=1:length(files)
        load(files{i})
        if sum(contains(metrics, 'Absolute error'))
            j = find(ismember(metrics, 'Absolute error')) ;
            param_errors(:, i) = plot_data(:, j) ;
        end
        if sum(contains(metrics, 'MAE'))
            j = find(ismember(metrics, 'MAE')) ;
            param_MAE(:, i) = plot_data(:, j) ;
        end
        if sum(contains(metrics, 'SD'))
            j = find(ismember(metrics, 'SD')) ;
            param_SD(:, i) = plot_data(:, j) ;
        end
        if sum(contains(metrics, 'Validation error'))
            j = find(ismember(metrics, 'Validation error')) ;
            valid_errors(:, i) = plot_data(:, j) ;
        end
        if sum(contains(metrics, 'Validation SD'))
            j = find(ismember(metrics, 'Validation SD')) ;
            valid_SD(:, i) = plot_data(:, j) ;
        end
        if sum(contains(metrics, 'Validation accuracy'))
            j = find(ismember(metrics, 'Validation accuracy')) ;
            valid_acc(:, i) = 1 - plot_data(:, j) ;
        end

    end

    % Test log-normalized parameter absolute errors
    disp("Comparing parameter errors for 2 protocols:")
    [h,p,ci,stats] = ttest2(param_errors(:,1), param_errors(:,2))
    
    % Test mean log-normalized parameter errors
    disp("Comparing mean parameter errors for 2 protocols:")
    [h,p,ci,stats] = ttest2(param_MAE(:,1), param_MAE(:,2))    
    
    % Test log-normalized parameter SD
    disp("Comparing parameter SDs for 2 protocols:")
    [h,p,ci,stats] = ttest2(param_SD(:,1), param_SD(:,2))      
    
    % Test validation errors
    disp("Comparing validation errors for 2 protocols:")
    [h,p,ci,stats] = ttest2(valid_errors(:,1), valid_errors(:,2))   
  
    % Test validation SDs
    disp("Comparing validation SDs for 2 protocols:")
    [h,p,ci,stats] = ttest2(valid_SD(:,1), valid_SD(:,2))   

    % Test validation accuracies
    disp("Comparing validation accuracies for 2 protocols:")
    [h,p,ci,stats] = ttest2(valid_SD(:,1), valid_SD(:,2)) 

end

%%

% [p,tbl,stats] = kruskalwallis(param_MAE) ;

% [p,tbl,stats] = kruskalwallis(param_SD) ;
% 
% [p,tbl,stats] = kruskalwallis(threshold_MAE) ;