% Plot protocol error & SD 
clear
close all

% Parameters
realdata = false ;
folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select experiment DIRECTORY/IES') ;
n = 10 ; % runs per experiment

base = 'Analysis/ProtocolRankings' ; % folder for storing results
% Get current date & time for filename
format shortg
date_time = char(strjoin(string(clock), '_')) ;

metrics = {'Absolute error', 'SD'} ;
plot_data = zeros(length(folders), length(metrics)) ; % store data
plot_error = zeros(length(folders), length(metrics)) ; % store error bar data
experiment_names = cell(length(folders), 1) ; 
load GA/x_names.mat
data = zeros(length(metrics), length(folders), n) ;

%% Get metrics for each experiment
for i=1:length(folders)
    folder = folders{i} ;
    if ispc % Windows
        experiment = strsplit(folders{i}, '\') ;
    else
        experiment = strsplit(folders{i}, '/') ;    
    end
    experiment_name = experiment{end} ;
    experiment_names{i} = experiment_name ;
    
    % Get stats file
    load([folder, '/stats_data.mat'], 'log_norm_params', 'log_abs_errors')
    
    if sum(contains(metrics, 'Absolute error'))
        j = find(ismember(metrics, 'Absolute error')) ;
        data(j, i, :) = mean(log_abs_errors, 2) ; % avg across runs
    end
    
    if sum(contains(metrics, 'SD'))
        j = find(ismember(metrics, 'SD')) ;
        data(j, i, :) = std(log_norm_params, 0, 2) ;      
    end
    
    if sum(contains(metrics, 'Range'))
        j = find(ismember(metrics, 'Range')) ;
        data(j, i, :) = range(log_norm_params, 2) ;
    end
    
    if sum(contains(metrics, 'Validation error'))
        j = find(ismember(metrics, 'Validation error')) ;
        data(j, i, :) = validation_errors.sum(1:end-1) ;
    end
    
    if sum(contains(metrics, 'Fitting error'))
        j = find(ismember(metrics, 'Fitting error')) ;
        data(j, i, :) = fitting_errors.average(1:end-1) ;
    end
    
end

%% Plots
if sum(contains(metrics, 'Absolute error'))
    figure
    hold on
    j = find(ismember(metrics, 'Absolute error')) ;
    plot_data = reshape(data(j, :, :), length(folders), n) ;
    boxplot(plot_data')
    xticks(1:length(experiment_names)) 
    xticklabels(experiment_names)
    xtickangle(90)
    ylabel('Average absolute error')
    title('Absolute error')
    savefig([base, '/', date_time, '_MAE.fig'])
end

if sum(contains(metrics, 'SD'))
    figure
    hold on
    j = find(ismember(metrics, 'SD')) ;
    plot_data = reshape(data(j, :, :), length(folders), n) ;
    boxplot(plot_data')
    xticks(1:length(experiment_names)) 
    xticklabels(experiment_names)
    xtickangle(90)
    ylabel('Average SD')
    title('SD')
    savefig([base, '/', date_time, '_SD.fig'])
end

if sum(contains(metrics, 'Range'))
    figure
    hold on
    j = find(ismember(metrics, 'Range')) ;
    plot_data = reshape(data(j, :, :), length(folders), n) ;
    boxplot(plot_data')
    xticks(1:length(experiment_names)) 
    xticklabels(experiment_names)
    xtickangle(90)
    ylabel('Average range')  
    title('Range')
    savefig([base, '/', date_time, '_range.fig'])
end

if sum(contains(metrics, 'Validation error'))
    figure
    hold on
    j = find(ismember(metrics, 'Validation error')) ;
    plot_data = reshape(data(j, :, :), length(folders), n) ;
    boxplot(plot_data')
    xticks(1:length(experiment_names)) 
    xticklabels(experiment_names)
    xtickangle(90)
    ylabel('Average validation error')  
    title('Validation error')
    savefig([base, '/', date_time, '_validationError.fig'])
end

if sum(contains(metrics, 'Fitting error'))
    figure
    hold on
    j = find(ismember(metrics, 'Fitting error')) ;
    plot_data = reshape(data(j, :, :), length(folders), n) ;
    boxplot(plot_data')
    xticks(1:length(experiment_names)) 
    xticklabels(experiment_names)
    xtickangle(90)
    ylabel('Average fitting error')  
    title('Fitting error')
    savefig([base, '/', date_time, '_fittingError.fig'])
end


% Save data
save([base, '/', date_time, '_data.mat'], ...
    'metrics', 'data', 'experiment_names')

