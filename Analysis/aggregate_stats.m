clear
close all

% Calculate calibration error, spread, validation metrics over all cells
% For single protocol - run section
realdata = false ;
protocol_name = 'p18-5-7_APCaT_norm' ; % will use as part of output file name
base = 'Analysis/AggregateStats' ; % folder for storing results
metrics = {'Absolute error', 'MAE', 'SD', 'Validation error', 'Validation accuracy', 'Validation SD'} ;

load GA/x_names.mat

folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select experiment DIRECTORY/IES') ;
n = 10 ; % runs per experiment

% Get current date & time for filename
format shortg
date_time = char(strjoin(string(clock), '_')) ;

set(groot, 'defaultFigureRenderer', 'painters') % for saving as svg

plot_data = zeros(n * length(folders), length(metrics)) ; % store data
cell_numbers = zeros(n*length(folders), 1) ; % store cell order

% Get metrics for each experiment
idx_data = 1 ;
for i=1:length(folders)
    folder = folders{i} ;
    load([folder, '/Run_0/Details.mat'], 'cell_number')
    cell_numbers(idx_data:idx_data+n-1) = ones(n, 1)*cell_number ;
    
    % Get stats file
    load([folder, '/stats_data.mat'], 'log_norm_params', 'log_abs_errors')

    if sum(contains(metrics, 'Absolute error'))
        j = find(ismember(metrics, 'Absolute error')) ;
        plot_data(idx_data:idx_data+n-1, j) = mean(log_abs_errors, 2) ; 
    end

    if sum(contains(metrics, 'MAE'))
        j = find(ismember(metrics, 'MAE')) ;
        plot_data(idx_data:idx_data+n-1, j) = mean(mean(log_abs_errors, 2)) ; 
    end
    
    if sum(contains(metrics, 'SD'))
        j = find(ismember(metrics, 'SD')) ;
        plot_data(idx_data:idx_data+n-1, j) = std(log_norm_params, 0, 2) ;        
    end
    
    if sum(contains(metrics, 'Range'))
        j = find(ismember(metrics, 'Range')) ;
        plot_data(idx_data:idx_data+n-1, j) = range(log_norm_params, 2) ;
    end
    
    if sum(contains(metrics, 'Validation error')) % square error
        % Get true validation threshold or calculations - IKr threshold
        load([folder, '/thresholds_IKr.mat'], 'predictedThreshold') % 1x10 thresholds
        actual_IKr_thresholds = readtable('Pseudodataset/saved_data/IKr block thresholds.xlsx') ;
        actual_IKr = actual_IKr_thresholds.Threshold(cell_number) ;

        % Calculations
        IKr_errors = abs(predictedThreshold - actual_IKr) ;

        j = find(ismember(metrics, 'Validation error')) ;
        plot_data(idx_data:idx_data+n-1, j) = IKr_errors' ;

    end

    if sum(contains(metrics, 'Validation accuracy'))
        % Get true validation threshold or calculations - IKr threshold
        load([folder, '/thresholds_IKr.mat'], 'predictedThreshold') % 1x10 thresholds
        actual_IKr_thresholds = readtable('Pseudodataset/saved_data/IKr block thresholds.xlsx') ;
        actual_IKr = actual_IKr_thresholds.Threshold(cell_number) ;

        % Calculations
        IKr_errors = abs(((1-predictedThreshold) - (1-actual_IKr))) / (1-actual_IKr) ;

        j = find(ismember(metrics, 'Validation accuracy')) ;
        plot_data(idx_data:idx_data+n-1, j) = IKr_errors' ;
    end

    if sum(contains(metrics, 'Validation SD'))
        % Get true validation threshold or calculations - IKr threshold
        load([folder, '/thresholds_IKr.mat'], 'predictedThreshold') % 1x10 thresholds

        % Calculations
        IKr_std = std(predictedThreshold, 0, 2) ;        

        j = find(ismember(metrics, 'Validation SD')) ;
        plot_data(idx_data:idx_data+n-1, j) = IKr_std' ;
    end
    
    if sum(contains(metrics, 'Fitting error'))
        % Open fitting error file
        j = find(ismember(metrics, 'Fitting error')) ;
        plot_data(idx_data:idx_data+n-1, j) = fitting_errors.average(1:end-1) ;
    end

    idx_data = idx_data + n ;
    
end

% % Plot
% figure
% hold on
% boxplot(plot_data)
% xticks(1:length(metrics))
% xticklabels(metrics)
% xtickangle(45)
% savefig([base, '/', protocol_name, '_MAE.fig'])
% print('-dsvg', [base, '/', protocol_name, '_MAE.svg'])
% print('-dpng', [base, '/', protocol_name, '_MAE.png'])

% Save data
save([base, '/', protocol_name, '_data.mat'], ...
    'metrics', 'plot_data', 'protocol_name', 'cell_numbers')
tbl = array2table([plot_data, cell_numbers], 'VariableNames',[metrics, 'Cell Number']) ;
writetable(tbl, [base, '/', protocol_name, '.xlsx'])

%% Plot comparing aggregate stats for multiple protocols
clear
close all

% Use files with same metrics
files = uigetfile_n_dir([pwd, '/Analysis/AggregateStats'], 'Select stats files') ;
base = 'Analysis/AggregateStats/4cells_numconds' ; % output file name
n_cells = 4 ;
n_runs = 10 ;
cells_vector = reshape(repmat(1:n_cells, n_runs, 1), 1, []) ;
cells_matrix = repmat(cells_vector, 1, length(files))' ;
load(files{1}, 'metrics')

set(groot, 'defaultFigureRenderer', 'painters') % for saving as svg
for j=1:length(metrics)
    metric_name = metrics{j} ;
    all_data = zeros(length(files), n_cells*n_runs) ;
    all_protocols = cell(1, length(files)) ;
    all_cells = zeros(n_cells*n_runs, length(files)) ;

    figure
    hold on
    title(metric_name)
    colormap jet
    
    for i=1:length(files)
        load(files{i}, 'plot_data', 'protocol_name', 'cell_numbers')
        all_data(i, :) = plot_data(:, j) ;
        all_protocols{i} = protocol_name ;
        all_cells(:, i) = cell_numbers ;

        swarmchart(i*ones(1, length(plot_data(:,j))), plot_data(:,j), ...
            [], cell_numbers.*2, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
    end

    x = repmat(1:length(files), n_cells*n_runs, 1) ;
    boxplot(all_data')
    xticks(1:length(files))
    xticklabels(all_protocols)
    xtickangle(45)
    savefig([base, '_', metric_name, '.fig'])
    print('-dsvg', [base, '_', metric_name, '.svg'])
    print('-dpng', [base, '_', metric_name, '.png'])

    figure
    hold on
    title(metric_name)
    boxplot(all_data')
    xticks(1:length(files))
    xticklabels(all_protocols)
    xtickangle(45)
    savefig([base, '_', metric_name, '_boxonly.fig'])
    print('-dsvg', [base, '_', metric_name, '_boxonly.svg'])
    print('-dpng', [base, '_', metric_name, '_boxonly.png'])
    close
    
end

