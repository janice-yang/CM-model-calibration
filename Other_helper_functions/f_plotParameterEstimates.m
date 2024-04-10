function [] = f_plotParameterEstimates(folders, objective, logfactor, logplot, lim, data_model, fit_model)
%%%%%% CHANGE FUNCTION NAME
% Plot distribution of GA-optimized parameter estimates vs actual parameters
% 
% Input:
%   folders = (cell array) directory/ies containing run results
%   objective = (string) 'single' or 'multi' indicating single- or multi-objetive GA
%   logfactor = (int) scaling factor for parameters
%   logplot = (bool) whether to plot log normalized params
%   data_model = (string) which model generated the data
%       'experiment' - normalize parameters to baseline values
%       'Kernik19'
%       'Paci18'
%       'Tomek19'
%   fit_model = (string) which model the data were fitted to
%       'Kernik19'
%       'Paci18'
% Output:
%   Plots and saves parameter distribution
%%%%%%

%%% Catch weird inputs
if strcmp(objective, 'multi') && (length(folders) ~= 1)
    error("Currently only support 1 run for multi-objective GA.")
end
if ~(sum(strcmp(objective, {'single', 'multi'})))
    error("GA objective options: 'single', 'multi'")
end
if ~sum(strcmp(data_model, {'experiment', 'Kernik19', 'Paci18', 'Tomek19'}))
    error("Data model options: 'experiment', 'Kernik19', 'Paci18', 'Tomek19'")
end
if ~sum(strcmp(fit_model, {'Kernik19', 'Paci18'}))
    error("Fitted model options: 'Kernik19', 'Paci18'")
end

% Set output file path
if ispc % Windows
    foldersplit = strsplit(folders{1}, '\') ;
else
    foldersplit = strsplit(folders{1}, '/') ;    
end
base = strjoin(foldersplit(1:end-1), '/') ; % experiment folder with all runs
fp_out = [base, '/log_optimparams_distribution'] ;

% Load experiment info & initialize results array
load([folders{1}, '/Details.mat'], 'cell_number', 'protocol_number')
sequence = regexprep(num2str(protocol_number), '\s*' , '-') ;
results = {} ;

% Conductance names
load GA/x_names_Paci.mat names
Paci_names = names ;
load GA/x_names.mat names
Kernik_names = names ;
load GA/x_names_Tomek.mat names
Tomek_names = names ;

if strcmp(fit_model, data_model) || strcmp(data_model, 'experiment')
    if strcmp(fit_model, 'Kernik19')
        names = Kernik_names ;
        i_K = 1:length(Kernik_names) ;
    elseif strcmp(fit_model, 'Paci18')
        names = Paci_names ;
        i_P = 1:length(Paci_names) ;
    else
        error("Model options: 'Kernik19', 'Paci18'")
    end
else  % fitted model doesn't match model/experiment that produced data
    if ~(sum(strcmp(data_model, {'Kernik19', 'Paci18', 'Tomek19'})) && sum(strcmp(fit_model, {'Kernik19', 'Paci18'})))
        error("Currently only support Kernik19, Paci18, and Tomek19 for mixing data & fitted models.")
    elseif strcmp(data_model, 'Paci18')
        [names, i_K, i_P] = intersect(Kernik_names, Paci_names) ;
    elseif strcmp(data_model, 'Tomek19')
        [names, i_K, i_T] = intersect(Kernik_names, Tomek_names) ;        
    end
end

% Load ground truth parameter values
if strcmp(data_model, 'experiment')
    ground_truth = ones(1, length(names)) ;
elseif strcmp(data_model, 'Kernik19')
    ground_truth = readmatrix('Pseudodataset/saved_data/ground_truth_conductances.xlsx') ; 
    ground_truth = ground_truth(cell_number, i_K) ;
elseif strcmp(data_model, 'Paci18')
    ground_truth = readmatrix('Pseudodataset/saved_data/Paci18_ground_truth_conductances_remaining.xlsx') ; 
    ground_truth = ground_truth(cell_number, i_P) ; 
elseif strcmp(data_model, 'Tomek19')
    ground_truth = readmatrix('Pseudodataset/saved_data/Tomek_ground_truth_conductances.xlsx') ; 
    ground_truth = ground_truth(cell_number, i_T) ; 
end

% Load & normalize parameter estimates
if strcmp(objective, 'single')
    for i=1:length(folders)
        load([folders{i}, '/minoptimparams.mat'], 'minoptimparams')
        if strcmp(fit_model, 'Kernik19')
            minoptimparams = minoptimparams(i_K) ;
        elseif strcmp(fit_model, 'Paci18')
            minoptimparams = minoptimparams(i_P) ;
        end
        results{i} = logfactor.^minoptimparams ;
    end
    
    normalized_params = zeros(length(results), length(names)) ;
    for i=1:length(results)
        normalized_params(i, :) = results{i} ./ ground_truth ;    
    end
elseif strcmp(objective, 'multi')
    % pareto front individuals for multi-objective GA
        folder = folders{1} ;
    load([folder, '/minoptimparams.mat'], 'optimparams')
    if strcmp(fit_model, 'Kernik19')
        optimparams = optimparams(:,i_K) ;
    elseif strcmp(fit_model, 'Paci18')
        optimparams = optimparams(:,i_P) ;
    end
    % Calculate scaled params + normalize to ground truth
    optimparams = logfactor.^optimparams ;
    normalized_params = optimparams ./ ground_truth ;
else
    error("Enter GA objective as 'single' or 'multi'.")
end

% Plot
figure
hold on
if logplot
    if length(folders) == 1 % single run
        plot(log2(normalized_params), 'bo', "LineWidth", 2)
    elseif length(folders) > 1 % multiple rng seed runs
        UnivarScatter(log2(normalized_params), 'Label', names, 'Width', 0.5, 'Whiskers', 'lines', ...
            'WhiskersWidthRatio', 1.25);
    end
    ylabel("log(predicted/ground truth)")
    ylim([-lim lim])
    yticks(-lim:lim)
    yline(0, 'g:', 'LineWidth', 2, 'HandleVisibility', 'off')
else
    if length(folders) == 1 % single run
        plot(normalized_params, 'bo', "LineWidth", 2)
    elseif length(folders) > 1 % multiple rng seed runs
        UnivarScatter(normalized_params, 'Label', names, 'Width', 0.5, 'Whiskers', 'lines', ...
            'WhiskersWidthRatio', 1.25);
    end
    ylabel("predicted / ground truth")
    ylim([0 2^lim])
    yline(1, 'g:', 'LineWidth', 2, 'HandleVisibility', 'off')
end
title(['Optimized parameters distribution for protocol ' sequence])
xticks(1:length(names))
xticklabels(names)
xtickangle(45)

savefig([fp_out, '.fig'])

end
   