clear
close all

load Pseudodataset/saved_data/gaKernik_population_unfiltered

T = 310.0 ;         % (K)
Nao = 151.0 ;       % (mM)
Ko = 5.4 ;          % (mM)
Cao = 1.8 ;         % (mM)

%% Run spontaneous beating for all cells
for i=1:length(gaKernik_cells)
   setUpPacingProtocol(gaKernik_cells{i}, 0, 300, 999, 1) ;
   setEnvironment(gaKernik_cells{i}, T, Nao, Cao, Ko) ;
   odeSolver(gaKernik_cells{i}) ;
end
disp('Finished running simulations of all cells!')

%% Calculate features
sbr_all = zeros(1, length(gaKernik_cells)) ;
for i=1:length(gaKernik_cells)
    t = gaKernik_cells{i}.state.t ;
    V = gaKernik_cells{i}.state.Y(:, 1) ;
    [~, ~, features] = APextract_custom(t, V, []) ;
    gaKernik_cells{i}.protocol.spontfreq = features(4) ; % store beating frequency
    sbr_all(i) = gaKernik_cells{i}.protocol.spontfreq ;
end

% Plot frequency distribution
figure
histogram(sbr_all, 50)
title("Beating frequency distribution")
xlabel("Frequency (Hz)")
ylabel("Cell count")
% xline(0.3, 'r--')
% xline(1.0, 'r--')
print('-dpng', 'Pseudodataset/saved_data/gaKernik_pop_beating_frequency')

%% Filter by spontaneous beating rate - sbr_filter function
min_sbr = 0.3 ; % Hz
max_sbr = 1.0 ; % Hz
filtered_cells = sbr_filter(gaKernik_cells, max_sbr, min_sbr) ;
save Pseudodataset/saved_data/gaKernik_cells_filtered filtered_cells min_sbr max_sbr

%% Randomly select 10 cells
clear
load Pseudodataset/saved_data/gaKernik_cells_filtered
rng(0) % set seed for reproducibility
n = 10 ; % number of cells to select
selected_idx = randperm(length(filtered_cells), n) ;
selected_cells = cell(1, length(selected_idx)) ;
for i=1:n
   selected_cells{i} = filtered_cells{selected_idx(i)} ; 
end

% Save conductances of selected cells
columns = {'Cell number', ...
    'g_Na', 'g_f', 'p_CaL', 'g_to', 'g_Ks', 'g_Kr', 'g_K1', 'g_PCa', 'g_b_Na', 'g_b_Ca', 'VmaxUp', 'g_irel_max', 'kNaCa', 'PNaK', 'V_leak', 'g_CaT'} ;
data = zeros(length(selected_cells), length(columns)) ;
for i=1:length(selected_cells)
   data(i, 1) = i ;
   data(i, 2:end) = selected_cells{i}.conductances.scaling ;
end
% data(:, 1) = int16(data(:, 1)) ;

writematrix(data(:, 2:end), 'Pseudodataset/saved_data/ground_truth_conductances.xlsx')
tbl = array2table(data, 'VariableNames', columns) ;
writetable(tbl, 'Pseudodataset/saved_data/ground_truth_conductances_labeled.xlsx')

