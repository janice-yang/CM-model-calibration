% iPSC Genetic Algorithm
% Optimization of Model Current Parameters to Experimental Data
% Multiple runs in parallel - only changing rng seed
% clear
% close all

%load GA/runNum

% Protocol numbers:
%{ 
1 - Cao High, PCL = 800ms	
2 - Cao High, spontaneous
3 - Cao Low, PCL = 800ms	
4 - Cao Low, spontaneous	
5 - ICaL Block = 25%, PCL = 800ms	
6 - ICaL Block = 25%, spontaneous	
7 - ICaL Block = 50%, PCL = 800ms	
8 - ICaL Block = 50%, spontaneous	
9 - IKr Block = 15%, PCL = 800ms	
10- IKr Block = 15%, spontaneous	
11- IKr Block = 30%, PCL = 800ms	
12- IKr Block = 30%, spontaneous	
13- Ko High, PCL = 800ms	
14- Ko High, spontaneous	
15- Ko Low, PCL = 800ms	
16- Ko Low, spontaneous	
17- PCL = 600ms	
18- PCL = 800ms	
19- Spontaneous	
%}

addpath('Classes_all')
addpath('dydts_all')
addpath('Other_helper_functions')
addpath('GA')
% addpath('Setup_and_parameters')
% Initial conditions determined using f_getInit_Y.m

% Parameters to determine with GA - names here for convenience, but values
% inside gaK19.m class, we just feed in the scale factors here.

% % Set in command line
% popsize = 20 ;
% n_gen = 5 ;
% n_stallgen = 5 ;
% seeds = [1, 2] ; % rng seeds (array of non-negative integers)
% cell_number = 1 ;
% protocol_number = [18,5,7] ;
% isNormalized = true; % whether or not data will be normalized

% Open parallel worker pool
% n_cores = 48 ; % set in command line
n_workers = length(seeds) ;
parpool(n_cores) ;

%% Parfor loop through rng seeds
% for w=seeds
parfor (w=1:length(seeds), n_workers)
    runNum = seeds(w) ;
    par_sga(runNum, cell_number, protocol_number, isNormalized, scaleCaT, popsize, n_gen, n_stallgen)
end
% Delete parallel pool
delete(gcp('nocreate'));

%% Single cell, multiple protocol sequences - not ready yet
% 
% if (length(cell_numbers) == 1 && length(protocol_seqs) > 1)
%     parfor w=1:length(protocol_seqs)
%         cell_number = cell_numbers(1) ;
%         protocol_number = protocol_seqs{w} ;
%         par_sga(w, cell_number, protocol_number, isNormalized, popsize, n_gen, n_stallgen)
%     end
%     % Move all events of the cell to first run folder
%     try
%     movefile('GA/Results/events*.mat', 'GA/Results/Run_1/Events');
%     catch E
%     end
% 
% %%-- Multiple cells, single protocol sequence
% elseif (length(cell_numbers) > 1 && length(protocol_seqs) == 1)
%     parfor w=1:length(cell_numbers)
%         cell_number = cell_numbers(w) ;
%         protocol_number = protocol_seqs{1} ;
%         par_sga(w, cell_number, protocol_number, isNormalized, popsize, n_gen, n_stallgen)
% 
%         try
%         movefile('GA/Results/events*.mat', ['GA/Results/', 'Run_',int2str(runNum),'/Events']);
%         catch E
%         end
%     end
% 
% %%-- Throw error for any other case
% else
%     error(['Parallel run currently only supports either: ', ...
%         '(1) Running multiple protocol sequences for 1 cell, or ', ...
%         '(2) Running one protocol sequence on multiple cells. ', ...
%         'For a single cell + protocol sequence run, use the batch scripts.'])
% end

