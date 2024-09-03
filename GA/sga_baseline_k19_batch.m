% iPSC Genetic Algorithm
%% FOR RUNNING ON A COMPUTE CLUSTER USING LSF
% Optimization of Model Current Parameters to Experimental Data

% clear
% close all

% load GA/initials
% load GA/runNum
% runNum = runNum + 1 ; % comment out this line to overwrite prior run
% seed = 3 ;
% runNum = seed ;
% rng(seed) ;

% save GA/runNum runNum
mkdir('GA/Results', ['Run_',int2str(runNum)]);
% 
% realdata = false ;
% cell_number = 4 ;
% protocol_number = [19] ;
% isNormalized = true; % whether or not data will be normalized
% datatype = 'APCaT' ;
% scaleCaT = false ; % whether or not to scale CaT to approximate AP range
% sigmaAP = 0.0 ; % SD of noise to be added to pseudodata (AP)
% sigmaCaT = 0.0 ; % SD of noise to be added to pseudodata (CaT)
% 
% realdata = true ;
% filename = '20220113_paced.xlsx' ; % Spreadsheet with data in 2-3 columns (Time, AP, CaT)
% sheetnames = {'1.8Ca 1Hz', '1.8Ca 1.25Hz'} ;
% protocol_number = [28,27] ;
% isNormalized = true ;
% scaleCaT = false ;
% datatype = 'APCaT' ;
% nbeats = [11,11] ;

save GA/curr_cell_protocol protocol_number isNormalized scaleCaT datatype
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

% Parameters to determine with GA - names here for convenience, but values
% inside gaK19.m class, we just feed in the scale factors here.

names = {...
'g_Na',...        1
'g_f',...         2
'p_CaL',...       3   
'g_to',...        4  
'g_Ks',...        5
'g_Kr',...        6
'g_K1',...        7
'g_PCa',...       8
'g_b_Na',...      9
'g_b_Ca',...      10
'VmaxUp',...      11
'g_irel_max',...  12
'kNaCa',...       13
'PNaK',...        14
'V_leak',...      15
'g_CaT',...       16
% 'taud_scale',...  17
% 'tauf_scale',...  18
} ;

save GA/x_names names
       
nvars = length(names) ;

% Fitness Function (includes evaluation of model)
% fitnessfcn = @sga_fitness_k19;
fitnessfcn = @(x)sga_fitness_k19(x, runNum, 1); 
outputfcn = @(options,state,flag)ga_output_k19(options, state, flag, runNum, 1) ;

%%Load experimental data:
if realdata
    experimental_dataset = f_getExperimentData(filename, sheetnames, protocol_number, isNormalized, datatype, nbeats) ;
else
    experimental_dataset = f_getPseudodata(cell_number, protocol_number, isNormalized, sigmaAP, sigmaCaT);
end
%%
% Population Size and Variation

popsize = 150 ; 
%variation = 1 ;

% lower and upper bounds
  lb = ones(nvars,1)*-2;   
  ub = ones(nvars,1)*2;   
  % To set boundaries for specific parameters:
  %{l
  param_index = find(strcmp(names, 'paramName'));
  ub(param_index)=0.25; %Scale_factor=1.78
  lb(param_index)=-0.1; %Scale_factor=0.8
  %}

%Initial population inside boundaries
%initpop = (ub-lb)'.*lhsdesign(popsize,nvars) + lb'; %hypercube sampling
initpop = 4*rand(popsize, nvars)-2; %uniform distribution from -2 to 2 --> multiplied to parameters as 2^x, or 0.5 to 2

options = optimoptions(@ga, ... 
    'SelectionFcn', {@selectiontournament, 4}, ... % was 8
    'UseParallel', true,...
    'CrossoverFcn', {@crossoverscattered},... 
    'CrossoverFraction', 0.85, ...
    'MutationFcn', {@mutationadaptfeasible}, ...
    'EliteCount', ceil(0.05*popsize), ...
    'InitialPopulationMatrix',initpop, ...
    'PopulationSize', popsize, ...
    'PlotFcn', {@gaplotscorediversity, @gaplotbestf, @gaplotbestindiv}, ...
    'OutputFcn', outputfcn, ... 
    'Display', 'iter',...
    'MaxGenerations', 20,... % 100
    'MaxStallGenerations', 5);   % if there is no change in the best fitness value for some number of generations, GA stops

if realdata
    save(['GA/Results/Run_',int2str(runNum), '/Details'], 'filename', 'sheetnames', 'protocol_number', 'isNormalized', 'datatype', 'scaleCaT', 'options', 'popsize', 'runNum');
else
    save(['GA/Results/Run_',int2str(runNum), '/Details'], 'cell_number','protocol_number', 'isNormalized', 'datatype', 'scaleCaT', 'options', 'popsize', 'runNum');
end

%% Run Genetic Algorithm
  [optimparams,bestfitness,exitflag,output,population,scores] = ...
    ga(fitnessfcn,nvars,[],[],[],[],lb,ub,[],options);

%% Evaluate results:
% run one more time with the best parameter set from that generation 
% (optimparams), and save the results
orig_conductance = ones(1,length(optimparams(1,:)));
[optimparams, index] = unique(optimparams,'rows');
bestfitness = bestfitness(index,:);

[value, ind] = min(bestfitness);
minoptimparams = optimparams(ind,:);
save(['GA/Results/', 'Run_',int2str(runNum), '/minoptimparams'], 'minoptimparams','optimparams');

%% Simulate protocol/protocol series for optimal parameters
for i = 1:size(minoptimparams,1)
params = minoptimparams(i,:);
x_conductance = (2.^params);
[t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance,names); % see how the optimal parameters do

% Extract for each protocol in series
figure 
p = 1 ;
for j=1:length(protocol_number)
    [keepT, V, CaT,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},stimtimes{j}, 1);
    
    if strcmp(datatype,'AP') || strcmp(datatype,'APCaT')
        % Align ends of exp and simulated extracted waveforms
        exp_T_V = experimental_dataset{j}.Time_AP{1} ; % - tpeak + 100;
        exp_V = experimental_dataset{j}.AP{1};
        [exp_T_V, exp_V, T_V, V] = f_alignWaveformEnds(exp_T_V, exp_V, keepT, V) ;
        V_sim = interp1(T_V,V,exp_T_V);
        %Calculate errors: --> update if using normalized data
        rel_error_AP = 100*(abs(exp_V - V_sim)./exp_V);
        compare_AP = [exp_V, V_sim, rel_error_AP];
        
        subplot(length(protocol_number),2,p)
        plot(exp_T_V,V_sim,'LineWidth',2)
        hold on
        plot(exp_T_V,exp_V,'Color','red','LineWidth',2)
        title(['Protocol ',int2str(protocol_number(j)),' Action Potential'])
        legend('Prediction','Observation')
        set(gca,'FontSize',12)
    end
    
    %% Adjust for AP or CaT only
    if strcmp(datatype,'CaT') || strcmp(datatype,'APCaT')
        exp_T_CaT = experimental_dataset{j}.Time_CaT{1} ; % - tpeak + 100;
        exp_CaT = experimental_dataset{j}.CaT{1};
        [exp_T_CaT, exp_CaT, T_CaT, CaT] = f_alignWaveformEnds(exp_T_CaT, exp_CaT, keepT, CaT) ;
        CaT_sim = interp1(T_CaT, CaT, exp_T_CaT);
        
        subplot(length(protocol_number),2,p+1)
        plot(exp_T_CaT, CaT_sim,'LineWidth',2)
        hold on
        plot(exp_T_CaT, exp_CaT, 'Color', 'red','LineWidth',2)
        title(['Protocol ',int2str(protocol_number(j)),' CaT'])
        set(gca,'FontSize',12)
        hold off
        
        rel_error_CaT = 100*(abs(exp_CaT - CaT_sim)./exp_CaT);
        compare_CaT = [exp_CaT, CaT_sim, rel_error_CaT];
    end
    
    % Increment subplot position
    p = p + 2 ; 
    
end

fig = gcf;
fig.Units = 'normalized';
fig.Position = [0.25 0.3 0.5 0.4];
print(['GA/Results/', 'Run_',int2str(runNum), '/Outputs_',int2str(i)],'-dpng','-r0')

end

movefile(['GA/Results/', 'Run_',int2str(runNum), 'k19*.mat'], ['GA/Results/', 'Run_',int2str(runNum)]);

%Heatmap optimparams:
figure
h = heatmap(optimparams);
h.Colormap = parula;
h.XData = names;
saveas(gcf, ['GA/Results/Run_',int2str(runNum),'/optimparams_heatmap.fig'])

try
movefile('GA/Results/events*.mat', ['GA/Results/', 'Run_',int2str(runNum),'/Events']);
catch E
end

% Save whole workspace?
% save(['GA/Results/Run_',int2str(runNum),'/final_workspace.mat'])
