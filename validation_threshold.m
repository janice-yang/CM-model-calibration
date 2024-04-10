clear
close all

realdata = false ;
folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select DIRECTORY/IES containing GA results') ;

logfactor = 2 ;
%%% CHANGE names FILE TO MATCH RUNS
load GA/x_names.mat names
% conds = readmatrix('Pseudodataset/saved_data/ground_truth_conductances.xlsx') ;

init_Y = [-75.7444536163477,0.338969592726776,0.000203113729306209,7.16928093750999,104.748824394112,0,0.000386686486786781,0.165948056753057,0.927145173320106,0.321984775889061,0.452222061313948,0.157787210225653,0.743251868606151,0.121059208476135,0.0292207648412020,0.00620538308203812,0.736108314102295,0.000264118925707198,0.000263623380304203,0.746780367359425,0.0122283684412607,0.000154416932298621,0.0123158737520428,0] ;

% Environment settings
T = 310.0 ;         % (K)
Nao = 151.0 ;       % (mM)
Ko = 5.4 ;          % (mM)
Cao = 1.8 ;         % (mM)
% Stimulation setup 
simulationTime = 300000 ; % ms
pcl = 800 ; % pacing CL (ms) - 0 if spontaneous
amplitudes = 60 ;
numPulses = simulationTime / pcl ;
precedingTime = pcl - 1 ;
pulseDurations = 1 ;

if ispc % Windows
    foldersplit = strsplit(folders{1}, '\') ;
else
    foldersplit = strsplit(folders{1}, '/') ;    
end
base = strjoin(foldersplit(1:end-1), '/') ; % experiment folder with all runs

%%% Validate estimated parameters using stepped I_Kr block
tol = 0.01 ;
tol_freq = 0.01 ; % frequency error tolerance for arrhythmia 
max_iter = 20 ;

predictedThreshold = zeros(1, length(folders)) ;

for i=1:length(folders)
    load([folders{i}, '/Details.mat'], 'cell_number')
    load([folders{i}, '/minoptimparams.mat'], 'minoptimparams')
    x_conductance = logfactor.^minoptimparams ;
    
    top = 1.0 ; % IKr scaling
    bottom = 0.0 ; % IKr scaling
    
    gap = top - bottom ;
    counter = 0 ;
    
    while gap > tol
        gap = top - bottom ;
        mid = bottom + (top - bottom) ./ 2 ;
        % Simulate mid
        k19 = gaKernik19(init_Y, 1) ; 
        setUpPacingProtocol(k19, amplitudes, numPulses, precedingTime, pulseDurations) ;
        setEnvironment(k19, T, Nao, Cao, Ko) ;

        ikrblock = ones(1,16);
        ikrblock(6) = mid ;
        setUpDrugApplication(k19, ikrblock, zeros(1,16), ones(1,16)*simulationTime)

        saveX_conductance(k19, x_conductance) ;
        scaleConductances(k19, x_conductance, names); 
        scaleParameters(k19, x_conductance, names) ;

        odeSolver(k19);

        t = k19.state.t ;
        V = k19.state.Y(:,1) ;
        Cai = k19.state.Y(:,3) ;

        % Get last 5 seconds + calculate features
        idx_start = find(t > t(end) - 5000, 1) ;
        t = t(idx_start:end) ;
        V = V(idx_start:end) ;
        Cai = Cai(idx_start:end) ;
        dV = diff(V) ;
        Vrange = max(V) - min(V) ;
        Vhalf = 0.5*Vrange + min(V) ;

        below_dices = [find(V < min(V) + 0.5*Vrange);length(t)+1] ;
        times_up = t(below_dices(diff(below_dices) > 5)) ;
        cyclelength = mean(diff(times_up)) ;
        cycle_std = std(diff(times_up)) ;
        APfrequency = 1000/cyclelength ;    

        if cycle_std > 100 || Vrange < 10 || sum(dV) == 0 || APfrequency > 1000/pcl + tol_freq
            bottom = mid ;
        else % no arrhythmia detected
            top = mid ;
        end   

        % Stop loop if exceeded max iterations
        counter = counter + 1 ;
        if counter > max_iter
            break
        end
    end
    
    predictedThreshold(i) = bottom ;
    disp(['Run ', int2str(i), ': ', num2str(bottom)])

end

%%% Calculate / plot threshold errors
save([base, '/thresholds_IKr.mat'], 'predictedThreshold')

%% Comparisons & error calculations
clear
close all
f_out = ['Analysis/Thresholds', '/cell1_4-19-5'] ;
cell_number = 1 ;
realdata = false ;

folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select DIRECTORY/IES from experiments to compare') ;
numruns = 10 ;

protocols = 1:length(folders) ;
xlabels = strsplit(num2str(protocols), ' ') ;
set(groot, 'defaultFigureRenderer', 'painters') % for saving as svg

if realdata % Experimental data thresholds
    actual_Ko = 0.7 ; % replace w/ measured threshold
    actual_IKr = 0.5 ; % replace w/ measured threshold
else % Pseudodata thresholds
    actual_Ko_thresholds = readtable('Pseudodataset/saved_data/Ko thresholds.xlsx') ;
    actual_Ko = actual_Ko_thresholds.Threshold(cell_number) ;
    actual_IKr_thresholds = readtable('Pseudodataset/saved_data/IKr block thresholds.xlsx') ;
    actual_IKr = actual_IKr_thresholds.Threshold(cell_number) ;
end

predicted_Ko = zeros(length(folders), numruns) ;
predicted_IKr = zeros(length(folders), numruns) ;
experiments = cell(1, length(folders)) ;
for i=1:length(folders)
    load([folders{i}, '/thresholds_IKr.mat'], 'predictedThreshold')
    predicted_IKr(i, :) = predictedThreshold ;
    if ispc % Windows
        foldersplit = strsplit(folders{i}, '\') ;
    else
        foldersplit = strsplit(folders{i}, '/') ;    
    end
    experiments{i} = foldersplit{end} ;
end

handle_IKr = figure ;
hold on
UnivarScatter(predicted_IKr', 'Label', xlabels, 'Width', 0.5, 'Whiskers', 'lines', ...
    'WhiskersWidthRatio', 1.25);
yline(actual_IKr, 'g-', 'LineWidth', 2, 'HandleVisibility', 'off')
xlim([0, length(folders)+1])
xticks(protocols)
xticklabels(protocols)
xlabel("Protocols")
ylim([0, 1])
ylabel("% I_{Kr} threshold")
legend(experiments, 'Location', 'bestoutside')
savefig([f_out, '/IKr_threshold.fig'])
print('-dpng', [f_out, '/IKr_threshold'])
print('-dsvg', [f_out, '/IKr_threshold'])

% Calculations
IKr_errors = (predicted_IKr - actual_IKr).^2 ;
IKr_MSE = mean(IKr_errors, 2) ;
IKr_std = std(predicted_IKr, 0, 2) ; 

save([f_out, '/IKr_calculations.mat'], 'IKr_errors', 'IKr_MSE', 'IKr_std', 'experiments')

