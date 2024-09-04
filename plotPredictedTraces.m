clear
close all

realdata = true ; 
% filename = '20220113_paced.xlsx' ; % experimental data file
% filename = 'DMG240.xlsx' ;
filename = 'TestDavidDataProcessed2.xlsx' ;
cell_number = 12 ;
protocol_number = [35] ;
% sheetnames = {'1.0Ca 1Hz','1.8Ca 1Hz','1.8Ca 1.25Hz'} ;
% sheetnames = {'1.8Ca 1Hz'}
% sheetnames = {'DMG240_70Na_1Hz'} ;
% sheetnames = {'DMG240_50pNifedipine50nM_1Hz'} ;
sheetnames = {'DMG240_50pDofetilide1nM_0.5Hz'}
nbeats = [15] ;
n_extract = 2 ;
isNormalized = true ;  
datatype = 'APCaT' ; 
logfactor = 2 ;

if ~(strcmp(datatype,'AP') || strcmp(datatype,'CaT') || strcmp(datatype,'APCaT'))
    error("Please enter a valid datatype - 'AP', 'CaT', or 'APCaT'")
end

folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select DIRECTORY/IES containing GA results') ;
if ispc % Windows
    foldersplit = strsplit(folders{1}, '\') ;
else
    foldersplit = strsplit(folders{1}, '/') ;    
end
base = strjoin(foldersplit(1:end-1), '/') ; % experiment folder with all runs
protocol_str = join(split(num2str(protocol_number), '  '), '-') ;
protocol_str = protocol_str{1} ;

% load([folders{1}, '/Details.mat'], 'protocol_number') ;
save GA/curr_cell_protocol.mat protocol_number isNormalized
load GA/x_names.mat names

% For storing experimental data
expT = cell(1, length(protocol_number)) ;
expV = cell(1, length(protocol_number)) ;
expCai = cell(1, length(protocol_number)) ;

set(groot, 'defaultFigureRenderer', 'painters') % for saving as svg

%%
extracted = figure ;
hold on

% Plot experimental
if realdata
    [experimental_dataset] = f_getExperimentData(filename, sheetnames, protocol_number, isNormalized, datatype, nbeats) ;

    % Plot
    for j=1:length(protocol_number)       
        figure(extracted)
        t = experimental_dataset{j}.Time_AP{1} ;
        V = experimental_dataset{j}.AP{1} ;
        Cai = experimental_dataset{j}.CaT{1} ;
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            subplot(2, length(protocol_number)+1, j) % AP
            hold on
            plot(t, V, 'r--', 'LineWidth', 2)
            xlabel("Time (ms)")
            ylabel("mV")
        end
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT
            hold on
            plot(t, Cai, 'r--', 'LineWidth', 2)
            xlabel("Time (ms)")
            ylabel("mM")
        end
        
        % Store data for fitness calculation
        expT{j} = t ;
        expV{j} = V ;
        expCai{j} = Cai ;
    end
else 
    % load([folders{1}, '/Details.mat'], 'cell_number', 'protocol_number') ;
    [experimental_dataset] = f_getPseudodata(cell_number, protocol_number, isNormalized, 0, 0) ;
    set(gcf, 'Position', [0 0 length(protocol_number)*400 400])
    for j=1:length(protocol_number)
        t = experimental_dataset{j}.Time_AP{1} ;
        V = experimental_dataset{j}.AP{1} ;
        Cai = experimental_dataset{j}.CaT{1} ;

        subplot(2, length(protocol_number)+1, j) % AP
        hold on
        plot(t, V, 'r--', 'LineWidth', 2)
        xlabel("Time (ms)")
        ylabel({"AP","mV"})
        subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT
        hold on
        plot(t, Cai, 'r--', 'LineWidth', 2)
        xlabel("Time (ms)")
        ylabel({"Ca_i","mM"})

        % Store data for fitness calculation
        expT{j} = t ;
        expV{j} = V ;
        expCai{j} = Cai ;
    end
end

% Predicted conductances
fitness_all = zeros(length(folders)+1, 3) ;
for i=1:length(folders)
    load([folders{i}, '/Details.mat'], 'runNum') ;
    load([folders{i}, '/minoptimparams.mat'], 'minoptimparams') ;
    x_conductance = logfactor.^minoptimparams ;
    [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names) ;
    y = zeros(length(protocol_number), 2) ; % store fitness
    single = figure ;
    for j=1:length(protocol_number)
        if sum(t_stim{j})
            if ~realdata
                [t, V, Cai,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},stimtimes{j},n_extract);
            else
                [t, V, Cai,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},[],n_extract);
            end
        figure(single)
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            subplot(2, length(protocol_number)+1, j) % AP
            hold on
            plot(t, V)
            plot(expT{j}, expV{j}, 'r--', 'LineWidth', 2)
        end
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT
            hold on
            plot(t, Cai)
            plot(expT{j}, expCai{j}, 'r--', 'LineWidth', 2)
        end
        
        figure(extracted)
        subplot(2, length(protocol_number)+1, j) % AP
        hold on
        plot(t, V)
        subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT
        hold on
        plot(t, Cai)
        
        % Fitness calculation for protocol j
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
        [exp_T_V, exp_V, T_V, V] = f_alignWaveformEnds(expT{j}, expV{j}, t, V) ;
        V_sim = interp1(T_V, V, exp_T_V);
        y(j,1) = sum((V_sim-exp_V).^2) / length(V_sim) ;
        end
        
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            [exp_T_CaT, exp_CaT, T_CaT, CaT] = f_alignWaveformEnds(expT{j}, expCai{j}, t, Cai) ;
            CaT_sim = interp1(T_CaT, CaT, exp_T_CaT);

            % Scale CaT using ratio from experimental APA/CaTA
            CaT_factor = range(exp_V) / range(exp_CaT) ;
            
            % Record new range(V_sim)/range(CaT_sim)
            AP_CaT_ratio = range(V_sim) / range(CaT_sim) ;
            protocol = protocol_number(j) ;
            tbl = table(AP_CaT_ratio, protocol, x_conductance) ; % change to actual conductances (scale baseline by x_conductance)
            writetable(tbl, [folders{i}, '/predicted_CaT_scale_factors_', protocol_str, '.txt'], 'WriteMode','Append')
            y(j,2)= sum((CaT_sim-exp_CaT).^2) / length(CaT_sim) ;
        end
        
        end
    
    end
    figure(single)
    if isNormalized
        savefig([base, '/predicted_traces_norm_', protocol_str, '_Run', int2str(runNum)]) ;
    else
        savefig([base, '/predicted_traces_raw_', protocol_str, '_Run', int2str(runNum)]) ;
    end
    close(single)
    % Average and total fitness score for run i
    % Save average, total, & individual fitness scores in run folder
    avg_protocols = mean(y, 2) ; % average (AP+CaT/2) per protocol
    avg_APCaT = mean(y, 1) ; % avg_AP, avgCaT across all protocols
    total_fitness = sum(y, 'all') ;
    save([folders{i}, '/fitness_predicted_traces_', protocol_str, '.mat'],'y','avg_protocols','avg_APCaT','total_fitness')
    fitness_all(runNum+1,1:2) = avg_APCaT ;
    fitness_all(runNum+1,3) = sum(avg_APCaT) ;
end

% Average fitness across all runs in base folder
fitness_all(end,:) = mean(fitness_all, 1) ;
tbl = array2table(fitness_all, 'VariableNames', {'AP', 'CaT', 'average'}) ;
if isNormalized
    writetable(tbl, [base, '/fitness_predicted_norm_', protocol_str, '.xlsx'])
else
    writetable(tbl, [base, '/fitness_predicted_raw_', protocol_str, '.xlsx'])
end

figure(extracted)
lSub = subplot(2, length(protocol_number)+1, length(protocol_number)+1); 
hold on
plot(1, 1, 'r--');
set(lSub, 'Visible', 'off');
legend(lSub, 'Experiment', 'Location', 'best')

% Save
if isNormalized
    print('-dpng', [base, '/predicted_traces_norm_', protocol_str])
    print('-dsvg', [base, '/predicted_traces_norm_', protocol_str])
    savefig([base, '/predicted_traces_norm_', protocol_str, '.fig'])
else
    print('-dpng', [base, '/predicted_traces_raw_', protocol_str])
    print('-dsvg', [base, '/predicted_traces_raw_', protocol_str])
    savefig([base, '/predicted_traces_raw_', protocol_str, '.fig'])
end