clear
close all

realdata = false ; 
filename = '20220113_paced.xlsx' ; % experimental data file
% sheetnames = {'1.0Ca 1Hz','1.8Ca 1Hz','1.8Ca 1.25Hz'} ;
sheetnames = {'1.8Ca 1.25Hz'} ;
isNormalized = false ;
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
load([folders{1}, '/Details.mat'], 'protocol_number') ;
save GA/curr_cell_protocol.mat protocol_number isNormalized
load GA/x_names.mat

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
    % For normalization through multiple protocols
    V_stim_all = []; % full AP data for normalization
    Cai_stim_all = []; % full CaT data for normalation
    ends = zeros(1, length(protocol_number)) ; % to keep track of protocol indices for normalization

    % Plot
    for j=1:length(protocol_number)
        data = readmatrix(['ExperimentalData/', filename], ...
            'Sheet', sheetnames{j}) ; 
        time = data(:, 1) ;
        stimV = data(:, end) ;
        stim_idx = ismembertol(data(:,2), stimV, 2.5e-7) ;
        stimtimes = time(stim_idx) ;
        
        % Noise processing
        if strcmp(datatype, 'APCaT') || strcmp(datatype, 'AP')
            tempV = data(:, 2) ;
            erodedSignal = imerode(tempV, ones(100, 1)) ; 
            tempV = tempV - erodedSignal ;
            tempV = medfilt1(tempV, 10) ;
        end
        if strcmp(datatype, 'APCaT') || strcmp(datatype, 'CaT')
            tempCaT = data(:, 3) ;
            erodedSignal = imerode(tempCaT, ones(100, 1)) ; 
            tempCaT = tempCaT - erodedSignal ;
            tempCaT = medfilt1(tempCaT, 10) ;
        end
        if strcmp(datatype, 'APCaT')
            [t, V, Cai,tinit,errorcode] = waveform_extract_new(time, tempV, tempCaT, stimtimes);
        elseif strcmp(datatype, 'AP')
            [t,V,errorcode] = APextract_new(time,tempV,stimtimes) ;
        else
            [t,Cai,errorcode] = CaTextract_custom(time,tempCaT,stimtimes) ;
        end
        
        expT{j} = t ;
        V_stim_all = [V_stim_all; V] ;
        Cai_stim_all = [Cai_stim_all; Cai] ;
        if j==1
            ends(j) = length(t) ;
        else
            ends(j) = length(t) + (ends(j-1)) ;
        end
    end
    if isNormalized
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            V_stim_all = (V_stim_all - min(V_stim_all)) ./ max(V_stim_all - min(V_stim_all)) ;
        end
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            Cai_stim_all = (Cai_stim_all - min(Cai_stim_all)) ./ max(Cai_stim_all - min(Cai_stim_all)) ;
        end
    end
    
    for j=1:length(protocol_number)
        if j==1
            idx_start = 1 ;
            idx_end = ends(j) ;
        else
            idx_start = ends(j-1) + 1 ;
            idx_end = ends(j) ;
        end
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            expV{j} = V_stim_all(idx_start:idx_end);
        end
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            expCai{j} = Cai_stim_all(idx_start:idx_end) ;
        end
        
        figure(extracted)
        t = expT{j} ;
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            V = expV{j} ;
            subplot(2, length(protocol_number)+1, j) % AP
            hold on
            plot(t, V, 'r--', 'LineWidth', 2)
            xlabel("Time (ms)")
            ylabel("mV")
        end
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            Cai = expCai{j} ;
            subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT
            hold on
            plot(t, Cai, 'r--', 'LineWidth', 2)
            xlabel("Time (ms)")
            ylabel("mM")
        end
    end
else 
    load([folders{1}, '/Details.mat'], 'cell_number', 'protocol_number') ;
    [experimental_dataset] = f_getPseudodata(cell_number, protocol_number, isNormalized, 0, 0) ;
    set(gcf, 'Position', [0 0 length(protocol_number)*400 400])
%     conds = readmatrix('Pseudodataset/saved_data/ground_truth_conductances.xlsx') ;
%     x_conductance = conds(cell_number, :) ;  
%     [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names) ;
    for j=1:length(protocol_number)
        t = experimental_dataset{j}.Time_AP{1} ;
        V = experimental_dataset{j}.AP{1} ;
        Cai = experimental_dataset{j}.CaT{1} ;
%         [t, V, Cai,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},stimtimes{j});

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

        %% Store data for fitness calculation
        expT{j} = t ;
        expV{j} = V ;
        expCai{j} = Cai ;
    end
end
    
% % Plot baseline
% x_conductance = ones(1, length(names)) ;
% [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names) ;
% for j=1:length(protocol_number)
%     [t, V, Cai,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},stimtimes{j});
%     
%     subplot(2, length(protocol_number)+1, j) % AP
%     hold on
%     plot(t, V, 'g--', 'LineWidth', 2)
%     subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT
%     hold on
%     plot(t, Cai, 'g--', 'LineWidth', 2)
%     
% end

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
                [t, V, Cai,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},stimtimes{j});
            else
                [t, V, Cai,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},[]);
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
        
        %% Fitness calculation for protocol j
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
        [exp_T_V, exp_V, T_V, V] = f_alignWaveformEnds(expT{j}, expV{j}, t, V) ;
        V_sim = interp1(T_V, V, exp_T_V);
        y(j,1) = sum((V_sim-exp_V).^2) / length(V_sim) ;
        end
        
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            [exp_T_CaT, exp_CaT, T_CaT, CaT] = f_alignWaveformEnds(expT{j}, expCai{j}, t, Cai) ;
            CaT_sim = interp1(T_CaT, CaT, exp_T_CaT);
            % Average of 600ms PCL, 800ms PCL, and spontaneous factors
            CaT_factor = (248317.2776 + 275865.6087 + 258347.0044) / 3 ;
            CaT_sim = CaT_sim*CaT_factor ;
            exp_CaT = exp_CaT*CaT_factor ;
            
            % Record new range(V_sim)/range(CaT_sim)
            AP_CaT_ratio = range(V_sim) / range(CaT_sim) ;
            protocol = protocol_number(j) ;
            tbl = table(AP_CaT_ratio, protocol, x_conductance) ; % change to actual conductances (scale baseline by x_conductance)
            writetable(tbl, [folders{i}, '/predicted_waveforms_CaT_scale_factors.txt'], 'WriteMode','Append')
            y(j,2)= sum((CaT_sim-exp_CaT).^2) / length(CaT_sim) ;
        end
        
        end
    
    end
    figure(single)
    if isNormalized
        savefig([base, '/predicted_waveforms_norm_Run', int2str(runNum)]) ;
        print('-dsvg', [base, '/predicted_waveforms_norm_Run', int2str(runNum)])
        print('-dpng', [base, '/predicted_waveforms_norm_Run', int2str(runNum)])
    else
        savefig([base, '/predicted_waveforms_raw_Run', int2str(runNum)]) ;
        print('-dsvg', [base, '/predicted_waveforms_raw_Run', int2str(runNum)])
        print('-dpng', [base, '/predicted_waveforms_raw_Run', int2str(runNum)])
    end
    close(single)
    % Average and total fitness score for run i
    % Save average, total, & individual fitness scores in run folder
    avg_protocols = mean(y, 2) ; % average (AP+CaT/2) per protocol
    avg_APCaT = mean(y, 1) ; % avg_AP, avgCaT across all protocols
    total_fitness = sum(y, 'all') ;
    save([folders{i}, '/fitness_simulated.mat'],'y','avg_protocols','avg_APCaT','total_fitness')
    fitness_all(runNum+1,1:2) = avg_APCaT ;
    fitness_all(runNum+1,3) = sum(avg_APCaT) ;
end

%% Average fitness across all runs in base folder
fitness_all(end,:) = mean(fitness_all, 1) ;
tbl = array2table(fitness_all, 'VariableNames', {'AP', 'CaT', 'average'}) ;
if isNormalized
    writetable(tbl, [base, '/fitness_simulated_all_norm.xlsx'])
else
    writetable(tbl, [base, '/fitness_simulated_all_raw.xlsx'])
end

figure(extracted)
lSub = subplot(2, length(protocol_number)+1, length(protocol_number)+1); 
hold on
plot(1, 1, 'r--');
% plot(1, 1, 'g--');
set(lSub, 'Visible', 'off');
% legend(lSub, 'Experiment', 'Baseline model', 'Location', 'best')
legend(lSub, 'Experiment', 'Location', 'best')

%% Save
if isNormalized
    print('-dpng', [base, '/predicted_waveforms_norm'])
    print('-dsvg', [base, '/predicted_waveforms_norm'])
    savefig([base, '/predicted_waveforms_norm.fig'])
else
    print('-dpng', [base, '/predicted_waveforms_raw'])
    print('-dsvg', [base, '/predicted_waveforms_raw'])
    savefig([base, '/predicted_waveforms_raw.fig'])
end
