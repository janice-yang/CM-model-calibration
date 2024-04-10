clear
close all

realdata = false ; 
filename = '20220113_paced.xlsx' ; % experimental data file
sheetnames = {'1.0Ca 1Hz','1.8Ca 1.25Hz'} ;
isNormalized = false ;
if isNormalized
    scaleCaT = false ;
else
    scaleCaT = true ;
end
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
save GA/curr_cell_protocol protocol_number isNormalized scaleCaT
load GA/x_names.mat

% For storing experimental data
expT = cell(1, length(protocol_number)) ;
expV = cell(1, length(protocol_number)) ;
expCai = cell(1, length(protocol_number)) ;

%% Experimental data
extracted = figure ;
set(gcf, 'Position', [0 0 length(protocol_number)*400 400])
hold on
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
    disp('Experimental data loaded')
else 
    load([folders{1}, '/Details.mat'], 'cell_number', 'protocol_number') ;
    [experimental_dataset] = f_getPseudodata(cell_number, protocol_number, isNormalized, 0, 0) ;
    conds = readmatrix('Pseudodataset/saved_data/ground_truth_conductances.xlsx') ;
    x_conductance = conds(cell_number, :) ;  
    [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names) ;
    for j=1:length(protocol_number)
        t = experimental_dataset{j}.Time_AP{1} ;
        V = experimental_dataset{j}.AP{1} ;
        Cai = experimental_dataset{j}.CaT{1} ;
        [t, V, Cai,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},stimtimes{j});
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
    disp('Pseudodata loaded')
end

% Predicted conductances
fitness_all = zeros(length(folders)+1, 4) ;
for i=1:length(folders)
    single = figure ;
    load([folders{i}, '/Details.mat'], 'runNum') ;
    load([folders{i}, '/minoptimparams.mat'], 'minoptimparams') ;
    x_conductance = logfactor.^minoptimparams ;
    [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names) ;
    y = zeros(length(protocol_number), 2) ; % store fitness
    for j=1:length(protocol_number)
        if sum(t_stim{j})
        [t, V, Cai,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},stimtimes{j});
        
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
            if scaleCaT
                % Average of 600ms PCL, 800ms PCL, and spontaneous factors
                CaT_factor = (248317.2776 + 275865.6087 + 258347.0044) / 3 ;
                CaT_sim = CaT_sim*CaT_factor ;
                exp_CaT = exp_CaT*CaT_factor ;

                % Record new range(V_sim)/range(CaT_sim)
                AP_CaT_ratio = range(V_sim) / range(CaT_sim) ;
                protocol = protocol_number(j) ;
                tbl = table(AP_CaT_ratio, protocol, x_conductance) ; % change to actual conductances (scale baseline by x_conductance)
                writetable(tbl, [folders{i}, '/fitting_CaT_scale_factors.txt'], 'WriteMode','Append')
            end
            y(j,2)= sum((CaT_sim-exp_CaT).^2) / length(CaT_sim) ;
        end
        else
            figure(single)
            subplot(2, length(protocol_number)+1, j) % AP, protocol i
            hold on
            plot(expT{j}, expV{j}, 'r-', 'LineWidth', 1)
            yline(0, 'b-', 'LineWidth', 1)
            title(['Protocol ', int2str(protocol_number(j))])
            xlim([0 max(expT{j})])
            xlabel("Time (ms)")
            ylabel({'AP', 'mV'})

            subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT protocol i
            hold on
            plot(expT{j}, expCai{j}, 'r-', 'LineWidth', 1)
            yline(0, 'b-', 'LineWidth', 1)
            xlim([0 max(expT{j})])
            xlabel('Time (ms)')
            ylabel({'[Ca^{2+}]_i', 'mM'})

            figure(extracted)
            subplot(2, length(protocol_number)+1, j) % AP, protocol j
            hold on
%             plot(expT_extract, expV_extract, 'r--', 'LineWidth', 2)
            yline(0)
            title(['Protocol ', int2str(protocol_number(j))])
            xlabel("Time (ms)")
            ylabel({'AP', 'mV'})

            subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT protocol i
            hold on
%             plot(expT_extract, expCai_extract, 'r--', 'LineWidth', 2)
            yline(0)
            xlabel('Time (ms)')
            ylabel({'[Ca^{2+}]_i', 'mM'})

            y(j,1) = 1e10 ;
            y(j,2) = 1e10 ;
            warning(['Simulation error during protocol ',num2str(protocol_number(j)),...
                ' in run ', int2str(runNum)])
        end
            
    end
    
    figure(single)
    if isNormalized
        savefig([base, '/predicted_waveforms_norm_Run', int2str(runNum)]) ;
    else
        savefig([base, '/predicted_waveforms_raw_Run', int2str(runNum)]) ;
    end
    close(single)
    
    % Average and total fitness score for run i
    % Save average, total, & individual fitness scores in run folder
    avg_fitness = mean(y, 2) ; % average (AP+CaT/2)
    avg_APCaT = mean(y, 1) ; % avg_AP, avgCaT across all protocols
    total_fitness = sum(y, 'all') ;
    save([folders{i}, '/fitness_fitting_protocol.mat'],'y','avg_fitness','avg_APCaT','total_fitness')
    fitness_all(runNum+1,1:2) = avg_APCaT ;
    fitness_all(runNum+1,3) = sum(avg_APCaT) ;
    fitness_all(runNum+1,4) = runNum ;
    disp(['Run ' int2str(runNum) ' finished.'])
end

% Sort & Average fitness across all runs in base folder
[~, idx] = sort(fitness_all(1:end-1,3), 'ascend') ;
fitness_all(1:end-1, :) = fitness_all(idx, :) ;
fitness_all(end,:) = mean(fitness_all, 1) ;
tbl = array2table(fitness_all, 'VariableNames', {'AP', 'CaT', 'average', 'runNum'}) ;
if isNormalized
    writetable(tbl, [base, '/ranked_all_norm.xlsx'], 'Sheet', 'fitting protocol')
else
    writetable(tbl, [base, '/ranked_all_raw.xlsx'], 'Sheet', 'fitting protocol')
end
figure(extracted)
lSub = subplot(2, length(protocol_number)+1, length(protocol_number)+1);
hold on
plot(1, 1, 'r--');
set(lSub, 'Visible', 'off');
% legend(lSub, 'Experiment', 'Baseline model', 'Location', 'best')
legend(lSub, 'Experiment', 'Location', 'best')
if isNormalized
    savefig([base, '/predicted_waveforms_norm.fig'])
else
    savefig([base, '/predicted_waveforms_raw.fig'])
end
disp('Fitting protocol plots/calculations done.')

%% Validation run
fitness_all = zeros(length(folders)+1, 4) ; % reset for validation
if realdata
%     exp_file = uigetfile_n_dir([pwd, '/ExperimentalData/'], 'Select file containing VALIDATION data') ;
    exp_file = {'F:\My Drive\Sobie lab\ivs-tox\ExperimentalData\20220113_1.0Ca_2Hz_valid.xlsx'} ;
    sheetnames = {'Sheet1'} ; % each sheet = 1 validation protocol
    protocol_number = 29 ;

    exp_file = exp_file{1} ;
    datamatrix = readmatrix(exp_file) ;
    stimV = datamatrix(:, end) ;
else
    load([folders{1}, '/Details.mat'], 'cell_number')
    protocol_number = [12] ; % validation protocol(s)
end
sequence_str = replace(num2str(protocol_number), '  ', '-') ;
sequence = {sequence_str} ;
if ispc % Windows
    foldersplit = strsplit(folders{1}, '\') ;
else
    foldersplit = strsplit(folders{1}, '/') ;    
end
base = strjoin(foldersplit(1:end-1), '/') ; % experiment folder with all runs

save GA/curr_cell_protocol protocol_number isNormalized scaleCaT
load 'GA/x_names.mat' names

colors = repmat('krgbcm', 1, 300) ;
if realdata
    all_params = zeros(length(folders), length(names)) ;
    for i=1:length(folders)
        load([folders{i}, '/minoptimparams.mat'], 'minoptimparams')
        all_params(i, :) = minoptimparams ;
    end

    logfactor = 2 ; 
    x_conductances = (logfactor.^all_params);
    for i=1:length(folders)
        folder = folders{i} ; 
        load([folder, '/Details.mat'], 'runNum')
        x_conductance = x_conductances(i, :) ;
        [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names); 
        y = zeros(length(protocol_number), 2) ;

        for j=1:length(protocol_number)
            if sum(t_stim{j})
                % Predicted for protocol i
                simT = t_stim{j} ;
                simV = V_stim{j} ;
                simCai = Cai_stim{j} ;
                sim_stimtimes = stimtimes{j} ;

                % Experimental for protocol i
                expdata = readmatrix(exp_file, 'Sheet', sheetnames{j}) ;
                expT = expdata(:,1) ;
                stim_idx = ismembertol(expdata(:,2), stimV, 2.5e-7) ;
                exp_stimtimes = expT(stim_idx) ;
                if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
                    expV = expdata(:, 2) ;
                    % Subtract and smooth baseline
                    erodedSignal = imerode(expV, ones(100, 1)) ; 
                    expV = expV - erodedSignal ;
                    expV = medfilt1(expV, 10) ;
                    if isNormalized
                        expV = (expV - min(expV)) ./ max(expV - min(expV)) ;
                    end
                end

                if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
                    expCai = expdata(:, 3) ;
                    erodedSignal = imerode(expCai, ones(100, 1)) ; 
                    expCai = expCai - erodedSignal ;
                    expCai = medfilt1(expCai, 10) ;
                    if isNormalized
                        expCai = (expCai - min(expCai)) ./ max(expCai - min(expCai)) ;
                    end
                end
                
                % Align by first AP upstroke
                if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
                    % Experimental
                    Vderiv = [expV(2:end);expV(end)] - expV ;
                    [~,dexmax] = maxk(Vderiv,5) ;
                    expT_init = expT(min(dexmax)) ;
                    expT = expT - expT_init ;
                    exp_stimtimes = exp_stimtimes - expT_init ;
                    % Predicted
                    Vderiv = [simV(2:end);simV(end)] - simV ;
                    [~,dexmax] = maxk(Vderiv,5) ;
                    simT_init = simT(min(dexmax)) ;
                    simT = simT - simT_init ;
                else
                    % Experimental
                    Cderiv = [expCai(2:end);expCai(end)] - expCai ;
                    [~,dexmax] = maxk(Cderiv,5) ;
                    expT_init = expT(min(dexmax)) ;
                    expT = expT - expT_init ;
                    exp_stimtimes = exp_stimtimes - expT_init ;
                    % Predicted
                    Cderiv = [simCai(2:end);simCai(end)] - simCai ;
                    [~,dexmax] = maxk(Cderiv,5) ;
                    simT_init = simT(min(dexmax)) ;
                    simT = simT - simT_init ;                    
                end

                % Extracted traces
                if strcmp(datatype, 'APCaT')
                    [expT_extract, expV_extract, expCai_extract,tinit,errorcode] = waveform_extract_new(expT,expV,expCai,exp_stimtimes);
                elseif strcmp(datatype, 'AP')
                    [expT_extract,expV_extract,errorcode] = APextract_new(expT,expV,exp_stimtimes) ;
                else
                    [expT_extract,expV_extract,errorcode] = CaTextract_custom(expT,expCai,exp_stimtimes) ;
                end
                [simT_extract, simV_extract, simCai_extract,tinit,errorcode] = waveform_extract_new(simT,simV,simCai,sim_stimtimes);
                
                % Fitness calculation
                if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
                    [expT_V, expV, simT_V, simV] = f_alignWaveformEnds(expT_extract, expV_extract, simT_extract, simV_extract) ;
                    V_sim = interp1(simT_V, simV, expT_V);
                    y(j,1) =sum((V_sim - expV).^2)/length(V_sim);
                end
                if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
                    [expT_Cai, expCai, simT_Cai, simCai] = f_alignWaveformEnds(expT_extract, expCai_extract, simT_extract, simCai_extract) ;
                    CaT_sim = interp1(simT_Cai, simCai, expT_Cai);
                    y(j,2)= sum((CaT_sim - expCai).^2)/length(CaT_sim);
                end

            else % If Model Crashes           
                y(j,1) = 1e10 ;
                y(j,2) = 1e10 ;
                warning(['Simulation error during protocol ',num2str(protocol_number(j)),...
                    ' in run ', int2str(runNum)])%                  close all
%                  error(['Simulation error during protocol ',num2str(protocol_number(j))])
            end
        end
        
        % Fitness
        avg_fitness = mean(y, 2) ; % average (AP+CaT/2)
        avg_APCaT = mean(y, 1) ; % avg_AP, avgCaT across all protocols
        total_fitness = sum(y, 'all') ;
        save([folders{i}, '/fitness_val_protocol.mat'],'y','avg_fitness','avg_APCaT','total_fitness','protocol_number')
        fitness_all(runNum+1,1:2) = avg_APCaT ;
        fitness_all(runNum+1,3) = sum(avg_APCaT) ;
        fitness_all(runNum+1,4) = runNum ;

        disp(['Run ' int2str(runNum) ' validation done.'])
    end

else     % Pseudodata
    load Pseudodataset/ranges ranges
    all_params = zeros(length(folders), length(names)) ;
    for i=1:length(folders)
        load([folders{i}, '/minoptimparams.mat'], 'minoptimparams')
        all_params(i, :) = minoptimparams ;
    end

    logfactor = 2 ; 
    x_conductances = (logfactor.^all_params);
    for i=1:length(folders)
        folder = folders{i} ; 
        load([folder, '/Details.mat'], 'runNum')
        x_conductance = x_conductances(i, :) ;
        [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names); 
        % Simulate for Kernik baseline conductances
%         [base_t_stim, base_V_stim, base_Cai_stim, base_stimtimes] = ga_simulation_k19(ones(1, length(names)), names); 
        y = zeros(length(protocol_number), 2) ;
        
        for j=1:length(protocol_number)
            % Experimental for protocol j
            expT = readmatrix('Pseudodataset/saved_data/pseudodataset.xlsx',...
                'Sheet',['Cell ',int2str(cell_number)],'Range','A3:A50003');
            expVCaT = readmatrix('Pseudodataset/saved_data/pseudodataset.xlsx',...
                'Sheet',['Cell ',int2str(cell_number)],'Range',ranges{protocol_number(j)});
            expV = expVCaT(:, 1) ;
            expCai = expVCaT(:, 2) ;
            if isNormalized
                expV = (expV - min(expV)) ./ max(expV - min(expV)) ;
                expCai = (expCai - min(expCai)) ./ max(expCai - min(expCai)) ;
            end

            % Align by first AP upstroke
            % % Experimental
            Vderiv = [expV(2:end);expV(end)] - expV ;
            [~,dexmax] = maxk(Vderiv,5) ;
            expT_init = expT(min(dexmax)) ;
            expT = expT - expT_init ;
            
            % Extracted
            if ismember(protocol_number(j), [1,3,5,7,9,11,13,15,18,27,30]) % 800ms PCL
                exp_stimtimes = [2896,3696,4496] ;
            elseif protocol_number(j) == 17
                exp_stimtimes = [3496, 4096, 4696] ;
            elseif ismember(protocol_number(j), [21,23,25]) % 1000ms PCL
                exp_stimtimes = [2496, 3496, 4496] ;
            elseif ismember(protocol_number(j), [26,29]) % 500ms PCL
                exp_stimtimes = [3499.5, 3999.5, 4499.5] ;
            elseif ismember(protocol_number(j), [28,31]) % 2000ms PCL
                exp_stimtimes = [496, 2496, 4496] ;
            else % spontaneous beating
                exp_stimtimes = [] ;
            end
            [expT_extract, expV_extract, expCai_extract,tinit,errorcode] = waveform_extract_new(expT,expV,expCai,exp_stimtimes);
            
            if (sum(t_stim{j}))
                % Predicted for protocol i
                simT = t_stim{j} ;
                simV = V_stim{j} ;
                simCai = Cai_stim{j} ;
                sim_stimtimes = stimtimes{j} ;

                % % Predicted
                Vderiv = [simV(2:end);simV(end)] - simV ;
                [~,dexmax] = maxk(Vderiv,5) ;
                simT_init = simT(min(dexmax)) ;
                simT = simT - simT_init ;

                % Fitness calculation
                [expT_V, expV, simT_V, simV] = f_alignWaveformEnds(expT, expV, simT, simV) ;
                [expT_Cai, expCai, simT_Cai, simCai] = f_alignWaveformEnds(expT, expCai, simT, simCai) ;
                V_sim = interp1(simT_V, simV, expT_V);
                CaT_sim = interp1(simT_Cai, simCai, expT_Cai);
                
                if scaleCaT
                    % Scale CaT using factor from baseline Kernik APA/CaTA ratio
                    CaT_factor = (248317.2776 + 275865.6087 + 258347.0044) / 3 ;
                    CaT_sim = CaT_sim*CaT_factor ;
                    expCai = expCai*CaT_factor ;
                end
                
                y(j,1) = sum((V_sim - expV).^2)/length(V_sim);
                y(j,2)= sum((CaT_sim - expCai).^2)/length(CaT_sim);

            else % If Model Crashes - stop run                
                y(j,1) = 1e10 ;
                y(j,2) = 1e10 ;
%                  close all
                 warning(['Simulation error during protocol ',num2str(protocol_number(j))])
            end
        end
        
        % Fitness
        avg_fitness = mean(y, 2) ; % average (AP+CaT/2)
        avg_APCaT = mean(y, 1) ; % avg_AP, avgCaT across all protocols
        total_fitness = sum(y, 'all') ;
        save([folders{i}, '/fitness_val_protocol.mat'],'y','avg_fitness','avg_APCaT','total_fitness')
        fitness_all(runNum+1,1:2) = avg_APCaT ;
        fitness_all(runNum+1,3) = sum(avg_APCaT) ;
        fitness_all(runNum+1,4) = runNum ;
        
        disp(['Run ' int2str(runNum) 'validation done.'])
    end
end

% Sort & Average fitness across all runs in base folder
[~, idx] = sort(fitness_all(1:end-1,3), 'ascend') ;
fitness_all(1:end-1, :) = fitness_all(idx, :) ;
fitness_all(end,:) = mean(fitness_all, 1) ;
tbl = array2table(fitness_all, 'VariableNames', {'AP', 'CaT', 'sum', 'runNum'}) ;

disp('Validation plots/calculations done.')
