function [] = f_validation_run(folders, cell_number, protocol_number, method, isNormalized)

% Validate fitted parameters on unseen experiments
%{
Inputs:
    folders: cell array of directories containing individual run results
    cell_number: experimental cell ID
    protocol_number: array of validation experiments (can be multiple)
    method ('single', 'median', or 'mean'): which parameter(s) to use for
        validation (single = from individual runs, median/mean = summary from
        all runs)

Protocols for validation:
Ko low --> normal --> high, PCL 800ms: 15,18,13
Ko low --> normal --> high, spontaneous: 16,19,14
IKr block 0% --> 15% --> 30%, PCL 800ms: 18,9,11
IKr block 0% --> 15% --> 30%, spontaneous: 19,10,12
%}

scaleCaT = false ;
sequence_str = replace(num2str(protocol_number), '  ', '-') ;
sequence = {sequence_str} ;
if ispc % Windows
    foldersplit = strsplit(folders{1}, '\') ;
else
    foldersplit = strsplit(folders{1}, '/') ;    
end
base = strjoin(foldersplit(1:end-1), '/') ; % experiment folder with all runs

save GA/curr_cell_protocol cell_number protocol_number isNormalized scaleCaT

load 'GA/x_names.mat' names
load Pseudodataset/ranges ranges
all_params = zeros(length(folders), length(names)) ;
for i=1:length(folders)
    load([folders{i}, '/minoptimparams.mat'], 'minoptimparams')
    all_params(i, :) = minoptimparams ;
end

logfactor = 2 ; 
x_conductances = (logfactor.^all_params);

% % For features calculation
colnames = strsplit(num2str(protocol_number)) ; % e.g. {'1','2','3'}
rownames = cell(1,length(folders)+4) ;
freq_all = zeros(length(folders)+4, length(protocol_number)) ;
APD90_all = zeros(length(folders)+4, length(protocol_number)) ;
CaTA_all = zeros(length(folders)+4, length(protocol_number)) ;

% Simulate for each set of predicted conductances
if strcmp(method, 'single')
    colors = repmat('krgbmc', 1, 300) ;
    extracted = figure ;
    set(gcf, 'Position', [623,287,1047,793]) ;
    for i=1:length(folders)
        folder = folders{i} ; 
        x_conductance = x_conductances(i, :) ;
        [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names); 
        y = zeros(1, 2) ;

        % Plot
        single = figure ;
        set(gcf, 'Position', [50 50 1000 500])
        hold on
        for j=1:length(protocol_number)
            if (sum(t_stim{j}))
                % Predicted for protocol i
                simT = t_stim{j} ;
                simV = V_stim{j} ;
                simCai = Cai_stim{j} ;
                sim_stimtimes = stimtimes{j} ;
                
                % % Features
                try
                    % Frequency
                    [~, ~, features] = APextract_custom(simT,simV,sim_stimtimes) ;
                    freq_all(i,j) = features(4) ;
                    % APD90 and CaTA
                    [keepT, v, ca,tinit,errorcode] = waveform_extract_new(simT,simV,simCai,sim_stimtimes);
                    [outputs,outputlabels] = calculate_features(v,ca,keepT) ;
                    APD90_all(i,j) = outputs(7) ;
                    CaTA_all(i,j) = outputs(9) ;
                catch E
                    freq_all(i,j) = 0 ;
                    APD90_all(i,j) = 0 ;
                    CaTA_all(i,j) = 0 ;
                end
                
                % Experimental for protocol i
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
                % % Predicted
                Vderiv = [simV(2:end);simV(end)] - simV ;
                [~,dexmax] = maxk(Vderiv,5) ;
                simT_init = simT(min(dexmax)) ;
                simT = simT - simT_init ;

                % Plot for protocol i
                figure(single)
                subplot(2, length(protocol_number)+1, j) % AP, protocol i
                hold on
                plot(expT, expV, 'r-', 'LineWidth', 1)
                plot(simT, simV, 'b-', 'LineWidth', 1)
                xlim([0 max(expT)])
                title(['Protocol ', int2str(protocol_number(j))])
                xlabel("Time (ms)")
                ylabel({'AP', 'mV'})

                subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT protocol i
                hold on
                plot(expT, expCai, 'r-', 'LineWidth', 1)
                plot(simT, simCai, 'b-', 'LineWidth', 1)
                xlim([0 max(expT)])
                xlabel('Time (ms)')
                ylabel({'[Ca^{2+}]_i', 'mM'})
                
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
                [simT_extract, simV_extract, simCai_extract,tinit,errorcode] = waveform_extract_new(simT,simV,simCai,sim_stimtimes);
                
                figure(extracted)
                subplot(2, length(protocol_number)+1, j) % AP, protocol j
                hold on
                plot(expT_extract, expV_extract, 'r--', 'LineWidth', 2)
                plot(simT_extract, simV_extract, 'LineWidth', 1)
                title(['Protocol ', int2str(protocol_number(j))])
                xlabel("Time (ms)")
                ylabel({'AP', 'mV'})
                
                subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT protocol i
                hold on
                plot(expT_extract, expCai_extract, 'r--', 'LineWidth', 2)
                plot(simT_extract, simCai_extract, 'LineWidth', 1)
                xlabel('Time (ms)')
                ylabel({'[Ca^{2+}]_i', 'mM'})

                % Fitness calculation
                [expT_V, expV, simT_V, simV] = f_alignWaveformEnds(expT, expV, simT, simV) ;
                [expT_Cai, expCai, simT_Cai, simCai] = f_alignWaveformEnds(expT, expCai, simT, simCai) ;
                V_sim = interp1(simT_V, simV, expT_V);
                CaT_sim = interp1(simT_Cai, simCai, expT_Cai);
                
                % Scale CaT using factor from baseline Kernik APA/CaTA ratio
                CaT_factor = (248317.2776 + 275865.6087 + 258347.0044) / 3 ;
                CaT_sim = CaT_sim*CaT_factor ;
                expCai = expCai*CaT_factor ;
                
                y(1) = y(1) + (sum((V_sim - expV).^2)/length(V_sim));
                y(2)= y(2) + (sum((CaT_sim - expCai).^2)/length(CaT_sim));

            else % If Model Crashes - stop run
                 error(['Simulation error during protocol ',num2str(protocol_number(j))])
            end
        end
        figure(single)
        % Legend: https://stackoverflow.com/questions/44349293/matlab-place-legend-outside-the-plot
        lSub = subplot(2, length(protocol_number)+1, length(protocol_number)+1); 
        hold on
        plot(1, 1, 'r-');
        plot(1, 1, 'b-');
        plot(1, 1, 'g--');
        set(lSub, 'Visible', 'off');
        legend(lSub, 'Experiment', 'Predicted', 'Location', 'best')
        
        % Fitness
        AP_fitness = y(1) ;
        CaT_fitness = y(2) ;
        sum_fitness = sum(y) ;
        tbl = table(sequence, AP_fitness, CaT_fitness, sum_fitness) ;

        figure(single)
        if isNormalized
            savefig([folder, '/validation_norm_', sequence_str]) ;
            writetable(tbl, [folder, '/validation_fitness_norm.xlsx'], 'WriteMode','Append')
        else
            savefig([folder, '/validation_raw_', sequence_str]) ;
            writetable(tbl, [folder, '/validation_fitness_raw.xlsx'], 'WriteMode','Append')
        end
        close(single)
        rownames{i} = ['Run ' int2str(i-1)] ;
        
    end
    figure(extracted)
    % Legend: https://stackoverflow.com/questions/44349293/matlab-place-legend-outside-the-plot
    lSub = subplot(2, length(protocol_number)+1, length(protocol_number)+1); 
    hold on
    plot(1, 1, 'r--', 'LineWidth', 2);
    set(lSub, 'Visible', 'off');
    legend(lSub, 'Experiment', 'Location', 'best')

    if isNormalized
        savefig([base, '/extracted_norm_', sequence_str]) ;
    else
        savefig([base, '/extracted_raw_', sequence_str]) ;
    end
    
    % % Features
    % Mean of all runs
    freq_all(end-3,:) = mean(freq_all(1:length(folders),:), 1) ;
    APD90_all(end-3,:) = mean(APD90_all(1:length(folders),:), 1) ;
    CaTA_all(end-3,:) = mean(CaTA_all(1:length(folders),:), 1) ;
    rownames{end-3} = 'Mean' ;
    % Median
    freq_all(end-2,:) = median(freq_all(1:length(folders),:), 1) ;
    APD90_all(end-2,:) = median(APD90_all(1:length(folders),:), 1) ;
    CaTA_all(end-2,:) = median(CaTA_all(1:length(folders),:), 1) ;
    rownames{end-2} = 'Median' ;
    % SD
    freq_all(end-1,:) = std(freq_all(1:length(folders),:), 0, 1) ;
    APD90_all(end-1,:) = std(APD90_all(1:length(folders),:), 0, 1) ;
    CaTA_all(end-1,:) = std(CaTA_all(1:length(folders),:), 0, 1) ;
    rownames{end-1} = 'SD' ;
    % Range
    freq_all(end,:) = range(freq_all(1:length(folders),:), 1) ;
    APD90_all(end,:) = range(APD90_all(1:length(folders),:), 1) ;
    CaTA_all(end,:) = range(CaTA_all(1:length(folders),:), 1) ;
    rownames{end} = 'Range' ;
    
    tbl_freq = array2table(freq_all, 'VariableNames', colnames, 'RowNames', rownames) ;
    tbl_APD90 = array2table(APD90_all, 'VariableNames', colnames, 'RowNames', rownames) ;
    tbl_CaTA = array2table(CaTA_all, 'VariableNames', colnames, 'RowNames', rownames) ;
    % Save data
%     writetable(tbl_freq, [base,'/features.xlsx'], 'Sheet', 'frequency')
%     writetable(tbl_APD90, [base,'/features.xlsx'], 'Sheet', 'APD90')
%     writetable(tbl_CaTA, [base,'/features.xlsx'], 'Sheet', 'CaTA')
    
else
    if strcmp(method, 'median')
        x_conductance = median(x_conductances, 1) ;
    elseif strcmp(method, 'mean')
        x_conductance = mean(x_conductances, 1) ;
    else
        error("Please enter one of the following methods for using fitted parameters: single, median, or mean.")
    end
    
    % Simulate for predicted conductances
    [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names); 
    % Simulate for Kernik baseline conductances
    [base_t_stim, base_V_stim, base_Cai_stim, base_stimtimes] = ga_simulation_k19(ones(1, length(names)), names); 
    y = zeros(1, 2) ;
    
    figure
    set(gcf, 'Position', [50 50 1000 500])
    hold on
    for i=1:length(protocol_number)
        if (sum(t_stim{i}) && sum(base_t_stim{i}))
            % Predicted for protocol i
            simT = t_stim{i} ;
            simV = V_stim{i} ;
            simCai = Cai_stim{i} ;
            % Experimental for protocol i
            expT = readmatrix('Pseudodataset/saved_data/pseudodataset.xlsx',...
                'Sheet',['Cell ',int2str(cell_number)],'Range','A3:A50003');
            expVCaT = readmatrix('Pseudodataset/saved_data/pseudodataset.xlsx',...
                'Sheet',['Cell ',int2str(cell_number)],'Range',ranges{protocol_number(i)});
            expV = expVCaT(:, 1) ;
            expCai = expVCaT(:, 2) ;
            % Baseline for protocol i
            baseT = base_t_stim{i} ;
            baseV = base_V_stim{i} ;
            baseCai = base_Cai_stim{i} ;

            % Match max timepoints
            if (max(expT) <= max(simT)) && (max(expT) <= max(baseT))
                % Adjust predicted and baseline results to match experimental T
                simV = simV(simT <= max(expT)) ;
                simCai = simCai(simT <= max(expT)) ;
                simT = simT(simT <= max(expT)) ;
                baseV = baseV(baseT <= max(expT)) ;
                baseCai = baseCai(baseT <= max(expT)) ;
                baseT = baseT(baseT <= max(expT)) ;
            elseif (max(simT) <= max(expT)) && (max(simT) <= max(baseT))
                % Adjust experimental and baseline results to match predicted T
                expV = expV(expT <= max(simT)) ;
                expCai = expCai(expT <= max(simT)) ;
                expT = expT(expT <= max(simT)) ;
                baseV = baseV(baseT <= max(simT)) ;
                baseCai = baseCai(baseT <= max(simT)) ;
                baseT = baseT(baseT <= max(simT)) ;
            else % baseline simulation is shortest
                % Adjust experimental and predicted results to match baseline T
                expV = expV(expT <= max(baseT)) ;
                expCai = expCai(expT <= max(baseT)) ;
                expT = expT(expT <= max(baseT)) ;
                simV = simV(simT <= max(baseT)) ;
                simCai = simCai(simT <= max(baseT)) ;
                simT = simT(simT <= max(baseT)) ;
            end

            % Align by first AP upstroke
            % % Experimental
            Vderiv = [expV(2:end);expV(end)] - expV ;
            [~,dexmax] = maxk(Vderiv,5) ;
            expT_init = expT(min(dexmax)) ;
            expT = expT - expT_init ;
            % % Predicted
            Vderiv = [simV(2:end);simV(end)] - simV ;
            [~,dexmax] = maxk(Vderiv,5) ;
            simT_init = simT(min(dexmax)) ;
            simT = simT - simT_init ;
            % % Baseline
            Vderiv = [baseV(2:end);baseV(end)] - baseV ;
            [~,dexmax] = maxk(Vderiv,5) ;
            baseT_init = baseT(min(dexmax)) ;
            baseT = baseT - baseT_init ;

            % Plot for protocol i
            subplot(2, length(protocol_number)+1, i) % AP, protocol i
            hold on
            plot(expT, expV, 'r-', 'LineWidth', 1)
            plot(simT, simV, 'b-', 'LineWidth', 1)
            plot(baseT, baseV, 'g--', 'LineWidth', 1)
            xlim([0 max(expT)])
            title(['Protocol ', int2str(protocol_number(i))])
            xlabel("Time (ms)")
            ylabel({'AP', 'mV'})

            subplot(2, length(protocol_number)+1, i+length(protocol_number)+1) % CaT protocol i
            hold on
            plot(expT, expCai, 'r-', 'LineWidth', 1)
            plot(simT, simCai, 'b-', 'LineWidth', 1)
            plot(baseT, baseCai, 'g--', 'LineWidth', 1)   
            xlim([0 max(expT)])
            xlabel('Time (ms)')
            ylabel({'[Ca^{2+}]_i', 'mM'})

            % Fitness calculation
            [expT_V, expV, simT_V, simV] = f_alignWaveformEnds(expT, expV, simT, simV) ;
            [expT_Cai, expCai, simT_Cai, simCai] = f_alignWaveformEnds(expT, expCai, simT, simCai) ;
            V_sim = interp1(simT_V, simV, expT_V);
            CaT_sim = interp1(simT_Cai, simCai, expT_Cai);
            
            % Scale CaT using factor from baseline Kernik APA/CaTA ratio
            CaT_factor = (248317.2776 + 275865.6087 + 258347.0044) / 3 ;
            CaT_sim = CaT_sim*CaT_factor ;
            expCai = expCai*CaT_factor ;
            
            y(1) = y(1) + (sum((V_sim - expV).^2)/length(V_sim));
            y(2)= y(2) + (sum((CaT_sim - expCai).^2)/length(CaT_sim));

        else % If Model Crashes - stop run
%              close all
             error(['Simulation error during protocol ',num2str(protocol_number(i))])
        end
    end
    % Legend: https://stackoverflow.com/questions/44349293/matlab-place-legend-outside-the-plot
    lSub = subplot(2, length(protocol_number)+1, length(protocol_number)+1); 
    hold on
    plot(1, 1, 'r-');
    plot(1, 1, 'b-');
    plot(1, 1, 'g--');
    set(lSub, 'Visible', 'off');
    legend(lSub, 'Experiment', 'Predicted', 'Baseline model', 'Location', 'best')

    % Fitness
    AP_fitness = y(1) ;
    CaT_fitness = y(2) ;
    sum_fitness = sum(y) ;
    tbl = table(sequence, AP_fitness, CaT_fitness, sum_fitness) ;

    savefig([base, '/validation_', sequence_str]) ;
    writetable(tbl, [base, '/validation_fitness.xlsx'], 'WriteMode','Append')

end
