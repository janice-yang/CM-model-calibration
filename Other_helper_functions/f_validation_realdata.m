function [] = f_validation_realdata(folder, isNormalized, exp_file, sheetnames, datatype, stimV, protocol_number)
% Validate fitted parameters on unseen experiments
%{
Inputs:
    folders: cell array of directories containing individual run results
    isNormalized: boolean indicating whether data is normalized 
    exp_file: validation protocol experimental data file name (Excel file)
    sheetnames: cell array indicating which sheets in exp_file to read
    datatype: string indicating whether experimental data includes AP, CaT, or both
    stimV: array of stimulation times if paced (sorry about the name)
    protocol_number: only single number for now
    method ('single', 'median', or 'mean'): which parameter(s) to use for
        validation (single = from individual runs, median/mean = summary from
        all runs)

Protocols for validation:
Ko low --> normal --> high, PCL 800ms: 15,18,13
Ko low --> normal --> high, spontaneous: 16,19,14
IKr block 0% --> 15% --> 30%, PCL 800ms: 18,9,11
IKr block 0% --> 15% --> 30%, spontaneous: 19,10,12
Experimental 1nM dofetilide (40% IKr block): 22 (spont.)
Experimental 1nM dofetilide (40% IKr block): 23 (1 Hz pacing)
%}

if ~(strcmp(datatype, 'AP') || strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT'))
    error("Please enter a valid datatype - 'AP', 'CaT', or 'APCaT'")
end
    
if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT') % always normalize if CaT data is included
    isNormalized = true ;
end 
scaleCaT = false ;
sequence_str = replace(num2str(protocol_number), '  ', '-') ;
sequence = {sequence_str} ;

save GA/curr_cell_protocol protocol_number isNormalized scaleCaT
load 'GA/x_names.mat' names

logfactor = 2 ; 

% Simulate for each set of predicted conductances
colors = repmat('krgbmc', 1, 300) ;
load([folder, '/minoptimparams.mat'], 'minoptimparams')
x_conductance = logfactor.^minoptimparams ;
[t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names); 
y = zeros(1, 2) ;

figure
set(gcf, 'Position', [50 50 1000 500])
hold on
for j=1:length(protocol_number)
    if sum(t_stim{j})
        % Predicted for protocol i
        simT = t_stim{j} ;
        simV = V_stim{j} ;
        simCai = Cai_stim{j} ;

        % Get last 5 seconds
        idx_start = find(simT > simT(end) - 5000, 1) ;
        simT = simT(idx_start:end) ;
        simT = simT - simT(1) ; 
        simV = simV(idx_start:end) ;
        simCai = simCai(idx_start:end) ;
        
        % Experimental for protocol i
        expdata = readmatrix(exp_file, 'Sheet', sheetnames{j}) ;
        expT = expdata(:,1) ;
        stim_idx = ismembertol(expdata(:,1), stimV, 2.5e-7) ;
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
            
        % Match max timepoints
        if max(expT) <= max(simT)
            % Adjust predicted results to match experimental T
            simT = simT(simT <= max(expT)) ;                
            if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
                simV = simV(simT <= max(expT)) ;
            end
            if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
                simCai = simCai(simT <= max(expT)) ;
            end
        else % simulation is shorter than experiment
            expT = expT(expT <= max(simT)) ;
            exp_stimtimes = exp_stimtimes(exp_stimtimes <= max(simT)) ;
            if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
                expV = expV(expT <= max(simT)) ;
            end
            if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
                expCai = expCai(expT <= max(simT)) ;
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

        % Plot for protocol i
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            subplot(2, length(protocol_number)+1, j) % AP, protocol i
            hold on
            plot(expT, expV, 'r-', 'LineWidth', 2)
            plot(simT, simV, 'b-', 'LineWidth', 2)
            xlim([0 max(expT)])
            title(['Protocol ', int2str(protocol_number(j))])
            xlabel("Time (ms)")
            ylabel('mV')
        end
        
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            subplot(2, length(protocol_number)+1, j+length(protocol_number)+1) % CaT protocol i
            hold on
            plot(expT, expCai, 'r-', 'LineWidth', 2)
            plot(simT, simCai, 'b-', 'LineWidth', 2)
            xlim([0 max(expT)])
            xlabel('Time (ms)')
            ylabel('mM')
        end

        % Fitness calculation
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            [expT_V, expV, simT_V, simV] = f_alignWaveformEnds(expT, expV, simT, simV) ;
            V_sim = interp1(simT_V, simV, expT_V);
            y(1) = y(1) + (sum((V_sim - expV).^2)/length(V_sim));
        end
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            [expT_Cai, expCai, simT_Cai, simCai] = f_alignWaveformEnds(expT, expCai, simT, simCai) ;
            CaT_sim = interp1(simT_Cai, simCai, expT_Cai);
            y(2)= y(2) + (sum((CaT_sim - expCai).^2)/length(CaT_sim));
        end

    else % If Model Crashes - stop run
         close all
         error(['Simulation error during protocol ',num2str(protocol_number(j))])
    end
end

% Fitness
AP_fitness = y(1) ;
CaT_fitness = y(2) ;
sum_fitness = sum(y) ;
tbl = table(sequence, AP_fitness, CaT_fitness, sum_fitness) ;

writetable(tbl, [folder, '/validation_fitness.xlsx'], 'WriteMode','Append')

end
