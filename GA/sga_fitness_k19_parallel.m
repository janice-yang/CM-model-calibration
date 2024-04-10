% Fitness Function FOR PARALLEL COMPUTING

function [fitness] = sga_fitness_k19_parallel(x, runNum) % Takes inputs from main script initial population and returns error values as output
%% Standard Opening

format compact


%% User Input for Variable Parameters

% Pacing
% pace = 1.0;         %Frequency (Hz) 
load 'GA/curr_cell_protocol.mat' cell_number protocol_number scaleCaT


%% SET MODEL PARAMETERS
  
  %==========Maximum Conductances & Currents==========%
  % As always, adjusted for capacitance

logfactor = 2 ; 
x_conductance = (logfactor.^x);

%% ODE SOLVER
%%%%%%%%%%%%%%%%%%%%%%%

load 'GA/x_names.mat' names
load experimental_dataset.mat experimental_dataset

[t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, names); % returns time, V and Cai vectors from simulation with Kernik19 model. If spontaneous, stimtimes = []
y = zeros(1, 2) ;
for i=1:length(protocol_number)
    if sum(t_stim{i})
        [keepT, V, CaT,tinit,errorcode] = waveform_extract_new(t_stim{i}, V_stim{i},Cai_stim{i},stimtimes{i});
        
        % Stop if no AP was extracted, write culprit conductances to file
        if (~isvector(keepT) || ~isvector(V) || ~isvector(CaT))
            T = array2table(x_conductance, "VariableNames", names) ;
            writetable(T, ['GA/noAP_errors/', 'protocol_',int2str(protocol_number(i)),'.xlsx'], 'WriteMode', "append", 'WriteRowNames', true);
            break
        end
        
        % Align ends of exp and simulated extracted waveforms
        exp_T_V = experimental_dataset{i}.Time_AP{1} ;
        exp_V = experimental_dataset{i}.AP{1};
        exp_T_CaT = experimental_dataset{i}.Time_CaT{1} ;
        exp_CaT = experimental_dataset{i}.CaT{1};
        
        [exp_T_V, exp_V, T_V, V] = f_alignWaveformEnds(exp_T_V, exp_V, keepT, V) ;
        [exp_T_CaT, exp_CaT, T_CaT, CaT] = f_alignWaveformEnds(exp_T_CaT, exp_CaT, keepT, CaT) ;
        
        %% Fitness Calculation
        
        V_sim = interp1(T_V, V, exp_T_V);
        CaT_sim = interp1(T_CaT, CaT, exp_T_CaT);
        
        if scaleCaT
            % Scale CaT using ratio from experimental APA/CaTA
            CaT_factor = range(exp_V) / range(exp_CaT) ;
            % Adjust CaTs using scale factor
            CaT_sim = CaT_sim*CaT_factor ;
            exp_CaT = exp_CaT*CaT_factor ;
            
            % Record new range(V_sim)/range(CaT_sim)
            AP_CaT_ratio = range(V_sim) / range(CaT_sim) ;
            protocol = protocol_number(i) ;
            tbl = table(AP_CaT_ratio, protocol, x_conductance) ; % change to actual conductances (scale baseline by x_conductance)
            writetable(tbl, ['GA/Results/Run_' int2str(runNum), '/CaT_scale_factors.txt'], 'WriteMode','Append')
        end

        %Normalize CaT_sim:
        %CaT_sim=(CaT_sim-min(CaT_sim))/(max(CaT_sim)-min(CaT_sim));

        y(1) = y(1) + sum((V_sim-exp_V).^2);
        y(2)= y(2) + sum((CaT_sim-exp_CaT).^2);


    else % If Model Crashes - if not on first protocol, add 1e10?
         y(1) = 1e10;
         y(2) = 1e10;
    end
end

if ~isreal(y(1))
     y(1) = 1e10;
     y(2) = 1e10;
end

if isnan(y)
     y(1) = 1e10;
     y(2) = 1e10;
end

fitness  = sum(y);
return
