function [state,options,optchanged] = ga_output_k19(options,state,flag,runNum) 
persistent iteration_number

if isempty(iteration_number)
  iteration_number = 0;
end
iteration_number = iteration_number + 1 ;

currentpop = state.Population ;
currentscore = state.Score ;

[~,mindex] = min(currentscore);
bestparams = currentpop(mindex,:);

x_conductance = (2.^bestparams);

optchanged = false ; 

% load runNum runNum 
% load initials initials
save(['GA/Results/', 'Run_', int2str(runNum), 'k19output_',int2str(iteration_number)] ...
  ,'currentpop','currentscore','state','bestparams') ;

% %% Plot results: remove when GA is working
 load x_names names
 load experimental_dataset experimental_dataset
 load 'GA/curr_cell_protocol.mat' protocol_number datatype
 % Plot simulation vs experiment for each set of optimal parameters
 for i = 1:size(bestparams,1) 
     [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance(i,:), names);
     p = 1 ; % for subplot positioning
     figure
     hold on
     handle = gcf ;
     for j=1:length(protocol_number)
        if sum(t_stim{j})
            [keepT, V, CaT,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},stimtimes{j});

            % Stop if no AP was extracted, write culprit conductances to file
            if (~isvector(keepT) || ~isvector(V) || ~isvector(CaT))
            T = array2table(x_conductance(i,:), "VariableNames", names) ;
            writetable(T, ['GA/noAP_errors/', 'protocol_',int2str(protocol_number(j)),'.xlsx'], 'WriteMode', "append", 'WriteRowNames', true);
            continue
            else
            % Align ends of exp and simulated extracted waveforms
            if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            exp_T_V = experimental_dataset{j}.Time_AP{1} ;
            exp_V = experimental_dataset{j}.AP{1};
            [exp_T_V, exp_V, T_V, V] = f_alignWaveformEnds(exp_T_V, exp_V, keepT, V) ;

            V_sim = interp1(T_V, V, exp_T_V);
            subplot(length(protocol_number),2,p)
            hold on
            plot(exp_T_V,V_sim)
            plot(exp_T_V,exp_V)
            ylabel('mV')
            title(['Protocol ', num2str(protocol_number(j)), ' AP'])  
            end
            
            if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
                exp_T_CaT = experimental_dataset{j}.Time_CaT{1} ;
                exp_CaT = experimental_dataset{j}.CaT{1};
                [exp_T_CaT, exp_CaT, T_CaT, CaT] = f_alignWaveformEnds(exp_T_CaT, exp_CaT, keepT, CaT) ;
                CaT_sim = interp1(T_CaT, CaT, exp_T_CaT);
                subplot(length(protocol_number),2,p+1)
                hold on
                plot(exp_T_CaT,CaT_sim)
                plot(exp_T_CaT,exp_CaT)
                ylabel('CaT')
                title(['Protocol ', num2str(protocol_number(j)), ' CaT'])
            end

            legend({'prediction', 'observation'}, 'Location', 'best')
            p = p + 2 ; 
            end
        end
     end
     print('-djpeg',['GA/Results/', 'Run_',int2str(runNum),'/k19plot_',int2str(iteration_number),'_',int2str(i)]) ;
     close(handle)
 end
return

