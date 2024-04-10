%{
ToRORd19 is a concrete class for a cardiomyocyte based on the model published by
Tomek et al. in 2019. It inherits from the superclass AbsCM which is an abstract
class representing cardiomyocytes based on any model.
%}

classdef ToRORd19 < AbsCM
    properties
        % defined in @f_ToRORd19_loadBaselineProperties
        universals; celltype; name ; ODEModel ; isCellAlive ; cellID ; geometry ; state ; YNames ; YUnits ;
        currentNames ; currentUnits ; parameterNames ; parameterUnits ;
        parameters ; conductanceNames ; conductances ; protocol ;
        x_conductances ; y_allevents ; x_allevents ; event_number ; isPaced ;
    end
    
    methods
        function tord19cm = ToRORd19(init_Y, celltype, cellID)
            [universals, name, model, isCellAlive, geometry, state, YNames, YUnits,...
                currentNames, currentUnits, parameterNames, parameterUnits,...
                parameters, conductanceNames, conductances, protocol, isPaced]...
                = f_ToRORd19_loadBaselineProperties(init_Y) ;
            tord19cm.universals = universals ;
            tord19cm.celltype = celltype ;
            tord19cm.name = name ;
            tord19cm.ODEModel = model ;
            tord19cm.isCellAlive = isCellAlive ;
            tord19cm.cellID = cellID ;
            tord19cm.geometry = geometry ;
            tord19cm.state = state ;
            tord19cm.YNames = YNames ;
            tord19cm.YUnits = YUnits ;
            tord19cm.currentNames = currentNames ;
            tord19cm.currentUnits = currentUnits ;
            tord19cm.parameterNames = parameterNames ;
            tord19cm.parameterUnits = parameterUnits ;
            tord19cm.parameters = parameters ;
            tord19cm.conductanceNames = conductanceNames ;
            tord19cm.conductances = conductances ;
            tord19cm.protocol = protocol ;
            tord19cm.isPaced = isPaced ;
        end
        
        function saveX_conductance(tord19cm, x_conductance)
           tord19cm.x_conductances = x_conductance ; 
        end
        
        function scaleParameters(tord19cm, scaleFactors, names)
            [~, indA, indB] = intersect(names, tord19cm.parameterNames, 'stable') ;
            tord19cm.parameters.scaling(indB) = tord19cm.parameters.scaling(indB) .* scaleFactors(indA);
        end
        
        function scaleConductances(tord19cm, scaleFactors, names)
            [~, indA, indB] = intersect(names, tord19cm.conductanceNames, 'stable'); 
            tord19cm.conductances.scaling(indB) = tord19cm.conductances.scaling(indB) .* scaleFactors(indA);
        end
        
        function setEnvironment(tord19cm, T, nao, cao, ko)
            tord19cm.universals.T = T ;
            tord19cm.universals.nao = nao ;
            tord19cm.universals.cao = cao ;
            tord19cm.universals.ko = ko ;
        end
        
        
        function setUpPacingProtocol(tord19cm, amplitudes, numPulses, precedingTime, pulseDurations)
            [tord19cm.protocol.intervalTimes, tord19cm.protocol.stimulus]...
                = f_setUpPacingProtocol(amplitudes, numPulses, precedingTime, pulseDurations) ;
            
            tord19cm.protocol.amplitudes = amplitudes ;
            tord19cm.protocol.numPulses = numPulses ;
            tord19cm.protocol.precedingTime = precedingTime ;
            tord19cm.protocol.pulseDuration = pulseDurations ;
            
            tord19cm.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations) ; % for reference
            tord19cm.protocol.phaseLengths(1) = 0.2 + (numPulses(1) - 1) * precedingTime(1) + numPulses(1) * pulseDurations(1) ; % first pulse doesn't have full preceding time
            tord19cm.protocol.phaseLengths(end) = tord19cm.protocol.phaseLengths(end) + 1 ; % add extra second at the end to capture last AP
            tord19cm.protocol.frequencies = 1000./(precedingTime + pulseDurations) ; % for reference
            
            tord19cm.protocol.totalTime = tord19cm.protocol.intervalTimes(end,2) ;
            tord19cm.protocol.nPhases = length(amplitudes) ;
            tord19cm.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations)/1000 ;
            tord19cm.protocol.intervalTimes = tord19cm.protocol.intervalTimes * 1000; % change to ms
            
        end
        
        % FUNCTION INFO: setUpDrugApplication
        %{
        This function takes in a ToRORd19 object as well as information about the
        drug protocol and sets the 'drugEffects' field of conductances to a
        vector of scaling factors that are the effects the drugs have on the
        different channel conductances in the order they are listed in the
        'conductanceNames' property.

        Drugs that are not applied should be denoted as 1, antagonists
        should be a number between 0 - 1 signifying the fraction of initial
        conductance that will remain after drug application on that target, and
        agonists should be a number greater than 1. The 'conductances' field
        'applicationTimes' is set to contain the onTimes and offTimes
        corresponding to the time each drug is added and "magic wand" removed.
        Any times can be used for drugs that are not added (effect = 1).
        
        The methodology used to ensure maximum precision in the pacing protocol
        was not repeated for drug application as drugs are applied for a much
        longer period of time than a current stimulus so precision to the
        fraction of a second is not essential.
        %}
        function setUpDrugApplication(tord19cm, drugEffects, onTimes, offTimes)
            tord19cm.conductances.applicationTimes = [onTimes; offTimes] ;
            tord19cm.conductances.drugEffects = drugEffects ;
        end
        
        function odeSolver(tord19cm)
            warning off
            warning('tord19cm odesolver starting') % create a baseline warning
            
            if isempty(tord19cm.protocol.intervalTimes)
                error('Error: Please set up pacing protocol first')
            else
                options = odeset('MaxStep', 1, 'InitialStep', 2e-2, 'Event', @f_tinyStepCatcher_t19) ; % because t in ms? original: odeset('MaxStep', 1e-3, 'InitialStep', 2e-5) ; % solver settings
                
                % Allocate space initially to save time. Extra zeros will be
                % removed later to save space.
                initialrows = ceil(tord19cm.protocol.totalTime*2000) ;
%                 initialrows = 25000 ; % enough for 9 beats
                
                tord19cm.state.Y = zeros(initialrows, length(tord19cm.state.init_Y)) ;
                tord19cm.state.t = zeros(initialrows, 1) ; % first row at T0, no need to update
                
                tord19cm.state.Y(1,:) = tord19cm.state.init_Y ; % first row, input values (generally steady state)
                input_values = tord19cm.state.init_Y ; % initialize, changes in for-loop
                index = 2 ; % initialize, changes in for-loop
                
                stimIndices = find(tord19cm.protocol.stimulus > 0);
%                 i_stimToKeep = stimIndices(end-9) ; % last 10 beats
                i_stimToKeep = 0;%stimIndices(end - end + 1); % whole simulation
                
                % convert properties stored as structs to cells to pass to
                % ode15s
                universals_arr = cell2mat(struct2cell(tord19cm.universals));
                geometry_arr = cell2mat(struct2cell(tord19cm.geometry));
                parameters_arr = tord19cm.parameters.baseline.*tord19cm.parameters.scaling ;
                conductances_arr = tord19cm.conductances.baseline.*tord19cm.conductances.scaling ;
                drugEffects_arr = tord19cm.conductances.drugEffects ; 
                applicationTimes_arr = tord19cm.conductances.applicationTimes;
                               
                timer = 0;
                for i = 1:length(tord19cm.protocol.stimulus)
                    tic;
                    if tord19cm.isCellAlive == 1
                        input_values(end) = tord19cm.protocol.stimulus(i) ; % last input value is the stimulus amplitude
                        [t_current, Y_current, t_event, y_event, ~] = ode15s(...
                            tord19cm.ODEModel,...
                            tord19cm.protocol.intervalTimes(i,:),...
                            input_values,...
                            options,...
                            tord19cm.celltype,...
                            universals_arr,...
                            geometry_arr,...
                            parameters_arr,...
                            conductances_arr,...
                            drugEffects_arr,...
                            applicationTimes_arr);
                        
                        current_length = length(t_current)-1 ; % length of this interval in the stimulus protocol
                        
                        %% check for weird short intervals
                        if current_length == 0 % another way the model can break with bad parameters ---> try bundle this into the existing check
                            %disp('Die, bad cell, die!!! CL0')
                            tord19cm.isCellAlive = 0;
                            cell_state = tord19cm.state;
                            cell_parameters = tord19cm.x_conductances;
                            save(['GA/Results/events', num2str(tord19cm.cellID), '_CL0.mat'], 'cell_state', 'cell_parameters');
                            tord19cm.state.t = zeros(10,1); % resultsless
                            tord19cm.state.Y = zeros(10,length(tord19cm.YNames)); % resultsless
                            break
                        end
                        
                        %% check if there was an event
                        if (t_event) % if there was an event
                            numEvents = length(t_event);
                            if isempty(tord19cm.event_number)
                                tord19cm.event_number = 1 ;
                            end
                            tord19cm.y_allevents(:,tord19cm.event_number:tord19cm.event_number + numEvents-1) = y_event' ;
                            tord19cm.x_allevents(:,tord19cm.event_number) = tord19cm.x_conductances' ;
                            y_allEvents = tord19cm.y_allevents;
                            x_allEvents = tord19cm.x_allevents;
                            save(['GA/Results/events', num2str(tord19cm.cellID), '_', num2str(tord19cm.event_number), '.mat'], 'y_allEvents', 'x_allEvents');
                            %disp('Die, bad cell, die!!! E')
                            tord19cm.isCellAlive = 0;
                            tord19cm.state.t = zeros(10,1); % resultsless
                            tord19cm.state.Y = zeros(10,length(tord19cm.YNames)); % resultsless
                            tord19cm.event_number = tord19cm.event_number + numEvents ;
                            break
                        end
                        
                    end
                    
                    % cells may have been killed off above so we check again for
                    % survival
                    if tord19cm.isCellAlive == 1
                        if i >= i_stimToKeep
                            tord19cm.state.t(index:index+current_length-1) = t_current(2:end) ; % save outcome of this run
                            tord19cm.state.Y(index:index+current_length-1,:) = Y_current(2:end,:) ; % save outcome of this run
                            index = index + current_length; % update the index only if storing all or if last 3 beats
                            input_values = tord19cm.state.Y(index - 1, :); % set input values for next interval to end of current interval
                        else
                            input_values = Y_current(end, :) ; % set input values for next interval to end of current interval
                        end
                    end
                    timer = timer + toc ; 
                    
                    if timer > 7000 %40
                        cell_state = tord19cm.state;
                        cell_parameters = tord19cm.x_conductances;
                        save(['GA/Results/events', num2str(tord19cm.cellID), '_timer.mat'], 'cell_state', 'cell_parameters', 'timer');
                        %disp('Die, bad cell, die!!! T')
                        tord19cm.isCellAlive = 0;
                        tord19cm.state.t = zeros(10,1); % resultsless
                        tord19cm.state.Y = zeros(10,length(tord19cm.YNames)); % resultsless
                        break
                    end   
                end
                %disp(timer)
                
                if ~isempty(tord19cm.state.t) && tord19cm.isCellAlive == 1
                % remove extra zeros to save space
                zeroInds = find(~tord19cm.state.t) ; % returns indices of all zeroes
                firstExtraZero = zeroInds(2) ; % first zero is T0
                tord19cm.state.t = tord19cm.state.t(1:firstExtraZero-1) ;
                tord19cm.state.Y = tord19cm.state.Y(1:firstExtraZero-1,:) ;
                
                tord19cm.state.t = tord19cm.state.t(2:end);
                tord19cm.state.Y = tord19cm.state.Y(2:end, :);

                if sum(tord19cm.protocol.amplitudes) % if not spontaneous
                tord19cm.protocol.last3stimTimes = tord19cm.protocol.intervalTimes([end-5, end-3, end-1], 1) - tord19cm.state.t(1);
                else 
                tord19cm.protocol.last3stimTimes = [];
                end
                tord19cm.state.t = tord19cm.state.t - tord19cm.state.t(1); % set t0 to 0
%                 tord19cm.protocol.last3stimTimes = tord19cm.protocol.intervalTimes(stimIndices(end-2:end),1) - tord19cm.state.t(1);
%                 tord19cm.state.t = tord19cm.state.t - tord19cm.state.t(1); % set t0 to 0
                end
            end
            %disp('tord19cm odesolver done!')
        end
        
        function getCurrents(tord19cm)
            %disp('tord19cm getcurrents starting')
            if isempty(tord19cm.state.t)
                error('Error: Please run the odeSolver on this ToRORd19 CM first')
            else
                tord19cm.state.currents = zeros(length(tord19cm.state.t), length(tord19cm.currentNames));
                numLoops = length(tord19cm.state.t) ;
                
                % convert properties stored as structs to cells to pass to
                % ode15s
                universals_arr = cell2mat(struct2cell(tord19cm.universals));
                geometry_arr = cell2mat(struct2cell(tord19cm.geometry));
                parameters_arr = tord19cm.parameters.baseline.*tord19cm.parameters.scaling ;
                conductances_arr = tord19cm.conductances.baseline.*tord19cm.conductances.scaling ;
                drugEffects_arr = cell2mat(struct2cell(tord19cm.conductances.drugEffects));
                applicationTimes_arr = cell2mat(struct2cell(tord19cm.conductances.applicationTimes));
                
                for i= 1:numLoops
                    [~, ij_data]    = f_ToRORd19(tord19cm.state.t(i), tord19cm.state.Y(i,:),...
                        tord19cm.celltype,...
                        universals_arr,...
                        geometry_arr,...
                        parameters_arr,...
                        conductances_arr,...
                        drugEffects_arr,...
                        applicationTimes_arr); % ~ says only return data, not dY output
                    tord19cm.state.currents(i,1:end-1) = ij_data;
                end
                %tord19cm.state.currents(:,end) = sum(tord19cm.state.currents(:,1:end),2); *** what would we want to include here? sum for Itot; the 2 at the end is so we sum the rows not the columns
                tord19cm.state.currents(:,end) = tord19cm.state.currents(:,8) + tord19cm.state.currents(:,9) ; % INaCa = INaCa_i + INaCa_ss
                
                % plot all the currents!
                %{
                figure
                lastStim = find(tord19cm.state.t == tord19cm.protocol.intervalTimes(end-1,1));
                for i = 1:length(tord19cm.currentNames)
                    subplot(5,8,i)
                    plot(tord19cm.state.t(lastStim:end)-tord19cm.state.t(lastStim),tord19cm.state.currents(lastStim:end,i))
                    title(tord19cm.currentNames(i))
                end
                %}
            end
            %disp('tord19cm getcurrents done!')
        end
    end
    
end