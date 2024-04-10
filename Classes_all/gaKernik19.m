%{
Kernik19 is a concrete class for a cardiomyocyte based on the model published by
Kernik et al. in 2019. It inherits from the superclass AbsCM which is an abstract
class representing cardiomyocytes based on any model.
%}
classdef gaKernik19 < AbsCM
    properties
        name = 'gaKernik19 CM' ;
        ODEModel = @f_gaKernik19 ;
        
        geometry = struct('V_tot', 3960, ... µm^3
            'Vc_tenT', 16404, ... µm^3
            'VSR_tenT', 1094, ... µm^3
            'V_tot_tenT', [], ... µm^3
            'V_SR', [], ... µm^3
            'Vc', [], ... µm^3
            'Cm', 60) ; % pF, we keep capacitance in geometry for convenience
        
        % For informational queries & keeping order consistent across versions
        YNames = {'Vm', 'Ca_SR', 'Cai', 'Nai', 'Ki', 'Ca_ligand', 'd', 'f1', 'fCa', 'Xr1', 'Xr2', 'Xs', 'h', 'j', 'm', 'Xf', 's', 'r', 'dCaT', 'fCaT', 'R', 'O', 'I', 'stim'} ;
        YUnits = {'mV', 'mM',    'mM',  'mM',  'mM', 'mM',        '-', '-',  '-',   '-',   '-',   '-',  '-', '-', '-', '-',  '-', '-', '-',    '-',    '-', '-', '-', 'A?'} ;
        
        currentNames = {'i_K1', 'i_to', 'i_Kr', 'i_Ks', 'i_CaL', 'i_NaK', 'i_Na', 'i_NaCa', 'i_PCa', 'i_f', 'i_b_Na', 'i_b_Ca', 'i_rel', 'i_up', 'i_leak', 'i_stim', 'i_CaT'};
        %         currentUnits = {'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'mV', 'mV', 'A', 'A'} ;
        
        state = struct('init_Y', [], 'Y', [], 't', [], 'currents', []) ;
        
        % For informational queries & keeping order consistent across versions
        paramNames = {'xK11', 'xK12', 'xK13', 'xK14', 'xK15',... i_K1
            'Xr1_1', 'Xr1_2', 'Xr1_5', 'Xr1_6', 'Xr2_1', 'Xr2_2', 'Xr2_5', 'Xr2_6', 'tauXr1_const', 'tauXr2_const',... i_Kr
            'ks1', 'ks2', 'ks5', 'ks6', 'tauks_const',... i_Ks
            'r1', 'r2', 'r5', 'r6', 's1', 's2', 's5', 's6', 'tau_r_const', 'tau_s_const',... i_to
            'd1', 'd2', 'd5', 'd6', 'f1', 'f2', 'f5', 'f6', 'taud_const', 'tauf_const',... i_CaL
            'm1', 'm2', 'm5', 'm6', 'h1', 'h2', 'h5', 'h6', 'j1', 'j2', 'tau_m_const', 'tau_h_const', 'tau_j_const',... i_Na
            'xF1', 'xF2', 'xF5', 'xF6', 'xF_const',... i_f
            'taud_scale', 'tauf_scale',... tau scalings
            } ;
        %?paramUnits = {} ;
        parameters = struct( 'baseline', [0.477994972217041, 27.2427558793487, 4.92502331781412, 8.72223760006882, 56.6361974998244,... i_K1
            0.00574885237435, 13.6234926362576, 0.047630571181836, -7.06808742965549, 0.012456640526827, -25.9944581644377, 37.3426331501041, 22.0919642353902, 50, 0,... i_Kr
            0.00116558448, 66726.8386758936, 0.28045890825, -18.86697157291, 4.74115e-06,... i_Ks
            0.0553614181713, 11.6842023429669, 3.9891810803775, -11.0471393012032, 0.0003442309443, -17.6344722898096, 186.760536909695, 8.1809338733227, 0.6967584211715, 11.2244577239469,... i_to
            12.966294189722, 7.079145964711, 0.044909415507, -6.909880369242, 0.000512589826, -49.50571203387, 1931.21122351432, 5.730027499699, 1.65824694683, 100.462559171103,... i_CaL
            108.045846384818, 13.107015733941, 0.002326914367, -7.917726289513, 0.003626598864, -19.839358860026, 9663.29497711474, 7.395503564613, 0.000512257182, -66.583755502652, 0.031977580384, 0.167331502516, 0.951088724962,... i_Na
            5.7897e-7, -14.5897121702, 20086.6502378844, 10.20235284528, 23.94529134653,... i_f
            1, 1,... tau scalings
            ],...
            'scaling', ones(1,60)) ;
        
        % TO DO: create properties for three-tiered parameters - commonly
        % varied, locally varied, constant
        
        conductanceNames = {        'g_Na',         'g_f',   'p_CaL',          'g_to',          'g_Ks',  'g_Kr',   'g_K1',            'g_PCa', 'g_b_Na', 'g_b_Ca',   'VmaxUp',  'g_irel_max', 'kNaCa', 'PNaK',   'V_leak', 'g_CaT'} ;
        conductanceUnits = {        'nS/pF',        'nS/pF', 'nS/pF',          'nS/pF',         'nS/pF', 'nS/pF',  'g_K1unit',        'A/F',   'nS/pF',  'nS/pF',    'mM/s',    '/ms',        'A/F',   'A/F',    '/ms',    'nS/pF'} ;
        conductances = struct( 'baseline', [9.720613409241, 0.0435,  0.308027691379,   0.1178333333333, 0.0077,  0.218025, 0.133785777797606, 0.2625,  4.35e-04, 3.6704e-04, 1.105e-04, 12.5,         1100,    2.4761,   1.6e-06,  0.185],...
            'scaling', ones(1,16),...
            'drugEffects', ones(1,16),...
            'applicationTimes', zeros(2,16)) ;
        
        protocol = struct('intervalTimes', [] , 'stimulus', [], 'amplitudes', [],...
            'numPulses', [], 'precedingTime', [], 'pulseDuration', [], 'totalTime', [],...
            'phaseLengths', [], 'frequencies', []) ;
        
        %for ga:
        universals = struct('F', 96485.3415, 'R', 8.314472, 'T', 310.0,...
            'Nao', 151.0, 'Ko', 5.4, 'Cao', 1.8) ;
        isCellAlive = 1;
        cellID;
        x_conductances = [];
        y_allevents = [];
        x_allevents = [];
        event_number = [];
    end
    
    methods
        
        % CONSTRUCTOR:
        %{
        creates an instance of Paci18 with initial state variables equal to
        init_Y.
        %}
        function k19CM = gaKernik19(init_Y, cellID)
            % some options for init_Y:
            %{
            [-75.5966016388547,0.335086796732326,0.000219191642424964,7.16928091250999,104.748824394112,0,0.000394925342652924,0.170990105585540,0.877798946134089,0.309767485715433,0.450577185148519,0.153788281650949,0.739543607812429,0.124515982574505,0.0297549962926414,0.00640338504912616,0.746802810614006,0.000267597833344161,0.000270195573471577,0.756032904368393,0.0113120363433751,0.000165045105312396,0.0142153622323012,0]
            %}
            k19CM.state.init_Y = init_Y ;
            k19CM.cellID = cellID ;
            
            % set geometry from ten Tusscher values
            k19CM.geometry.V_tot_tenT = k19CM.geometry.Vc_tenT + k19CM.geometry.VSR_tenT ; % V_total data from Hwang et al., V_c and V_SR  proportionally scaled from Ten Tusscher 2004 values
            k19CM.geometry.Vc = k19CM.geometry.V_tot * (k19CM.geometry.Vc_tenT/k19CM.geometry.V_tot_tenT) ; % = 3712.4 µm^3 (93.7% total volume)
            k19CM.geometry.V_SR = k19CM.geometry.V_tot * (k19CM.geometry.VSR_tenT/k19CM.geometry.V_tot_tenT) ;% = 247.6 µm^3 (6.3% total volume)
        end
        
        function saveX_conductance(k19CM, x_conductance)
            k19CM.x_conductances = x_conductance ;
        end
        
        % FUNCTION INFO: scaleCell
        %{
        This function takes two sets of scale factors to modify the parameters
        and conductances and saves it in the relevant struct. The default value
        for the scaling fields is a vector of 1s, so this function need only be
        used when alternative scaling is required.
        %}
        function scaleCell(k19CM, conductanceScaleFactors, parameterScaleFactors)
            k19CM.conductances.scaling = conductanceScaleFactors ;
            k19CM.parameters.scaling = parameterScaleFactors ;
        end
        
        function scaleConductances(k19CM, conductanceScaleFactors, names)
            [~, indA, indB] = intersect(names, k19CM.conductanceNames, 'stable') ;
            k19CM.conductances.scaling(indB) = k19CM.conductances.scaling(indB).* conductanceScaleFactors(indA) ;
        end
        
        function scaleParameters(k19CM, parameterScaleFactors, names)
            [~, indA, indB] = intersect(names, k19CM.paramNames, 'stable') ;
            k19CM.parameters.scaling(indB) = k19CM.parameters.scaling(indB).* parameterScaleFactors(indA) ;
        end
        
        function setEnvironment(k19CM, T, nao, cao, ko)
            k19CM.universals.T = T ;
            k19CM.universals.Nao = nao ;
            k19CM.universals.Cao = cao ;
            k19CM.universals.Ko = ko ;
        end
        
        % FUNCTION INFO: setUpPacingProtocol
        %{
        This function takes in a Paci18 object as well as information about the
        pacing protocol and sets the protocol property's intervalTimes field to
        a two-column matrix containing the start and end times of each interval
        in the protocol, where an interval is defined as either the time
        preceding the stimulus, or the time during the stimulus.
        The stimulus field is set to the amplitude of the stimulus
        applied to the k19CM in each of these intervals, zero before and some
        other value during. This methodology is used because of the inconsistent
        time steps used by MATLAB's ODE solvers. While the inconsistent time
        steps speed up the process overall, it may overshoot the exact moment
        each stimulus begins, which may be consequential given the short
        stimulus durations.
        
        This function employs an external function called f_setUpPacingProtocol
        that should be included in the same folder as the Paci18 class.

        It is important when using this function that the amplitudes,
        numPulses, precedingTime and pulseDuration vectors are all the same
        length, otherwise the function will throw an error.
 
        Amplitudes should be denoted as a number in (A/F) that will be
        multiplied by the membrane capacitance of the cell (F) to produce a
        stimulus (A). For phases where stimulus is not desired, amplitude should
        be set to zero.
        
        A phase is defined as a period of time (ms) where the protocol remains
        consistent, i.e. with a stimulus of fixed amplitude being applied at a
        fixed rate.
        
        precedingTime & pulseDuration are denoted in milliseconds. Note that the
        the intervals in the intervalTimes matrix will be in seconds.
        
        frequencies are calculated in Hertz (/s).
        
        *** Note about space:
        When designing pacing protocols, remember that in its current form,
        this program stores high resolution data from the whole protocol.
        This can be fairly space intensive for long protocols which may
        need to be broken into a few steps by saving the Y & t-matrices and
        using the last row as the init_Y for the following run.
        Alternatively, to capture, for example, only the last few seconds
        of a protocol that has been run to steady state, modifications can
        be made in the loop in the odeSolver.
        %}
        function setUpPacingProtocol(k19CM, amplitudes, numPulses, precedingTime, pulseDurations)
            [k19CM.protocol.intervalTimes, k19CM.protocol.stimulus]...
                = f_setUpPacingProtocol(amplitudes, numPulses, precedingTime, pulseDurations) ;
            
            
            k19CM.protocol.amplitudes = amplitudes ;
            k19CM.protocol.numPulses = numPulses ;
            k19CM.protocol.precedingTime = precedingTime ;
            k19CM.protocol.pulseDuration = pulseDurations ;
            
            k19CM.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations) ; % for reference
            k19CM.protocol.phaseLengths(1) = 0.2 + (numPulses(1) - 1) * precedingTime(1) + numPulses(1) * pulseDurations(1) ; % first pulse doesn't have full preceding time
            k19CM.protocol.phaseLengths(end) = k19CM.protocol.phaseLengths(end) + 1 ; % add extra second at the end to capture last AP
            k19CM.protocol.frequencies = 1000./(precedingTime + pulseDurations) ; % for reference
            
            k19CM.protocol.totalTime = k19CM.protocol.intervalTimes(end,2) ;
            k19CM.protocol.nPhases = length(amplitudes) ;
            k19CM.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations)/1000 ;
            
            k19CM.protocol.intervalTimes = k19CM.protocol.intervalTimes*1000 ; % Kernik19 model uses ms not s
            k19CM.protocol.stimtimes = k19CM.protocol.intervalTimes(k19CM.protocol.stimulus > 0,1) ;
        end
        
        % FUNCTION INFO: setUpDrugApplication
        %{
        This function takes in a Paci18 object as well as information about the
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
        function setUpDrugApplication(k19CM, drugEffects, onTimes, offTimes)
            k19CM.conductances.applicationTimes = [onTimes; offTimes] ;
            k19CM.conductances.drugEffects = drugEffects ;
        end
        
        % FUNCTION INFO: odeSolver
        %{
        This function should be used on a Paci18 cell after it has been
        instantiated and its pacing protocol and drug protocol have been set up.
        The function will run through the simulation of the Paci18 CM's behavior
        according to these protocols and will fill out the t and Y fields of the
        state property to reflect how the state variables change over time.
        %}
        function odeSolver(k19CM)
            warning off
            warning(['k19cm odesolver starting for cell ', num2str(k19CM.cellID)]) % create a baseline warning
            
            if isempty(k19CM.protocol.intervalTimes)
                error('Error: Please set up pacing protocol first')
            else
                options = odeset('MaxStep', 1, 'InitialStep', 2e-2, 'Event', @f_tinyStepCatcher_k19) ; % solver settings % originally 1e-3, 2e-5 (s -> ms)
                
                % Allocate space initially to save time (~1220 rows per second
                % is usually enough). Extra zeros will be removed later to save
                % space.
                initialrows = 12500; % enough for 4 beats... for whole simulation: ceil(k19CM.protocol.totalTime*2000) ;
%                 initialrows = ceil(k19CM.protocol.totalTime*2000) ;
                k19CM.state.Y = zeros(initialrows, length(k19CM.state.init_Y)) ;
                k19CM.state.t = zeros(initialrows, 1) ; % first row at T0, no need to update
                
                k19CM.state.Y(1,:) = k19CM.state.init_Y ; % first row, input values (generally steady state)
                input_values = k19CM.state.init_Y ; % initialize, changes in for-loop
                index = 2 ; % initialize, changes in for-loop
                
                i_stimToKeep = length(k19CM.protocol.stimulus) - 9 ; % last 4 beats
%                 i_stimToKeep = 1 ; % whole simulation
                
                % convert structs to arrays for passing to ode15s:
                universals_arr = cell2mat(struct2cell(k19CM.universals));
                geometry_arr = cell2mat(struct2cell(k19CM.geometry));
                parameters_arr = k19CM.parameters.baseline .* k19CM.parameters.scaling;
                conductances_arr = k19CM.conductances.baseline .* k19CM.conductances.scaling;
                drugEffects_arr = k19CM.conductances.drugEffects ; 
                applicationTimes_arr = k19CM.conductances.applicationTimes;
                
                timer = 0;
                for i = 1:length(k19CM.protocol.stimulus)
                    tic;
                    if k19CM.isCellAlive == 1
                        input_values(end) = k19CM.protocol.stimulus(i) ; % last input value is the stimulus amplitude
                        
                        [t_current, Y_current, t_event, y_event, ~] = ode15s(...
                            k19CM.ODEModel,...
                            k19CM.protocol.intervalTimes(i,:),...
                            input_values,...
                            options,...
                            universals_arr,...
                            geometry_arr,...
                            parameters_arr,...
                            conductances_arr,...
                            drugEffects_arr,...
                            applicationTimes_arr);
                        
                        current_length = length(t_current)-1 ; % length of this interval in the stimulus protocol
                        
                        % check for short cycle lengths:
                        if current_length == 0
                            k19CM.isCellAlive = 0;
                            cell_state = k19CM.state;
                            cell_parameters = k19CM.x_conductances;
                            save(['GA/Results/events', num2str(k19CM.cellID), '_CL0.mat'], 'cell_state', 'cell_parameters');
                            k19CM.state.t = zeros(10, 1); %resultsless
                            k19CM.state.Y = zeros(10, length(k19CM.YNames)); %resultsless
                            break
                        end
                        
                        % check if there was an event
                        if (t_event) %if there was an event
                            numEvents = length(t_event) ;
                            if isempty(k19CM.event_number)
                                k19CM.event_number = 1;
                            end
                            k19CM.y_allevents(:, k19CM.event_number:k19CM.event_number + numEvents - 1) = y_event' ;
                            k19CM.x_allevents(:, k19CM.event_number) = k19CM.x_conductances' ;
                            y_allEvents = k19CM.y_allevents ;
                            x_allEvents = k19CM.x_allevents ;
                            save(['GA/Results/events', num2str(k19CM.cellID), '_', num2str(k19CM.event_number), '.mat'], 'y_allEvents', 'x_allEvents') ;
                            k19CM.isCellAlive = 0 ;
                            k19CM.state.t = zeros(10, 1); %resultsless
                            k19CM.state.Y = zeros(10, length(k19CM.YNames)); %resultsless
                            k19CM.event_number = k19CM.event_number + numEvents ;
                            break
                        end
                    end
                    
                    % cells may have been killed off in the above loops so
                    % we check again for survival
                    if k19CM.isCellAlive == 1
                        if i >= i_stimToKeep
                            
                            k19CM.state.t(index:index+current_length-1) = t_current(2:end) ; % save outcome of this run
                            k19CM.state.Y(index:index+current_length-1,:) = Y_current(2:end,:) ; % save outcome of this run
                            index = index + current_length;
                            input_values = k19CM.state.Y(index - 1, :); % set input values for next interval to end of current interval
                        else
                            input_values = Y_current(end, :) ;
                        end
                    end
                    timer = timer + toc ;
                    %%
                    if timer > 180 % seconds; adjust this to change maximum simulation time
                        cell_state = k19CM.state;
                        cell_parameters = k19CM.x_conductances ;
                        save(['GA/Results/events', num2str(k19CM.cellID), '_timer.mat'], 'cell_state', 'cell_parameters', 'timer') ;
                        k19CM.isCellAlive = 0 ;
                        k19CM.state.t = zeros(10, 1); %resultsless
                        k19CM.state.Y = zeros(10, length(k19CM.YNames)); %resultsless
                        break
                    end
                end
                %disp(timer)
                
                if ~isempty(k19CM.state.t) && k19CM.isCellAlive == 1
                    % remove extra zeros to save space
                    zeroInds = find(~k19CM.state.t) ; % returns indices of all zeroes
                    firstExtraZero = zeroInds(2) ; % first zero is T0
                    k19CM.state.t = k19CM.state.t(1:firstExtraZero-1) ;
                    k19CM.state.Y = k19CM.state.Y(1:firstExtraZero-1,:) ;
                    
                    k19CM.state.t = k19CM.state.t(2:end) ;
                    k19CM.state.Y = k19CM.state.Y(2:end,:) ;
                    
                    if sum(k19CM.protocol.amplitudes) % if not spontaneous
                    k19CM.protocol.last3stimTimes = k19CM.protocol.intervalTimes([end-5, end-3, end-1], 1) - k19CM.state.t(1);
                    else 
                    k19CM.protocol.last3stimTimes = [];
                    end
                    k19CM.state.t = k19CM.state.t - k19CM.state.t(1); % set t0 to 0
                end
            end
        end
        
        % FUNCTION INFO: getCurrents
        %{
        This function should be used on a Kernik19 cell after the odeSolver has
        been run. It will use the ODE model for the cell to find  & store the
        values of the different currents , and E_ions at the time points
        collected by the ODE solver.
        %}
        function getCurrents(k19CM)
            if isempty(k19CM.state.t)
                error('Error: Please run the odeSolver on this Paci18 CM first')
            else
                k19CM.state.currents = zeros(length(k19CM.state.t), length(k19CM.currentNames));
                numLoops = length(k19CM.state.t) ;
                
                % convert structs to arrays for passing to ode15s:
                universals_arr = cell2mat(struct2cell(k19CM.universals));
                geometry_arr = cell2mat(struct2cell(k19CM.geometry));
                parameters_arr = k19CM.parameters.baseline .* k19CM.parameters.scaling;
                conductances_arr = k19CM.conductances.baseline .* k19CM.conductances.scaling;
                drugEffects_arr = k19CM.conductances.drugEffects ; 
                applicationTimes_arr = k19CM.conductances.applicationTimes;
                
                for i= 1:numLoops
                    [~, data]    = f_gaKernik19(k19CM.state.t(i), k19CM.state.Y(i,:),...
                        universals_arr,...
                            geometry_arr,...
                            parameters_arr,...
                            conductances_arr,...
                            drugEffects_arr,...
                            applicationTimes_arr);
                    k19CM.state.currents(i,:) = data;
                    
                end
                %k19CM.state.currents(:,end) = sum(k19CM.state.currents(:,1:16), 2); % sum for Itot; the 2 at the end is so we sum the rows not the columns
                %'k19 Currents Done!'
            end
        end
    end
end
