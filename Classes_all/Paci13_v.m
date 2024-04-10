%{
Paci13_v is a concrete class for a ventricular cardiomyocyte based on the model
published by Paci et al. in 2013. It inherits from the superclass AbsCM which is
an abstract class representing cardiomyocytes based on any model.
%}
classdef Paci13_v < AbsCM
    properties % I want to restrict access to this class
        name = 'Paci13v CM' ; 
        ODEModel = @f_Paci13_v ;
        
        geometry = struct('V_SR', 583.73, ... �m^3
            'Vc', 8800.0, ... �m^3
            'Cm', 9.87109e-11) ; % F, we keep capacitance in geometry for convenience
        
        % For informational queries & keeping order consistent across versions
        
        YNames = {'h', 'j', 'm', 'd', 'fCa', 'f1', 'f2', 'r', 'q', 'Xr1', 'Xr2', 'Xs', 'Xf', 'g', 'V', 'Nai', 'Cai', 'Ca_S_R', 'stim'} ;
        YUnits = {'-', '-', '-', '-', '-',   '-',  '-',  '-', '-', '-',   '-',   '-',  '-',  '-', 'V', 'mM',  'mM',  'mM',     'A'} ;
        
        currentNames = {'INa', 'If', 'ICaL', 'Ito', 'IKs', 'IKr', 'IK1', 'INaCa', 'INaK', 'IPCa', 'IbNa', 'IbCa', 'Istim', 'E_K', 'E_Na', 'Itot'};
        currentUnits = {'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'mV', 'mV', 'A', 'A'} ;
        
        state = struct('init_Y', [], 'Y', [], 't', [], 'currents', []) ;
        
        % the choice of what is in conductances versus parameters below is
        % the legacy of prior versions of the model. I intend to discuss
        % this with Dr. Sobie to divide them more thoughtfully, based on
        % sets that we would want to vary simultaneously. At present, both
        % structs contain a 'scaling' field. This may change if it is
        % determined that only one set should be scaled.
        conductanceNames = {'GNa', 'Gf',  'GCaL',       'Gto', 'GKs', 'GKr', 'GK1', 'GpCa', 'GbNa', 'GbCa', 'Vmaxup', 'G_RyR', 'KNaCA', 'PNaK', 'Vleak', 'G_SERCA'} ;
        conductanceUnits = {'S/F', 'S/F', 'm^3 /(Fxs)', 'S/F', 'S/F', 'S/F', 'S/F', 'A/F',  'S/F',  'S/F',  'mM/s',   '-',     'A/F',   'A/F',  '1/s', '-'} ;
        
        conductances = struct('baseline', [3671.2302, 30.10312, 8.64E-05, 29.9038, 2.041, 29.8667, 28.1492, 0.4125, 0.9, 0.69264, 0.56064, 1, 4900, 1.841424, 0.00044444, 1],...
            'scaling', ones(1,16),...
            'drugEffects', ones(1,16),...
            'applicationTimes', zeros(2,16)) ;
        
        parameterNames = {'arel', 'brel', 'crel', 'Bufc', 'Bufsr', 'Kbufc', 'Kbufsr', 'Kup', 'KpCa', 'L0', 'Pkna', 'Ksat', 'KmCa', 'KmNai', 'alpha', 'gamma', 'KmNa', 'Kmk'} ;
        parameterUnits = {'mM/s', 'mM',   'mM/s', 'mM',   'mM',    'mM',    'mM',     'mM',  'mM',   '-',  '-',    '-',    'mM',   'mM',    '-',     '-',     'mM',   'mM'} ;
        
        parameters = struct('baseline', [16.464, 0.25, 8.232, 0.25, 10, 0.001, 0.3, 0.00025, 0.0005, 0.025, 0.03, 0.1, 1.38, 87.5, 2.8571432, 0.35, 40, 1],...
            'scaling', ones(1,18))
        
        protocol = struct('intervalTimes', [] , 'stimulus', [], 'amplitudes', [],...
            'numPulses', [], 'precedingTime', [], 'pulseDuration', [], 'totalTime', [],...
            'phaseLengths', [], 'frequencies', []) ;
    end
    
    methods
        
        % CONSTRUCTOR:
        %{
        creates an instance of Paci13_v with initial state variables equal to
        init_Y.
        %}
        function p13vCM = Paci13_v(init_Y)
            % some options for init_Y: --> save others as convenient
            %{
            Published in Paci et al. 2013:
            Y = [0.75, 0.75, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0.1, 1, -70e-3, 10, 0.0002, 0.3, 0];
            
            Spontaneous beating for 1600s:
            Yspont1600 = [0.332847821564613,0.0293131941172596,0.157776996637716,0.000170600470008095,0.996053809629183,0.539449140622959,0.965935296023459,0.00870364766566940,0.704692092161293,0.696795347884297,0.394780776090026,0.0549507879446617,0.0881092584899003,0.999952129071328,-0.0667641338617266,11.2442465250989,6.57896709839818e-05,0.170871363004459,0] ;
            %}
            p13vCM.state.init_Y = init_Y ;
        end
        
        % FUNCTION INFO: scaleCell
        %{
        This function takes a set of scale factors for the conductances and
        a set of scale factors for the parameters and saves them in the
        corresponding struct. The default value for both scaling fields is
        a vector of 1s, so this function need only be used when alternative
        scaling is required.
        %}
        function scaleCell(p13vCM, conductanceScaleFactors, parameterScaleFactors)
            p13vCM.parameters.scaling = parameterScaleFactors ;
            p13vCM.conductances.scaling = conductanceScaleFactors ;
        end
        
        % FUNCTION INFO: setUpPacingProtocol
        %{
        This function takes in a Paci13_v object as well as information about the
        pacing protocol and sets the protocol property's intervalTimes field to
        a two-column matrix containing the start and end times of each interval
        in the protocol, where an interval is defined as either the time
        preceding the stimulus, or the time during the stimulus.
        The stimulus field is set to the amplitude of the stimulus
        applied to the p18CM in each of these intervals, zero before and some
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
        function setUpPacingProtocol(p13vCM, amplitudes, numPulses, precedingTime, pulseDurations)
            [p13vCM.protocol.intervalTimes, p13vCM.protocol.stimulus]...
                = f_setUpPacingProtocol(amplitudes, numPulses, precedingTime, pulseDurations) ;
            
            p13vCM.protocol.amplitudes = amplitudes ;
            p13vCM.protocol.numPulses = numPulses ;
            p13vCM.protocol.precedingTime = precedingTime ;
            p13vCM.protocol.pulseDuration = pulseDurations ;
            
            p13vCM.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations) ; % for reference
            p13vCM.protocol.phaseLengths(1) = 0.2 + (numPulses(1) - 1) * precedingTime(1) + numPulses(1) * pulseDurations(1) ; % first pulse doesn't have full preceding time
            p13vCM.protocol.phaseLengths(end) = p13vCM.protocol.phaseLengths(end) + 1 ; % add extra second at the end to capture last AP
            p13vCM.protocol.frequencies = 1000./(precedingTime + pulseDurations) ; % for reference
            
            p13vCM.protocol.totalTime = p13vCM.protocol.intervalTimes(end,2) ;
            p13vCM.protocol.nPhases = length(amplitudes) ;
            p13vCM.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations)/1000 ;
        end
        
        % FUNCTION INFO: setUpDrugApplication
        %{
        This function takes in a Paci13_v object as well as information about the
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
        function setUpDrugApplication(p13vCM, drugEffects, onTimes, offTimes)
            if sum(onTimes > offTimes) >= 1 % if any onTimes are come after their correspondng off times throw an error
                error("Error: onTimes must be less than their corresponding offTimes")
            else
                p13vCM.conductances.applicationTimes = [onTimes; offTimes] ;
                p13vCM.conductances.drugEffects = drugEffects ;
            end
        end
        
        % FUNCTION INFO: odeSolver
        %{
        This function should be used on a Paci13_v cell after it has been
        instantiated and its pacing protocol and drug protocol have been set up.
        The function will run through the simulation of the Paci13_v CM's behavior
        according to these protocols and will fill out the t and Y fields of the
        state property to reflect how the state variables change over time.
        %}
        function odeSolver(p13vCM)
            if isempty(p13vCM.protocol.intervalTimes)
                error('Error: Please set up pacing protocol first')
            else
                options = odeset('MaxStep', 1e-3, 'InitialStep', 2e-5) ; % solver settings
                
                % Allocate space initially to save time (~1220 rows per second
                % is usually enough). Extra zeros will be removed later to save
                % space.
                initialrows = ceil(p13vCM.protocol.totalTime*2000) ;
                p13vCM.state.Y = zeros(initialrows, length(p13vCM.state.init_Y)) ;
                p13vCM.state.t = zeros(initialrows, 1) ; % first row at T0, no need to update
                
                p13vCM.state.Y(1,:) = p13vCM.state.init_Y ; % first row, input values (generally steady state)
                input_values = p13vCM.state.init_Y ; % initialize, changes in for-loop
                index = 2 ; % initialize, changes in for-loop
                
                for i = 1:length(p13vCM.protocol.stimulus)
                    input_values(end) = p13vCM.protocol.stimulus(i) ; % last input value is the stimulus amplitude
                    [t_current, Y_current] = ode15s(p13vCM.ODEModel, p13vCM.protocol.intervalTimes(i,:), input_values, options, p13vCM.geometry, p13vCM.conductances, p13vCM.parameters);
                    current_length = length(t_current)-1 ; % length of this interval in the stimulus protocol
                    
                    p13vCM.state.t(index:index+current_length-1) = t_current(2:end) ; % save outcome of this run
                    p13vCM.state.Y(index:index+current_length-1,:) = Y_current(2:end,:) ; % save outcome of this run
                    index = index + current_length;

                    input_values = p13vCM.state.Y(index - 1, :) ; % set input values for next interval to end of current interval
                end
                
                % remove extra zeros to save space
                zeroInds = find(~p13vCM.state.t) ; % returns indices of all zeroes
                firstExtraZero = zeroInds(2) ; % first zero is T0
                p13vCM.state.t = p13vCM.state.t(1:firstExtraZero-1) ;
                p13vCM.state.Y = p13vCM.state.Y(1:firstExtraZero-1,:) ;
                'p13v ODE Done!'
            end
        end
        
        % FUNCTION INFO: getCurrents
        %{
        This function should be used on a Paci18 cell after the odeSolver has
        been run. It will use the ODE model for the cell to find  & store the
        values of the different currents , and E_ions at the time points
        collected by the ODE solver.
        %}
        function getCurrents(p13vCM)
            if isempty(p13vCM.state.t)
                error('Error: Please run the odeSolver on this Paci18 CM first')
            else
                p13vCM.state.currents = zeros(length(p13vCM.state.t), length(p13vCM.currentNames));
                numLoops = length(p13vCM.state.t) ;
                for i= 1:numLoops
                    [~, data]    = f_Paci13_v(p13vCM.state.t(i), p13vCM.state.Y(i,:), p13vCM.geometry, p13vCM.conductances, p13vCM.parameters); % ~ says only return data, not dY output
                    p13vCM.state.currents(i,1:end-1) = data;
                   % {'p13v currents:', i/numLoops*100} % display in command window for a progress tracker, runs 0 to 100
                end
                p13vCM.state.currents(:,end) = sum(p13vCM.state.currents(:,1:13), 2); % sum for Itot; the 2 at the end is so we sum the rows not the columns
                'p13v Currents Done!'
            end
        end
        
        % FUNCTION INFO: makePlots
        %{
        This function should be used on a Paci13_v cell after the odeSolver and
        getCurrents functions have been run. It will produce plots associated
        with the Paci13_v CM. Customize according to the protocol you choose and
        the figures you are trying to generate. This function currently includes
        potentially useful code blocks for specific types of figures which can
        be modified or commented out if not needed.
        %}
        function makePlots(p13vCM)
            if isempty(p13vCM.state.t) % check that simulation has been run
                error('Error: Please run the odeSolver on this Paci13_v CM first')
                %             elseif isempty(p13vCM.state.currents)
                %                 error('Error: Please run the getCurrents function on this Paci13_v CM first')
            else
 
                %{
                % ------------------------------------------------------------------
                % the following block of code divides Y, t, and currents data into
                % pacing phases, useful for comparing different parts of a protocol.
                % Required for 'SUPER ZOOMED OVERLAY' plot.
%                 if p13vCM.protocol.nPhases > 1
%                     cutoffs = cell(1, p13vCM.protocol.nPhases);
%                     cutoff_i = zeros(1,p13vCM.protocol.nPhases);
%                     for i = 1:p13vCM.protocol.nPhases
%                         cutoffs{i} = find(p13vCM.state.t == sum(p13vCM.protocol.phaseLengths(1:i)));
%                         cutoff_i(i) = cutoffs{i}(1) ;
%                     end
%                     cutoff_i = [1, cutoff_i]; % indices where each phase switches over
%                     
%                     tPhases = cell(p13vCM.protocol.nPhases,1); % t for each of the phases
%                     YPhases = cell(p13vCM.protocol.nPhases,1); % Y for each of the phases
%                     currentPhases = cell(p13vCM.protocol.nPhases,1); % currents for each of the phases
%                     
%                     for i = 1:p13vCM.protocol.nPhases
%                         tPhases{i} = p13vCM.state.t(cutoff_i(i):cutoff_i(i+1))...
%                             -p13vCM.state.t(cutoff_i(i)) ; % set so time starts at zero for each phase
%                         YPhases{i} = p13vCM.state.Y(cutoff_i(i):cutoff_i(i+1),:);
%                         currentPhases{i} = p13vCM.state.currents(cutoff_i(i):cutoff_i(i+1),:);
%                     end
%                     
%                 else
%                     tPhases = p13vCM.state.t ;
%                     YPhases = p13vCM.state.Y ;
%                 end
%                 
%                 % ------------------------------------------------------------------
%                 % SIMPLE FIGURES - for individual currents or state variables, any
%                 % range of time
%                 
%                 % dataForPlot - choose from the following according to your needs:
%                 %{
%                 x:
%                 p18CM.state.t
%                 tPhases{i}  ***
%                 y:
%                 p18CM.state.Y
%                 p18CM.state.currents
%                 YPhases{i}  ***
%                 currentPhases{i}  ***
%             
%                 *** where i is the phase number of interest
%                 %}
%}
                dataForPlot_x = p13vCM.state.t ;
                dataForPlot_y = p13vCM.state.Y ;
                n = length(dataForPlot_x) ;
                lower = ceil(1/n * n); % selects first index for plot
                upper = floor(n/n * n); % selects last index for plot
                % max range is (1/n * n) to (n/n * n) but can narrow it down
                % e.g. (0.1 * n) to (0.9 * n) plots 10th percentile index to
                % 90th percentile index of data
                
                % Repeat & modify the below figure units as required.
                % If necessary, reset dataForPlot_x, dataForPlot_y, n, lower, and
                % upper between figures
               
                % YNames = {'h', 'j', 'm', 'd', 'fCa', (1-5)
                %           'f1', 'f2', 'r', 'q', 'Xr1', (6-10)
                %           'Xr2', 'Xs', 'Xf', 'g', 'V', (7-15)
                %           'Nai', 'Cai', 'Ca_S_R', 'stim'} (16-19)

                
                figure
                % consult YNames or currentNames to select correct column of
                % dataForPlot_y, ylabel, and title
                hold on
            plot(dataForPlot_x(lower:upper),dataForPlot_y(lower:upper,15),'b', 'LineWidth', 1.1)
            plot(dataForPlot_x(lower:upper),dataForPlot_y(lower:upper,19)/1000,'r', 'LineWidth', 1.1)
                xlabel('Time (s)')
                ylabel('V (V)')
                title('p13v Action potential with different pacing rates')
                
%                 figure
%                 % consult YNames or currentNames to select correct column of
%                 % dataForPlot_y, ylabel, and title
%                 plot(dataForPlot_x(lower:upper),dataForPlot_y(lower:upper,18),'b')
%                 xlabel('Time (s)')
%                 ylabel('[Ca]_SR (mM)')
%                 title('SR-Calcium with different pacing rates')
                
%{
                 % -----------------------------------------------------------------
%                 % SUPER ZOOMED OVERLAY - use when you want to compare a single
%                 % action potential at different pacing rates for example
%                 plotend = min(p13vCM.protocol.phaseLengths) - 0.3; % where
%                 % to zoom to, make sure this is within all phases you want to
%                 % test. Possibility to reset this for each phase with a loop if
%                 % you want a different time point for each phase
%                 
%                 plotstart = plotend - 2; % this plot will cover this amount of
%                 % time. It helps to zoom in gradually & tweak as you go, and
%                 % make sure you don't overshoot the beginning of any of your
%                 % .phases
%                 
%                 % sometimes you'll want to nudge the trace of a spontaneously
%                 % beating CM, such as when you're comparing traces from different
%                 % pacing frequencies. Set spontoffset to do so by trial and error.
%                 spontOffset = 0;
%                 dataForPlot_y = YPhases ; % choose YPhases or currentPhases
%                 namesForPlot_y = p13vCM.YNames ; % choose YNames or currentNames
%                 unitsForPlot_y = p13vCM.YUnits ; % choose YUnits or currentUnits
%                 % choose columnsToPlot by consulting YNames or currentNames e.g.
%                 % [1, 3, 4]. We will generate a separate figure for each.
%                 columnsToPlot = 15 ;
%                 
%                 phase1Spont = (p13vCM.protocol.amplitudes(1) == 0) ; % boolean = 0 or 1 depending on if first phase is spontaneous
%                 
%                 % find indices in each phase for the specific peak you want to plot
%                 plotlengthindex = cell(p13vCM.protocol.nPhases, 1) ;
%                 plotlengthindex_i = zeros(2, p13vCM.protocol.nPhases);
%                 if phase1Spont % special handling if first phase is spontaneous
%                     plotlengthindex{1} = find((tPhases{1} >= plotstart-spontOffset)&(tPhases{1} <= plotend-spontOffset)) ;
%                 end
%                 for i = 1+phase1Spont:p13vCM.protocol.nPhases
%                     plotlengthindex{i} = find((tPhases{i} >= plotstart)&(tPhases{i} <= plotend)) ;
%                 end
%                 for i = 1:p13vCM.protocol.nPhases
%                     plotlengthindex_i(:,i) = [plotlengthindex{i}(1);plotlengthindex{i}(end)] ;
%                 end
%                 
%                 % plot the plot(s)!
%                 for j = 1:length(columnsToPlot)
%                     figure
%                     hold on
%                     legendstr = cell(1,p13vCM.protocol.nPhases) ;
%                     if phase1Spont % special handling if first phase is spontaneous
%                         legendstr{1} = ['spontaneous,\newlineshifted by ', num2str(spontOffset),'s'] ;
%                         
%                         plot(tPhases{1}(plotlengthindex_i(1,1):plotlengthindex_i(2,1))+spontOffset,...
%                             dataForPlot_y{1}(plotlengthindex_i(1,1):plotlengthindex_i(2,1),columnsToPlot(j)),...
%                             'LineWidth', 1.5)
%                     end
%                     
%                     for i = 1+phase1Spont:p13vCM.protocol.nPhases
%                         plot(tPhases{i}(plotlengthindex_i(1,i):plotlengthindex_i(2,i)),...
%                             dataForPlot_y{i}(plotlengthindex_i(1,i):plotlengthindex_i(2,i),columnsToPlot(j)),...
%                             'LineWidth', 1.5)
%                         legendstr{i} = [num2str(p13vCM.protocol.frequencies(i)), ' Hz'];
%                     end
%                     
%                     xlabel('Time (s)')
%                     ylabel([namesForPlot_y{columnsToPlot(j)},' (', unitsForPlot_y{columnsToPlot(j)},')'])
%                     title([namesForPlot_y{columnsToPlot(j)}, ' at steady state at different rates'])
%                     legend(legendstr)
%                 end
%}
             end
         end
        
        % FUNCTION INFO: getInfo
        %{
        This function returns information about the Paci13_v CM it takes as input.
        My intention is that the Paci13_v properties be kept private unless
        obtained through this getter but for now I am having difficulty managing
        access, so the class is public.
        %}
        function info = getInfo(p13vCM)
            info.geometry = p13vCM.geometry ;
            info.paramNames = p13vCM.paramNames ;
            info.paramUnits = p13vCM.paramUnits ;
            info.YNames = p13vCM.YNames ;
            info.YUnits = p13vCM.YUnits ;
            info.currentNames = p13vCM.currentNames ;
            info.state = p13vCM.state;
        end
    end
    
end
