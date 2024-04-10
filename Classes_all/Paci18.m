%{
Paci18 is a concrete class for a cardiomyocyte based on the model published by
Paci et al. in 2018. It inherits from the superclass AbsCM which is an abstract
class representing cardiomyocytes based on any model.
%}
classdef Paci18 < AbsCM
    properties % I want to restrict access to this class
        name = 'Paci18 CM' ; 
        ODEModel = @f_Paci18 ;
        
        geometry = struct('V_SR', 583.73, ... 痠^3
            'Vc', 8800.0, ... 痠^3
            'Cm', 9.87109e-11) ; % F, we keep capacitance in geometry for convenience
        
        % For informational queries & keeping order consistent across versions
        YNames = {'Vm', 'Ca_S_R', 'Cai', 'g', 'd', 'f1', 'f2', 'fCa', 'Xr1', 'Xr2', 'Xs', 'h', 'j', 'm', 'Xf', 'q', 'r', 'Nai', 'm_L', 'h_L', 'RyRa', 'RyRo', 'RyRc', 'stim'} ;
        YUnits = {'V',  'mM',     'mM',  '-', '-', '-',  '-',  '-',   '-',   '-',   '-',  '-', '-', '-', '-',  '-', '-', 'mM',  '-',   '-',   '-',    '-',    '-',    'A?'} ;
        currentNames = {'INa', 'If', 'ICaL', 'Ito', 'IKs', 'IKr', 'IK1', 'INaCa', 'INaK', 'IpCa', 'IbNa', 'IbCa', 'Irel', 'Iup', 'Ileak', 'Istim', 'E_K', 'E_Na', 'INaL', 'Itot'};
        currentUnits = {'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'mV', 'mV', 'A', 'A'} ; 
        
        state = struct('init_Y', [], 'Y', [], 't', [], 'currents', []) ;
        
        % For informational queries & keeping order consistent across versions
        paramNames = {'RyRa1', 'RyRa2', 'RyRahalf', 'RyRohalf', 'RyRchalf', 'Kup', 'alpha'} ;
        paramUnits = {'然',    '然',    '然',       '然',       '然',       'mM',  '-'} ;
        
        parameters = struct( 'baseline', [0.05354, 0.0488, 0.02427, 0.01042, 0.00144, 3.1928e-4, 2.5371],...
            'scaling', ones(1,7)) ;
        
        % TO DO: create properties for three-tiered parameters - commonly
        % varied, locally varied, constant
        
        % TO DO: increase sophistication & coverage of drug protocol
        
        conductanceNames = {'g_Na','g_f', 'g_CaL',     'g_to', 'g_Ks', 'g_Kr', 'g_K1', 'g_PCa', 'g_b_Na', 'g_b_Ca', 'VmaxUp', 'g_irel_max', 'kNaCa', 'PNaK', 'V_leak', 'GNaLmax', 'g_serca'} ;
        conductanceUnits = {'S/F', 'S/F', 'm^3/(F.s)', 'S/F',  'S/F',  'S/F',  'S/F',  'A/F',   'S/F',    'S/F',    'mM/s',   'mM/s',       'A/F',   'A/F',  '/s',     'S/F',     '-'} ;
        conductances = struct( 'baseline', [3671.2302, 30.10312, 8.635702e-5, 29.9038, 2.041, 29.8667, 28.1492, 0.4125, 0.95, 0.727272, 0.5113, 62.5434, 3917.0463, 2.6351, 4.7279e-4, 2.3*7.5, 1],...
            'scaling', ones(1,17),...
            'drugEffects', ones(1,17),...
            'applicationTimes', zeros(2,17)) ;
        
        protocol = struct('intervalTimes', [] , 'stimulus', [], 'amplitudes', [],...
            'numPulses', [], 'precedingTime', [], 'pulseDuration', [], 'totalTime', [],... 
            'phaseLengths', [], 'frequencies', []) ;
    end
    
    methods
        
        % CONSTRUCTOR:
        %{
        creates an instance of Paci18 with initial state variables equal to
        init_Y.
        %}
        function p18CM = Paci18(init_Y)
            % some options for init_Y:
            %{
            Published in Paci et al. 2018:
            Y = [-0.070 0.32 0.0002 0 0 1 1 1 0 1 0 0.75 0.75 0 0.1 1 0 9.2 0 0.75 0.3 0.9 0.1 0];
        
            SS at 800s spontaneous beating:
            Yspont = [-0.0749228904740065 0.0936532528714175 3.79675694306440e-05 0 8.25220533963093e-05 0.741143500777858 0.999983958619179 0.997742015033076 0.266113517200784 0.434907203275640 0.0314334976383401 0.745356534740988 0.0760523580322096 0.0995891726023512 0.0249102482276486 0.841714924246004 0.00558005376429710 8.64821066193476 0.00225383437957339 0.0811507312565017 0.0387066722172937 0.0260449185736275 0.0785849084330126 0];
        
            SS at 1600s @ 1Hz Pacing:
            Y1Hz = [-0.0729913843090347 0.0708338064144874 1.43599887702576e-05 0 0.000108074980807618 0.952981801622850 0.999977361637709 0.999140156035565 0.0224706562993318 0.425731226641029 0.0349693801962028 0.748585039181593 0.190134709490443 0.110909372549615 0.0738817114617895 0.824952784655803 0.00614766851990035 6.70121650951907 0.00323518754688945 0.128881879258908 0.0371809524458325 1.48498886591972e-05 0.914405999740954 0];
        
            SS at 800s @ 0.75Hz Pacing:
            Y0pt75Hz = [-0.0655767307105754 0.0570603513259588 1.18827126579712e-05 0 0.000312323580542382 0.994066949323603 0.999866825022176 0.999286148680894 0.00419112964847023 0.390002974058969 0.0535535727025202 0.550146287415602 0.240034242204755 0.168629280685409 0.0857386539061896 0.733732403170870 0.00907517028499561 5.76272608086791 0.0131540039527741 0.0811077552198685 0.0262038617397998 0.000239253301063336 0.989415842418175 0];
            %}
            p18CM.state.init_Y = init_Y ;
        end
        
        % FUNCTION INFO: scaleCell
        %{
        This function takes two sets of scale factors to modify the parameters
        and conductances and saves it in the relevant struct. The default value 
        for the scaling fields is a vector of 1s, so this function need only be 
        used when alternative scaling is required.
        %}
        function scaleCell(p18CM, conductanceScaleFactors, parameterScaleFactors)
            p18CM.conductances.scaling = conductanceScaleFactors ;
            p18CM.parameters.scaling = parameterScaleFactors ;
        end
        
        % FUNCTION INFO: setUpPacingProtocol
        %{
        This function takes in a Paci18 object as well as information about the
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
        function setUpPacingProtocol(p18CM, amplitudes, numPulses, precedingTime, pulseDurations)
            [p18CM.protocol.intervalTimes, p18CM.protocol.stimulus]...
                = f_setUpPacingProtocol(amplitudes, numPulses, precedingTime, pulseDurations) ;
            
            p18CM.protocol.amplitudes = amplitudes ;
            p18CM.protocol.numPulses = numPulses ; 
            p18CM.protocol.precedingTime = precedingTime ; 
            p18CM.protocol.pulseDuration = pulseDurations ;
            
            p18CM.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations) ; % for reference
            p18CM.protocol.phaseLengths(1) = 0.2 + (numPulses(1) - 1) * precedingTime(1) + numPulses(1) * pulseDurations(1) ; % first pulse doesn't have full preceding time
            p18CM.protocol.phaseLengths(end) = p18CM.protocol.phaseLengths(end) + 1 ; % add extra second at the end to capture last AP
            p18CM.protocol.frequencies = 1000./(precedingTime + pulseDurations) ; % for reference
            
            p18CM.protocol.totalTime = p18CM.protocol.intervalTimes(end,2) ;
            p18CM.protocol.nPhases = length(amplitudes) ;
            p18CM.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations)/1000 ;
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
        function setUpDrugApplication(p18CM, drugEffects, onTimes, offTimes)
            p18CM.conductances.applicationTimes = [onTimes; offTimes] ;
            p18CM.conductances.drugEffects = drugEffects ;
        end
        
        % FUNCTION INFO: odeSolver
        %{
        This function should be used on a Paci18 cell after it has been
        instantiated and its pacing protocol and drug protocol have been set up.
        The function will run through the simulation of the Paci18 CM's behavior
        according to these protocols and will fill out the t and Y fields of the
        state property to reflect how the state variables change over time.
        %}
        function odeSolver(p18CM)
            if isempty(p18CM.protocol.intervalTimes)
                error('Error: Please set up pacing protocol first')
            else
                options = odeset('MaxStep', 1e-3, 'InitialStep', 2e-5) ; % solver settings
                
                % Allocate space initially to save time (~1220 rows per second
                % is usually enough). Extra zeros will be removed later to save
                % space.
                initialrows = ceil(p18CM.protocol.totalTime*2000) ;
                p18CM.state.Y = zeros(initialrows, length(p18CM.state.init_Y)) ;
                p18CM.state.t = zeros(initialrows, 1) ; % first row at T0, no need to update
                
                p18CM.state.Y(1,:) = p18CM.state.init_Y ; % first row, input values (generally steady state)
                input_values = p18CM.state.init_Y ; % initialize, changes in for-loop
                index = 2 ; % initialize, changes in for-loop
                
                for i = 1:length(p18CM.protocol.stimulus)
                    input_values(end) = p18CM.protocol.stimulus(i) ; % last input value is the stimulus amplitude
                    
                    [t_current, Y_current] = ode15s(p18CM.ODEModel, p18CM.protocol.intervalTimes(i,:), input_values, options, p18CM.geometry, p18CM.conductances, p18CM.parameters);
                    current_length = length(t_current)-1 ; % length of this interval in the stimulus protocol
                    
                    p18CM.state.t(index:index+current_length-1) = t_current(2:end) ; % save outcome of this run
                    p18CM.state.Y(index:index+current_length-1,:) = Y_current(2:end,:) ; % save outcome of this run
                    index = index + current_length;

                    input_values = p18CM.state.Y(index - 1, :) ; % set input values for next interval to end of current interval
                end
                
                % remove extra zeros to save space
                zeroInds = find(~p18CM.state.t) ; % returns indices of all zeroes
                firstExtraZero = zeroInds(2) ; % first zero is T0
                p18CM.state.t = p18CM.state.t(1:firstExtraZero-1) ;
                p18CM.state.Y = p18CM.state.Y(1:firstExtraZero-1,:) ;
                'p18 ODE Done!'
            end
        end
        
        % FUNCTION INFO: getCurrents
        %{
        This function should be used on a Paci18 cell after the odeSolver has
        been run. It will use the ODE model for the cell to find  & store the
        values of the different currents , and E_ions at the time points
        collected by the ODE solver.
        %}
        function getCurrents(p18CM)
            if isempty(p18CM.state.t)
                error('Error: Please run the odeSolver on this Paci18 CM first')
            else
                p18CM.state.currents = zeros(length(p18CM.state.t), length(p18CM.currentNames));
                numLoops = length(p18CM.state.t) ;
                for i= 1:numLoops
                    [~, data]    = f_Paci18(p18CM.state.t(i), p18CM.state.Y(i,:), p18CM.geometry, p18CM.conductances, p18CM.parameters); % ~ says only return data, not dY output
                    p18CM.state.currents(i,1:end-1) = data;
           %         {'p18 currents:', i/numLoops*100} % display in command window for a progress tracker, runs 0 to 100
                    
                end
                p18CM.state.currents(:,end) = sum([p18CM.state.currents(:,1:16), p18CM.state.currents(:,19)], 2); % sum for Itot; the 2 at the end is so we sum the rows not the columns
                'p18 Currents Done!'
            end
        end
        
        % FUNCTION INFO: makePlots
        %{
        This function should be used on a Paci18 cell after the odeSolver and 
        getCurrents functions have been run. It will produce plots associated 
        with the Paci18 CM. Customize according to the protocol you choose and 
        the figures you are trying to generate. This function currently includes
        potentially useful code blocks for specific types of figures which can
        be modified or commented out if not needed.
        %}
        function makePlots(p18CM)
             if isempty(p18CM.state.t) % check that simulation has been run
                error('Error: Please run the odeSolver on this Paci18 CM first')
            elseif isempty(p18CM.state.currents)
                error('Error: Please run the getCurrents function on this Paci18 CM first')
            else 
                
            % ------------------------------------------------------------------     
            % the following block of code divides Y, t, and currents data into 
            % pacing phases, useful for comparing different parts of a protocol.
            % Required for 'SUPER ZOOMED OVERLAY' plot. 
            
%             cutoffs = cell(1, p18CM.protocol.nPhases); 
%             cutoff_i = zeros(1,p18CM.protocol.nPhases);
%             for i = 1:p18CM.protocol.nPhases
%                 cutoffs{i} = find(p18CM.state.t == sum(p18CM.protocol.phaseLengths(1:i)));
%                 cutoff_i(i) = cutoffs{i}(1) ;
%             end
%             cutoff_i = [1, cutoff_i]; % indices where each phase switches over
%             
%             tPhases = cell(p18CM.protocol.nPhases,1); % t for each of the phases
%             YPhases = cell(p18CM.protocol.nPhases,1); % Y for each of the phases
%             currentPhases = cell(p18CM.protocol.nPhases,1); % currents for each of the phases
%             
%             for i = 1:p18CM.protocol.nPhases
%                 tPhases{i} = p18CM.state.t(cutoff_i(i):cutoff_i(i+1))...
%                     -p18CM.state.t(cutoff_i(i)) ; % set so time starts at zero for each phase
%                 YPhases{i} = p18CM.state.Y(cutoff_i(i):cutoff_i(i+1),:);
%                 currentPhases{i} = p18CM.state.currents(cutoff_i(i):cutoff_i(i+1),:);
%             end
            
            % ------------------------------------------------------------------            
            % SIMPLE FIGURES - for individual currents or state variables, any 
            % range of time
            
            % dataForPlot - choose from the following according to your needs:
            %{
            x:
            p18CM.state.t
            tPhases{i}  ***
            y:
            p18CM.state.Y
            p18CM.state.currents
            YPhases{i}  ***
            currentPhases{i}  ***
            
            *** where i is the phase number of interest
            %} 
            dataForPlot_x = p18CM.state.t ;
            dataForPlot_y = p18CM.state.Y ;
            n = length(dataForPlot_x) ;
            lower = ceil(1/n * n); % selects first index for plot
            upper = floor(n/n * n); % selects last index for plot         
                % max range is (1/n * n) to (n/n * n) but can narrow it down 
                % e.g. (0.1 * n) to (0.9 * n) plots 10th percentile index to  
                % 90th percentile index of data
            
            % Repeat & modify the below figure units as required.
            % If necessary, reset dataForPlot_x, dataForPlot_y, n, lower, and 
            % upper between figures
            
            figure
            % consult YNames or currentNames to select correct column of 
            % dataForPlot_y, ylabel, and title
            hold on
            plot(dataForPlot_x(lower:upper),dataForPlot_y(lower:upper,1),'b')
            plot(dataForPlot_x(lower:upper),dataForPlot_y(lower:upper,24)/1000,'r')
            xlabel('Time (s)')
            ylabel('Vm (V)')
            title('p18 Action potential with different pacing rates + stim')
            
%             figure
%             % consult YNames or currentNames to select correct column of 
%             % dataForPlot_y, ylabel, and title
%             plot(dataForPlot_x(lower:upper),dataForPlot_y(lower:upper,3),'b')
%             xlabel('Time (s)')
%             ylabel('[Ca]_i (mM)')
%             title('Ca-transient with different pacing rates')
%             
%             figure
%             % consult YNames or currentNames to select correct column of 
%             % dataForPlot_y, ylabel, and title
%             plot(dataForPlot_x(lower:upper),dataForPlot_y(lower:upper,2),'b')
%             xlabel('Time (s)')
%             ylabel('[Ca]_SR (mM)')
%             title('SR-Calcium with different pacing rates')
            
%             % -----------------------------------------------------------------
%             % SUPER ZOOMED OVERLAY - use when you want to compare a single
%             % AP/CaT at different pacing rates for example
%             plotend = min(p18CM.protocol.phaseLengths) - 0; % where 
%             % to zoom to, make sure this is within all phases you want to  
%             % test. Possibility to reset this for each phase with a loop if 
%             % you want a different time point for each phase
%             
%             plotstart = plotend - 10; % this plot will cover this amount of time 
%             % It helps to zoom in gradually & tweak as you go, and make sure you  
%             % don't overshoot the beginning of any of your phases
%             
%             % sometimes you'll want to nudge the trace of a spontaneously 
%             % beating CM, such as when you're comparing traces from different 
%             % pacing frequencies. Set spontoffset to do so by trial and error.
%             spontOffset = 0;
%             dataForPlot_y = YPhases ; % choose YPhases or currentPhases
%             namesForPlot_y = p18CM.YNames ; % choose YNames or currentNames
%             unitsForPlot_y = p18CM.YUnits ; % choose YUnits or currentUnits
%             % choose columnsToPlot by consulting YNames or currentNames e.g. 
%             % [1, 3, 4]. We will generate a separate figure for each.
%             columnsToPlot = 1 ; 
%             
%             phase1Spont = (p18CM.protocol.amplitudes(1) == 0) ; % boolean = 0 or 1 depending on if first phase is spontaneous
%             
%             % find indices in each phase for the specific peak you want to plot
%             plotlengthindex = cell(p18CM.protocol.nPhases, 1) ;
%             plotlengthindex_i = zeros(2, p18CM.protocol.nPhases);
%             if phase1Spont % special handling if first phase is spontaneous
%             plotlengthindex{1} = find((tPhases{1} >= plotstart-spontOffset)&(tPhases{1} <= plotend-spontOffset)) ;
%             end
%             for i = 1+phase1Spont:p18CM.protocol.nPhases
%                 plotlengthindex{i} = find((tPhases{i} >= plotstart)&(tPhases{i} <= plotend)) ;
%             end
%             for i = 1:p18CM.protocol.nPhases
%                 plotlengthindex_i(:,i) = [plotlengthindex{i}(1);plotlengthindex{i}(end)] ;
%             end
%             
%             % plot the plot(s)!
%             for j = 1:length(columnsToPlot)
%                 figure
%                 hold on
%                 legendstr = cell(1,p18CM.protocol.nPhases) ;
%                 if phase1Spont % special handling if first phase is spontaneous
%                 legendstr{1} = ['spontaneous,\newlineshifted by ', num2str(spontOffset),'s'] ;
%                 
%                 plot(tPhases{1}(plotlengthindex_i(1,1):plotlengthindex_i(2,1))+spontOffset,...
%                     dataForPlot_y{1}(plotlengthindex_i(1,1):plotlengthindex_i(2,1),columnsToPlot(j)),...
%                     'LineWidth', 1.5)
%                 end
%                 
%                 for i = 1+phase1Spont:p18CM.protocol.nPhases
%                     plot(tPhases{i}(plotlengthindex_i(1,i):plotlengthindex_i(2,i)),...
%                         dataForPlot_y{i}(plotlengthindex_i(1,i):plotlengthindex_i(2,i),columnsToPlot(j)),...
%                         'LineWidth', 1.5)
%                     legendstr{i} = [num2str(p18CM.protocol.frequencies(i)), ' Hz'];
%                 end
%                 
%                 xlabel('Time (s)')
%                 ylabel([namesForPlot_y{columnsToPlot(j)},' (', unitsForPlot_y{columnsToPlot(j)},')'])
%                 title([namesForPlot_y{columnsToPlot(j)}, ' at different rates'])
%                 legend(legendstr)
%             end
            end
        end

        % FUNCTION INFO: getInfo
        %{
        This function returns information about the Paci18 CM it takes as input.
        My intention is that the Paci18 properties be kept private unless
        obtained through this getter but for now I am having difficulty managing
        access, so the class is public.
        %}
        function info = getInfo(p18CM)
            info.geometry = p18CM.geometry ;
            info.paramNames = p18CM.paramNames ;
            info.paramUnits = p18CM.paramUnits ;
            info.YNames = p18CM.YNames ;
            info.YUnits = p18CM.YUnits ;
            info.currentNames = p18CM.currentNames ;
            info.state = p18CM.state;
        end
    end

end