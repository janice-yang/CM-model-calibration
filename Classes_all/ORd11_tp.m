%{
ORd11 is a concrete class for a cardiomyocyte based on the model published by
Luo and Rudy in 1991. It inherits from the superclass AbsCM which is an abstract
class representing cardiomyocytes based on any model.
%}
classdef ORd11_tp < AbsCM
    properties 
        name = 'ORd11 CM' ; 
        celltype = 0 ;
        ODEModel = @f_ORd11_tp ;
        geometry = struct(... 
            'L', 0.01, ...      % cell length (cm)
            'r', 0.0011, ...    % radius (cm)
            'R_CG', 2, ...      % ratio between capacitive and geometric membrane areas 
            'A_geo', [], ...    % geometric area (cm^2)
            'A_cap', [], ...    % capacitive area (cm^2)
            'v_cell', [], ...   % volume of cell (µL)
            'v_myo', [], ...    % volume of myoplasmic compartment (µL)
            'v_nsr', [], ...    % volume of network SR compartment (µL)
            'v_jsr', [], ...    % volume of junctional SR compartment (µL)
            'v_ss', [], ...     % volume of subspace compartment (µL)
            'Cm', 1, ...        % capacitance (µF)
            'z_Na', 1, ...      % charge of Na
            'z_K', 1, ...       % charge of K
            'z_Ca', 2) ;        % charge of Ca
     
        % For informational queries & keeping order consistent across versions
        YNames = {...
            'Vm', ... 1 
            'Na_i', 'Na_ss', 'K_i', 'K_ss', 'Ca_i', 'Ca_ss', 'Ca_nsr', 'Ca_jsr', ... 8 
            'm', 'h_f', 'h_s', 'j', 'h_CaMKs', 'j_CaMK', 'm_l', 'h_l', 'h_lCaMK', ... 9
            'd', 'f_f', 'f_s', 'f_Caf', 'f_Cas', 'j_Ca', 'f_CaMKf', 'f_CaCaMKf', 'n' ... 9
            'a', 'i_f', 'i_s', 'a_CaMK', 'i_CaMKf', 'i_CaMKs', ... 6
            'x_rf', 'x_rs', ... 2
            'x_s1', 'x_s2', ... 2
            'x_K1', ... 1
            'J_relNP', 'J_relCaMK', ... 2
            'CaMK_trap', ... 1
            'I_stim' ... 1
            }; 

        YUnits = {...
            'mV', ... 1
            'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', ... 8
            '-', '-', '-', '-', '-', '-', '-', '-', '-', ... 9
            '-', '-', '-', '-', '-', '-', '-', '-', '-', ... 9
            '-', '-', '-', '-', '-', '-', ... 6
            '-', '-', ... 2
            '-', '-', ... 2
            '-', ... 1
            'mM ms^(-1)', 'mM ms^(-1)', ... 2
            '-', ... 1
            'A' ... 1
            };
        
        currentNames = {'I_Na', 'I_to', 'I_CaL', 'I_CaNa', 'I_CaK', ... 5
            'I_Kr', 'I_Ks', 'I_K1', 'I_NaCa', 'I_NaK', ... 5
            'I_Nab', 'I_Cab', 'I_Kb', 'I_pCa', ... 4
            'I_stim' ... 1
            };
        
        currentUnits = {'A', 'A', 'A', 'A', 'A', ... 5
            'A', 'A', 'A', 'A', 'A', ... 5
            'A', 'A', 'A', 'A', ... 4
            'A' ... 1
            }; 
        
        state = struct('init_Y', [], 'Y', [], 't', [], 'currents', []) ;
        
   
        
        paramNames = {'Gbar_Naf', 'Gbar_Nal', 'Gbar_to', 'Gbar_Kr', 'Gbar_Ks', 'Gbar_K1', 'Gbar_NaCa', 'Gbar_Kb', 'Gbar_pCa', ... Maximal conductances 
            'PR_NaK', ...
            'P_Ca', 'P_CaNa', 'P_CaK', 'P_CaCaMK', 'P_CaNaCaMK', 'P_CaKCaMK', 'P_Nab', 'P_Cab', ... Permeabilities
            'gamma_Nai', 'gamma_Nao', 'gamma_Ki', 'gamma_Ko', 'gamma_Cai', 'gamma_Cao', ... Activity coefficients
            'MgADP', 'MgATP', 'H', 'SigmaP', 'cmdnbar', 'trpnbar', 'bsrbar', 'bslbar', 'csqnbar', ... Average concentrations
            'K_mCaM', 'K_mCaMK_active', 'K_mCass', 'K_mCaact', 'K_mcmdn', 'K_mtrpn', 'K_mbsr', 'K_mbsl', 'K_mcsqn', ... Half-saturation values
            'tau_diffNa', 'tau_diffK', 'tau_diffCa', 'tau_tr', ... Time-constants
            'A_hf', 'A_hs', 'A_ff', 'A_fs', 'A_fCaMKf', 'A_fCaMKs', ... Percentages
            'alpha_CaMK', 'beta_CaMK', 'CaMK_o', ... CaMK parameters
            'k_Na1', 'k_Na2', 'k_Na3', 'k_asymm', 'omega_Na', 'omega_Ca', 'omega_NaCa', 'k_Caon', 'k_Caoff', 'q_Na', 'q_Ca', ... I_NaCa parameters
            'k_1plus', 'k_1neg', 'k_2plus', 'k_2neg', 'k_3plus', 'k_3neg', 'k_4plus', 'k_4neg', ... I_NaK parameters
            'Delta', ...  
            'K_oNai', 'K_oNao', 'K_Ki', 'K_Ko', 'K_MgATP', 'K_HP', 'K_NaP', 'K_KP', ... I_NaK dissociation constants 
            'beta_tau', 'beta_tauCaMK', 'alpha_rel', 'alpha_relCaMK', ... J_rel parameters
            'DeltaKbar_mPLB', 'DeltaJbar_upCaMK', ... J_up parameters
            };
        
        paramUnits = {'mS microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', 'microA microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', ... Maximal conductances 
            '-', ...
            'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', ... Permeabilities
            '-', '-', '-', '-', '-', '-', ... Activity coefficients
            'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', ... Average concentrations
            'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', ... Half-saturation values
            'ms', 'ms', 'ms', 'ms', ... Time-constants
            '-', '-', '-', '-', '-', '-', ... Percentages
            'ms^(-1)', 'ms^(-1)', '-', ... CaMK parameters
            'mM', 'mM', 'mM', '-', 'Hz', 'Hz', 'Hz', 'mM ms^(-1)', 'Hz', '-', '-', ... I_NaCa parameters
            'Hz', 'mM^(-1)', 'Hz', 'Hz', 'Hz', 'Hz mM^(-1)', 'Hz', 'Hz', ... I_NaK parameters
            '-', ...  
            'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', ... I_NaK dissociation constants 
            'ms', 'ms', 'ms', 'ms', ... J_rel parameters
            'mM', 'flux units?', ... J_up parameters}; 
            };
        parameters = struct(...
            'baseline', [],...
            'scaling', []);
        
        % conductances and permeabilities
        conductanceNames = {'Gbar_Naf', 'Gbar_Nal', 'Gbar_to', 'Gbar_Kr', 'Gbar_Ks', 'Gbar_K1', 'Gbar_NaCa', 'Gbar_Kb', 'Gbar_pCa', 'G_RyR', 'G_SERCA', 'G_NaK', 'G_up_epi',... Maximal conductances 
            'PR_NaK', ...
            'P_Ca', 'P_CaNa', 'P_CaK', 'P_CaCaMK', 'P_CaNaCaMK', 'P_CaKCaMK', 'P_Nab', 'P_Cab'};
        conductanceUnits = {'mS microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', 'microA microF^(-1)', 'mS microF^(-1)', 'mS microF^(-1)', '-', '-', ... Maximal conductances 
            '-', ...
            'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)', 'cm s^(-1)'};
        conductances = struct(... 
            'baseline', [],...
            'scaling', [],...
            'drugEffects', [],...
            'applicationTimes', []);
        
        % ***note: protocol is set up separately from the constructor to reduce
        % the number of arguments that the constructor must take, simplifying
        % practical use of this program
        protocol = struct('intervalTimes', [] , 'stimulus', [], 'amplitudes', [],...
            'numPulses', [], 'precedingTime', [], 'pulseDuration', [], 'totalTime', [],... 
            'phaseLengths', [], 'frequencies', []) ;
    end
    
    methods
        
        % CONSTRUCTOR:
        %{
        creates an instance of ORd11 with initial state variables equal to
        init_Y.
        %}
        function ORd11CM = ORd11_tp(init_Y, celltype) % celltype = 0 for endo, 1 for epi
            ORd11CM.celltype = celltype;
            allParameters = f_ORd11_parameters(celltype) ; % a function that returns a 
            % vector of parameters including conductances and permeabilities 
 
            % populate baseline conductances and parameters
            ORd11CM.conductances.baseline = allParameters(1:22) ;
            ORd11CM.parameters.baseline = allParameters(23:end) ;
            
            % initiate no-scaling, no-drug condition
            ORd11CM.parameters.scaling = ones(1,length(ORd11CM.parameters.baseline)) ;
            nConductances = length(ORd11CM.conductances.baseline);
            ORd11CM.conductances.scaling = ones(1,nConductances) ;
            ORd11CM.conductances.drugEffects = ones(1,nConductances) ;
            ORd11CM.conductances.applicationTimes = zeros(2,nConductances) ;
            
            % set initial state variables:
            ORd11CM.state.init_Y = init_Y ;
            
            % calculate surface areas and volumes:
            ORd11CM.geometry.A_geo = 2 * pi * ORd11CM.geometry.r^2 + 2 * pi * ORd11CM.geometry.r * ORd11CM.geometry.L;
            ORd11CM.geometry.A_cap = ORd11CM.geometry.R_CG * ORd11CM.geometry.A_geo;
            ORd11CM.geometry.v_cell = pi * ORd11CM.geometry.r^2 * ORd11CM.geometry.L;
            ORd11CM.geometry.v_myo  = 0.68 * ORd11CM.geometry.v_cell;
            ORd11CM.geometry.v_nsr  = 0.0552 * ORd11CM.geometry.v_cell;
            ORd11CM.geometry.v_jsr  = 0.0048 * ORd11CM.geometry.v_cell;
            ORd11CM.geometry.v_ss   = 0.02 * ORd11CM.geometry.v_cell;           
        end
        
        % FUNCTION INFO: scaleCell
        %{
        This function takes a set of scale factors to modify the parameters
        and saves it in the parameters struct. The default value for the
        scaling field is a vector of 1s, so this function need only be used
        when alternative scaling is required. I may modify this function to take
        in multiple sets of scaling factors to scale different parts of the
        model, e.g. conductances.
        %}
        function scaleCell(ORd11CM, conductanceScaleFactors, parameterScaleFactors)
            ORd11CM.conductances.scaling = conductanceScaleFactors ;
            ORd11CM.parameters.scaling = parameterScaleFactors ;
        end
        
        % FUNCTION INFO: setUpPacingProtocol
        %{
        This function takes in a ORd11 object as well as information about the
        pacing protocol and sets the pacing field intervalTimes to a two-column
        matrix containing the start and end times of each interval in the
        protocol, where an interval is defined as either the time preceding the
        stimulus, or the time during the stimulus.
        The pacing field stimulus is set to the amplitude of the stimulus
        applied to the ORd11CM in each of these intervals, zero before and some 
        other value during. This methodology is used because of the inconsistent 
        time steps used by MATLAB's ODE solvers. While the inconsistent time 
        steps speed up the process overall, it may overshoot the exact moment 
        each stimulus begins, which may be consequential given the short 
        stimulus durations.
        
        This function employs an external function called f_setUpPacingProtocol
        that should be included in the same folder as the ORd11 class.

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
        function setUpPacingProtocol(ORd11CM, amplitudes, numPulses, precedingTime, pulseDurations)
            [ORd11CM.protocol.intervalTimes, ORd11CM.protocol.stimulus]...
                = f_setUpPacingProtocol(amplitudes, numPulses, precedingTime, pulseDurations) ;
            
            ORd11CM.protocol.amplitudes = amplitudes ;
            ORd11CM.protocol.numPulses = numPulses ; 
            ORd11CM.protocol.precedingTime = precedingTime ; 
            ORd11CM.protocol.pulseDuration = pulseDurations ;
            
            ORd11CM.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations) ; % for reference
            ORd11CM.protocol.phaseLengths(1) = 0.2 + (numPulses(1) - 1) * precedingTime(1) + numPulses(1) * pulseDurations(1) ; % first pulse doesn't have full preceding time
            ORd11CM.protocol.phaseLengths(end) = ORd11CM.protocol.phaseLengths(end) + 1 ; % add extra second at the end to capture last AP
            ORd11CM.protocol.frequencies = 1000./(precedingTime + pulseDurations) ; % for reference
            
            ORd11CM.protocol.totalTime = ORd11CM.protocol.intervalTimes(end,2) ;
            ORd11CM.protocol.nPhases = length(amplitudes) ;
            ORd11CM.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations)/1000 ;
            
%             ORd11CM.protocol.intervalTimes = ORd11CM.protocol.intervalTimes ; % Tomek19 model uses ms not s
            ORd11CM.protocol.stimtimes = ORd11CM.protocol.intervalTimes(ORd11CM.protocol.stimulus > 0,1) ;            
        end
        
        % FUNCTION INFO: setUpDrugApplication
        %{
        This function takes in a ORd11 object as well as information about the
        drug protocol and sets the 'drugs' field 'drugEffects' to a vector of
        scaling factors that are the effects the drugs have on the different
        channel conductances in the order listed in the 'drugTargets' property.
        Drugs that are not applied should be denoted as 1, while drugs that are
        applied should be a number between 0 - 1 signifying the fraction of
        initial conductance that will remain after drug application on that
        target. The 'drugs' field 'applicationTimes' is set  as a vector of the
        corresponding time each drug is added. Any number can be used for drugs
        that are not added.
        The methodology used to ensure maximum precision in the pacing protocol 
        was not repeated for drug application as drugs are applied for a much 
        longer period of time than a current stimulus so precision to the 
        fraction of a second is not essential.
        %}
        function setUpDrugApplication(ORd11CM, drugEffects, onTimes, offTimes)
            ORd11CM.conductances.applicationTimes = [onTimes; offTimes] ;
            ORd11CM.conductances.drugEffects = drugEffects ;
        end
        
        % FUNCTION INFO: odeSolver
        %{
        This function should be used on a ORd11 cell after it has been
        instantiated and its pacing protocol and drug protocol have been set up.
        The function will run through the simulation of the ORd11 CM's behavior
        according to these protocols and will fill out the t and Y fields of the
        state property to reflect how the state variables change over time.
        %}
        function odeSolver(ORd11CM)
            if isempty(ORd11CM.protocol.intervalTimes)
                error('Error: Please set up pacing protocol first')
            else
                options = odeset('MaxStep', 1e-3, 'InitialStep', 2e-5) ; % solver settings
                
                % Allocate space initially to save time (~1220 rows per second
                % is usually enough). Extra zeros will be removed later to save
                % space.
                initialrows = ceil(ORd11CM.protocol.totalTime*1250) ;
                ORd11CM.state.Y = zeros(initialrows, length(ORd11CM.state.init_Y)) ;
                ORd11CM.state.t = zeros(initialrows, 1) ; % first row at T0, no need to update
                
                ORd11CM.state.Y(1,:) = ORd11CM.state.init_Y ; % first row, input values (generally steady state)
                input_values = ORd11CM.state.init_Y ; % initialize, changes in for-loop
                index = 2 ; % initialize, changes in for-loop
                
                for i = 1:length(ORd11CM.protocol.stimulus)
                    input_values(end) = ORd11CM.protocol.stimulus(i) ; % last input value is the stimulus amplitude
                    [t_current, Y_current] = ode15s(...
                        ORd11CM.ODEModel,...
                        ORd11CM.protocol.intervalTimes(i,:),...
                        input_values,...
                        options,...
                        ORd11CM.geometry,...
                        ORd11CM.parameters,...
                        ORd11CM.conductances,...
                        ORd11CM.celltype);
                    current_length = length(t_current)-1 ; % length of this interval in the stimulus protocol
                    
                    ORd11CM.state.t(index:index+current_length-1) = t_current(2:end) ; % save outcome of this run
                    ORd11CM.state.Y(index:index+current_length-1,:) = Y_current(2:end,:) ; % save outcome of this run
                    index = index + current_length;

                    input_values = ORd11CM.state.Y(index - 1, :) ; % set input values for next interval to end of current interval
                end
                
                % remove extra zeros to save space
                zeroInds = find(~ORd11CM.state.t) ; % returns indices of all zeroes
                firstExtraZero = zeroInds(2) ; % first zero is T0
                ORd11CM.state.t = ORd11CM.state.t(1:firstExtraZero-1) ;
                ORd11CM.state.Y = ORd11CM.state.Y(1:firstExtraZero-1,:) ;
                'ORd11_endo ODE Done!'
            end
        end
        
        % FUNCTION INFO: getCurrents
        %{
        This function should be used on a ORd11 cell after the odeSolver has
        been run. It will use the ODE model for the cell to find  & store the
        values of the different currents , and E_ions at the time points
        collected by the ODE solver.
        %}
        function getCurrents(ORd11CM)
            if isempty(ORd11CM.state.t)
                error('Error: Please run the odeSolver on this ORd11 CM first')
            else
                ORd11CM.state.currents = zeros(length(ORd11CM.state.t), length(ORd11CM.currentNames));
                numLoops = length(ORd11CM.state.t) ;
                for i= 1:numLoops
                    [~, data]    = f_ORd11_tp(ORd11CM.state.t(i), ORd11CM.state.Y(i,:), ORd11CM.geometry, ORd11CM.parameters, ORd11CM.conductances, ORd11CM_endo.celltype); % ~ says only return data, not dY output
                    ORd11CM.state.currents(i,1:end) = data;
                    %{'currents:', i/numLoops*100} % display in command window for a progress tracker, runs 0 to 100
                    
                end
                ORd11CM.state.currents(:,end) = sum(ORd11CM.state.currents(:,1:end),2); % sum for Itot; the 2 at the end is so we sum the rows not the columns
                'ORd11_endo Currents Done!'
            end
        end
        
        % FUNCTION INFO: makePlots
        %{
        This function should be used on a ORd11 cell after the odeSolver and 
        getCurrents functions have been run. It will produce plots associated 
        with the ORd11 CM. Customize according to the protocol you choose and 
        the figures you are trying to generate. This function currently includes
        potentially useful code blocks for specific types of figures which can
        be modified or commented out if not needed.
        %}
         function makePlots(ORd11CM)
%              if isempty(ORd11CM.state.t) % check that simulation has been run
%                 error('Error: Please run the odeSolver on this ORd11 CM first')
%             elseif isempty(ORd11CM.state.currents)
%                 error('Error: Please run the getCurrents function on this ORd11 CM first')
%             else 
%                 
%             % ------------------------------------------------------------------     
%             % the following block of code divides Y, t, and currents data into 
%             % pacing phases, useful for comparing different parts of a protocol.
%             % Required for 'SUPER ZOOMED OVERLAY' plot. 
%             
% %             cutoffs = cell(1, ORd11CM.protocol.pacing.nPhases); 
% %             cutoff_i = zeros(1,ORd11CM.protocol.pacing.nPhases);
% %             for i = 1:ORd11CM.protocol.pacing.nPhases
% %                 cutoffs{i} = find(ORd11CM.state.t == sum(ORd11CM.protocol.pacing.phaseLengths(1:i)));
% %                 cutoff_i(i) = cutoffs{i}(1) ;
% %             end
% %             cutoff_i = [1, cutoff_i]; % indices where each phase switches over
% %             
% %             tPhases = cell(ORd11CM.protocol.pacing.nPhases,1); % t for each of the phases
% %             YPhases = cell(ORd11CM.protocol.pacing.nPhases,1); % Y for each of the phases
% %             currentPhases = cell(ORd11CM.protocol.pacing.nPhases,1); % currents for each of the phases
% %             
% %             for i = 1:ORd11CM.protocol.pacing.nPhases
% %                 tPhases{i} = ORd11CM.state.t(cutoff_i(i):cutoff_i(i+1))...
% %                     -ORd11CM.state.t(cutoff_i(i)) ; % set so time starts at zero for each phase
% %                 YPhases{i} = ORd11CM.state.Y(cutoff_i(i):cutoff_i(i+1),:);
% %                 currentPhases{i} = ORd11CM.state.currents(cutoff_i(i):cutoff_i(i+1),:);
% %             end
%             
%             % ------------------------------------------------------------------            
%             % SIMPLE FIGURES - for individual currents or state variables, any 
%             % range of time
%             
%             % dataForPlot - choose from the following according to your needs:
%             %{
%             x:
%             ORd11CM.state.t
%             tPhases{i}  ***
%             y:
%             ORd11CM.state.Y
%             ORd11CM.state.currents
%             YPhases{i}  ***
%             currentPhases{i}  ***
%             
%             *** where i is the phase number of interest
%             %} 
%             dataForPlot_x = ORd11CM.state.t ;
%             dataForPlot_y = ORd11CM.state.Y ;
%             n = length(dataForPlot_x) ;
%             lower = ceil(1/n * n); % selects first index for plot
%             upper = floor(n/n * n); % selects last index for plot         
%                 % max range is (1/n * n) to (n/n * n) but can narrow it down 
%                 % e.g. (0.1 * n) to (0.9 * n) plots 10th percentile index to  
%                 % 90th percentile index of data
%             
%             % Repeat & modify the below figure units as required.
%             % If necessary, reset dataForPlot_x, dataForPlot_y, n, lower, and 
%             % upper between figures
%             
%             figure(1)
%             % consult YNames or currentNames to select correct column of 
%             % dataForPlot_y, ylabel, and title
%             hold on
%             plot(dataForPlot_x(lower:upper),dataForPlot_y(lower:upper,1),'b')
%             plot(dataForPlot_x(lower:upper),dataForPlot_y(lower:upper,end),'r')
%             xlabel('Time (s)')
%             ylabel('Vm (mV)')
%             title('Action potential with different pacing rates + stim')
            
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
%             plotend = min(ORd11CM.protocol.pacing.phaseLengths) - 0; % where 
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
%             namesForPlot_y = ORd11CM.YNames ; % choose YNames or currentNames
%             unitsForPlot_y = ORd11CM.YUnits ; % choose YUnits or currentUnits
%             % choose columnsToPlot by consulting YNames or currentNames e.g. 
%             % [1, 3, 4]. We will generate a separate figure for each.
%             columnsToPlot = 1 ; 
%             
%             phase1Spont = (ORd11CM.protocol.pacing.amplitudes(1) == 0) ; % boolean = 0 or 1 depending on if first phase is spontaneous
%             
%             % find indices in each phase for the specific peak you want to plot
%             plotlengthindex = cell(ORd11CM.protocol.pacing.nPhases, 1) ;
%             plotlengthindex_i = zeros(2, ORd11CM.protocol.pacing.nPhases);
%             if phase1Spont % special handling if first phase is spontaneous
%             plotlengthindex{1} = find((tPhases{1} >= plotstart-spontOffset)&(tPhases{1} <= plotend-spontOffset)) ;
%             end
%             for i = 1+phase1Spont:ORd11CM.protocol.pacing.nPhases
%                 plotlengthindex{i} = find((tPhases{i} >= plotstart)&(tPhases{i} <= plotend)) ;
%             end
%             for i = 1:ORd11CM.protocol.pacing.nPhases
%                 plotlengthindex_i(:,i) = [plotlengthindex{i}(1);plotlengthindex{i}(end)] ;
%             end
%             
%             % plot the plot(s)!
%             for j = 1:length(columnsToPlot)
%                 figure
%                 hold on
%                 legendstr = cell(1,ORd11CM.protocol.pacing.nPhases) ;
%                 if phase1Spont % special handling if first phase is spontaneous
%                 legendstr{1} = ['spontaneous,\newlineshifted by ', num2str(spontOffset),'s'] ;
%                 
%                 plot(tPhases{1}(plotlengthindex_i(1,1):plotlengthindex_i(2,1))+spontOffset,...
%                     dataForPlot_y{1}(plotlengthindex_i(1,1):plotlengthindex_i(2,1),columnsToPlot(j)),...
%                     'LineWidth', 1.5)
%                 end
%                 
%                 for i = 1+phase1Spont:ORd11CM.protocol.pacing.nPhases
%                     plot(tPhases{i}(plotlengthindex_i(1,i):plotlengthindex_i(2,i)),...
%                         dataForPlot_y{i}(plotlengthindex_i(1,i):plotlengthindex_i(2,i),columnsToPlot(j)),...
%                         'LineWidth', 1.5)
%                     legendstr{i} = [num2str(ORd11CM.protocol.pacing.frequencies(i)), ' Hz'];
%                 end
%                 
%                 xlabel('Time (s)')
%                 ylabel([namesForPlot_y{columnsToPlot(j)},' (', unitsForPlot_y{columnsToPlot(j)},')'])
%                 title([namesForPlot_y{columnsToPlot(j)}, ' at different rates'])
%                 legend(legendstr)
%             end
            %end
        end

        % FUNCTION INFO: getInfo
        %{
        This function returns information about the ORd11 CM it takes as input.
        My intention is that the ORd11 properties be kept private unless
        obtained through this getter but for now I am having difficulty managing
        access, so the class is public.
        %}
        function info = getInfo(ORd11CM)
            info.geometry = ORd11CM.geometry ;
            info.paramNames = ORd11CM.paramNames ;
            info.paramUnits = ORd11CM.paramUnits ;
            info.YNames = ORd11CM.YNames ;
            info.YUnits = ORd11CM.YUnits ;
            info.currentNames = ORd11CM.currentNames ;
            info.state = ORd11CM.state;
        end
    end

end
