function [t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance, x_names)

% Cell ID = barcode for individual cells in simulation (not experimental cell number)
persistent cellID

if isempty(cellID)
    cellID = 0;
end
cellID = cellID + 1;

load 'GA/curr_cell_protocol.mat' protocol_number isNormalized

%% Prepare cell - combination multiple protocol
% Store sequence
t_stim = {};
V_stim = {};
Cai_stim = {};
stimtimes = {};
V_stim_all = []; % full AP data for normalization
Cai_stim_all = []; % full CaT data for normalation
for i=1:length(protocol_number)
    % Initiate cell, set environment & protocol
    if i == 1 % first run
        [k19] = f_prepareCell(protocol_number(i), cellID) ; % default parameters
    else
        [k19] = f_prepareCell(protocol_number(i), cellID, init) ; % final parameters from previous protocol
    end
    saveX_conductance(k19, x_conductance) ;
    scaleConductances(k19, x_conductance, x_names); 
    scaleParameters(k19, x_conductance, x_names) ;
    %% Run simulation           
    odeSolver(k19);

    %% Extract necessary outputs
    prot_Cai_stim = k19.state.Y(:, 3); % mM
    prot_V_stim = k19.state.Y(:, 1); % mV
    prot_t_stim = k19.state.t ; % ms
    % Concatenate for multiple protocols
    Cai_stim{end + 1} = prot_Cai_stim ;
    V_stim{end + 1} = prot_V_stim ;
    t_stim{end + 1} = prot_t_stim ;
    V_stim_all = [V_stim_all; prot_V_stim] ; % full AP data for normalization
    Cai_stim_all = [Cai_stim_all; prot_Cai_stim] ; % full CaT data for normalization
    % Add extraction of final parameter values for re-initialization
    init = k19.state.Y(end, :) ;
    
    if k19.isCellAlive == 1
        prot_stimtimes = k19.protocol.last3stimTimes ; % ms
        stimtimes{end + 1} = prot_stimtimes ;
    else
        prot_stimtimes = [0, 0, 0]; % ms
        stimtimes{end + 1} = prot_stimtimes ;
    end
    
end

%% Data normalization - 0 to 1, where 0 and 1 are min and max within experiment (respectively)
if isNormalized
    norm_V_stim = (V_stim_all - min(V_stim_all)) / max(V_stim_all - min(V_stim_all)) ;
    norm_Cai_stim = (Cai_stim_all - min(Cai_stim_all)) / max(Cai_stim_all - min(Cai_stim_all)) ;

    for i=1:length(protocol_number) % redistribute back to cell (separate data from protocols
        V_stim{i} = norm_V_stim(1:length(t_stim{i})) ;
        Cai_stim{i} = norm_Cai_stim(1:length(t_stim{i})) ;

        if i~=length(protocol_number) % if we're not on last protocol
            % reset normalized vectors to start from next protocol
            idx_start = length(t_stim{i}) + 1 ;
            norm_V_stim = norm_V_stim(idx_start:end) ; 
            norm_Cai_stim = norm_Cai_stim(idx_start:end) ;
        end
    end
end
%%
return



