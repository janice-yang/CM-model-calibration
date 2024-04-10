function [new_exp_t, new_exp_wave, new_sim_t, new_sim_wave] = f_alignWaveformEnds(exp_t, exp_wave, sim_t, sim_wave)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aligns extracted waveforms (AP or CaT) from experiment and simulation

% Inputs (all column vectors)
%   exp_t: time vector from experimental data
%   exp_wave: waveform from experimental data
%   sim_t: time vector extracted from simulation
%   sim_wave:: waveform extracted from simulation

% Outputs: new vectors, aligned in time

% NOTE: simulation waveforms must start/end at same or slightly
% before/after experimental, for interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up data as matrices
new_exp = [exp_t, exp_wave] ;
new_sim = [sim_t, sim_wave] ;

% Front end: align to later start (if not already aligned)
if exp_t(1) ~= sim_t(1)
    if exp_t(1) > sim_t(1) % exp starts later
        min_sim = find(new_sim(:, 1) >= exp_t(1), 1) ;
        new_sim = new_sim(min_sim-1:end, :) ;
    else % sim starts later
        new_exp = new_exp(new_exp(:, 1) >= sim_t(1), :) ;
    end
end
    
% Tail end: chop at earlier end 
if exp_t(end) ~= sim_t(end)
    if exp_t(end) > sim_t(end) % exp ends later
        new_exp = new_exp(new_exp(:, 1) <= sim_t(end), :) ;
    else % sim ends later
        max_sim = find(new_sim(:, 1) <= exp_t(end), 1, 'last') ;
        new_sim = new_sim(1:max_sim+1, :) ;
    end
end

new_exp_t = new_exp(:, 1) ;
new_exp_wave = new_exp(:, 2) ;
new_sim_t = new_sim(:, 1) ;
new_sim_wave = new_sim(:, 2) ;

end