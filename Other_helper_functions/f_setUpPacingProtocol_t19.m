function [intervalTimes, iStimulus] = f_setUpPacingProtocol_t19(amplitudes, numPulses, precedingTime, pulseDurations) % precedingTime & pulseDuration in ms
if (length(amplitudes) == length(numPulses) && length(amplitudes) == length(precedingTime))
    nPhases = length(amplitudes) ; % number of phases (a phase is a consistent set of conditions: amplitude, numPulses, precedingTime)
    protocol_nPulses = sum(numPulses) ; % total number of pulses
    
    intervalTimes = zeros(2*protocol_nPulses + 1, 2) ;
    iStimulus = zeros(2*protocol_nPulses+1, 1) ; % value of stimulus in each interval
    
    pulseIndex = 1 ;
    
    for i = 1:nPhases
        for j = 1:numPulses(i)
            
            if pulseIndex == 1 % can't use recursive method for first row
                intervalTimes(1, 1) = 0 ;
                intervalTimes(1, 2) = 0.001 ; % want a short delay at the beginning -> not enough for a spontaneous AP to squeeze in.
            else
                intervalTimes(2*pulseIndex-1, 1) = intervalTimes(2*(pulseIndex-1), 2) ; % first in quadrant
                intervalTimes(2*pulseIndex-1, 2) = intervalTimes(2*pulseIndex-1, 1) + precedingTime(i)/1000 ; % second in quadrant
            end
            intervalTimes(2*pulseIndex, 1) = intervalTimes(2*pulseIndex-1, 2); % third in quadrant
            intervalTimes(2*pulseIndex, 2) = intervalTimes(2*pulseIndex, 1) + pulseDurations(i)/1000 ; % fourth in quadrant
            
            iStimulus(2*pulseIndex) = amplitudes(i) ; %/Cm ;   % switch on stimulus for 'during' intervals
            pulseIndex = pulseIndex + 1 ;
        end
    end
    intervalTimes(2*pulseIndex-1, 1) = intervalTimes(2*(pulseIndex-1), 2) ; % last row
    intervalTimes(2*pulseIndex-1, 2) = intervalTimes(2*pulseIndex-1, 1) + precedingTime(end)/1000 ; % last row allows for last beat
else
    error('Error: amplitude, numPulses, and precedingTime vectors are not of equal length')
end
end



