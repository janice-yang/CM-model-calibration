function f_primeCM(CM, freq, seconds)

isV = zeros(1,length(CM.YNames)) ; 
for i = 1:length(isV)
    isV(i) = isequal(CM.YNames{i},'Vm') || isequal(CM.YNames{i},'V') || isequal(CM.YNames{i},'V0');
end
v_index = find(isV) ;
if length(v_index) ~= 1
    error('Error: Prime manually as membrane voltage is not named Vm or V or V0')
end


setUpPacingProtocol(CM, 20, ceil(seconds*freq), (1000/freq)-2, 2) ;
odeSolver(CM) ;

figure
hold on
plot(CM.state.t(end-10000:end), CM.state.Y(end-10000:end,v_index(1)))

if abs(CM.state.init_Y(v_index)) < 1
plot(CM.state.t(end-10000:end), CM.state.Y(end-10000:end,end)/1000) ;
elseif abs(CM.state.init_Y(v_index)) >= 1
plot(CM.state.t(end-10000:end), CM.state.Y(end-10000:end,end)) ;
end

last_AP_t = find(CM.state.t >= CM.state.t(end)-2) ; 
last_AP_Y = CM.state.Y(last_AP_t(1):end,:) ; 
last_AP_v = CM.state.Y(last_AP_t(1):end,v_index(1)) ;
[~,I]= min(last_AP_v) ;
new_init = last_AP_Y(I,:) ;

CM.state.init_Y = new_init ;

% clear priming steps
CM.state.Y = [] ; 
CM.state.t = [] ; 
CM.protocol.intervalTimes = [] ;
CM.protocol.stimulus = [] ;
CM.protocol.amplitudes = [] ;            
CM.protocol.numPulses = [] ; 
CM.protocol.precedingTime = [] ;
CM.protocol.pulseDuration = [] ; 
CM.protocol.totalTime = [] ; 
CM.protocol.phaseLengths = [] ; 
CM.protocol.frequencies = [] ;
"Primed!"
end