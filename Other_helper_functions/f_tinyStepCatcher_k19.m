function [value,isterminal,direction] = f_tinyStepCatcher_k19(~, ~, ~, ~, ~, ~, ~, ~)
 timeout = 1.5;
 time = toc<timeout;
 
 %warn = 1;
 [~,msgID] = lastwarn;
 warn = ~sum(strcmp(msgID, {'MATLAB:ode15s:IntegrationTolNotMet', 'MATLAB:singularMatrix', 'MATLAB:nearlySingularMatrix'}));
 warning('-');
 warning('off');
 
 if time == 0
     %disp('timeout!')
 end
 if warn == 0
     %disp('warnout!')
 end
 
value = time * warn;
isterminal = 1;
direction = 0;
end