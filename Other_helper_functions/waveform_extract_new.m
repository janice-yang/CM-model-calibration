
function [keepT,keepV,keepCai,tinit,errorcode] = waveform_extract_new(t,V,Cai,stimtimes,n) 

%%%% When there's an error, the function still needs to return waveforms
%%%% When error, return last 2 seconds of the waveform (defined by duration_default)
%%%% Also return the type of error 
errorcode = 0 ;

tinit = 0 ;

i_keep = 10 ; % number of points before upstroke to keep

duration_default = 2000 ;
if (max(t) > duration_default)
  [~,min_t_dex] = min(t - duration_default) ;
  keepT = t(min_t_dex:end) ;
  keepV = V(min_t_dex:end) ;
  keepCai = Cai(min_t_dex:end) ;
else
  keepT = t ;
  keepV = V ;
  keepCai = Cai ;
end


% % For spontaneously beating cells
if isempty(stimtimes)
    
    [peakV,~] = max(V) ;
    restV = min(V) ;
    Vamplitude = peakV - restV ;
    
    below_dices = [find(V < restV + 0.5*Vamplitude);length(t)+1] ;
    above_dices = [find(V > restV + 0.5*Vamplitude);length(t)+1] ;
    
    intervals1 = find(diff(below_dices) > 5) ;
    dices_up = below_dices(intervals1) ;
    times_up = t(below_dices(intervals1)) ;
    intervals2 = find(diff(above_dices) > 5) ;
    dices_down = above_dices(intervals2) ;
    times_down = t(above_dices(intervals2)) ;
    
    cyclelength = mean(diff(times_up)) ;
    APfrequency = 1000/cyclelength ;

    if ( isempty(dices_up) || isempty(dices_down) )
        errorcode = 1 ;
        return
    end
    
    if ( length(dices_up) < n && length(dices_down) < n ) % check if there are enough beats to extract
        keepT = t ;
        keepV = V ;
        keepCai = Cai ;
        errorcode = 2 ;
        return
    end
    
    lastup = (max(dices_up) > max(dices_down)) ;
    if (lastup && length(dices_down) >= n+1) % 2 for single AP, 3 for double
      dexrange = dices_down(end):dices_up(end) ;
      [~,dex_end] = min(V(dexrange)) ;
      dex_end = dex_end + min(dexrange) ;
      
      % dexrange = dices_down(end-1):dices_up(end-1) ;
      dexrange = dices_down(end-n):dices_up(end-n) ;
      [~,dex_begin] = min(V(dexrange)) ;
      dex_begin = dex_begin + min(dexrange) ;

      keepT = t(dex_begin-i_keep:dex_end) ;
      keepV = V(dex_begin-i_keep:dex_end) ;
      keepCai = Cai(dex_begin-i_keep:dex_end) ;
      errorcode = 0 ;
    elseif (~lastup && length(dices_down) >= n+2) % 3 for single AP, 4 for double
      dexrange = dices_down(end-1):dices_up(end) ;
      [~,dex_end] = min(V(dexrange)) ;
      dex_end = dex_end + min(dexrange) ;
      
      % dexrange = dices_down(end-2):dices_up(end-1) ;
      dexrange = dices_down(end-(n+1)):dices_up(end-n) ;
      [~,dex_begin] = min(V(dexrange)) ;
      dex_begin = dex_begin + min(dexrange) ;

      keepT = t(dex_begin-i_keep:dex_end) ;
      keepV = V(dex_begin-i_keep:dex_end) ;
      keepCai = Cai(dex_begin-i_keep:dex_end) ;
      errorcode = 0 ;
    end
    
% % For stimulated cells
else  
%   [~,dex_t] = min(abs(t-stimtimes(end))) ;   % last AP 
%   keepT = t(dex_t:end) ;
%   keepV = V(dex_t:end) ;
%   keepCai = Cai(dex_t:end) ;
%   stimtimes = stimtimes / 1000 ; % seconds
%   [~,dex_t] = min(abs((t+295000)-stimtimes(end-1))) ;    % 2nd to last AP
%   [~,dex_t_end] = min(abs((t+295000)-stimtimes(end))) ; % end for single AP extraction   
  % [~,dex_t] = min(abs(t-stimtimes(end-1))) ;    % 2nd to last AP
  [~,dex_t] = min(abs(t-stimtimes(end-n))) ;    % n-th to last AP
  [~,dex_t_end] = min(abs(t-stimtimes(end))) ; % end for single AP extraction   
  keepT = t(dex_t-i_keep:dex_t_end-i_keep) ;
  keepV = V(dex_t-i_keep:dex_t_end-i_keep) ;
  keepCai = Cai(dex_t-i_keep:dex_t_end-i_keep) ;
end % IF, difference between stimulated & spontaneously beating cells

%%%%%%%%%%%%%%%%%%%%
% % AP duration calculation common to both stimulated & spontaneous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vderiv = [keepV(2:end);keepV(end)] - keepV ;
% [~,dexmax] = max(Vderiv) ;
[~,dexmax] = maxk(Vderiv,5) ;
dexmax = min(dexmax) ;
tinit = keepT(dexmax) ;
% tinit = keepT(1) ;
    
keepT = keepT - tinit ;

return
