
function [keepT,keepV,errorcode] = APextract_new(t,V,stimtimes) 

%%%% When there's an error, the function still needs to return waveforms
%%%% When error, return last 2 seconds of the waveform (defined by durtaion_default)
%%%% Also return the type of error 
errorcode = 0 ;

duration_default = 2000 ;
if (max(t) > duration_default)
  [~,min_t_dex] = min(t - duration_default) ;
  keepT = t(min_t_dex:end) ;
  keepV = V(min_t_dex:end) ;
else
  keepT = t ;
  keepV = V ;
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
    
    if ( length(dices_up) == 1 && length(dices_down) == 1 )
        keepT = t ;
        keepV = V ;
        errorcode = 2 ;
        return
    end
    
    lastup = (max(dices_up) > max(dices_down)) ;

    if (lastup && length(dices_down) >= 2)
      dexrange = dices_down(end):dices_up(end) ;
      [~,dex_end] = min(V(dexrange)) ;
      dex_end = dex_end + min(dexrange) ;
      
      dexrange = dices_down(end-1):dices_up(end-1) ;
      [~,dex_begin] = min(V(dexrange)) ;
      dex_begin = dex_begin + min(dexrange) ;

      keepT = t(dex_begin:dex_end) ;
      keepV = V(dex_begin:dex_end) ;
      errorcode = 0 ;
    elseif (~lastup && length(dices_down) >= 3)
      dexrange = dices_down(end-1):dices_up(end) ;
      [~,dex_end] = min(V(dexrange)) ;
      dex_end = dex_end + min(dexrange) ;
      
      dexrange = dices_down(end-2):dices_up(end-1) ;
      [~,dex_begin] = min(V(dexrange)) ;
      dex_begin = dex_begin + min(dexrange) ;

      keepT = t(dex_begin:dex_end) ;
      keepV = V(dex_begin:dex_end) ;
      errorcode = 0 ;
    end

    
% % For stimulated cells
else  
  [~,dex_t] = min(abs(t-stimtimes(end))) ;    
  keepT = t(dex_t:end) ;
  keepV = V(dex_t:end) ;
end % IF, difference between stimulated & spontaneously beating cells

%%%%%%%%%%%%%%%%%%%%
% % AP duration calculation common to both stimulated & spontaneous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vderiv = [keepV(2:end);keepV(end)] - keepV ;
[peakV,peakdex] = max(keepV) ;
[restV,restdex] = min(keepV) ;
tpeak = keepT(peakdex) ;

[dVdtmax,dexmax] = max(Vderiv) ;
tinit = keepT(dexmax) ;
    
    
apd30dex = find(keepT > tpeak & keepV < (restV + 0.7*(peakV - restV))) ;
if isempty(apd30dex), apd30dex = length(keepT) ; , end
intervals30 = find(diff(apd30dex)>5) ;

if isempty(intervals30)
    APD30 = keepT(apd30dex(1)) - tinit ;
else
    APD30 = keepT(apd30dex(intervals30(end)+1))-tinit ;
end

apd50dex = find(keepT > tpeak & keepV < (restV + 0.5*(peakV - restV))) ;
if isempty(apd50dex), apd50dex = length(keepT) ; , end
intervals50 = find(diff(apd50dex)>5) ;

if isempty(intervals50)
    APD50 = keepT(apd50dex(1)) - tinit ;
else
    APD50 = keepT(apd50dex(intervals50(end)+1))-tinit ;
end

apd90dex = find(keepT > tpeak & keepV < (restV + 0.1*(peakV - restV))) ;
if isempty(apd90dex), apd90dex = length(keepT) ;  end
intervals90 = find(diff(apd90dex)>5) ;

if isempty(intervals90)
    APD90 = keepT(apd90dex(1)) - tinit ;
else
    APD90 = keepT(apd90dex(intervals90(end)+1))-tinit ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(stimtimes)
  APfeature = [APD30,APD50,APD90,APfrequency] ;
else
  APfeature = [APD30,APD50,APD90] ;
end

keepT = keepT - tinit ;

return