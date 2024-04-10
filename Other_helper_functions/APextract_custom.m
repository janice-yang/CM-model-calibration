
function [keepT, keepV, APfeature] = APextract_custom(t,V,stimtimes) 

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
        APfeature = zeros(1,4) ; %1e10*ones(1,4) ;
        keepT = [] ;
        keepV = [] ;
        return
    end
    
    if ( length(dices_up) == 1 && length(dices_down) == 1 )
        keepT = t ;
        keepV = V ;
    end
    
    if (length(dices_down) > 1 && dices_down(end) > dices_up(end))
        tempdices = dices_down(end-1):dices_up(end) ;
        [~,localmindex] = min(V(tempdices)) ;
        keepT = t(min(tempdices)+localmindex-1:end) ;
        keepV = V(min(tempdices)+localmindex-1:end) ;
    end
    
    if (length(dices_up) > 1 && dices_up(end) > dices_down(end)) 
        tempdices = dices_down(end):dices_up(end) ;
        [~,localmindex] = min(V(tempdices)) ;
        dex_end = min(tempdices) + localmindex - 1 ;
        if (length(dices_down) > 1)
            dex_begin = dices_down(end-1) ;
        else
            dex_begin = 1 ;
        end
        keepT = t(dex_begin:dex_end) ;
        keepV = V(dex_begin:dex_end) ;
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
if isempty(apd30dex), apd30dex = length(keepT) ;  end
intervals30 = find(diff(apd30dex)>5) ;

if isempty(intervals30)
    APD30 = keepT(apd30dex(1)) - tinit ;
else
    APD30 = keepT(apd30dex(intervals30(end)+1))-tinit ;
end

apd50dex = find(keepT > tpeak & keepV < (restV + 0.5*(peakV - restV))) ;
if isempty(apd50dex), apd50dex = length(keepT) ;  end
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

if APD90 < 0 || APD50 < 0 || APD30 < 0
    APD30 = 0;
    APD50 = 0;
    APD90 = 0;
    APfrequency = 0;
end

if isempty(stimtimes)
  APfeature = [APD30,APD50,APD90,APfrequency] ;
else
  APfeature = [APD30,APD50,APD90, NaN] ;
end

% keepT = keepT-keepT(1); % start at 0
keepT = keepT-tinit ; % upstroke at 0
return