
function [keepT, keepCa, CaTfeature] = CaTextract_custom(t,Cai,stimtimes) 

% % For spontaneously beating cells
if isempty(stimtimes)
    
    [peakCa,~] = max(Cai) ;
    restCa = min(Cai) ;
 
    CaAmplitude = peakCa - restCa ;
    
    below_dices = [find(Cai < restCa + 0.5*CaAmplitude);length(t)+1] ;
    above_dices = [find(Cai > restCa + 0.5*CaAmplitude);length(t)+1] ;
    
    intervals1 = find(diff(below_dices) > 5) ;
    dices_up = below_dices(intervals1) ;
    times_up = t(below_dices(intervals1)) ;
    intervals2 = find(diff(above_dices) > 5) ;
    dices_down = above_dices(intervals2) ;
    times_down = t(above_dices(intervals2)) ;
    
    
    if ( isempty(dices_up) || isempty(dices_down) )
        CaTfeature = zeros(1,4) ; %1e10*ones(1,4) ;
        keepT = [] ;
        keepCa = [] ;
        return
    end
    
    if ( length(dices_up) == 1 && length(dices_down) == 1 )
        keepT = t ;
        keepCa = Cai ;
    end
    
    if (length(dices_down) > 1 && dices_down(end) > dices_up(end))
        tempdices = dices_down(end-1):dices_up(end) ;
        [~,localmindex] = min(Cai(tempdices)) ;
        keepT = t(min(tempdices)+localmindex-1:end) ;
        keepCa = Cai(min(tempdices)+localmindex-1:end) ;
    end
    
    if (length(dices_up) > 1 && dices_up(end) > dices_down(end)) ;
        tempdices = dices_down(end):dices_up(end) ;
        [~,localmindex] = min(Cai(tempdices)) ;
        dex_end = min(tempdices) + localmindex - 1 ;
        if (length(dices_down) > 1)
            dex_begin = dices_down(end-1) ;
        else
            dex_begin = 1 ;
        end
        keepT = t(dex_begin:dex_end) ;
        keepCa = Cai(dex_begin:dex_end) ;
    end
    
% % For stimulated cells
else  
  [~,dex_t] = min(abs(t-stimtimes(end))) ;    
  keepT = t(dex_t:end) ;
  keepCa = Cai(dex_t:end) ;
end % IF, difference between stimulated & spontaneously beating cells

%%%%%%%%%%%%%%%%%%%%
% % AP duration calculation common to both stimulated & spontaneous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Caderiv = [keepCa(2:end);keepCa(end)] - keepCa ;
[peakCa,peakdex] = max(keepCa) ;
[restCa,restdex] = min(keepCa) ;
tpeak = keepT(peakdex) ;

[dCadtmax,dexmax] = max(Caderiv) ;
tinit = keepT(dexmax) ;
    
time2peak = tpeak - tinit ;
    
apd30dex = find(keepT > tpeak & keepCa < (restCa + 0.7*(peakCa - restCa))) ;
if isempty(apd30dex), apd30dex = length(keepT) ; , end
intervals30 = find(diff(apd30dex)>5) ;

if isempty(intervals30)
    CaD30 = keepT(apd30dex(1)) - tinit ;
else
    CaD30 = keepT(apd30dex(intervals30(end)+1))-tinit ;
end

apd50dex = find(keepT > tpeak & keepCa < (restCa + 0.5*(peakCa - restCa))) ;
if isempty(apd50dex), apd50dex = length(keepT) ; , end
intervals50 = find(diff(apd50dex)>5) ;

if isempty(intervals50)
    CaD50 = keepT(apd50dex(1)) - tinit ;
else
    CaD50 = keepT(apd50dex(intervals50(end)+1))-tinit ;
end

apd90dex = find(keepT > tpeak & keepCa < (restCa + 0.1*(peakCa - restCa))) ;
if isempty(apd90dex), apd90dex = length(keepT) ; , end
intervals90 = find(diff(apd90dex)>5) ;

if isempty(intervals90)
    CaD90 = keepT(apd90dex(1)) - tinit ;
else
    CaD90 = keepT(apd90dex(intervals90(end)+1))-tinit ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CaTfeature = [CaD30,CaD50,CaD90,time2peak] ;
% keepT = keepT-keepT(1); % start at 0
keepT = keepT-tinit ; % upstroke at 0
return