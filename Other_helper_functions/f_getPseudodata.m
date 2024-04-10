function [experimental_dataset] = f_getPseudodata(cell_num, protocol_num, isNormalized, sigmaAP, sigmaCaT)
% Protocol numbers:
%{ 
1 - Cao High, PCL = 800ms	
2 - Cao High, spontaneous
3 - Cao Low, PCL = 800ms	
4 - Cao Low, spontaneous	
5 - ICaL Block = 25%, PCL = 800ms	
6 - ICaL Block = 25%, spontaneous	
7 - ICaL Block = 50%, PCL = 800ms	
8 - ICaL Block = 50%, spontaneous	
9 - IKr Block = 15%, PCL = 800ms	
10- IKr Block = 15%, spontaneous	
11- IKr Block = 30%, PCL = 800ms	
12- IKr Block = 30%, spontaneous	
13- Ko High, PCL = 800ms	
14- Ko High, spontaneous	
15- Ko Low, PCL = 800ms	
16- Ko Low, spontaneous	
17- PCL = 600ms	
18- PCL = 800ms	
19- Spontaneous	
%}
try
load experimental_dataset experimental_dataset
% If cell+protocols match, keep experimental_dataset; otherwise replace experimental_dataset 
if isequal(experimental_dataset{end}.protocol_num, protocol_num) && experimental_dataset{end}.cell_num == cell_num && experimental_dataset{end}.isNormalized == isNormalized && experimental_dataset{end}.sigmaAP == sigmaAP && experimental_dataset{end}.sigmaCaT == sigmaCaT 
    disp('protocol/s already loaded')
    return
else
    clearvars experimental_dataset
end
catch E_nofile 
experimental_dataset = {};
end

disp('loading new protocol/s...')
load Pseudodataset/ranges ranges

%% Read in multiple protocols
t_all = [] ;
V_stim_all = []; % full AP data for normalization
Cai_stim_all = []; % full CaT data for normalation
ends = zeros(1, length(protocol_num)) ; % for keep track of protocol indices for normalization
for i=1:length(protocol_num)
    expT = readmatrix('Pseudodataset/saved_data/pseudodataset.xlsx',...
        'Sheet',['Cell ',int2str(cell_num)],'Range','A3:A50003');
    expVCaT = readmatrix('Pseudodataset/saved_data/pseudodataset.xlsx',...
        'Sheet',['Cell ',int2str(cell_num)],'Range',ranges{protocol_num(i)});
    
    if ismember(protocol_num(i), [1,3,5,7,9,11,13,15,18,27,30]) % 800ms PCL
        stimtimes = [2896,3696,4496] ;
    elseif protocol_num(i) == 17
        stimtimes = [3496, 4096, 4696] ;
    elseif ismember(protocol_num(i), [21,23,25]) % 1000ms PCL
        stimtimes = [2496, 3496, 4496] ;
    elseif ismember(protocol_num(i), [26,29]) % 500ms PCL
        stimtimes = [3499.5, 3999.5, 4499.5] ;
    elseif ismember(protocol_num(i), [28,31]) % 2000ms PCL
        stimtimes = [496, 2496, 4496] ;
    else % spontaneous beating
        stimtimes = [] ;
    end
    [keepT, v, ca,tinit,errorcode] = waveform_extract_new(expT, expVCaT(:, 1),expVCaT(:, 2),stimtimes);
    
    header = {'Cell', 'Protocol', 'Time_AP', 'AP', 'Time_CaT', 'CaT'};
    lengthToExtract = length(keepT);
    experimental_data = {ones(lengthToExtract, 1)*cell_num, ones(lengthToExtract, 1)*protocol_num(i), keepT, v, keepT, ca};
    experimental_data = cell2table(experimental_data, 'VariableNames', header) ;
    experimental_dataset{i} = experimental_data ; 
    % Keep track of length of experiment from each protocol
    if i==1
        ends(i) = length(experimental_data.AP{1}) ;
    else
        ends(i) = length(experimental_data.AP{1}) + (ends(i-1)) ;
    end
    t_all = [t_all; keepT] ;
    V_stim_all = [V_stim_all; experimental_data.AP{1}] ; % full AP data for normalization
    Cai_stim_all = [Cai_stim_all; experimental_data.CaT{1}] ; % full CaT data for normalization
    
end

% Add Gaussian noise
rng(0)
Vscale = normrnd(0, sigmaAP, size(V_stim_all)) ;   % Generate noise scalings
% dV = [0; diff(V_stim_all)./diff(t_all)] ;
Vamp = max(V_stim_all) - min(V_stim_all) ;
dV = [diff(V_stim_all); 1] ;
V_stim_all = V_stim_all + Vamp.*Vscale ;
% V_stim_all = V_stim_all + scale.*(max(dV) - dV + 1) ;
% dCai = [0; diff(Cai_stim_all)./diff(t_all)] ;
CaTscale = normrnd(0, sigmaCaT, size(Cai_stim_all)) ;   % Generate noise scalings
CaiAmp = max(Cai_stim_all) - min(Cai_stim_all) ;
dCai = [diff(Cai_stim_all); 1] ;
Cai_stim_all = Cai_stim_all + CaiAmp.*CaTscale ;

if isNormalized
    % Normalize 0 to 1
    norm_V_stim = (V_stim_all - min(V_stim_all)) ./ max(V_stim_all - min(V_stim_all)) ;
    norm_Cai_stim = (Cai_stim_all - min(Cai_stim_all)) ./ max(Cai_stim_all - min(Cai_stim_all)) ;

    for i=1:length(protocol_num) % redistribute back to cell (separate data from protocols
        if i==1
            idx_start = 1 ;
            idx_end = ends(i) ;
        else
            idx_start = ends(i-1) + 1 ;
            idx_end = ends(i) ;
        end
        experimental_dataset{i}.AP{1} = norm_V_stim(idx_start:idx_end);
        experimental_dataset{i}.CaT{1} = norm_Cai_stim(idx_start:idx_end) ;
    end
else
    for i=1:length(protocol_num) % redistribute back to cell (separate data from protocols
        if i==1
            idx_start = 1 ;
            idx_end = ends(i) ;
        else
            idx_start = ends(i-1) + 1 ;
            idx_end = ends(i) ;
        end
        experimental_dataset{i}.AP{1} = V_stim_all(idx_start:idx_end);
        experimental_dataset{i}.CaT{1} = Cai_stim_all(idx_start:idx_end) ;
    end
end

% Add metadata  at end
metadata = table(protocol_num, cell_num, isNormalized, sigmaAP, sigmaCaT) ;
experimental_dataset{end + 1} = metadata ;
% experimental_dataset = cell2table(experimental_data, 'VariableNames', header);
save experimental_dataset experimental_dataset
