function [experimental_dataset] = f_getExperimentData(filename, sheetnames, protocol_num, isNormalized, datatype, nbeats)
%{
Extract real experimental data for GA fitting.
Inputs:
    filename (str) - name of file in ExperimentalData folder (not full path)
    sheetnames (cell of strs) - name(s) of sheet(s) to read
    protocol_num (array) - integer representing experimental protocol (see below)
    isNormalized (bool) - true if data should be normalized, false if not
    datatype (str) - AP only ('AP'), CaT only ('CaT'), or both ('APCaT')
    nbeats (array) - number of beats (individual APs and CaTs) in each protocol
Output:
    experimental_dataset - cell array of tables, to be used in GA
%}

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
...
%}
try
    load experimental_dataset experimental_dataset
    % If file+protocols match, keep experimental_dataset; otherwise replace experimental_dataset 
    if width(experimental_dataset{end}) ~= 5 && isequal(experimental_dataset{end}.protocol_num, protocol_num) && experimental_dataset{end}.filename == filename && experimental_dataset{end}.isNormalized == isNormalized && experimental_dataset{end}.datatype == datatype
        disp('protocol/s already loaded')
        return
    else
        clearvars experimental_dataset
    end
catch E_nofile 
    experimental_dataset = {};
end

disp('loading new dataset/s...')

%% Read in multiple protocols
n_filter = 5 ;
n_erode = 600 ;
n_diff = 1 ;
t_all = [] ;
V_stim_all = []; % full AP data for normalization
Cai_stim_all = []; % full CaT data for normalation
ends = zeros(1, length(protocol_num)) ; % for keep track of protocol indices for normalization
for i=1:length(protocol_num)
    datamatrix = readmatrix(['ExperimentalData/', filename], ...
        'Sheet', sheetnames{i});
    expT = datamatrix(:, 1) ;
    % stimtimes = datamatrix(:, end) ;
    % stimtimes = stimtimes(~isnan(stimtimes)) ;

    if strcmp(datatype, 'APCaT')
        expV = datamatrix(:, 2) ;
        expCaT = datamatrix(:, 3) ;
        
        % Noise processing
        erodedSignal = imerode(expV, ones(n_erode, 1)) ; 
        expV = expV - erodedSignal ;
        expV = medfilt1(expV, n_filter) ;
        erodedSignal = imerode(expCaT, ones(n_erode, 1)) ; 
        expCaT = expCaT - erodedSignal ;
        expCaT = medfilt1(expCaT, n_filter) ;

        % Find upstroke times 
        d2V = diff(expV, n_diff) ;
        [~, stim_idx] = maxk(d2V, nbeats(i)) ; 
        stim_idx = sort(stim_idx, 'ascend') ;
        stimtimes = expT(stim_idx) ;

        % Average last 5 APs and CaTs
        n = 5 ;
        last_idx = stim_idx(end-n:end) ;
        avgV = [0;0] ;
        avgCaT = [0;0] ;
        for ii=1:length(last_idx)-1
            addV = expV(last_idx(ii):last_idx(ii+1)) ;
            if length(addV) > length(avgV)
                avgV(numel(addV)) = 0 ;
            elseif length(addV) < length(avgV)
                addV(numel(avgV)) = 0 ;
            end
            avgV = avgV + addV ;

            addCaT = expCaT(last_idx(ii):last_idx(ii+1)) ;
            if length(addCaT) > length(avgCaT)
                avgCaT(numel(addCaT)) = 0 ;
            elseif length(addCaT) < length(avgCaT)
                addCaT(numel(avgCaT)) = 0 ;
            end
            avgCaT = avgCaT + addCaT ;
        end
        v = avgV ./ n ;
        ca = avgCaT ./ n ;
        expT = expT(last_idx(end)-length(v):last_idx(end)-1) ;

        keepT = expT - expT(1) ;
        % [keepT, v, ca,tinit,errorcode] = waveform_extract_new(expT, expV,expCaT,stimtimes);
        
        header = {'Filename', 'Protocol', 'Time_AP', 'AP', 'Time_CaT', 'CaT'};
        lengthToExtract = length(keepT);
        experimental_data = {repmat(filename,lengthToExtract,1), ones(lengthToExtract, 1)*protocol_num(i), keepT, v, keepT, ca};
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
    elseif strcmp(datatype, 'AP')
        expV = datamatrix(:, 2) ;

        % Noise processing
        erodedSignal = imerode(expV, ones(n_erode, 1)) ; 
        expV = expV - erodedSignal ;
        expV = medfilt1(expV, n_filter) ;
        % erodedSignal = imerode(expCaT, ones(n_erode, 1)) ; 
        % expCaT = expCaT - erodedSignal ;
        % expCaT = medfilt1(expCaT, 10) ;

        % Find upstroke times 
        d2V = diff(expV, n_diff) ;
        [~, stim_idx] = maxk(d2V, nbeats(i)) ; 
        stim_idx = sort(stim_idx, 'ascend') ;
        stimtimes = expT(stim_idx) ;

        % Average last 5 APs
        n = 5 ;
        last_idx = stim_idx(end-n:end) ;
        avgV = [0;0] ;
        % avgCaT = [0;0] ;
        for ii=1:length(last_idx)-1
            addV = expV(last_idx(ii):last_idx(ii+1)) ;
            if length(addV) > length(avgV)
                avgV(numel(addV)) = 0 ;
            elseif length(addV) < length(avgV)
                addV(numel(avgV)) = 0 ;
            end
            avgV = avgV + addV ;

            % addCaT = expCaT(last_idx(ii):last_idx(ii+1)) ;
            % if length(addCaT) > length(avgCaT)
            %     avgCaT(numel(addCaT)) = 0 ;
            % elseif length(addCaT) < length(avgCaT)
            %     addCaT(numel(avgCaT)) = 0 ;
            % end
            % avgCaT = avgCaT + addCaT ;
        end
        v = avgV ./ n ;
        % ca = avgCaT ./ n ;
        expT = expT(last_idx(end)-length(v):last_idx(end)-1) ;

        keepT = expT - expT(1) ;
        % [keepT, v, ~] = APextract_custom(expT,expV,stimtimes);
        
        header = {'Filename', 'Protocol', 'Time_AP', 'AP'};
        lengthToExtract = length(keepT);
        experimental_data = {repmat(filename,lengthToExtract,1), ones(lengthToExtract, 1)*protocol_num(i), keepT, v};
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

    elseif strcmp(datatype, 'CaT')
        expCaT = datamatrix(:, 2) ;

        % Noise processing
        % erodedSignal = imerode(expV, ones(n_erode, 1)) ; 
        % expV = expV - erodedSignal ;
        % expV = medfilt1(expV, 10) ;
        erodedSignal = imerode(expCaT, ones(n_erode, 1)) ; 
        expCaT = expCaT - erodedSignal ;
        expCaT = medfilt1(expCaT, n_filter) ;

        % Find upstroke times 
        d2CaT = diff(expCaT, n_diff) ;
        [~, stim_idx] = maxk(d2CaT, nbeats(i)) ; 
        stim_idx = sort(stim_idx, 'ascend') ;
        stimtimes = expT(stim_idx) ;

        % Average last 5 APs
        n = 5 ;
        last_idx = stim_idx(end-n:end) ;
        % avgV = [0;0] ;
        avgCaT = [0;0] ;
        for ii=1:length(last_idx)-1
            % addV = expV(last_idx(ii):last_idx(ii+1)) ;
            % if length(addV) > length(avgV)
            %     avgV(numel(addV)) = 0 ;
            % elseif length(addV) < length(avgV)
            %     addV(numel(avgV)) = 0 ;
            % end
            % avgV = avgV + addV ;

            addCaT = expCaT(last_idx(ii):last_idx(ii+1)) ;
            if length(addCaT) > length(avgCaT)
                avgCaT(numel(addCaT)) = 0 ;
            elseif length(addCaT) < length(avgCaT)
                addCaT(numel(avgCaT)) = 0 ;
            end
            avgCaT = avgCaT + addCaT ;
        end
        % v = avgV ./ n ;
        ca = avgCaT ./ n ;
        expT = expT(last_idx(end)-length(ca):last_idx(end)-1) ;

        keepT = expT - expT(1) ;
        % [keepT, ca, ~] = CaTextract_custom(expT,expCaT,stimtimes);
        
        header = {'Filename', 'Protocol', 'Time_CaT', 'CaT'};
        lengthToExtract = length(keepT);
        experimental_data = {repmat(filename,lengthToExtract,1), ones(lengthToExtract, 1)*protocol_num(i), keepT, ca};
        experimental_data = cell2table(experimental_data, 'VariableNames', header) ;
        experimental_dataset{i} = experimental_data ; 
        % Keep track of length of experiment from each protocol
        if i==1
            ends(i) = length(experimental_data.CaT{1}) ;
        else
            ends(i) = length(experimental_data.CaT{1}) + (ends(i-1)) ;
        end
        t_all = [t_all; keepT] ;
        Cai_stim_all = [Cai_stim_all; experimental_data.CaT{1}] ; % full CaT data for normalization

    else
        error("Please enter 'AP', 'CaT', or 'APCaT' as datatype.")
    end

end

% % Add Gaussian noise
% rng(0)
% Vscale = normrnd(0, sigmaAP, size(V_stim_all)) ;   % Generate noise scalings
% % dV = [0; diff(V_stim_all)./diff(t_all)] ;
% Vamp = max(V_stim_all) - min(V_stim_all) ;
% dV = [diff(V_stim_all); 1] ;
% V_stim_all = V_stim_all + Vamp.*Vscale ;
% % V_stim_all = V_stim_all + scale.*(max(dV) - dV + 1) ;
% % dCai = [0; diff(Cai_stim_all)./diff(t_all)] ;
% CaTscale = normrnd(0, sigmaCaT, size(Cai_stim_all)) ;   % Generate noise scalings
% CaiAmp = max(Cai_stim_all) - min(Cai_stim_all) ;
% dCai = [diff(Cai_stim_all); 1] ;
% Cai_stim_all = Cai_stim_all + CaiAmp.*CaTscale ;

if isNormalized
    % Normalize 0 to 1
    if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
        norm_V_stim = (V_stim_all - min(V_stim_all)) ./ max(V_stim_all - min(V_stim_all)) ;
    end
    if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
        norm_Cai_stim = (Cai_stim_all - min(Cai_stim_all)) ./ max(Cai_stim_all - min(Cai_stim_all)) ;
    end
    
    for i=1:length(protocol_num) % redistribute back to cell (separate data from protocols
        if i==1
            idx_start = 1 ;
            idx_end = ends(i) ;
        else
            idx_start = ends(i-1) + 1 ;
            idx_end = ends(i) ;
        end
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            experimental_dataset{i}.AP{1} = norm_V_stim(idx_start:idx_end);
        end
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            experimental_dataset{i}.CaT{1} = norm_Cai_stim(idx_start:idx_end) ;
        end
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
        if strcmp(datatype, 'AP') || strcmp(datatype, 'APCaT')
            experimental_dataset{i}.AP{1} = V_stim_all(idx_start:idx_end);
        end
        if strcmp(datatype, 'CaT') || strcmp(datatype, 'APCaT')
            experimental_dataset{i}.CaT{1} = Cai_stim_all(idx_start:idx_end) ;
        end
    end
end

% Add metadata  at end
metadata = table({filename}, protocol_num, isNormalized, {datatype}) ;
experimental_dataset{end + 1} = metadata ;
% experimental_dataset = cell2table(experimental_data, 'VariableNames', header);
save experimental_dataset experimental_dataset
end
