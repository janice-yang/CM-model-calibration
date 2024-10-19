clear
close all

% Extracting voltage & Ca imaging data for f_getExperimentalData
output_file = 'DMG242.xlsx' ;
% sheetname = 'DMG242PF_Dofetilide50nM_0.5Hz' ;
sheetname = 'Nifedipine300nM_1Hz' ;
headers = {'t', 'AP', 'CaT'} ;                      
t_end = 11000 ;                                     % ms
data = zeros(1000000, length(headers)) ; 

% Get AP data
file = uigetfile_n_dir([pwd, '/ExperimentalData'] , 'Select file(s) containing AP data (ch2):') ;
AP = imread(file{1}) ;
AP_med = mean(AP, 2) ;
data(1:length(AP_med), 1) = linspace(0, t_end, length(AP_med)) ;    % time
data(1:length(AP_med), 2) = AP_med ;                                % AP


% Get CaT data
file = uigetfile_n_dir([pwd, '/ExperimentalData'], 'Select file(s) containing CaT data (ch1):') ;
CaT = imread(file{1}) ;
CaT_med = mean(CaT, 2) ;
data(1:length(CaT_med), 3) = CaT_med ;                              % CaT

% Remove trailing zeros
i_end = find(data(:,1), 1, 'last') ;
data = data(1:i_end, :) ;

writematrix(data, ['ExperimentalData/', output_file], 'Sheet', sheetname) ;

% Manually fill in headers & stimV times :(

figure
subplot(2,2,1)
image(AP', 'CDataMapping', 'scaled')
yticks([])
yticklabels([])
ylabel('28 um')
xticks([])
xticklabels([])
xlabel('Time')
subplot(2,2,2)
image(CaT', 'CDataMapping', 'scaled')
yticks([])
yticklabels([])
ylabel('28 um')
xticks([])
xticklabels([])
xlabel('Time')
subplot(2,2,3)
plot(data(:,1), data(:,2), 'k-')
ylabel('Fluor intensity')
xlabel('Time (ms)')
title('AP')
xlim([0 t_end])
subplot(2,2,4)
plot(data(:,1), data(:,3), 'k-')
ylabel('Fluor intensity')
xlabel('Time (ms)')
title('CaT')
xlim([0 t_end])
sgtitle(sheetname)

% tbl = table(data, 'VariableNames', headers) ;
% writetable(tbl, output_file, 'Sheet', sheetname) ;

%% Mean filter test (post-erosion and last 5 beats tests)
n_filters = 5 ;
% n_erodes = 600:100:1000 ;
% n_filter = 5 ;
n_erode = 1000 ;
p = 1 ; % Subplot positioning

nbeats = 11 ;
n_diff = 1 ;

figure
set(gcf, 'Position', [0 0 500 length(n_filters)*500])
hold on
for i=1:length(n_filters)
    expT = data(:, 1) ;
    expV = data(:, 2) ;
    expCaT = data(:, 3) ;
    n_filter = n_filters(i) ;

    % Find upstroke times 
    d2V = diff(expV, n_diff) ;
    [~, stim_idx] = maxk(d2V, nbeats) ; 
    stim_idx = sort(stim_idx, 'ascend') ;
    stimtimes = expT(stim_idx) ;
    % Last 5
    n = 5 ;
    start_idx = stim_idx(end-5) ;
    expT = expT(start_idx:end) ;
    expT = expT - expT(1) ;
    expV = expV(start_idx:end) ;
    expCaT = expCaT(start_idx:end) ; 


    % Noise processing
    erodedSignal = imerode(expV, ones(n_erode, 1)) ; 
    expV = expV - erodedSignal ;
    % expV = medfilt1(expV, n_filter) ;
    erodedSignal = imerode(expCaT, ones(n_erode, 1)) ; 
    expCaT = expCaT - erodedSignal ;
    % expCaT = medfilt1(expCaT, n_filter) ;

    % Trace averaging
    % Find upstroke times again
    d2V = diff(expV, n_diff) ;
    [~, stim_idx] = maxk(d2V, 5) ; 
    stim_idx = sort(stim_idx, 'ascend') ;
    stimtimes = expT(stim_idx) ;

    % Average last 5 APs and CaTs
    n = 5 ;
    last_idx = stim_idx ;
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
    % expT = expT(last_idx(end)-length(v):last_idx(end)-1) ;
    expT = linspace(0, expT(last_idx(end)) - expT(last_idx(end-1)), length(v)) ;
    % expT = linspace(0, 1000, length(v)) ;

    keepT = expT - expT(1) ;

    expT = keepT ;
    expV = v ;
    expCaT = ca ;

    % % Noise processing
    % erodedSignal = imerode(expV, ones(n_erode, 1)) ; 
    % expV = expV - erodedSignal ;
    expV = medfilt1(expV, n_filter) ;
    % erodedSignal = imerode(expCaT, ones(n_erode, 1)) ; 
    % expCaT = expCaT - erodedSignal ;
    expCaT = medfilt1(expCaT, n_filter) ;

    subplot(length(n_filters), 2, p)
    hold on
    % plot(data(:,1), data(:, 2), ...
    plot(expT, expV)
    xlabel('Time (ms)')
    title(['AP, filter size ' num2str(n_filter)])
    % title('AP')
    
    subplot(length(n_filters), 2, p+1)
    hold on
    % plot(data(:,1), data(:, 3), ...
    plot(expT, expCaT)        
    xlabel('Time (ms)')
    title(['CaT, filter size ' num2str(n_filter)])
    % title('CaT')
    p = p + 2 ;
end

%% Normalization

normV = (expV - min(expV)) ./ max(expV - min(expV)) ;
normCaT = (expCaT - min(expCaT)) ./ max(expCaT - min(expCaT)) ;