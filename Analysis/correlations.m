%% Getting compiled parameter stats
realdata = false ;
protocol_name = 'p5_APCaT_norm' ; % for output file name
base = 'Analysis/AggregateStats' ; % folder for storing results

load GA/x_names.mat

folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select experiment DIRECTORY/IES') ;
n = 10 ; % runs per experiment

% Get current date & time for filename
format shortg
date_time = char(strjoin(string(clock), '_')) ;

set(groot, 'defaultFigureRenderer', 'painters') % for saving as svg

all_errors = zeros(n*length(folders), length(names)) ; % store data
MAEs = zeros(length(folders), length(names)) ;
SDs = zeros(length(folders), length(names)) ;
ranges = zeros(length(folders), length(names)) ;

cell_numbers = zeros(n*length(folders), 1) ; % store cell order

% Get metrics for each experiment
idx_data = 1 ;
for i=1:length(folders)
    folder = folders{i} ;
    load([folder, '/Run_0/Details.mat'], 'cell_number')
    cell_numbers(idx_data:idx_data+n-1) = ones(n, 1)*cell_number ;
    
    % Get stats file
    load([folder, '/stats_data.mat'], 'log_norm_params', 'log_abs_errors')

    % Average error for each parameter
    all_errors(idx_data:idx_data+n-1, :) = log_abs_errors ;
    MAEs(i, :) = mean(log_abs_errors, 1) ; 
    SDs(i, :) = std(log_norm_params, 0, 1) ;
    ranges(i, :) = range(log_norm_params, 1) ;

    idx_data = idx_data + n ;
    
end

% Save data
save([base, '/', protocol_name, '_data_parameters.mat'], ...
    'all_errors', 'MAEs', 'SDs', 'ranges', 'protocol_name', 'cell_numbers')

%% Correlation: IKr block threshold accuracy vs num protocol conditions?
clear

IKr_files = uigetfile_n_dir([pwd, '/Analysis/AggregateStats'], 'Select file(s) containing current contributions data') ;
n_conditions = [2, 1, 3, 2, 2, 1, 1] ;
n_cells = 4 ; 
n_runs = 10 ;
metric_name = 'SD' ;

IKr_thresh_acc = zeros(n_cells*n_runs, length(n_conditions)) ;
labels = cell(1, length(IKr_files)) ;
colors = repmat('krgbmc', 1, 300) ;

for i = 1:length(IKr_files)
    % Get IKr_thresh_acc for each protocol
    load(IKr_files{i})
    idx = find(strcmp(metrics, metric_name)) ;
    IKr_thresh_acc(:, i) = plot_data(:, idx) ; 

    foldersplit = strsplit(IKr_files{i}, '/') ;
    labels{i} = foldersplit{end} ;
end

[rho, pval] = corrcoef(n_conditions, mean(IKr_thresh_acc)) ;

figure
hold on
for i=1:length(n_conditions)
    plot(ones(1,length(IKr_thresh_acc(:,i)))*n_conditions(i), IKr_thresh_acc(:,i), ...
        [colors(i), 'o']) % also the trend line?
end
legend(labels)
xlim([0, max(n_conditions)+1])
ylim([0, 1])
xlabel("Number of conditions")
if strcmp(metric_name, 'Validation error')
    ylabel("IKr threshold error")
elseif strcmp(metric_name, 'Absolute error')
    ylabel("Calibration error")
elseif strcmp(metric_name, 'SD')
    ylabel("Calibration SD")
end
title(['Rho = ', num2str(rho(1,2), 3), ', p = ', num2str(pval(1,2), 3)])

%% Predicted IKr block vs true IKr block for all 4 cells
clear
close all

metric_name = 'Validation error' ;
protocol_name = 'p4-19-5_APCaT_norm' ;
cell_numbers = [12,15,16,4] ;
n_runs = 10 ; % number of GA runs per experiment

IKr_files = uigetfile_n_dir([pwd, '/Analysis'], 'Select file(s) containing predicted IKr thresholds data') ;

true_thresholds = zeros(1, length(cell_numbers)) ;
predicted_thresholds = zeros(n_runs, length(cell_numbers)) ;
labels = strsplit(num2str(cell_numbers)) ;

% Get true thresholds
actual_IKr_thresholds = readtable('Pseudodataset/saved_data/IKr block thresholds.xlsx') ;
for i=1:length(cell_numbers)
    true_thresholds(i) = actual_IKr_thresholds.Threshold(cell_numbers(i)) ;
end

% Get predicted thresholds 
for i=1:length(IKr_files)
    load(IKr_files{i}, 'predictedThreshold')
    predicted_thresholds(:,i) = predictedThreshold ;
end

figure
hold on
plot(true_thresholds, true_thresholds, 'k--')
colors = repmat('krgbmc', 1, 300) ;
for i=1:length(labels)
    plot(true_thresholds(i), predicted_thresholds(:,i), [colors(i), 'o'])
end
xlabel("True I_{Kr} block threshold")
ylabel("Predicted I_{Kr} block threshold")
ylim([0, 1])

% Flatten data to calculate correlations/trend line
y = reshape(predicted_thresholds, n_runs*length(cell_numbers), 1) ;
x = ones(length(y),1) ;
counter = 1 ;
for i=1:length(cell_numbers)
    x(counter:counter+n_runs-1) = true_thresholds(i)*ones(n_runs, 1) ;
end

[rho, pval] = corrcoef(x, y) ;
title(['Rho = ', num2str(rho(1,2), 3), ', p = ', num2str(pval(1,2), 3)])

Fit = polyfit(x,y,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
plot(x, polyval(Fit,x), 'r--')

%% Parameter sensitivity APD90/CaTA vs parameter error from calibrations

clear
close all

addpath('sensitivity simulations/')

feature_idx = 1 ; % 1 = APD90; 4 = CaTA
cells = [12, 15, 16, 4] ;
protocol_name = 'p4-19-5_norm' ;
SD_file = uigetfile_n_dir([pwd, '/Analysis/AggregateStats'], 'Select file containing stats data') ;

% Sensitivity analysis coefficients
% load('sensitivity simulations/kernik 2019/MLR/APD90 and CaTA/sigma0.2_pop1000.mat')
load('sensitivity simulations/kernik 2019/APD90etc_PLS_1000.mat')
load GA/x_names.mat names

X = log(population) ;
Y = [log(APD90),Vrest,Vpeak, CaTA] ;
Y(isnan(Y)) = 0 ;
% Y = [log(APD90),Vrest,Vpeak, CaTA] ;
[n_samples,n_parameters] = size(X) ;
[n_samples_y,n_outputs] = size(Y) ;

if (n_samples ~= n_samples_y)
  disp('WARNING:  Input & Output matrices appear to be incompatible')
end

logtransform = [1,0,0,0] ;

outputlabels = { ...
    'APD90', ...
    'V_r_e_s_t', ...
    'V_p_e_a_k', ...
    'CaTA', ...
} ;

X_Z = zscore(X) ;
Y_Z = zscore(Y) ;

B = (X_Z'*X_Z)^(-1)*X_Z'*Y_Z ;

Yhat_Z = X_Z*B ;
Yhat = Yhat_Z.*(ones(n_samples,1)*std(Y)) + ones(n_samples,1)*mean(Y) ;

SSYT = sum((Y-ones(n_samples,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(n_samples,1)*mean(Y)).^2);
R2each = SSYR./SSYT; 

%%% PLS
nfactor = length(names); 
[Tpls,P,W,Wstar,U,b,C,B_pls,...
        Bpls_star,Xori_rec,Yori_rec,...
        R2_X,R2_Y]=PLS_nipals(X_Z,Y_Z,nfactor) ;

labels = strsplit(num2str(cells)) ;
dx = 0.01 ;  % for plot text
dy = 0.01 ;  % for plot text

load(SD_file{1})
x_errors_all = zeros(1, length(cells)*length(names)) ;
y_sum_all = zeros(1, length(cells)*length(names)) ;
names_all = cell(1, length(cells)*length(names)) ;

% Get parameter_SDs for each cell 
counter = 1 ;
for i=1:length(cells)
    x_error = MAEs(i, :) ;
    y = abs(B_pls(:, feature_idx)) ; 

    x_errors_all(counter:counter+length(names)-1) = x_error ;
    y_sum_all(counter:counter+length(names)-1) = y' ;
    names_all(counter:counter+length(names)-1) = names ;
    counter = counter + length(names) ;
end

[rho, pval] = corrcoef(x_errors_all, y_sum_all) ;

colors = repmat('krgbmc', 1, 300) ;
figure
hold on
counter = 1 ;
for i = 1:length(cells)
    plot(x_errors_all(counter:counter+length(names)-1), ...
        y_sum_all(counter:counter+length(names)-1), ...
        [colors(i), 'o']) % also the trend line?
    text(x_errors_all(counter:counter+length(names)-1) + dx, ...
        y_sum_all(:,counter:counter+length(names)-1) + dy, ...
        names_all(counter:counter+length(names)-1))
    counter = counter + length(names) ;
end

xlabel("Parameter error")
legend(labels, 'Location', 'best')
title(['All cells: Rho = ', num2str(rho(1,2), 3), ', p = ', num2str(pval(1,2), 3)])

if feature_idx == 4
    ylabel("|CaTA regression coefficient|")
    savefig(['Analysis/Correlations/allcells_', protocol_name, 'CaTAsens_paramError_pop1000.fig'])
    save(['Analysis/Correlations/allcells_', protocol_name, 'CaTAsens_paramError_pop1000.mat'],...
        'x_errors_all', 'y_sum_all', 'cells', 'protocol_name', 'rho', 'pval')
elseif feature_idx == 1 
    ylabel("|APD_{90} regression coefficient|")
    savefig(['Analysis/Correlations/allcells_', protocol_name, 'APD90sens_paramError_pop1000.fig'])
    save(['Analysis/Correlations/allcells_', protocol_name, 'APD90sens_paramError_pop1000.mat'],...
        'x_errors_all', 'y_sum_all', 'cells', 'protocol_name', 'rho', 'pval')
end


%% Parameter sensitivity APD90 vs parameter SD from calibrations
clear
close all

addpath('sensitivity simulations/')

feature_idx = 4 ; % 1 = APD90; 4 = CaTA
cells = [12, 15, 16, 4] ;
protocol_name = 'p4-19-5_norm' ;
SD_file = uigetfile_n_dir([pwd, '/Analysis/AggregateStats'], 'Select file containing stats data') ;

% Sensitivity analysis coefficients
% load('sensitivity simulations/kernik 2019/MLR/APD90 and CaTA/sigma0.2_pop1000.mat')
load('sensitivity simulations/kernik 2019/APD90etc_PLS_1000.mat')
load GA/x_names.mat names

X = log(population) ;
Y = [log(APD90),Vrest,Vpeak, CaTA] ;
Y(isnan(Y)) = 0 ;
% Y = [log(APD90),Vrest,Vpeak, CaTA] ;
[n_samples,n_parameters] = size(X) ;
[n_samples_y,n_outputs] = size(Y) ;

if (n_samples ~= n_samples_y)
  disp('WARNING:  Input & Output matrices appear to be incompatible')
end

logtransform = [1,0,0,0] ;

outputlabels = { ...
    'APD90', ...
    'V_r_e_s_t', ...
    'V_p_e_a_k', ...
    'CaTA', ...
} ;

X_Z = zscore(X) ;
Y_Z = zscore(Y) ;

B = (X_Z'*X_Z)^(-1)*X_Z'*Y_Z ;

Yhat_Z = X_Z*B ;
Yhat = Yhat_Z.*(ones(n_samples,1)*std(Y)) + ones(n_samples,1)*mean(Y) ;

SSYT = sum((Y-ones(n_samples,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(n_samples,1)*mean(Y)).^2);
R2each = SSYR./SSYT; 

%%% PLS
nfactor = length(names); 
[Tpls,P,W,Wstar,U,b,C,B_pls,...
        Bpls_star,Xori_rec,Yori_rec,...
        R2_X,R2_Y]=PLS_nipals(X_Z,Y_Z,nfactor) ;

labels = strsplit(num2str(cells)) ;
dx = 0.01 ;  % for plot text
dy = 0.01 ;  % for plot text

load(SD_file{1})
x_SD_all = zeros(1, length(cells)*length(names)) ;
y_sum_all = zeros(1, length(cells)*length(names)) ;
names_all = cell(1, length(cells)*length(names)) ;

% Get parameter_SDs for each cell 
counter = 1 ;
for i=1:length(cells)
    x_SD = SDs(i, :) ;
    y = abs(B_pls(:, feature_idx)) ; 

    x_SD_all(counter:counter+length(names)-1) = x_SD ;
    y_sum_all(counter:counter+length(names)-1) = y' ;
    names_all(counter:counter+length(names)-1) = names ;
    counter = counter + length(names) ;
end

[rho, pval] = corrcoef(x_SD_all, y_sum_all) ;

colors = repmat('krgbmc', 1, 300) ;
figure
hold on
counter = 1 ;
for i = 1:length(cells)
    plot(x_SD_all(counter:counter+length(names)-1), ...
        y_sum_all(counter:counter+length(names)-1), ...
        [colors(i), 'o']) % also the trend line?
    text(x_SD_all(counter:counter+length(names)-1) + dx, ...
        y_sum_all(:,counter:counter+length(names)-1) + dy, ...
        names_all(counter:counter+length(names)-1))
    counter = counter + length(names) ;
end

xlabel("Parameter SD")
legend(labels, 'Location', 'best')
title(['All cells: Rho = ', num2str(rho(1,2), 3), ', p = ', num2str(pval(1,2), 3)])

if feature_idx == 4
    ylabel("|CaTA regression coefficient|")
    savefig(['Analysis/Correlations/allcells_', protocol_name, 'CaTAsens_paramSD_pop1000.fig'])
    save(['Analysis/Correlations/allcells_', protocol_name, 'CaTAsens_paramSD_pop1000.mat'],...
        'x_SD_all', 'y_sum_all', 'cells', 'protocol_name', 'rho', 'pval')
elseif feature_idx == 1 
    ylabel("|APD_{90} regression coefficient|")
    savefig(['Analysis/Correlations/allcells_', protocol_name, 'APD90sens_paramSD_pop1000.fig'])
    save(['Analysis/Correlations/allcells_', protocol_name, 'APD90sens_paramSD_pop1000.mat'],...
        'x_SD_all', 'y_sum_all', 'cells', 'protocol_name', 'rho', 'pval')
end

%% Parameter sensitivity IKr threshold vs parameter SD from calibrations

clear
close all

addpath('sensitivity simulations/')

idx = 1 ; % 1 = APD90; 4 = CaTA
cells = [12, 15, 16, 4] ;
protocol_name = 'p4-19-5_norm' ;
SD_file = uigetfile_n_dir([pwd, '/Analysis/AggregateStats'], 'Select file containing stats data') ;

% Sensitivity analysis coefficients
% load('sensitivity simulations/kernik 2019/MLR/APD90 and CaTA/sigma0.2_pop1000.mat')
load('sensitivity simulations/kernik 2019/IKr_thresholds_1000.mat')
load GA/x_names.mat names

X = log(population) ;
Y = thresholds ;
% Y = [log(APD90),Vrest,Vpeak, CaTA] ;
Y(isnan(Y)) = 0 ;
% Y = [log(APD90),Vrest,Vpeak, CaTA] ;
[n_samples,n_parameters] = size(X) ;
[n_samples_y,n_outputs] = size(Y) ;

if (n_samples ~= n_samples_y)
  disp('WARNING:  Input & Output matrices appear to be incompatible')
end

logtransform = 0 ;
% logtransform = [1,0,0,0] ;

outputlabels = {'I_{Kr} threshold'} ;
% outputlabels = { ...
%     'APD90', ...
%     'V_r_e_s_t', ...
%     'V_p_e_a_k', ...
%     'CaTA', ...
% } ;

X_Z = zscore(X) ;
Y_Z = zscore(Y) ;

B = (X_Z'*X_Z)^(-1)*X_Z'*Y_Z ;

Yhat_Z = X_Z*B ;
Yhat = Yhat_Z.*(ones(n_samples,1)*std(Y)) + ones(n_samples,1)*mean(Y) ;

SSYT = sum((Y-ones(n_samples,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(n_samples,1)*mean(Y)).^2);
R2each = SSYR./SSYT; 

%%% PLS
nfactor = length(names); 
[Tpls,P,W,Wstar,U,b,C,B_pls,...
        Bpls_star,Xori_rec,Yori_rec,...
        R2_X,R2_Y]=PLS_nipals(X_Z,Y_Z,nfactor) ;

labels = strsplit(num2str(cells)) ;
dx = 0.01 ;  % for plot text
dy = 0.01 ;  % for plot text

load(SD_file{1})
x_SD_all = zeros(1, length(cells)*length(names)) ;
y_sum_all = zeros(1, length(cells)*length(names)) ;
names_all = cell(1, length(cells)*length(names)) ;

% Get parameter_SDs for each cell 
counter = 1 ;
for i=1:length(cells)
    x_SD = SDs(i, :) ;
    y = abs(B_pls(:, idx)) ; 

    x_SD_all(counter:counter+length(names)-1) = x_SD ;
    y_sum_all(counter:counter+length(names)-1) = y' ;
    names_all(counter:counter+length(names)-1) = names ;
    counter = counter + length(names) ;
end

[rho, pval] = corrcoef(x_SD_all, y_sum_all) ;

colors = repmat('krgbmc', 1, 300) ;
figure
hold on
counter = 1 ;
for i = 1:length(cells)
    plot(x_SD_all(counter:counter+length(names)-1), ...
        y_sum_all(counter:counter+length(names)-1), ...
        [colors(i), 'o']) % also the trend line?
    text(x_SD_all(counter:counter+length(names)-1) + dx, ...
        y_sum_all(:,counter:counter+length(names)-1) + dy, ...
        names_all(counter:counter+length(names)-1))
    counter = counter + length(names) ;
end

xlabel("Parameter SD")
ylabel("|I_{Kr} threshold regression coefficient|")
legend(labels, 'Location', 'best')
title(['All cells: Rho = ', num2str(rho(1,2), 3), ', p = ', num2str(pval(1,2), 3)])
savefig(['Analysis/Correlations/allcells_', protocol_name, 'IKrsens_paramSD_pop1000.fig'])
save(['Analysis/Correlations/allcells_', protocol_name, 'IKrsens_paramSD_pop1000.mat'],...
    'x_SD_all', 'y_sum_all', 'cells', 'protocol_name', 'rho', 'pval')

%% Parameter sensitivity IKr threshold vs parameter error from calibrations
clear
close all

addpath('sensitivity simulations/')

idx = 1 ; % 1 = APD90; 4 = CaTA
cells = [12, 15, 16, 4] ;
protocol_name = 'p4-19-5_norm' ;
SD_file = uigetfile_n_dir([pwd, '/Analysis/AggregateStats'], 'Select file containing stats data') ;

% Sensitivity analysis coefficients
% load('sensitivity simulations/kernik 2019/MLR/APD90 and CaTA/sigma0.2_pop1000.mat')
load('sensitivity simulations/kernik 2019/IKr_thresholds_1000.mat')
load GA/x_names.mat names

X = log(population) ;
Y = thresholds ;
% Y = [log(APD90),Vrest,Vpeak, CaTA] ;
Y(isnan(Y)) = 0 ;
% Y = [log(APD90),Vrest,Vpeak, CaTA] ;
[n_samples,n_parameters] = size(X) ;
[n_samples_y,n_outputs] = size(Y) ;

if (n_samples ~= n_samples_y)
  disp('WARNING:  Input & Output matrices appear to be incompatible')
end

logtransform = 0 ;
% logtransform = [1,0,0,0] ;

outputlabels = {'I_{Kr} threshold'} ;
% outputlabels = { ...
%     'APD90', ...
%     'V_r_e_s_t', ...
%     'V_p_e_a_k', ...
%     'CaTA', ...
% } ;

X_Z = zscore(X) ;
Y_Z = zscore(Y) ;

B = (X_Z'*X_Z)^(-1)*X_Z'*Y_Z ;

Yhat_Z = X_Z*B ;
Yhat = Yhat_Z.*(ones(n_samples,1)*std(Y)) + ones(n_samples,1)*mean(Y) ;

SSYT = sum((Y-ones(n_samples,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(n_samples,1)*mean(Y)).^2);
R2each = SSYR./SSYT; 

%%% PLS
nfactor = length(names); 
[Tpls,P,W,Wstar,U,b,C,B_pls,...
        Bpls_star,Xori_rec,Yori_rec,...
        R2_X,R2_Y]=PLS_nipals(X_Z,Y_Z,nfactor) ;

labels = strsplit(num2str(cells)) ;
dx = 0.01 ;  % for plot text
dy = 0.01 ;  % for plot text

load(SD_file{1})
x_errors_all = zeros(1, length(cells)*length(names)) ;
y_sum_all = zeros(1, length(cells)*length(names)) ;
names_all = cell(1, length(cells)*length(names)) ;

% Get parameter_SDs for each cell 
counter = 1 ;
for i=1:length(cells)
    x_error = MAEs(i, :) ;
    y = abs(B_pls(:, idx)) ; 

    x_errors_all(counter:counter+length(names)-1) = x_error ;
    y_sum_all(counter:counter+length(names)-1) = y' ;
    names_all(counter:counter+length(names)-1) = names ;
    counter = counter + length(names) ;
end

[rho, pval] = corrcoef(x_errors_all, y_sum_all) ;

colors = repmat('krgbmc', 1, 300) ;
figure
hold on
counter = 1 ;
for i = 1:length(cells)
    plot(x_errors_all(counter:counter+length(names)-1), ...
        y_sum_all(counter:counter+length(names)-1), ...
        [colors(i), 'o']) % also the trend line?
    text(x_errors_all(counter:counter+length(names)-1) + dx, ...
        y_sum_all(:,counter:counter+length(names)-1) + dy, ...
        names_all(counter:counter+length(names)-1))
    counter = counter + length(names) ;
end

xlabel("Parameter error")
ylabel("|I_{Kr} threshold regression coefficient|")
legend(labels, 'Location', 'best')
title(['All cells: Rho = ', num2str(rho(1,2), 3), ', p = ', num2str(pval(1,2), 3)])
savefig(['Analysis/Correlations/allcells_', protocol_name, 'IKrsens_paramError_pop1000.fig'])
save(['Analysis/Correlations/allcells_', protocol_name, 'IKrsens_paramError_pop1000.mat'],...
    'x_errors_all', 'y_sum_all', 'cells', 'protocol_name', 'rho', 'pval')


