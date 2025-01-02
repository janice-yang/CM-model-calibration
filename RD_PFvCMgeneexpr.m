clear
close all

% %% Main conductances to be changed:
% % OBS: The numbers multipling the conductance were calculated from RNA-Seq data
% GNa = optimparams(1)*2.37 ; % S/F
% GK1 = optimparams(2)*0.82 ; % S/F
% GKr = optimparams(3)*1.18 ; % S/F
% GKs = optimparams(4)*1.02 ; % S/F
% Gto = optimparams(5)*0.85 ; % S/F
% GCaL = optimparams(6)*0.96 ; % m^3 /(Fxs)
% KNaCA = optimparams(7)*1.22 ; % A/F
% PNaK = optimparams(8)*1.22 ; % A/F
% GpCa = optimparams(9)*1.06 ; % A/F
% GbNa = optimparams(10) ; % S/F
% GbCa = optimparams(11) ; % S/F
% G_RyR = optimparams(12)*1.13 ;
% Vmaxup = optimparams(13)*1.08;
% Vleak = optimparams(14) ; % (1/s)
% Gf = optimparams(15)*1.78; % S/F

% log scaling multipliers
logfactor = 2 ;
GNa = 2.37 ; % S/F
GK1 = 0.82 ; % S/F
GKr = 1.18 ; % S/F
GKs = 1.02 ; % S/F
Gto = 0.85 ; % S/F
GCaL = 0.96 ; % m^3 /(Fxs)
KNaCA = 1.22 ; % A/F
PNaK = 1.22 ; % A/F
GpCa = 1.06 ; % A/F
GbNa = 1 ; % S/F
GbCa = 1 ; % S/F
G_RyR = 1.13 ;
Vmaxup = 1.08 ;
Vleak = 1 ; % (1/s)
Gf = 1.78 ; % S/F

parameternames = {'GNa', 'Gf', 'NaK', 'GKr', 'RyR', 'NCX', 'SERCA', 'GpCa', 'GKs', 'GCaL', 'Gto', 'GK1'} ;
scalefactors = [GNa, Gf, PNaK, GKr, G_RyR, KNaCA, Vmaxup, GpCa, GKs, GCaL, Gto, GK1] ;
logscalefactors = zeros(1, length(scalefactors)) ;

for i=1:length(scalefactors)
    logscalefactors(i) = log2(scalefactors(i)) ;
end

bar(logscalefactors)
xticks(1:length(parameternames))
xticklabels(parameternames)
xtickangle(45)
ylim([-0.5 1])
title(['Log' num2str(logfactor) ' PF/VM'])

%% Parameter estimates PF vs VM - run by run
clear

logfactor = 2 ;
Kernik_baseline_names = {'GNa', 'Gf', 'GCaL', 'Gto', 'GKs', 'GKr', 'GK1', 'GpCa', 'GbNa', 'GbCa', 'GUp', 'Grel', 'GNCX', 'GNaK', 'Gleak', 'GCaT'} ;
Kernik_baseline_conds = [9.720613409241, 0.0435, 0.308027691379, 0.1178333333333, 0.0077, 0.218025, 0.133785777797606, 0.2625, 4.35e-04, 3.6704e-04, 1.105e-04, 12.5, 1100, 2.4761, 1.6e-06, 0.185] ;

PF_folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select DIRECTORY/IES containing PF results') ;
VM_folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select DIRECTORY/IES containing VM results') ;

PF_VM = zeros(length(PF_folders), length(Kernik_baseline_conds)) ;
for i=1:length(PF_folders)
    load([PF_folders{i}, '/minoptimparams.mat'], 'minoptimparams')
    PF_scaling = logfactor.^minoptimparams ;
    PF_conds = PF_scaling .* Kernik_baseline_conds ;

    load([VM_folders{i}, '/minoptimparams.mat'], 'minoptimparams')
    VM_scaling = logfactor.^minoptimparams ;
    VM_conds = VM_scaling .* Kernik_baseline_conds ;

    PF_VM(i, :) = PF_conds ./ VM_conds ;
end

figure
subplot(1,2,1)
hold on
ylabel('ln(PF/VM)')
title('Parameter estimates')
plot(1:length(Kernik_baseline_names), log(PF_VM), 'k*')
bar(median(log(PF_VM)))
xticks(1:length(Kernik_baseline_names))
xticklabels(Kernik_baseline_names)
xtickangle(45)
% ylim([-0.5 1])

subplot(1,2,2)
hold on
title('Gene expression')
ylabel('ln(PF/VM)')
GNa = 2.37 ; % S/F
GK1 = 0.82 ; % S/F
GKr = 1.18 ; % S/F
GKs = 1.02 ; % S/F
Gto = 0.85 ; % S/F
GCaL = 0.96 ; % m^3 /(Fxs)
KNaCA = 1.22 ; % A/F
PNaK = 1.22 ; % A/F
GpCa = 1.06 ; % A/F
GbNa = 1 ; % S/F
GbCa = 1 ; % S/F
G_RyR = 1.13 ;
Vmaxup = 1.08 ;
Vleak = 1 ; % (1/s)
Gf = 1.78 ; % S/F
DEG_conds = [GNa, Gf, GCaL, Gto, GKs, GKr, GK1, GpCa, Vmaxup, G_RyR, KNaCA, PNaK] ;
DEG_names = {'GNa', 'Gf', 'GCaL', 'Gto', 'GKs', 'GKr', 'GK1', 'GpCa', 'GUp', 'Grel', 'GNCX', 'GNaK'} ;
bar(log(DEG_conds))
xticks(1:length(DEG_names))
xticklabels(DEG_names)
xtickangle(45)
ylim([-0.5 1])

%% Parameter estimates PF vs VM - averaged over 10 runs
clear

logfactor = 2 ;
Kernik_baseline_names = {'GNa', 'Gf', 'GCaL', 'Gto', 'GKs', 'GKr', 'GK1', 'GpCa', 'GbNa', 'GbCa', 'GUp', 'Grel', 'GNCX', 'GNaK', 'Gleak', 'GCaT'} ;
Kernik_baseline_conds = [9.720613409241, 0.0435, 0.308027691379, 0.1178333333333, 0.0077, 0.218025, 0.133785777797606, 0.2625, 4.35e-04, 3.6704e-04, 1.105e-04, 12.5, 1100, 2.4761, 1.6e-06, 0.185] ;

PF_folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select DIRECTORY/IES containing PF results') ;
VM_folders = uigetfile_n_dir([pwd, '/GA/Results'], 'Select DIRECTORY/IES containing VM results') ;
%
PF_all = zeros(length(PF_folders), length(Kernik_baseline_conds)) ;
VM_all = zeros(length(VM_folders), length(Kernik_baseline_conds)) ;
for i=1:length(PF_folders)
    load([PF_folders{i}, '/minoptimparams.mat'], 'minoptimparams')
    PF_scaling = logfactor.^minoptimparams ;
    PF_all(i, :) = PF_scaling .* Kernik_baseline_conds ;

    load([VM_folders{i}, '/minoptimparams.mat'], 'minoptimparams')
    VM_scaling = logfactor.^minoptimparams ;
    VM_all(i, :) = VM_scaling .* Kernik_baseline_conds ;
end
mean_PF_VM = mean(PF_all) ./ mean(VM_all) ;
med_PF_VM = median(PF_all) ./ median(VM_all) ;

figure
subplot(3,1,1)
hold on
ylabel('log_2(PF/VM)')
title('Median parameter estimates')
bar(log2(med_PF_VM))
xticks(1:length(Kernik_baseline_names))
xticklabels(Kernik_baseline_names)
xtickangle(45)
% ylim([-0.5 1])

subplot(3,1,2)
hold on
ylabel('log_2(PF/VM)')
title('Mean parameter estimates')
bar(log2(mean_PF_VM))
xticks(1:length(Kernik_baseline_names))
xticklabels(Kernik_baseline_names)
xtickangle(45)
% ylim([-0.5 1])

subplot(3,1,3)
hold on
title('Gene expression')
ylabel('log_2(PF/VM)')
GNa = 2.37 ; % S/F
GK1 = 0.82 ; % S/F
GKr = 1.18 ; % S/F
GKs = 1.02 ; % S/F
Gto = 0.85 ; % S/F
GCaL = 0.96 ; % m^3 /(Fxs)
KNaCA = 1.22 ; % A/F
PNaK = 1.22 ; % A/F
GpCa = 1.06 ; % A/F
GbNa = 1 ; % S/F
GbCa = 1 ; % S/F
G_RyR = 1.13 ;
Vmaxup = 1.08 ;
Vleak = 1 ; % (1/s)
Gf = 1.78 ; % S/F
DEG_conds = [GNa, Gf, GCaL, Gto, GKs, GKr, GK1, GpCa, Vmaxup, G_RyR, KNaCA, PNaK] ;
DEG_names = {'GNa', 'Gf', 'GCaL', 'Gto', 'GKs', 'GKr', 'GK1', 'GpCa', 'GUp', 'Grel', 'GNCX', 'GNaK'} ;
bar(log2(DEG_conds))
xticks(1:length(DEG_names))
xticklabels(DEG_names)
xtickangle(45)
ylim([-0.5 1])