%addpath('gaToRORd19')

load experimental_dataset experimental_dataset

load x_names names
runNum=1;
      
nvars = length(names);

minoptimparams=bestparams;
pace=1;

for i=1:size(minoptimparams,1)
params=minoptimparams(i,:);
x_conductance = (10.^params);
[t_stim, V_stim, Cai_stim, stimtimes] = pvcm_function_tord19(x_conductance, pace,names); % see how the optimal parameters do
[T, V,CaT] = APextract_one(t_stim, V_stim,Cai_stim, stimtimes); % extract last AP only

%Vderiv = [experimental_dataset.mean_AP(2:end); experimental_dataset.mean_AP(end)] - experimental_dataset.mean_AP;
[peakV, peakdex] = max(experimental_dataset.mean_AP);
%[restV, restdex] = min(experimental_dataset.mean_AP);
tpeak = experimental_dataset.Time(peakdex);

exp_T = experimental_dataset.Time - tpeak + 100;
V_sim = interp1(T,V,exp_T);
CaT_sim = interp1(T, CaT, exp_T);


%Normalize CaT_sim for plot:
CaT_sim_norm=(CaT_sim-min(CaT_sim))/(max(CaT_sim)-min(CaT_sim));


figure % perhaps remove this figure once the GA is up and running
subplot(1,3,1)
plot(exp_T,V_sim,'LineWidth',2)
hold on
plot(exp_T,experimental_dataset.mean_AP,'Color','red','LineWidth',2)
title('Action Potential')
legend('Prediction','Observation')
set(gca,'FontSize',12)


subplot(1,3,2)
plot(exp_T, CaT_sim_norm,'LineWidth',2)
hold on
plot(exp_T, experimental_dataset.mean_CaT, 'Color', 'red','LineWidth',2)
title('CaT')
set(gca,'FontSize',12)
hold off

subplot(1,3,3)
plot(exp_T, CaT_sim,'LineWidth',2)
set(gca,'FontSize',12)
title('Amplitude')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 4];
print(['Results/Run_',int2str(runNum), '/Outputs_',int2str(i)],'-dpng','-r0')
end