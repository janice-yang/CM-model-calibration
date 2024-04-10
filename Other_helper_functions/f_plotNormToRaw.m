function [] = f_plotNormToRaw(folder, logfactor, legendLocation)
%%%%%%
% Plot raw AP & CaT (extracted) using minoptimparams from normalized GA run
%   Also plot corresponding experimental raw data for comparison
%   Save plot(s) to results folder
% 
% Input:
%   folder = (string) directory containing GA results on normalized data
%   logfactor = (int) scaling factor for parameters
%   legendLocation = (string) location of legend in plot
%       e.g. 'southwest', 'northeast'
%
% Output:
%   Plot of raw AP & CaT from predicted params + experimental data
%%%%%% 

% Conductance names
names = {...
'g_{Na}',...        1
'g_f',...         2
'p_{CaL}',...       3   
'g_{to}',...        4  
'g_{Ks}',...        5
'g_{Kr}',...        6
'g_{K1}',...        7
'g_{PCa}',...       8
'g_{b\_Na}',...      9
'g_{b\_Ca}',...      10
'VmaxUp',...      11
'g_{irel\_max}',...  12
'kNaCa',...       13
'PNaK',...        14
'V_{leak}',...      15
'g_{CaT}',...       16
} ;

folderFormat = folder(end) ; % to allow for '/' at end of folder path

if folderFormat == '/'
    load([folder, 'Details.mat'], 'cell_number', 'protocol_number')
    load([folder, '/minoptimparams.mat'], 'minoptimparams') % uncomment this line + comment next if loading from minoptimparams.mat works
    % load([folder, 'final_workspace.mat']) % workaround for minoptimparams loading
    close % close output figure from final_workspace.mat
else
    load([folder, '/Details.mat'], 'cell_number', 'protocol_number')
    load([folder, '/minoptimparams.mat'], 'minoptimparams') % uncomment this line + comment next if loading from minoptimparams.mat works
    % load([folder, '/final_workspace.mat']) 
    close 
end

isNormalized = false ;
save GA/curr_cell_protocol cell_number protocol_number isNormalized

experimental_dataset = f_getPseudodata(cell_number, protocol_number, isNormalized) ;

% Simulate raw AP and CaT from predicted params
disp('Simulating from predicted conductances...')

for i = 1:size(minoptimparams,1)
params = minoptimparams(i,:);
x_conductance = (logfactor.^params);
[t_stim, V_stim, Cai_stim, stimtimes] = ga_simulation_k19(x_conductance,names); % see how the optimal parameters do

% Extract for each protocol in series
figure 
p = 1 ;
for j=1:length(protocol_number)
    [keepT, V, CaT,tinit,errorcode] = waveform_extract_new(t_stim{j}, V_stim{j},Cai_stim{j},stimtimes{j});

    % Align ends of exp and simulated extracted waveforms
    exp_T_V = experimental_dataset{j}.Time_AP{1} ; % - tpeak + 100;
    exp_V = experimental_dataset{j}.AP{1};
    exp_T_CaT = experimental_dataset{j}.Time_CaT{1} ; % - tpeak + 100;
    exp_CaT = experimental_dataset{j}.CaT{1};

    [exp_T_V, exp_V, T_V, V] = f_alignWaveformEnds(exp_T_V, exp_V, keepT, V) ;
    [exp_T_CaT, exp_CaT, T_CaT, CaT] = f_alignWaveformEnds(exp_T_CaT, exp_CaT, keepT, CaT) ;

    V_sim = interp1(T_V,V,exp_T_V);
    CaT_sim = interp1(T_CaT, CaT, exp_T_CaT);
    
    subplot(length(protocol_number),2,p)
    plot(exp_T_V,V_sim,'LineWidth',2)
    hold on
    plot(exp_T_V,exp_V, 'Color', 'red','LineWidth',2)
    title(['Protocol ',int2str(protocol_number(j)),' Action Potential'])
    ylabel('mV')
    legend('Prediction','Observation', 'Location',legendLocation)
    set(gca,'FontSize',12)

    subplot(length(protocol_number),2,p+1)
    plot(exp_T_CaT, CaT_sim,'LineWidth',2)
    hold on
    plot(exp_T_CaT, exp_CaT, 'Color', 'red','LineWidth',2)
    title(['Protocol ',int2str(protocol_number(j)),' CaT'])
    ylabel('[Ca^{2+}]_i (mM)')
    set(gca,'FontSize',12)
    hold off
    
    % Increment subplot position
    p = p + 2 ; 
    
end

fig = gcf;
fig.Units = 'normalized';
fig.Position = [0.25 0.3 0.5 0.4];
if folderFormat == '/'
    print('-dpng', [folder, 'raw_output_minoptimparams_', num2str(i)]) ;
else
    print('-dpng', [folder, '/raw_output_minoptimparams_', num2str(i)]) ;
end

end

end