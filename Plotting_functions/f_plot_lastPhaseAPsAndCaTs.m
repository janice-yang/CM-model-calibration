function f_plot_lastPhaseAPsAndCaTs(cells, drugTarget, percentBlock) % cells is a cell array of CMs
numcells = length(cells)-1 ;

%figure
for i = 1:numcells
    CM = cells{i} ;
    
    % set time in all cells to seconds
    if CM.state.t(end) < 2000 % assuming we aren't running this for <= 2 seconds
        timeinseconds = CM.state.t ;
    elseif CM.state.t(end) >= 2000
        timeinseconds = CM.state.t/1000 ;
    end
    
    VmIndex = 0;
    CaiIndex = 0;
    lastPhaseLength = 0;
    firstPhaseLength = 0;
    secondPhaseLength = 0;
    %% find and plot AP & CaT
    if isequal(CM.name, 'Kernik19 CM')
        VmIndex = 1 ;
        CaiIndex = 3 ;
        lastPhaseLength = CM.protocol.phaseLengths(end);
        firstPhaseLength = CM.protocol.phaseLengths(1);
        secondPhaseLength = CM.protocol.phaseLengths(2);
        
    elseif isequal(CM.name, 'Paci13v CM')
        VmIndex = 15 ;
        CaiIndex = 17 ;
        lastPhaseLength = CM.protocol.phaseLengths(end);
        firstPhaseLength = CM.protocol.phaseLengths(1);
        secondPhaseLength = CM.protocol.phaseLengths(2);
        
    elseif isequal(CM.name, 'Paci18 CM')
        VmIndex = 1 ;
        CaiIndex = 3 ;
        lastPhaseLength = CM.protocol.phaseLengths(end);
        firstPhaseLength = CM.protocol.phaseLengths(1);
        secondPhaseLength = CM.protocol.phaseLengths(2);
        
    elseif isequal(CM.name, 'ORd11 CM')
        VmIndex = 1 ;
        CaiIndex = 6 ;
        lastPhaseLength = CM.protocol.pacing.phaseLengths(end);
        firstPhaseLength = CM.protocol.pacing.phaseLengths(1);
        secondPhaseLength = CM.protocol.pacing.phaseLengths(2);
    end
    
    % set Vm in all cells to millivolts
    if isequal(CM.YUnits{VmIndex}, 'mV')
        VminmV = CM.state.Y(:,VmIndex) ;
    elseif isequal(CM.YUnits{VmIndex}, 'V')
        VminmV = CM.state.Y(:,VmIndex)*1000 ;
    end
    
    drug_start = length(timeinseconds(timeinseconds < (timeinseconds(end)-lastPhaseLength-1.2)));
    drug_end = length(timeinseconds(timeinseconds < (timeinseconds(end)-1.2)));
    
    endTimeToPlot = timeinseconds(drug_start:drug_end)-timeinseconds(drug_start);
    %endIndicesToPlot = length(endTimeToPlot)-1;
    endVmToPlot = VminmV(drug_start:drug_end);
    endCaiToPlot = CM.state.Y(drug_start:drug_end,CaiIndex);
    
    norm_start = length(timeinseconds(timeinseconds < firstPhaseLength));
    norm_end = length(timeinseconds(timeinseconds < (firstPhaseLength+secondPhaseLength)));
    normTimeToPlot = timeinseconds(norm_start:norm_end)-timeinseconds(norm_start);
    %normIndicesToPlot = length(normTimeToPlot)-1;
    normVmToPlot = VminmV(norm_start:norm_end);
    normCaiToPlot = CM.state.Y(norm_start:norm_end,CaiIndex);

    %subplot(numcells, 2, 2*i-1)
    figure
    hold on
    plot(normTimeToPlot, normVmToPlot, 'k- ', 'LineWidth', 0.7)
    plot(endTimeToPlot, endVmToPlot, 'r- ', 'LineWidth', 0.7)
    ylim([-100 100]) ;
    title(strcat('AP for: ', CM.name, ' with ', num2str(100-percentBlock*100), '% block of ', drugTarget))
    legend('untreated', 'with drug')
    
    %subplot(numcells, 2, 2*i)
    figure
    hold on
    plot(normTimeToPlot, normCaiToPlot, 'k- ', 'LineWidth', 0.7)
    plot(endTimeToPlot, endCaiToPlot, 'r- ', 'LineWidth', 0.7)
    ylim([0 5e-4]) ;
    title(strcat('CaT for: ', CM.name, ' with ', num2str(100-percentBlock*100), '% block of ', drugTarget))
    legend('untreated', 'with drug')
   
end

end


















