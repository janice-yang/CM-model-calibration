function f_plotCaHandling(cells) % cells is a cell array of CMs
numcells = length(cells)-1 ;

%figure
for i = 1:numcells
    CM = cells{i} ; 
    
    if CM.state.t(end) < 2000 % assuming we aren't running this for <= 2 seconds
        timeinseconds = CM.state.t ;
    elseif CM.state.t(end) >= 2000
        timeinseconds = CM.state.t/1000 ;
    end
    
%% find and plot Cai
    isCai = zeros(1,length(CM.YNames)) ; 
for j = 1:length(isCai)
    isCai(j) = isequal(CM.YNames{j},'Cai') || isequal(CM.YNames{j},'Ca_i') ;
end
cai_index = find(isCai) ;
if length(cai_index) ~= 1
    error('Error: Generate plot manually as Cai is not named Cai or Ca_i')
end

%subplot(numcells, 2, 2*i-1)
figure
plot(timeinseconds, CM.state.Y(:,cai_index(1)), 'k- ', 'LineWidth', 0.7)
ylim([min(CM.state.Y(:,cai_index(1))) max(CM.state.Y(:,cai_index(1)))]) ; 
title(strcat('Cai plot for :', CM.name))

%% find and plot CaSR
isCaSR = zeros(1,length(CM.YNames)) ;

if ~isequal(CM.name,'ORd11 CM')
    for j = 1:length(isCaSR)
        isCaSR(j) = isequal(CM.YNames{j},'Ca_SR') || isequal(CM.YNames{j},'Ca_S_R') || isequal(CM.YNames{j},'Ca_nsr');
    end
    caSR_index = find(isCaSR) ;
    if length(caSR_index) ~= 1
        error('Error: Generate plot manually as CaSR is not named Ca_SR or Ca_S_R or Ca_nsr')
    end
    
%    subplot(numcells, 2, 2*i)
    figure
    plot(timeinseconds, CM.state.Y(:,caSR_index(1)), 'k- ', 'LineWidth', 0.7)
    ylim([min(CM.state.Y(:,caSR_index(1))) max(CM.state.Y(:,caSR_index(1)))]) ;
    title(strcat('CaSR plot for :', CM.name))
    
elseif isequal(CM.name,'ORd11 CM')
    for j = 1:length(isCaSR)
        isCaSR(j) = isequal(CM.YNames{j},'Ca_SR') || isequal(CM.YNames{j},'Ca_S_R') || isequal(CM.YNames{j},'Ca_nsr');
    end
    caSR_index = find(isCaSR) ;
    if length(caSR_index) ~= 1
        error('Error: Generate plot manually as CaSR is not named Ca_SR or Ca_S_R or Ca_nsr')
    end
    
  %  subplot(numcells, 2, 2*i)
    figure
    plot(timeinseconds, CM.state.Y(:,caSR_index(1)), 'k- ', 'LineWidth', 0.7)
    ylim([min(CM.state.Y(:,caSR_index(1))) max(CM.state.Y(:,caSR_index(1)))]) ;
    title(strcat('CaSR plot for :', CM.name))
end

end

end