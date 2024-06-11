clear
close all

% Select figure files
figures = uigetfile_n_dir([pwd, '/GA/Results'], 'Select FILES containing figures:') ;
main_idx = length(figures) ; % last figure will store all runs

for i=1:length(figures)
    openfig(figures{i}) ;
end

%% STOP! delete true trace (red dashed line) from all except the main figure
colors = repmat('kgbmc', 1, 300) ;
for j=1:length(figures)
    if j ~= main_idx    % Skip main figure
        figure(j)
        fig = gcf ;
        ax = fig.Children ;

        figure(main_idx)
        hold on

        % AP (axes are reversed in saved figure)
        lines = findobj(ax(2), 'Type', 'line') ;
        for i=1:length(lines)
            V = lines(i).YData ;
            t = lines(i).XData ;
            subplot(2,2,1)
            hold on
            plot(t,V, [colors(j), '-'], 'LineWidth', 1)
        end
    
        % CaT
        lines = findobj(ax(1), 'Type', 'line') ;
        for i=1:length(lines)
            Cai = lines(i).YData ;
            t = lines(i).XData ;
            subplot(2,2,3)
            hold on
            plot(t,Cai, [colors(j), '-'], 'LineWidth', 1)
        end
        
    end
end



