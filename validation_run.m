clear
close all

addpath('Other_helper_functions')
realdata = true ;
method = 'single' ; % single, median, or meanbjob

folders = uigetfile_n_dir([pwd, '/GA/Results/'], 'Select DIRECTORIES containing GA results') ;
isNormalized = false ;
if ispc % Windows
    foldersplit = strsplit(folders{1}, '\') ;
else
    foldersplit = strsplit(folders{1}, '/') ;    
end
base = strjoin(foldersplit(1:end-1), '/') ; % experiment folder with all runs

if realdata
    datatype = 'APCaT' ; % 'AP', 'CaT', or 'APCaT'
    exp_file = uigetfile_n_dir([pwd, '/ExperimentalData/'], 'Select file containing VALIDATION data') ;
    % sheetnames = {'Sheet1'} ; % each sheet = 1 validation protocol
    sheetnames = {'DMG240_50pDofetilide1nM_0.5Hz'} ; % each sheet = 1 validation protocol
    % protocol_numbers = 29 ;
    protocol_numbers = 35 ;
    
    exp_file = exp_file{1} ;
    datamatrix = readmatrix(exp_file) ;
    stimtimes = datamatrix(:, end) ;
%     expT = datamatrix(:, 1) ;
%     stim_idx = ismember(datamatrix(:,2), datamatrix(:, end)) ;
%     exp_stimtimes = expT(stim_idx) ;
else
    load([folders{1}, '/Details.mat'], 'cell_number')
    protocol_numbers = [12] ; % validation protocol(s)
end

set(groot, 'defaultFigureRenderer', 'painters') % for saving as svg

%%

if realdata
    for i=1:length(folders)
        for j=1:size(protocol_numbers, 1)
            f_validation_realdata(folders{i}, isNormalized, exp_file, sheetnames, ...
                datatype, stimtimes, protocol_numbers(j, :)) ;
        end
        savefig([folders{i}, '/validation_test_5s_', mat2str(protocol_numbers(j, :))])
        savefig([base, '/extracted_5s_', int2str(i)])
        close
        disp(num2str(i))
    end

else
    for i=1:size(protocol_numbers, 1)
        f_validation_run(folders, cell_number, protocol_numbers(i, :), method, isNormalized) ;
        disp(['Validation run finished for cell ', int2str(cell_number), ...
            ', protocol ' mat2str(protocol_numbers(i, :))])
    end
end



