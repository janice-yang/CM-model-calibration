function [bests,medians,iqrs,ngens,npops] = f_fitness_progression(folders)
%{
Analyze the progression 
Inputs:
    folders: cell array of folders containing GA results
Outputs: cell or numeric arrays storing progression of various metrics across generations
    bests: best fitness
    means: mean fitness
    vars: variance in fitness w/in population
    ngens: number of generations for each run
    npops: population size for each run
%}

bests = cell(1, length(folders)) ;
medians = cell(1, length(folders)) ;
iqrs = cell(1, length(folders)) ;
ngens = zeros(1, length(folders)) ;
npops = zeros(1, length(folders)) ;

for f=1:length(folders) % Loop through each run
    load([folders{f}, '/Details.mat'], 'options', 'popsize') ;
    npops(f) = popsize ;
    search = dir([folders{f}, '/Run_*k19output_*.mat']) ;
    matfiles = {search.name} ;
    ngens(f) = length(matfiles) - 2 ;
    
    data = zeros(3, length(matfiles) - 2) ; % best score, mean score, variance
    
    for g=1:length(matfiles) % Loop through each generation's data
        load([folders{f}, '/', matfiles{g}], 'currentscore', 'state')
        gen_num = state.Generation ;
        if ~(gen_num == 0)
            data(1, gen_num) = min(currentscore) ;
            data(2, gen_num) = median(currentscore) ;
            data(3, gen_num) = iqr(currentscore) ;
        end
    end
    
    bests{f} = data(1, :) ; 
    medians{f} = data(2, :) ;
    iqrs{f} = data(3, :) ;
end


end