function [small_pop] = sbr_filter(big_pop, max_sbr, min_sbr)
popsize = length(big_pop);
small_pop = {};
sp_dex = 1;
for i = 1:popsize
    if big_pop{i}.protocol.spontfreq < max_sbr && big_pop{i}.protocol.spontfreq > min_sbr
        small_pop{sp_dex} = big_pop{i};
        sp_dex = sp_dex + 1;
    end
end
end