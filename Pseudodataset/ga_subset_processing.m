clear
%%
load('Pseudodataset/saved_data/ga_subsets.mat')
subsets = {ss_caohigh, ss_caohigh_sp, ss_caolow, ss_caolow_sp, ss_ical25, ss_ical25_sp, ss_ical50, ss_ical50_sp, ss_ikr15, ss_ikr15_sp, ss_ikr30, ss_ikr30_sp, ss_kohigh, ss_kohigh_sp, ss_kolow, ss_kolow_sp, ss_pcl600, ss_pcl800, ss_spont} ;
%% get spontfreq
for i = 1:19
    for j = 1:10
        [~,~,features] = APextract_custom(subsets{i}{j}.state.t, subsets{i}{j}.state.Y(:,1), subsets{i}{j}.protocol.stimtimes);
        subsets{i}{j}.protocol.spontfreq = features(4);
    end
end
%% extract last 5s
spont = cell(1,10);
for j = 1:10
    spont{j} = extractlast5s(ss_spont{j});
end

pcl600 = cell(1,10);
for j = 1:10
    pcl600{j} = extractlast5s(ss_pcl600{j});
end

pcl800 = cell(1,10);
for j = 1:10
    pcl800{j} = extractlast5s(ss_pcl800{j});
end

caohigh = cell(1,10);
for j = 1:10
    caohigh{j} = extractlast5s(ss_caohigh{j});
end

caohigh_sp = cell(1,10);
for j = 1:10
    caohigh_sp{j} = extractlast5s(ss_caohigh_sp{j});
end

caolow = cell(1,10);
for j = 1:10
    caolow{j} = extractlast5s(ss_caolow{j});
end

caolow_sp = cell(1,10);
for j = 1:10
    caolow_sp{j} = extractlast5s(ss_caolow_sp{j});
end

kohigh = cell(1,10);
for j = 1:10
    kohigh{j} = extractlast5s(ss_kohigh{j});
end

kohigh_sp = cell(1,10);
for j = 1:10
    kohigh_sp{j} = extractlast5s(ss_kohigh_sp{j});
end

kolow = cell(1,10);
for j = 1:10
    kolow{j} = extractlast5s(ss_kolow{j});
end

kolow_sp = cell(1,10);
for j = 1:10
    kolow_sp{j} = extractlast5s(ss_kolow_sp{j});
end

ical25 = cell(1,10);
for j = 1:10
    ical25{j} = extractlast5s(ss_ical25{j});
end

ical25_sp = cell(1,10);
for j = 1:10
    ical25_sp{j} = extractlast5s(ss_ical25_sp{j});
end

ical50 = cell(1,10);
for j = 1:10
    ical50{j} = extractlast5s(ss_ical50{j});
end

ical50_sp = cell(1,10);
for j = 1:10
    ical50_sp{j} = extractlast5s(ss_ical50_sp{j});
end

ikr15 = cell(1,10);
for j = 1:10
    ikr15{j} = extractlast5s(ss_ikr15{j});
end

ikr15_sp = cell(1,10);
for j = 1:10
    ikr15_sp{j} = extractlast5s(ss_ikr15_sp{j});
end

ikr30 = cell(1,10);
for j = 1:10
    ikr30{j} = extractlast5s(ss_ikr30{j});
end

ikr30_sp = cell(1,10);
for j = 1:10
    ikr30_sp{j} = extractlast5s(ss_ikr30_sp{j});
end

data = {caohigh, caohigh_sp, caolow, caolow_sp, ical25,...5
    ical25_sp, ical50, ical50_sp, ikr15, ikr15_sp,...10
    ikr30, ikr30_sp, kohigh, kohigh_sp, kolow,...15
    kolow_sp, pcl600, pcl800, spont} ;...19

% save('Pseudodataset/saved_data/datapreinterp.mat', 'data', '-v7.3')

%% Get even timestep using interp1
time = 0:0.1:5001;
for i = 1:19
    for j = 1:10
        modV = interp1(data{i}{j}{2},data{i}{j}{3}, time);
        modCai = interp1(data{i}{j}{2},data{i}{j}{4}, time);
        data{i}{j}{2} = time';
        data{i}{j}{3} = modV';
        data{i}{j}{4} = modCai';
    end
end

% save('Pseudodataset/saved_data/datapostinterp.mat', 'data', '-v7.3')

%% to export
%
filename = 'Pseudodataset/saved_data/pseudodataset.xlsx';
for i = 1:10
    matrix = zeros(50011,39);
    matrix(:,1) = time';
    for j = 1:19
        matrix(:, 2*j) = data{j}{i}{3};
        matrix(:, 2*j+1) = data{j}{i}{4};
    end
    writematrix(matrix,filename,'Sheet',i, 'Range', 'A3:AM50013' )
end
%}
















