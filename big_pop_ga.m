
%% STARTING POPULATION:
clear

T = 310.0 ;         % (K) able to manipulate this?
Nao = 151.0 ;       % (mM) 145, 140 ? able to manipulate this?
Ko = 5.4 ;          % (mM) able to manipulate this?
Cao = 1.8 ;         % (mM) able to manipulate this?

popsize = 100;
load GA/x_names.mat names
numParams = length(names) ;
rng(0)

sigma = 0.1;
x_conductances = exp(randn(popsize,numParams)*sigma) ;
cells1 = cell(1,popsize);
for i = 1:popsize
    k19 = gaKernik19([-75.7444536163477,0.338969592726776,0.000203113729306209,7.16928093750999,104.748824394112,0,0.000386686486786781,0.165948056753057,0.927145173320106,0.321984775889061,0.452222061313948,0.157787210225653,0.743251868606151,0.121059208476135,0.0292207648412020,0.00620538308203812,0.736108314102295,0.000264118925707198,0.000263623380304203,0.746780367359425,0.0122283684412607,0.000154416932298621,0.0123158737520428,0], ...
       i) ; % requires "cellID" - just using cell # here
    scaleCell(k19, x_conductances(i,:), ones(1,58))
    saveX_conductance(k19, x_conductances(i,:))
    k19.name = [k19.name, '_', int2str(i)];
    cells1{i} = k19;
end
save big_pseudocells_0pt1.mat cells1

sigma = 0.2;
x_conductances = exp(randn(popsize,numParams)*sigma) ;
cells2 = cell(1,popsize);
for i = 1:popsize
    k19 = gaKernik19([-75.7444536163477,0.338969592726776,0.000203113729306209,7.16928093750999,104.748824394112,0,0.000386686486786781,0.165948056753057,0.927145173320106,0.321984775889061,0.452222061313948,0.157787210225653,0.743251868606151,0.121059208476135,0.0292207648412020,0.00620538308203812,0.736108314102295,0.000264118925707198,0.000263623380304203,0.746780367359425,0.0122283684412607,0.000154416932298621,0.0123158737520428,0], ...
       i) ; % requires "cellID" - just using cell # here
    scaleCell(k19, x_conductances(i,:), ones(1,58))
    saveX_conductance(k19, x_conductances(i,:))
    k19.name = [k19.name, '_', int2str(i)];
    cells2{i} = k19;
end
save big_pseudocells_0pt2.mat cells2

sigma = 0.3;
x_conductances = exp(randn(popsize,numParams)*sigma) ;
cells3 = cell(1,popsize);
for i = 1:popsize
    k19 = gaKernik19([-75.7444536163477,0.338969592726776,0.000203113729306209,7.16928093750999,104.748824394112,0,0.000386686486786781,0.165948056753057,0.927145173320106,0.321984775889061,0.452222061313948,0.157787210225653,0.743251868606151,0.121059208476135,0.0292207648412020,0.00620538308203812,0.736108314102295,0.000264118925707198,0.000263623380304203,0.746780367359425,0.0122283684412607,0.000154416932298621,0.0123158737520428,0], ...
       i) ; % requires "cellID" - just using cell # here
    scaleCell(k19, x_conductances(i,:), ones(1,58))
    saveX_conductance(k19, x_conductances(i,:))
    k19.name = [k19.name, '_', int2str(i)];
    cells3{i} = k19;
end
save big_pseudocells_0pt3.mat cells3
%%
for i = 1:popsize
    setEnvironment(cells1{i}, T, Nao, Cao, Ko) ;
    setUpPacingProtocol(cells1{i}, 0, 300, 999, 1)
    odeSolver(cells1{i}) ;
    
    setEnvironment(cells2{i}, T, Nao, Cao, Ko) ;
    setUpPacingProtocol(cells2{i}, 0, 300, 999, 1)
    odeSolver(cells2{i}) ;
    
    setEnvironment(cells3{i}, T, Nao, Cao, Ko) ;
    setUpPacingProtocol(cells3{i}, 0, 300, 999, 1)
    odeSolver(cells3{i}) ;
    
    disp(int2str(i));
end
save big_pseudocells_0pt1_run.mat cells1
save big_pseudocells_0pt2_run.mat cells2
save big_pseudocells_0pt3_run.mat cells3
%
freqs = zeros(3,popsize);
for i = 1:popsize
    [~,~,features1] = APextract_custom(cells1{i}.state.t, cells1{i}.state.Y(:,1), cells1{i}.protocol.stimtimes);
    [~,~,features2] = APextract_custom(cells2{i}.state.t, cells2{i}.state.Y(:,1), cells2{i}.protocol.stimtimes);
    [~,~,features3] = APextract_custom(cells3{i}.state.t, cells3{i}.state.Y(:,1), cells3{i}.protocol.stimtimes);
    freqs(:,i) = [features1(4); features2(4); features3(4)];
end
save beatfreqs.mat freqs
%%
h = figure;
freqs(isnan(freqs)) = 0;

subplot(3,1,1)
a = histogram(freqs(1,:), 'NumBins', 25);
subplot(3,1,2)
b = histogram(freqs(2,:), 'NumBins', 25);
subplot(3,1,3)
c = histogram(freqs(3,:), 'NumBins', 25);
%print('beatfreq','-dpng','-r0')