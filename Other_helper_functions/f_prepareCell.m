function [k19] = f_prepareCell(protocol_number, cellID, init)
T = 310.0 ;         % (K)
Nao = 151.0 ;       % (mM)
Ko = 5.4 ;          % (mM)
Cao = 1.8 ;         % (mM)

amplitudes = [60, 0, 60, 0, 60, 0, 60, 0, 60, 0,...
              60, 0, 60, 0, 60, 0, 60, 60, 0, ...
              0, 60, 0, 60, 0, 60, ... % 20-25 
              60, 60, 60, 60, 60, 60, ... % 26-31 
              60, 60, 60, 60, 60, 60, 60] ; % 32-38
numPulses = [125, 100, 125, 100, 125, 100, 125, 100, 125, 100,...
    125, 100, 125, 100, 125, 100, 167, 125, 100, ...
    100, 100, 100, 100, 100, 100, ...
    200, 125, 50, 200, 125, 50, ...
    100, 100, 100, 50, 50, 100, 100]; 
numPulses = numPulses .* 3 ;
precedingTime = [799, 999, 799, 999, 799, 999, 799, 999, 799, 999,...
    799, 999, 799, 999, 799, 999, 599, 799, 999, ...
    999, 999, 999, 999, 999, 999, ...
    499, 799, 1999, 499, 799, 1999, ...
    999, 999, 999, 1999, 1999, 999, 999];
pulseDurations = ones(1, 38); 

% Use default initial parameters if none supplied, otherwise use supplied params
if ~exist('init','var')
    k19 = gaKernik19([-75.7444536163477,0.338969592726776,0.000203113729306209,7.16928093750999,104.748824394112,0,0.000386686486786781,0.165948056753057,0.927145173320106,0.321984775889061,0.452222061313948,0.157787210225653,0.743251868606151,0.121059208476135,0.0292207648412020,0.00620538308203812,0.736108314102295,0.000264118925707198,0.000263623380304203,0.746780367359425,0.0122283684412607,0.000154416932298621,0.0123158737520428,0], cellID) ;
else
    k19 = gaKernik19(init, cellID) ; 
end

setUpPacingProtocol(k19, amplitudes(protocol_number), numPulses(protocol_number),...
    precedingTime(protocol_number), pulseDurations(protocol_number)) ;
setEnvironment(k19, T, Nao, Cao, Ko) ;  % Baseline environment

% Caohigh paced & spont
if protocol_number == 1 || protocol_number == 2
    setEnvironment(k19, T, Nao, 2.6, Ko) ;
end

% Caolow paced & spont
if protocol_number == 3 || protocol_number == 4
    setEnvironment(k19, T, Nao, 1.0, Ko) ;
end

% ICaL block 25% paced & spont
if protocol_number == 5 || protocol_number == 6
    icalblock = ones(1,16);
    icalblock(3) = 0.75;
    setUpDrugApplication(k19, icalblock, zeros(1,16), ones(1,16)*300000)
end

% ICaL block 50% paced & spont
if protocol_number == 7 || protocol_number == 8
    icalblock = ones(1,16);
    icalblock(3) = 0.5;
    setUpDrugApplication(k19, icalblock, zeros(1,16), ones(1,16)*300000)
end

% IKr block 15% paced & spont
if protocol_number == 9 || protocol_number == 10
    ikrblock = ones(1,16);
    ikrblock(6) = 0.85;
    setUpDrugApplication(k19, ikrblock, zeros(1,16), ones(1,16)*300000)
end

% IKr block 30% paced & spont
if protocol_number == 11 || protocol_number == 12
    ikrblock = ones(1,16);
    ikrblock(6) = 0.7;
    setUpDrugApplication(k19, ikrblock, zeros(1,16), ones(1,16)*300000)
end

% Kohigh paced & spont
if protocol_number == 13 || protocol_number == 14
    setEnvironment(k19, T, Nao, Cao, 5.8) ;
end

% Kolow paced & spont
if protocol_number == 15 || protocol_number == 16
    setEnvironment(k19, T, Nao, Cao, 5.0) ;
end

% % Nao 147.0 mM (experimental data)
% if protocol_number >= 20
%     Nao = 147.0 ;
%     setEnvironment(k19, T, Nao, Cao, Ko) ;
% end

if protocol_number == 22 || protocol_number == 23 % 40% IKr block (1nM)
    ikrblock = ones(1,16);
    ikrblock(6) = 0.6;
    setUpDrugApplication(k19, ikrblock, zeros(1,16), ones(1,16)*300000)
end

if protocol_number == 24 || protocol_number == 25 % Hypocalcemia
    setEnvironment(k19, T, Nao, 1.0, Ko) ;
end

% if protocol_number == 26 || protocol_number == 27 || protocol_number == 28
%     setEnvironment(k19, T, Nao, 1.8, Ko) ;
% end
if protocol_number == 29 || protocol_number == 30 || protocol_number == 31 %Hypocalcemia
    setEnvironment(k19, T, Nao, 1.0, Ko) ;
end

if protocol_number == 32 || protocol_number == 36 % Hyponatremia 70%
    setEnvironment(k19, T, 105.7, Cao, Ko) ; 
end
if protocol_number == 34 % 25% ICaL block
    icalblock = ones(1,16);
    icalblock(3) = 0.5;
    setUpDrugApplication(k19, icalblock, zeros(1,16), ones(1,16)*300000)    
end
if protocol_number == 35 % 40% IKr block?; trying 10% for 1nM
    ikrblock = ones(1,16);
    ikrblock(6) = 0.9;
    setUpDrugApplication(k19, ikrblock, zeros(1,16), ones(1,16)*300000)    
end

if protocol_number == 37 % 70% ICaL block?
    icalblock = ones(1,16);
    icalblock(3) = 0.3;
    setUpDrugApplication(k19, icalblock, zeros(1,16), ones(1,16)*300000) 
end
if protocol_number == 38 % 90% ICaL block?
    icalblock = ones(1,16);
    icalblock(3) = 0.1;
    setUpDrugApplication(k19, icalblock, zeros(1,16), ones(1,16)*300000)    
end

