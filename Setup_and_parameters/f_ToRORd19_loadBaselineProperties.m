% ToRORd19 baseline properties to be loaded by the constructor - keeps them from
% cluttering the main class space
function [universals,...
    name, model, isCellAlive,...
    geometry,...
    state, YNames, YUnits, currentNames, currentUnits,... 
    parameterNames, parameterUnits, parameters,...
    conductanceNames, conductances,...
    protocol,...
    event_number, y_allevents, x_allevents, x_conductances, isPaced]...
    = f_ToRORd19_loadBaselineProperties(initY)
%% universals
universals = struct('R', 8314.0, 'T', 310.0, 'F', 96485.0,...
            'nao', 140, 'cao', 1.8, 'ko', 5);
        
%% model
name = 'ToRORd19CM';
model = @f_ToRORd19 ; 
isCellAlive = 1 ; 
isPaced = 0 ; 

%% geometry
L = 0.01 ; %cm
rad = 0.0011 ;
vcell = 1000 * pi * rad * rad * L ;
Ageo = 2 * pi * rad * rad + 2 * pi * rad * L ; %full cell surface area
Acap = 2 * Ageo ;
vmyo = 0.68 * vcell ;
vnsr = 0.0552 * vcell ;
vjsr = 0.0048 * vcell ;
vss = 0.02 * vcell ;
geometry = struct('L',L, 'rad',rad, 'vcell',vcell, 'Ageo',Ageo, 'Acap',Acap,...
            'vmyo',vmyo, 'vnsr',vnsr, 'vjsr',vjsr, 'vss',vss) ;
        
%% state
state = struct('init_Y', initY, 'Y', [], 't', [], 'currents', []) ;
YNames = {'v', 'nai', 'nass', 'ki', 'kss',...5 
            'cai', 'cass', 'cansr', 'cajsr', 'm',...10 
            'hp', 'h', 'j', 'jp', 'mL',...15 
            'hL', 'hLp', 'a', 'iF', 'iS',...20 
            'ap', 'iFp', 'iSp', 'd', 'ff',... 25
            'fs', 'fcaf', 'fcas', 'jca', 'nca',...30 
            'nca_i', 'ffp', 'fcafp', 'xs1', 'xs2',...35 
            'Jrel_np', 'CaMKt', 'ikr_c0', 'ikr_c1', 'ikr_c2',...40 
            'ikr_o', 'ikr_i', 'Jrel_p', 'Istim'};% 44
YUnits = {};
    
currentNames = {'INa', 'INaL', 'Ito', 'ICaL', 'IKr',...5
    'IKs', 'IK1', 'INaCa_i', 'INaCa_ss', 'INaK',...10
    'IKb', 'INab', 'ICab', 'IpCa', 'Jdiff',...15
    'JdiffNa', 'JdiffK', 'Jup', 'Jleak', 'Jtr',...20
    'Jrel', 'CaMKa', 'Istim', 'fINap', 'fINaLp',...25
    'fICaLp', 'fJrelp', 'fJupp', 'cajsr', 'cansr',...30
    'PhiCaL_ss', 'v', 'ICaL_i', 'I_ClCa', 'I_Clbk',...35
    'ICaL_tot', 'INaCa'}; % 39
currentUnits = {};

%% parameters
ICaL_fractionSS = 0.8 ;
INaCa_fractionSS = 0.35 ;%0.35

parameterNames = {'ICaL_fractionSS', 'INaCa_fractionSS'} ;
parameterUnits = {'-', '-'} ;
parameters = struct('baseline', [ICaL_fractionSS, INaCa_fractionSS], 'scaling', ones(1, length(parameterNames))) ;

%% conductances
conductanceNames = {...
    'INa_Multiplier',...    INa
    'ICaL_Multiplier',...   ICaL
    'Ito_Multiplier',...    Ito
    'INaL_Multiplier',...   INaL
    'IKr_Multiplier',...    IKr
    'IKs_Multiplier',...    IKs
    'IK1_Multiplier',...    IK1
    'IKb_Multiplier',...    IKb
    'INaCa_Multiplier',...  INaCa
    'INaK_Multiplier',...   INaK
    'INab_Multiplier',...   INab
    'ICab_Multiplier',...   ICab
    'IpCa_Multiplier',...   IpCa
    'ICaCl_Multiplier',...  ICaCl
    'IClb_Multiplier',...   IClb
    'Jrel_Multiplier',...   Jrel
    'Jup_Multiplier'};    % Jup 
conductances = struct('baseline', ones(1, length(conductanceNames)), ...
    'scaling', ones(1, length(conductanceNames)),...
    'drugEffects', ones(1,length(conductanceNames)),...
    'applicationTimes', zeros(2,length(conductanceNames))) ;

%% protocol
protocol = struct('intervalTimes', [] , 'stimulus', [], 'amplitudes', [],...
            'numPulses', [], 'precedingTime', [], 'pulseDuration', [], 'totalTime', [],...
            'phaseLengths', [], 'frequencies', [], 'last3stimTimes', []) ;
        
%% ga related
event_number = []; % for autopsy
y_allevents = []; % for autopsy
x_allevents = []; % for autopsy
x_conductances = [];


end