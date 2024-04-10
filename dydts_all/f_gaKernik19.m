function [dY, currents] = f_gaKernik19(time, Y, universals, geometry, parameters, conductances, drugEffects, applicationTimes)
%% set universals
% physical constants & temperature
F = universals(1) ;          
R = universals(2) ; 
T = universals(3) ;
% ion concentrations
Nao = universals(4) ;
Ko = universals(5) ;
Cao = universals(6) ;
%--------------------------------------------------------------------------
%% Geometry
%--------------------------------------------------------------------------
V_SR = geometry(5) ;
Vc = geometry(6) ;
Cm  = geometry(7) ;
%--------------------------------------------------------------------------
%% Conductances
%--------------------------------------------------------------------------
c = conductances ;
for i = 1: length(c)
    if applicationTimes(1,i) <= time && applicationTimes(2,i) > time
        c(i) = c(i) .* drugEffects(i) ;
    end
end
g_Na = c(1) ;
g_f = c(2) ; 
p_CaL = c(3) ; 
g_to = c(4) ; 
g_Ks = c(5) ; 
g_Kr = c(6) ; 
g_K1 = c(7) ;
g_PCa = c(8) ;
g_b_Na = c(9) ;
g_b_Ca = c(10) ;
VmaxUp = c(11) ;
g_irel_max = c(12) ;
kNaCa = c(13) ;
PNaK = c(14) ;
V_leak = c(15) ;
g_CaT = c(16) ; 

%--------------------------------------------------------------------------
%% Parameters
%--------------------------------------------------------------------------
p = parameters ;
% i_K1:
xK11		= p(1) ;	
xK12		= p(2) ; 	
xK13        = p(3) ;   
xK14        = p(4) ;   	
xK15        = p(5) ;  	

% i_Kr:
Xr1_1		= p(6) ;
Xr1_2		= p(7) ;
Xr1_5		= p(8) ;
Xr1_6		= p(9) ;
Xr2_1		= p(10) ;
Xr2_2		= p(11) ;
Xr2_5		= p(12) ;
Xr2_6		= p(13) ;
tauXr1_const= p(14) ;
tauXr2_const= p(15) ;

% i_Ks:
ks1         = p(16) ; 
ks2         = p(17) ; 
ks5         = p(18) ; 
ks6         = p(19) ;
tauks_const = p(20) ;

% i_to:
r1          = p(21) ; 
r2          = p(22) ; 
r5          = p(23) ; 
r6          = p(24) ;
s1          = p(25) ; 
s2          = p(26) ; 
s5          = p(27) ; 
s6          = p(28) ;
tau_r_const = p(29) ;
tau_s_const = p(30) ;

% i_CaL:
d1 = p(31) ;
d2 = p(32) ;
d5 = p(33) ;
d6 = p(34) ;
f1 = p(35) ;
f2 = p(36) ;
f5 = p(37) ; 
f6 = p(38) ;
taud_const = p(39) ;
tauf_const = p(40) ;

% i_Na:
m1 = p(41) ;
m2 = p(42) ;
m5 = p(43) ;
m6 = p(44) ;
h1 = p(45) ;
h2 = p(46) ;
h5 = p(47) ;
h6 = p(48) ;
j1 = p(49) ;
j2 = p(50) ;
tau_m_const = p(51) ;
tau_h_const = p(52) ;
tau_j_const = p(53) ;

% i_f:
xF1 = p(54) ; 
xF2 = p(55) ; 
xF5 = p(56) ; 
xF6 = p(57) ; 
xF_const = p(58) ; 

%%% test scale constants
taud_scale = p(59) ;
tauf_scale = p(60) ;
 
%% Initial conditions
%--------------------------------------------------------------------------
%{
SS original
Y=[     -0.070  0.32    0.0002  0    0    1     1     1      0      1      0   0.75  0.75  0   0.1    1    0    9.2    0     0.75    0.3     0.9     0.1];
YNames = {'Vm', 'Ca_SR', 'Cai', 'g', 'd', 'f1', 'f2', 'fCa', 'Xr1', 'Xr2', 'Xs', 'h', 'j', 'm', 'Xf', 'q', 'r', 'Nai', 'm_L', 'h_L', 'RyRa', 'RyRo', 'RyRc'};
YUnits = {'V',   'mM',   'mM',  '-', '-', '-',  '-',  '-',   '-',   '-',   '-',  '-', '-', '-', '-',  '-', '-', 'mM',   '-',   '-',    '-',    '-',    '-'};

SS at 800 by sim_14 play tau RyR ROUND

Y = [-0.0749041438812569 0.0936432074356392 3.81472732752446e-05 0 8.27445474910059e-05 0.739375895533045 0.999983842545123 0.997731257249225 0.268791271739935 0.434813785079876 0.0314716535458102 0.744090874504749 0.0755303461442611 0.0996946953195255 0.0246909878620017 0.841491713475389 0.00558580115137648 8.64815457347901 0.00226188558003069 0.0806819453003550 0.0386969233870937 0.0277238422359321 0.0726455543836155];
%}
%--------------------------------------------------------------------------
% State variables
%--------------------------------------------------------------------------
%{
1: Vm     (V)     (in Membrane)
2: Ca_SR  (mM)    (in calcium_dynamics)
3: Cai    (mM)    (in calcium_dynamics)
4: g      (-)     (NOT USED)                                                            
5: d      (-)     (in i_CaL_d_gate)
6: f1     (-)     (in i_CaL_f1_gate)
7: f2     (-)     (in i_CaL_f2_gate)
8: fCa    (-)     (in i_CaL_fCa_gate)
9: Xr1    (-)     (in i_Kr_Xr1_gate)
10: Xr2   (-)     (in i_Kr_Xr2_gate)
11: Xs    (-)     (in i_Ks_Xs_gate)
12: h     (-)     (in i_Na_h_gate)
13: j     (-)     (in i_Na_j_gate)
14: m     (-)     (in i_Na_m_gate)
15: Xf    (-)     (in i_f_Xf_gate)
16: q     (-)     (in i_to_q_gate)
17: r     (-)     (in i_to_r_gate)
18: Nai   (mM)    (in sodium_dynamics)
19: m_L   (-)     (in i_NaL_m_gate)
20: h_L   (-)     (in i_NaL_h_gate)
21: RyRa  (-)     (in calcium_dynamics)
22: RyRo  (-)     (in calcium_dynamics)
23: RyRc  (-)     (in calcium_dynamics)
%}
%--------------------------------------------------------------------------
% Intracellular concentrations - all change according to model except Ki
%--------------------------------------------------------------------------
% Naio = 10     mM Y(18)
% Ki = 150.0;   % mM (in model_parameters)
% Cai  = 0.0002 mM Y(3)
% caSR = 0.3    mM Y(2)

% time (second)

%% Nernst potential
E_Ca = 0.5 * R * T / (F/1000) * log(Cao / Y(3)) ; % mV
E_Na = R * T / (F/1000) * log(Nao / Y(4)) ;  % mV
E_K  = R * T / (F/1000) * log(Ko / Y(5)) ;  % mV

%% i_Na 
% parameter-dependent values:
m3 = m5 * m1 ; 
m4 = 1 / ((1 / m2) + (1 / m6)) ;
h3 = h5 * h1 ; 
h4 = 1 / ((1 / h2) + (1 / h6)) ;
j5 = h5 ; 
j6 = h6 ;
j3 = j5 * j1 ; 
j4 = 1 / ((1 / j2) + (1 / j6)) ;

% 13: h (dimensionless) (inactivation in i_Na)
alpha_h = h1 .* exp(Y(1) ./ h2) ;
beta_h  = h3 .* exp(Y(1) ./ h4) ;
h_inf   = (alpha_h ./ (alpha_h + beta_h)) ;
tau_h   = (1 ./ (alpha_h + beta_h)) + tau_h_const ;
dY(13)  = (h_inf - Y(13)) ./ tau_h ;

% 14: j (dimensionless) (slow inactivation in i_Na)
alpha_j = j1 .* exp(Y(1) ./ j2) ;
beta_j  = j3 .* exp(Y(1) ./ j4) ;
j_inf   = alpha_j ./ (alpha_j + beta_j) ;
tau_j   = (1 ./ (alpha_j + beta_j)) + tau_j_const ;
dY(14)  = (j_inf - Y(14)) ./ tau_j ;

% 15: m (dimensionless) (activation in i_Na)
alpha_m = m1 .* exp(Y(1) ./ m2) ;
beta_m  = m3 .* exp(Y(1) ./ m4) ;
m_inf   = alpha_m ./ (alpha_m + beta_m) ;
tau_m   = (1 ./ (alpha_m + beta_m)) + tau_m_const ;
dY(15)  = (m_inf - Y(15)) ./ tau_m ;

%Current:
i_Na    = g_Na * Y(15)^3.0 * Y(13) * Y(14) * (Y(1) - E_Na) ;

%% i_f 
% parameter-dependent values:
xF3 = xF5 * xF1 ;  
xF4 = 1 / ((1 / xF2) + (1 / xF6)) ;

% 16: Xf (dimensionless) (inactivation in i_f)
alpha_Xf = xF1 .* exp(Y(1) ./ xF2) ;
beta_Xf = xF3 .* exp(Y(1) ./ xF4) ;
Xf_inf = alpha_Xf ./ (alpha_Xf + beta_Xf) ;
tau_Xf = (1 ./ (alpha_Xf + beta_Xf)) + xF_const ;
dY(16) = (Xf_inf - Y(16)) ./ tau_Xf ;

% Current:
NatoK_ratio = 0.491 ; % Verkerk et al. 2013
Na_frac = NatoK_ratio ./ (NatoK_ratio + 1) ;
i_fNa = Na_frac * g_f * Y(16) * (Y(1) - E_Na) ;
i_fK = (1 - Na_frac) * g_f * Y(16) * (Y(1) - E_K) ;
i_f = i_fNa + i_fK ;

%% i_CaL 
% parameter-dependent values:
d3 = d5 * d1 ;
d4 = 1 / ((1 / d2) + (1 / d6)) ;
f3 = f5 * f1 ;
f4 = 1 / ((1 / f2) + (1 / f6)) ;

% 7: d (dimensionless) (activation in i_CaL)
alpha_d = d1 .* exp(Y(1) ./ d2) ;
beta_d = d3 .* exp(Y(1) ./ d4) ;
d_inf = alpha_d ./ (alpha_d + beta_d) ;
tau_d = taud_scale .* ((1 ./ (alpha_d + beta_d)) + taud_const) ;
dY(7) = (d_inf - Y(7)) / tau_d ;

% 8: f (dimensionless) (inactivation  i_CaL)
alpha_f = f1 .* exp(Y(1) ./ f2) ;
beta_f = f3 .* exp(Y(1) ./ f4) ;
f_inf = alpha_f ./ (alpha_f + beta_f) ;
tau_f = tauf_scale .* ((1 ./ (alpha_f + beta_f)) + tauf_const) ;
dY(8) = (f_inf - Y(8)) / tau_f ;

% 9: fCa (dimensionless) (calcium-dependent inactivation in i_CaL)
% from Ten tusscher 2004
scale_Ical_Fca_Cadep = 1.2 ;
alpha_fCa = 1.0 / (1.0 + ((scale_Ical_Fca_Cadep .* Y(3)) / 0.000325)^8.0) ;
beta_fCa = 0.1 / (1.0 + exp((scale_Ical_Fca_Cadep .* Y(3) - 0.0005) / 0.0001)) ;
gamma_fCa = 0.2 / (1.0 + exp((scale_Ical_Fca_Cadep .* Y(3) - 0.00075) / 0.0008)) ;

fCa_inf = (alpha_fCa + beta_fCa + gamma_fCa + 0.23) / 1.46 ;
tau_fCa = 2 ; % ms
if fCa_inf > Y(9) && Y(1) > -60
    k_fca = 0 ;
else
    k_fca = 1 ;
end
dY(9) = k_fca .* (fCa_inf - Y(9)) / tau_fCa ;

% Current:
p_CaL_shannonCa = 5.4e-4 ;
p_CaL_shannonNa = 1.5e-8 ;
p_CaL_shannonK = 2.7e-7 ;
p_CaL_shannonTot = p_CaL_shannonCa + p_CaL_shannonNa + p_CaL_shannonK ;
p_CaL_shannonCap = p_CaL_shannonCa / p_CaL_shannonTot ;
p_CaL_shannonNap = p_CaL_shannonNa / p_CaL_shannonTot ;
p_CaL_shannonKp = p_CaL_shannonK / p_CaL_shannonTot ;

p_CaL_Ca = p_CaL_shannonCap * p_CaL ;
p_CaL_Na = p_CaL_shannonNap * p_CaL ;
p_CaL_K = p_CaL_shannonKp * p_CaL ;

ibarca = p_CaL_Ca * 4.0 * Y(1) * (F/1000)^2.0 / (R * T) * (0.341 * Y(3) * exp(2.0 * Y(1) * (F/1000) / (R * T)) - 0.341 * Cao) / (exp(2.0 * Y(1) * (F/1000) / (R * T)) - 1.0) ;
i_CaL_Ca = ibarca * Y(7) * Y(8) * Y(9) ;

ibarna = p_CaL_Na * Y(1) * (F/1000)^2.0 / (R * T) * (0.75 * Y(4) * exp(Y(1) * (F/1000) / (R * T)) - 0.75 * Nao) / (exp(Y(1) * (F/1000) / (R * T)) - 1.0) ;
i_CaL_Na = ibarna * Y(7) * Y(8) * Y(9) ;

ibark = p_CaL_K * Y(1) * (F/1000)^2.0 / (R * T) * (0.75 * Y(5) * exp(Y(1) * (F/1000) / (R * T)) - 0.75 * Ko) / (exp(Y(1) * (F/1000) / (R * T)) - 1.0) ;
i_CaL_K = ibark * Y(7) * Y(8) * Y(9) ;

i_CaL = i_CaL_Ca + i_CaL_Na + i_CaL_K ;

%% i_CaT
% SAN T-TYPE CA2+ model (Demir et al., Maltsev-Lakatta ), G_CaT determined by fit to Kurokawa IV:

% 19: dCaT (activation in i_CaT)
dcat_inf = 1 ./ (1 + exp(-(Y(1) + 26.3) ./ 6)) ;
tau_dcat = 1 ./ (1.068 * exp((Y(1) + 26.3) ./ 30) + 1.068 * exp(-(Y(1) + 26.3) ./ 30)) ;
dY(19) = (dcat_inf - Y(19)) / tau_dcat ;

% 20: fCaT (inactivation in i_CaT)
fcat_inf = 1 ./ (1 + exp((Y(1) + 61.7) ./ 5.6)) ;
tau_fcat = 1 ./ (0.0153 * exp(-(Y(1) + 61.7) ./ 83.3) + 0.015 * exp((Y(1) + 61.7) ./ 15.38)) ;
dY(20) = (fcat_inf - Y(20)) / tau_fcat ;

i_CaT = g_CaT * (Y(1) - E_Ca) * Y(19) * Y(20) ;

%% i_to 
% parameter-dependent values:
r3 = r5 * r1 ; 
r4 = 1 / ((1 / r2) + (1 / r6)) ;
s3 = s5 * s1 ; 
s4 = 1 / ((1 / s2) + (1 / s6)) ;

% 17: s (dimensionless) (inactivation in i_to)
alpha_s = s1 .* exp(Y(1) ./ s2) ;
beta_s = s3 .* exp(Y(1) ./ s4) ;
s_inf = alpha_s ./ (alpha_s + beta_s) ;
tau_s = (1 ./ (alpha_s + beta_s)) + tau_s_const ;
dY(17) = (s_inf - Y(17)) ./ tau_s ;

% 18: r (dimensionless) (activation in i_to)
alpha_r = r1 .* exp(Y(1) ./ r2) ;
beta_r = r3 .* exp(Y(1) ./ r4) ;
r_inf = alpha_r ./ (alpha_r + beta_r) ;
tau_r = (1 ./ (alpha_r + beta_r)) + tau_r_const ;
dY(18) = (r_inf - Y(18)) ./ tau_r ;

% Current:
i_to = g_to * (Y(1) - E_K) * Y(17) * Y(18) ;

%% i_Ks 
% parameter-dependent values:
ks3 = ks5 * ks1 ; 
ks4 = 1 / ((1 / ks2) + (1 / ks6)) ;

% 12: Xs (dimensionless) (activation in i_Ks)
alpha_Xs = ks1 .* exp(Y(1) ./ ks2) ;
beta_Xs = ks3 .* exp(Y(1) ./ ks4) ;
Xs_inf = alpha_Xs ./ (alpha_Xs + beta_Xs) ;
tau_Xs = (1 ./ (alpha_Xs + beta_Xs) ) + tauks_const ;
dY(12) = (Xs_inf - Y(12)) ./ tau_Xs ;

% Current:
i_Ks = g_Ks * (Y(1) - E_K) * (Y(12).^2) ;

%% i_Kr 
% parameter-dependent values:
Xr1_3 = Xr1_5 * Xr1_1 ; 
Xr2_3 = Xr2_5 * Xr2_1 ;
Xr1_4 = 1 / ((1 / Xr1_2) + (1 / Xr1_6)) ; 
Xr2_4 = 1 / ((1 / Xr2_2) + (1 / Xr2_6)) ;

% 10: Xr1 (dimensionless) (activation in i_Kr_Xr1)
alpha_Xr1 = Xr1_1 .* exp(Y(1) ./ Xr1_2) ;
beta_Xr1 = Xr1_3 .* exp(Y(1) ./ Xr1_4) ;
Xr1_inf = alpha_Xr1 ./ (alpha_Xr1 + beta_Xr1) ;
tau_Xr1 = (1 ./ (alpha_Xr1 + beta_Xr1)) + tauXr1_const ; 
dY(10) = (Xr1_inf - Y(10)) ./ tau_Xr1 ;

% 11: Xr2 (dimensionless) (inactivation in i_Kr_Xr2)
alpha_Xr2 = Xr2_1 .* exp(Y(1) ./ Xr2_2) ;
beta_Xr2 = Xr2_3 .* exp(Y(1) ./ Xr2_4) ;
Xr2_inf = alpha_Xr2 ./ (alpha_Xr2 + beta_Xr2) ;
tau_Xr2 = (1 ./ (alpha_Xr2 + beta_Xr2)) + tauXr2_const ;
dY(11) = (Xr2_inf - Y(11)) ./ tau_Xr2 ;

% Current:
i_Kr = g_Kr * (Y(1)-E_K) * Y(10) * Y(11) * sqrt(Ko/5.4) ;

%% i_K1 
alpha_xK1 = xK11 .* exp((Y(1) + xK13) ./ xK12) ;
beta_xK1 = exp((Y(1) + xK15) ./ xK14) ;
XK1_inf = alpha_xK1 ./ (alpha_xK1 + beta_xK1) ;

i_K1 = g_K1 * XK1_inf * (Y(1) - E_K) * sqrt(Ko/5.4) ;

%% i_NaCa 
% Ten Tusscher formulation
KmCa  = 1.38 ;      % Cai half-saturation constant millimolar (in i_NaCa)
KmNai = 87.5 ;      % Nai half-saturation constnat millimolar (in i_NaCa)
Ksat  = 0.1 ;       % saturation factor dimensionless (in i_NaCa)
gamma = 0.35 * 2 ;	% voltage dependence parameter dimensionless (in i_NaCa)
alpha = 2.5 * 1.1 ; % factor to enhance outward nature of inaca dimensionless (in i_NaCa)

i_NaCa = kNaCa * ((exp(gamma * Y(1) * (F / 1000) / (R * T)) * (Y(4)^3.0) * Cao) - (exp((gamma - 1.0) * Y(1) * (F / 1000) / (R * T)) * (Nao^3.0) * Y(3) * alpha)) / (((KmNai^3.0) + (Nao^3.0)) * (KmCa + Cao) * (1.0 + Ksat * exp((gamma - 1.0) * Y(1) * (F / 1000) / (R * T)))) ;

%% i_NaK 
% Ten Tusscher formulation
Km_K = 1.0 ;   % Ko half-saturation constant millimolar (in i_NaK)
Km_Na = 40.0 ;   %  Nai half-saturation constant millimolar (in i_NaK)

i_NaK = PNaK * ((Ko * Y(4)) / ((Ko + Km_K) * (Y(4) + Km_Na) * (1.0 + 0.1245 * exp(-0.1 * Y(1) * (F / 1000) / (R * T)) + 0.0353 * exp(-Y(1) * (F / 1000) / (R * T))))) ;

%% ------------------------- Sarcoplasmic reticulum ----------------------------
%% i_up (J_up/SR Uptake/SERCA) 
% Ten Tusscher formulation
Kup = 0.00025 * 0.702 ;   % millimolar (in calcium_dynamics)
i_up = VmaxUp / (1.0 + Kup^2.0 / Y(3)^2.0) ;

%% i_leak (J_leak/SR Leak) 
% Ten Tusscher formulation
i_leak = (Y(2) - Y(3)) * V_leak ;

%% i_rel (J_rel/SR Release/RYR) 
% re-fit parameters. scaled to account for differences in calcium concentration in
% cleft (cleft is used in shannon-bers model geometry, not in this model geometry)
koCa = 56320 * 11.43025 ;             % [mM^-2 1/ms]  
kiCa = 54 * 0.3425 ;                  % [1/mM/ms]
kom = 1.5 * 0.1429 ;                  % [1/ms]
kim = 0.001 * 0.5571 ;                % [1/ms]
ec50SR = 0.45 ;
MaxSR = 15 ;
MinSR = 1 ;

kCaSR = MaxSR - (MaxSR - MinSR) / (1 + (ec50SR / Y(2))^2.5) ;
koSRCa = koCa / kCaSR ;
kiSRCa = kiCa * kCaSR ;
RI = 1 - Y(21) - Y(22) - Y(23) ;

dY(21) = (kim * RI - kiSRCa * Y(3) * Y(21)) - (koSRCa * Y(3)^2 * Y(21) - kom * Y(22)) ;   % R
dY(22) = (koSRCa * Y(3)^2 * Y(21) - kom * Y(22)) - (kiSRCa * Y(3) * Y(22) - kim * Y(23)) ; % O
dY(23) = (kiSRCa * Y(3) * Y(22) - kim * Y(23)) - (kom * Y(23) - koSRCa * Y(3)^2 * RI) ;   % I

i_rel = g_irel_max * Y(22) * (Y(2) - Y(3)) * (V_SR / Vc) ; % g_irel_max previously called ks in original Kernik19, changed to match other models
%% -----------------------------------------------------------------------------

%% i_PCa (PMCA/Ca-ATPase/Ca Sarcolemmal Pump) 
% Ten Tusscher formulation
KPCa = 0.0005 ;   % millimolar (in i_PCa)
i_PCa = g_PCa * Y(3) / (Y(3) + KPCa) ;

%% i_b_Na & i_b_Ca (Background currents) 
% Ten Tusscher formulation
i_b_Na = g_b_Na * (Y(1) - E_Na) ;

i_b_Ca = g_b_Ca * (Y(1) - E_Ca) ;

%% CaSR, Cai, Nai, Ki (Ionic concentrations) 
% 2: CaSR (millimolar)
% rapid equilibrium approximation equations -- NOT as formulated in ten Tusscher 2004 text
Buf_SR = 10.0 * 1.2 ; % millimolar (in calcium_dynamics)
Kbuf_SR = 0.3 ; % millimolar (in calcium_dynamics)
Ca_SR_bufSR = 1 / (1.0 + Buf_SR * Kbuf_SR / (Y(2) + Kbuf_SR)^2.0) ;

dY(2) = Ca_SR_bufSR * Vc / V_SR * (i_up - (i_rel + i_leak)) ;

% 3: Cai (millimolar)
%rapid equilibrium approximation equations -- NOT as formulated in ten Tusscher 2004 text
Buf_C = 0.06 ; % millimolar (in calcium_dynamics)
Kbuf_C = 0.0006 ; % millimolar (in calcium_dynamics)
Cai_bufc = 1 / (1.0 + Buf_C * Kbuf_C / (Y(3) + Kbuf_C)^2.0) ;

dY(3) = Cai_bufc * (i_leak - i_up + i_rel - dY(6) - (i_CaL_Ca + i_CaT + i_b_Ca + i_PCa - 2*i_NaCa) * Cm / (2.0 * Vc * (F / 1000))) ;

% 4: Nai (millimolar) (in sodium_dynamics)
dY(4) = -Cm * (i_Na + i_b_Na + i_fNa + 3.0 * i_NaK + 3.0 * i_NaCa + i_CaL_Na) / ((F / 1000) * Vc) ;

% 5: Ki (millimolar) (in potatssium_dynamics)
dY(5) = -Cm * (i_K1 + i_to + i_Kr + i_Ks + i_fK - 2 .* i_NaK + i_CaL_K ) / ((F / 1000) * Vc) ;

%% Stimulation
i_stim = Y(end) ; 
dY(24) = 0 ;

%% Membrane potential
% I_stim:
% if (time >= i_stim_Start) && (time <= i_stim_End) && (mod(time-i_stim_Start-100, cyclelength)<i_stim_PulseDuration)
%     i_stim = stim_flag*i_stim_Amplitude;
% else
%     i_stim = 0.0;
% end

% %Voltage Clamp:
% voltageclamp = 0 ; % flag, want to protocolize this
% v_clamp_step=0;
% v_clamp_rest=-65;
% steplength=100;
% R_clamp = 0.02;
% if voltageclamp == 0 
%     v_clamp = Y(1) ; %set i_voltageclamp to 0    
% elseif voltageclamp == 1    % train of square pulse:
%     if mod(time, cyclelength) < cyclelength - steplength
%         v_clamp = v_clamp_rest ;
%     else
%         v_clamp = v_clamp_step ;
%     end
% end
%     
% i_voltageclamp = (v_clamp - Y(1)) / R_clamp ;
i_voltageclamp = 0;
dY(1) = -(i_K1 + i_to + i_Kr + i_Ks + i_CaL + i_CaT + i_NaK + i_Na + i_NaCa + i_PCa + i_f + i_b_Na + i_b_Ca - i_stim - i_voltageclamp) ;
currents = [i_K1, i_to, i_Kr, i_Ks, i_CaL, i_NaK, i_Na, i_NaCa, i_PCa, i_f, i_b_Na, i_b_Ca, i_rel, i_up, i_leak, i_stim, i_CaT];
dY = dY' ;
%===============================================================================
% End of file
%===============================================================================
