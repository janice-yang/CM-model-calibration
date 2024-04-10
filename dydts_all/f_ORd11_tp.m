function [dYdt,data] = f_ORd11_tp(time, Y, geometry, parameters, conductances, celltype) 

global R T F Nao Ko Cao 

%% From geometry structure

% Volume parameters 
L    = geometry.L;      %cm   cell length
r    = geometry.r;      %cm   radius 
R_CG = geometry.R_CG; 

% Areas (cm^2)
A_geo = geometry.A_geo; 
A_cap = geometry.A_cap; 

% Volumes (cm^3 = mL) 
v_cell = geometry.v_cell; 
v_myo  = geometry.v_myo; 
v_nsr  = geometry.v_nsr; 
v_jsr  = geometry.v_jsr; 
v_ss   = geometry.v_ss; 

% Membrane capacitance 
Cm   = geometry.Cm;       %microF

% Charges
z_Na = geometry.z_Na; 
z_K  = geometry.z_K; 
z_Ca = geometry.z_Ca; 

%% From parameters structure

pars = parameters.baseline .* parameters.scaling; 

c = conductances.baseline .* conductances.scaling; 

for i = 1:11
    if conductances.applicationTimes(1,i) <= time && conductances.applicationTimes(2,i) > time
        c(i) = c(i) .* conductances.drugEffects(i) ;
    end
end

%Maximal conductances (mS microF^(-1))
Gbar_Naf  = c(1); 
Gbar_Nal  = c(2); 
Gbar_to   = c(3);
Gbar_Kr   = c(4); 
Gbar_Ks   = c(5);
Gbar_K1   = c(6);
Gbar_NaCa = c(7); %(microA microF^(-1))
Gbar_Kb   = c(8);
Gbar_pCa  = c(9);
G_RyR     = c(10);
G_SERCA   = c(11);
G_NaK     = c(12);
G_up_epi  = c(13);

%Permeability ratio
PR_NaK = c(14);

%Permeabilities (cm s^(-1))
P_Ca       = c(15);
P_CaNa     = c(16);
P_CaK      = c(17);
P_CaCaMK   = c(18); 
P_CaNaCaMK = c(19); 
P_CaKCaMK  = c(20); 
P_Nab      = c(21); 
P_Cab      = c(22); 

%Activity coefficients (dimensionless)
gamma_Nai = pars(1); 
gamma_Nao = pars(2); 
gamma_Ki  = pars(3); 
gamma_Ko  = pars(4); 
gamma_Cai = pars(5); 
gamma_Cao = pars(6);

%Average concentrations (mM)
MgADP   = pars(7);
MgATP   = pars(8);
H       = pars(9); 
SigmaP  = pars(10);
cmdnbar = pars(11); 
trpnbar = pars(12); 
bsrbar  = pars(13); 
bslbar  = pars(14); 
csqnbar = pars(15); 

%Half-saturation values (mM)
K_mCaM         = pars(16); 
K_mCaMK_active = pars(17); %(dimensionless)
K_mCass        = pars(18); 
K_mCaact       = pars(19); 
K_mcmdn        = pars(20);
K_mtrpn        = pars(21); 
K_mbsr         = pars(22); 
K_mbsl         = pars(23); 
K_mcsqn        = pars(24);

%Time-scales (ms)
tau_diffNa = pars(25);
tau_diffK  = pars(26);
tau_diffCa = pars(27);
tau_tr     = pars(28); 

%Percentages (dimensionless)
A_hf     = pars(29); 
A_hs     = pars(30); 
A_ff     = pars(31); 
A_fs     = pars(32); 
A_fCaMKf = pars(33); 
A_fCaMKs = pars(34); 

%CaMK parameters
alpha_CaMK = pars(35);    %(ms^(-1))
beta_CaMK  = pars(36);    %(ms^(-1))
CaMKo     = pars(37);    %(dimensionless)  

%I_NaCa parameters
k_Na1      = pars(38);    %(mM)
k_Na2      = pars(39);    %(mM)
k_Na3      = pars(40);    %(mM)
k_asymm    = pars(41);    %(dimensionless)
omega_Na   = pars(42);    %(Hz)
omega_Ca   = pars(43);    %(Hz)
omega_NaCa = pars(44);    %(Hz)
k_Caon     = pars(45);    %(mM ms^(-1))
k_Caoff    = pars(46);    %(Hz)
q_Na       = pars(47);    %(dimensionless)
q_Ca       = pars(48);    %(dimensionless)

%I_NaK parameters
k_1plus = pars(49);       %(Hz)
k_1neg  = pars(50);       %(mM^(-1))
k_2plus = pars(51);       %(Hz)
k_2neg  = pars(52);       %(Hz)
k_3plus = pars(53);       %(Hz)
k_3neg  = pars(54);       %(Hz mM^(-2))
k_4plus = pars(55);       %(Hz)
k_4neg  = pars(56);       %(Hz)
Delta   = pars(57);       %(dimensionless)                                      

%I_NaK dissociation constants (mM)
KoNai  = pars(58);        
KoNao  = pars(59);
K_Ki    = pars(60); 
K_Ko    = pars(61);
K_MgATP = pars(62);
K_HP    = pars(63); 
K_NaP   = pars(64); 
K_KP    = pars(65); 

%J_rel parameters (ms)
beta_tau      = pars(66); 
beta_tauCaMK  = pars(67); 
alpha_rel     = pars(68); 
alpha_relCaMK = pars(69);

%J_up parameters 
DeltaKbar_mPLB   = pars(70);   %(ms)
DeltaJbar_upCaMK = pars(71);   %(flux units?)

%% State variables 

% Membrane voltage
Vm = Y(1); 

% Concentrations
Na_i   = Y(2); 
Na_ss  = Y(3); 
K_i    = Y(4);
K_ss   = Y(5); 
Ca_i   = Y(6); 
Ca_ss  = Y(7); 
Ca_nsr = Y(8); 
Ca_jsr = Y(9); 

% I_Na gating variables 
m       = Y(10);    %Activation for fast I_Na
h_f     = Y(11);    %Fast inactivation for fast I_Na
h_s     = Y(12);    %Slow inactivation for fast I_Na
j       = Y(13);    %Recovery from inactivation for fast I_Na
h_CaMKs = Y(14);    %Slow inactivation for CaMK phosphorylated fast I_Na
j_CaMK  = Y(15);    %Recovery from inactivation for CaMK phosphorylated fast I_Na 
m_l     = Y(16);    %Activation for late I_Na
h_l     = Y(17);    %Inactivation for late I_Na
h_lCaMK = Y(18);    %Inactivation for CaMK phosphorylated late I_Na

% I_to gating variables 
a       = Y(19);    %Activation for I_to 
i_f     = Y(20);    %Fast inactivation for I_to
i_s     = Y(21);    %Slow inactivation for I_to
a_CaMK  = Y(22);    %Activation for CaMK phosphorylated I_to
i_CaMKf = Y(23);    %Fast inactivation for CaMK phosphorylated I_to
i_CaMKs = Y(24);    %Fast inactivation for CaMK phosphorylated I_to

% I_LCa gating variables 
d         = Y(25);  %Activation for I_Cal
f_f       = Y(26);  %Fast voltage-dependent inactivation for I_CaL
f_s       = Y(27);  %Slow voltage-dependent inactivation for I_CaL
f_Caf     = Y(28);  %Fast Ca-dependent inactivation for I_CaL
f_Cas     = Y(29);  %Slow Ca-dependent inactivation for I_CaL
j_Ca      = Y(30);  %Recovery from Ca-dependent inactivation for I_CaL
f_CaMKf   = Y(31);  %Fast Ca-dependent inactivation for CaMK phosphorylated I_CaL
f_CaCaMKf = Y(32);  %Slow Ca-dependent inactivation for CaMK phosphorylated I_CaL
n         = Y(33);  %Fraction in Ca-dependent inactivation mode for I_CaL

% I_Kr gating variables
x_rf = Y(34);       %Fast activation/deactivation for I_Kr
x_rs = Y(35);       %Slow activation/deactivation for I_Kr

% I_Ks gating variables 
x_s1 = Y(36);       %Activation for I_Ks
x_s2 = Y(37);       %Deactivation for I_Ks

% I_K1 gating variables 
x_K1 = Y(38);       %Time-dependent inactivation for I_K1

% SR Ca release fluxes 
J_relNP = Y(39); 
J_relCaMK = Y(40); 

% Fraction of autonomous CaMK binding sites with trapped CaM
CaMK_trap = Y(41);  

% Stimulus current 
I_stim = Y(42); 

%% CaMK (page 16) 

% Fraction of CaMK binding sites bound to Ca/CaM
CaMK_bound = CaMKo * (1 - CaMK_trap) / (1 + K_mCaM / Ca_ss); 

% Fraction of active CaMK binding sites
CaMK_active = CaMK_bound + CaMK_trap; 

%% Nernst potentials (J / kC = mV) (page 6)

E_Na = ((R*1000) * T / F) * log(Nao / Na_i); 
E_K  = ((R*1000) * T / F) * log(Ko / K_i); 
E_Ks = ((R*1000) * T / F) * log( (Ko + PR_NaK * Nao) / (K_i + PR_NaK * Na_i)); 


%% I_Na - Na current (page 6)

% Steady-state values
m_inf      = 1 / (1 + exp( -(Vm + 39.57) / 9.871)); 
h_inf      = 1 / (1 + exp( (Vm + 82.9) / 6.086)); 
j_inf      = h_inf; 
h_CaMKinf  = 1 / (1 + exp( (Vm + 89.1) / 6.086)); 
j_CaMKinf  = j_inf; 
m_linf     = 1 / (1 + exp( -(Vm + 42.85) / 5.264)); 
h_linf     = 1 / (1 + exp( (Vm + 87.61) / 7.488)); 
h_lCaMKinf = 1 / (1 + exp( (Vm + 93.81) / 7.488)); 

% Time-scales (ms)
tau_m      = 1 / (6.765 * exp( (Vm + 11.64) / 34.77) + 8.552 * exp( -(Vm + 77.42) / 5.955)); 
tau_hf     = 1 / (1.432e-5 * exp( -(Vm + 1.196) / 6.285) + 6.149 * exp( (Vm + 0.5096) / 20.27)); 
tau_hs     = 1 / (0.009764 * exp( -(Vm + 17.95) / 28.05) + 0.3343 * exp( (Vm + 5.73) / 56.66)); 
tau_j      = 2.038 + 1 / (0.02136 * exp( -(Vm + 100.6) / 8.281) + 0.3052 * exp( (Vm + 0.9941) / 38.45)); 
tau_hCaMKs = 3*tau_hs; 
tau_jCaMK  = 1.46 * tau_j; 
tau_ml     = tau_m; 
tau_hl     = 200; 
tau_hLCaMK = 3 * tau_hl; 

% Inactivation gating variable
h       = A_hf * h_f + A_hs * h_s; 
h_CaMKf = h_f; 
h_CaMK  = A_hf * h_CaMKf + A_hs * h_CaMKs; 

% Fraction of channels of type I_Na phosphorylated by CaMK
phi_I_NaCaMK  = 1 / (1 + K_mCaMK_active / CaMK_active); 
phi_I_NalCaMK = 1 / (1 + K_mCaMK_active / CaMK_active); 

% Currents (mA)
I_Naf = Gbar_Naf * (Vm - E_Na) * m^3 * ( (1 - phi_I_NaCaMK) * h * j + phi_I_NaCaMK * h_CaMK * j_CaMK); 
I_Nal = Gbar_Nal * (Vm - E_Na) * m_l^3 * ( (1 - phi_I_NalCaMK) * h_l + phi_I_NalCaMK * h_lCaMK); 
I_Na  = I_Naf + I_Nal; 

%% I_to - Transient outward K current (page 8)

% Voltage-dependent percentage of fast and slow i and i_CaMK channels 
A_if = 1 / (1 + exp( (Vm - 213.6) / 151.2)); 
A_is = 1 - A_if; 
A_iCaMKf = A_if; 
A_iCaMKs = A_is; 

% I'm not sure what these values mean 
delta_CaMKdev = 1.354 + 1e-4 / (... 
    exp( (Vm - 167.4) / 15.89) + ...
    exp( -(Vm - 12.23) / 0.2154)); 
delta_CaMKrec = 1 - 0.5 / (1 + exp( (Vm + 70) / 20)); 

% Steady-state values 
a_inf = 1 / (1 + exp( -(Vm - 14.34) / 14.82)); 
i_inf = 1 / (1 + exp( (Vm + 43.94) / 5.711)); 
a_CaMKinf = 1 / (1 + exp( -(Vm - 24.34) / 14.82)); 
i_CaMKinf = i_inf; 

% Time-scales 
tau_a = 1.0515 / ( ... 
    (1 / (1.2089 * (1 + exp( -(Vm - 18.41) / 29.38)))) + ...
    (3.5 / (1 + exp( (Vm + 100) / 29.38))));
tau_if = 4.562 + (1 / ( ... 
    (0.3933 * exp( -(Vm + 100) / 100)) + (0.08004 * exp( (Vm + 50) / 16.59)))); 
tau_is = 23.62 + (1 / ( ... 
    (0.001416 * exp( -(Vm + 96.52) / 59.05)) + 1.7808e-8 * exp( (Vm + 114.1) / 8.079))); 
if celltype == 1 % adjust if celltype is epicardial (1)
    delta_epi = 1.0 - 0.95/(1.0 + exp((Vm + 70.0)/5.0));
    tau_if = tau_if * delta_epi;
    tau_is = tau_is * delta_epi;
end
tau_aCaMK = tau_a; 
tau_iCaMKf = tau_if * delta_CaMKdev * delta_CaMKrec; 
tau_iCaMKs = tau_is * delta_CaMKdev * delta_CaMKrec; 

% Inactivation of I_to 
i = A_if * i_f + A_is * i_s; 
i_CaMK = A_iCaMKf * i_CaMKf + A_iCaMKs * i_CaMKs; 

% Fraction of channels of type I_to phosphorylated by CaMK 
phi_I_toCaMK = 1 / (1 + K_mCaMK_active / CaMK_active); 

% Current (mA)
I_to = Gbar_to * (Vm - E_K) * ( (1 - phi_I_toCaMK) * a * i + phi_I_toCaMK * a_CaMK * i_CaMK); 

%% I_CaL - L-type Ca Current (page 9)

% Voltage-dependent percentages of fast and slow f_Ca and f_CaCaMK channels
A_fCaf     = 0.3 + 0.6 / (1 + exp( (Vm - 10)/10)); 
A_fCas     = 1 - A_fCaf; 
A_fCaCaMKf = A_fCaf;
A_fCaCaMKs = A_fCas; 

% Steady-state values 
d_inf       = 1 / (1 + exp( -(Vm + 3.94) / 4.23)); 
f_inf       = 1 / (1 + exp( (Vm + 19.58) / 3.696)); 
f_Cainf     = f_inf; 
j_Cainf     = f_Cainf; 
f_CaMKinf   = f_inf; 
f_CaCaMKinf = f_inf; 

% Time-scales (ms)
tau_d        = 0.6 + 1 / (exp( -0.05 * (Vm + 6)) + exp( 0.09 * (Vm + 14))); 
tau_ff       = 7 + 1 / ( 0.0045 * exp( -(Vm + 20) / 10) + 0.0045 * exp( (Vm + 20) / 10)); 
tau_fs       = 1000 + 1 / (0.000035 * exp( - (Vm + 5) / 4) + 0.000035 * exp( (Vm + 5) / 6)); 
tau_fCaf     = 7 + 1 / (0.04 * exp( -(Vm - 4) / 7) + 0.04 * exp( (Vm - 4) / 7)); 
tau_fCas     = 100 + 1 / (0.00012 * exp( -Vm / 3) + 0.00012 * exp( Vm / 7)); 
tau_jCa      = 75; 
tau_fCaMKf   = 2.5 * tau_ff; 
tau_fCaCaMKf = 2.5 * tau_fCaf;

% Inactivation for I_to 
f         = A_ff * f_f + A_fs * f_s; 
f_Ca      = A_fCaf * f_Caf + A_fCas * f_Cas; 
f_CaMKs   = f_s; 
f_CaMK    = A_fCaMKf * f_CaMKf + A_fCaMKs * f_CaMKs; 
f_CaCaMKs = f_Cas; 
f_CaCaMK  = A_fCaCaMKf * f_CaCaMKf + A_fCaCaMKs * f_CaCaMKs; 
%{
I think there's a typo in the paper because it originally had the
following:
f_CaCaMK = A_f_CaCaMKf * f_CACaMKf + A_fCaMKs * f_CaMKs; 
The second term doesn't use CaMK phosphorylated f_Ca but rather just f. 
%}

% Rate for gating variable n 
alpha_n = 1 / ( (1000 / (j_Ca * 1)) + (1 + K_mCass / Ca_ss)^4); 

% Maximum current through ion channel (microA / microF)
psi_Ca   = z_Ca^2 * (Vm * F^2 / ((R*1000) * T)) * (gamma_Cai * Ca_ss * exp( z_Ca * Vm * F / ((R*1000) * T)) ... 
    - gamma_Cao * Cao) / (exp( z_Ca * Vm * F / ((R*1000) * T)) - 1); 
psi_CaNa = z_Na^2 * (Vm * F^2 / ((R*1000) * T)) * (gamma_Nai * Na_ss * exp( z_Na * Vm * F / ((R*1000) * T)) ... 
    - gamma_Nao * Nao) / (exp( z_Na * Vm * F / ((R*1000) * T)) - 1); 
psi_CaK  = z_K^2 * (Vm * F^2 / ((R*1000) * T)) * (gamma_Ki * K_ss * exp( z_K * Vm * F / ((R*1000) * T)) ... 
    - gamma_Ko * Ko) / (exp( z_K * Vm * F / ((R*1000) * T)) - 1); 

Ibar_CaL      = P_Ca * psi_Ca; 
Ibar_CaLNa     = P_CaNa * psi_CaNa; 
Ibar_CaLK      = P_CaK * psi_CaK; 
Ibar_CaLCaCaMK  = P_CaCaMK * psi_Ca; 
Ibar_CaLNaCaMK = P_CaNaCaMK * psi_CaNa; 
Ibar_CaLKCaMK  = P_CaKCaMK * psi_CaK; 

% Fraction of channels of type I_CaL phosphorylated by CaMK 
phi_I_CaLCaMK = 1 / (1 + K_mCaMK_active / CaMK_active); 

% Currents (mA)
I_CaLCa  = Ibar_CaL * d * (1 - phi_I_CaLCaMK) * (f * (1 - n) + f_Ca * n * j_Ca) ...
    + Ibar_CaLCaCaMK * d * phi_I_CaLCaMK * (f_CaMK * (1 - n) + f_CaCaMK * n * j_Ca); 
I_CaLNa = Ibar_CaLNa * d * (1 - phi_I_CaLCaMK) * (f * (1 - n) + f_Ca * n * j_Ca) ...
    + Ibar_CaLNaCaMK * d * phi_I_CaLCaMK * (f_CaMK * (1 - n) + f_CaCaMK * n * j_Ca); 
I_CaLK  = Ibar_CaLK * d * (1 - phi_I_CaLCaMK) * (f * (1 - n) + f_Ca * n * j_Ca) ... 
    + Ibar_CaLKCaMK * d * phi_I_CaLCaMK * (f_CaMK * (1 - n) + f_CaCaMK * n * j_Ca); 

%% I_Kr - Rapid delayed rectifier K current (page 11)

% Voltage dependent percentage of fast and slow x_r channels 
A_xrf = 1 / (1 + exp( (Vm + 54.81) / 38.21)); 
A_xrs = 1 - A_xrf; 

% Steady-state values
x_rinf = 1 / (1 + exp( -(Vm + 8.337) / 6.789)); 

% Time-scales (ms)
tau_xrf = 12.98 + 1 / ( ... 
    0.3652 * exp( (Vm - 31.66) / 3.869) + ...
    4.123e-5 * exp( -(Vm - 47.78) / 20.38)); 
tau_xrs = 1.865 + 1 / ( ... 
    0.06629 * exp( (Vm - 34.7) / 7.355) + ...
    1.128e-5 * exp( -(Vm - 29.74) / 25.94)); 

% Instantaneous inactivation of I_Kr
R_Kr = 1 / ( (1 + exp( (Vm + 55) / 75)) * (1 + exp( (Vm - 10) / 30))); 

% Activation of I_Kr
x_r = A_xrf * x_rf + A_xrs * x_rs; 

% Current (mA)
I_Kr = Gbar_Kr * sqrt(Ko / 5.4) * x_r * R_Kr * (Vm - E_K); 

%% I_Ks - Slow delayed retifier K current (page 11)

% Steady-state values 
x_s1inf = 1 / (1 + exp( -(Vm + 11.6) / 8.932)); 
x_s2inf = x_s1inf; 

% Time-scales
tau_xs1 = 817.3 + 1 / ( ... 
    2.326e-4 * exp( (Vm + 48.28) / 17.8) + ...
    0.001292 * exp( -(Vm + 210) / 230)); 
tau_xs2 = 1 / ( 0.01 * exp( (Vm - 50) / 20) + 0.0193 * exp( -(Vm + 66.54) / 31)); 

% Current (mA)
I_Ks = Gbar_Ks * (1 + 0.6 / (1 + (3.8e-5 / Ca_i)^1.4)) * x_s1 * x_s2 * (Vm - E_Ks); 

%% I_K1 - Inward rectifier K current (page 12)

% Steady-state value
x_K1inf = 1 / (1 + exp( -(Vm + 2.5538 * Ko + 144.59) / (1.5692 * Ko + 3.8115))); 

% Time-scale (ms)
tau_xK1 = 122.2 / (exp( -(Vm + 127.2) / 20.36) + exp( (Vm + 236.8) / 69.33)); 

% Inactivation of I_K1
R_K1 = 1 / (1 + exp( (Vm + 105.8 - 2.6 * Ko) / 9.493)); 

% Current (mA)
I_K1 = Gbar_K1 * sqrt(Ko) * x_K1 * R_K1 * (Vm - E_K); 

%% I_NaCa - Na/Ca Exchanger Current (page 12) 
%{ 
This model is from Kang and Hilgemann 2004, which has an extensive setup. I
have made indications of what each of these components are with respect to
their original formulation, which was not given in the Supplemental
Materials of the O'Hara Rudy 2011 paper. 

Kang, T.M. and Hilgemann, D.W.. "Multiple transport modes of the cardiac
Na+/Ca2+ exchanger." Nature, 427, 544-548. 2004. 
%} 

h_Ca = exp(q_Ca * Vm * F / ((R*1000) * T)); %Voltage-dependent transport of Cao to binding site
h_Na = exp(q_Na * Vm * F / ((R*1000) * T)); %Voltage-dependent transport of Nao to 1Na binding site 

% I_NaCa_i----------------------------------------------------------------
h1 = 1 + Na_i / k_Na3 *(1 + h_Na); 
h2 = Na_i * h_Na / (k_Na3 * h1); 
h3 = 1 / h1; 
h4 = 1 + Na_i / k_Na1 * (1 + Na_i / k_Na2); 
h5 = Na_i^2 / (h4 * k_Na1 * k_Na2); 
h6 = 1 / h4; 
h7 = 1 + Nao / k_Na3 *(1 + 1 / h_Na); 
h8 = Nao / (k_Na3 * h_Na * h7); 
h9 = 1 / h7; 
h10 = k_asymm + 1 + Nao / k_Na1 * (1 + Nao / k_Na2); 
h11 = Nao^2 / (h10 * k_Na1 * k_Na2); 
h12 = 1 / h10; 

% Rate constants between transporter states (E_2)
k1       = h12 * Cao * k_Caon; 
k2       = k_Caoff; 
k3prime  = h9 * omega_Ca; 
k3dprime = h8 * omega_NaCa;
k3       = k3prime + k3dprime; 
k4prime  = h3 * omega_Ca / h_Ca;
k4dprime = h2 * omega_NaCa; 
k4       = k4prime + k4dprime; 
k5       = k_Caoff; 
k6       = h6 * Ca_i * k_Caon; 
k7       = h5 * h2 * omega_Na;
k8       = h8 * h11 * omega_Na;

% Used to calculate transporter states
x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3); 
x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8); 
x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3); 
x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8); 

% Four NCX transporter states 
E1 = x1 / (x1 + x2 + x3 + x4); 
E2 = x2 / (x1 + x2 + x3 + x4); 
E3 = x3 / (x1 + x2 + x3 + x4); 
E4 = x4 / (x1 + x2 + x3 + x4); 

% Ca-dependent allosteric current activation factor
allo_i = 1 / (1 + (K_mCaact / Ca_i)^2); 

% Flux of Na_i and Ca_i across membrane 
J_NaCaNa_i = 3 * (E4 * k7 - E1 * k8) + E3 * k4dprime - E2 * k3dprime; 
J_NaCaCa_i = E2 * k2 - E1 * k1; 

% Current 
I_NaCa_i = Gbar_NaCa * 0.8 * allo_i * (z_Na * J_NaCaNa_i + z_Ca * J_NaCaCa_i); 

% I_NaCa_ss---------------------------------------------------------------
h1 = 1 + Na_ss / k_Na3 *(1 + h_Na); 
h2 = Na_ss * h_Na / (k_Na3 * h1); 
h3 = 1 / h1; 
h4 = 1 + Na_ss / k_Na1 * (1 + Na_ss / k_Na2); 
h5 = Na_ss^2 / (h4 * k_Na1 * k_Na2); 
h6 = 1 / h4; 
h7 = 1 + Nao / k_Na3 *(1 + 1 / h_Na); 
h8 = Nao / (k_Na3 * h_Na * h7); 
h9 = 1 / h7; 
h10 = k_asymm + 1 + Nao / k_Na1 * (1 + Nao / k_Na2); 
h11 = Nao^2 / (h10 * k_Na1 * k_Na2); 
h12 = 1 / h10; 

% Rate constants between transporter states (E_'2)
k1       = h12 * Cao * k_Caon; 
k2       = k_Caoff; 
k3prime  = h9 * omega_Ca; 
k3dprime = h8 * omega_NaCa;
k3       = k3prime + k3dprime; 
k4prime  = h3 * omega_Ca / h_Ca;
k4dprime = h2 * omega_NaCa; 
k4       = k4prime + k4dprime; 
k5       = k_Caoff; 
k6       = h6 * Ca_ss * k_Caon; 
k7       = h5 * h2 * omega_Na;
k8       = h8 * h11 * omega_Na;

% Used to calculate transporter states
x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3); 
x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8); 
x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3); 
x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8); 

% Four NCX trasporter states 
E1 = x1 / (x1 + x2 + x3 + x4); %3 bound Na outside
E2 = x2 / (x1 + x2 + x3 + x4); %1 bound Ca (with possibly 1 bound Na) outside
E3 = x3 / (x1 + x2 + x3 + x4); %1 bound Ca (with possibly 1 bound Na) inside
E4 = x4 / (x1 + x2 + x3 + x4); %3 bound Na inside 

% Ca-dependent allosteric current activation factor
allo_ss = 1 / (1 + (K_mCaact / Ca_ss)^2); 

% Flux of Na_ss and Ca_ss across membrane 
J_NaCaNa_ss = 3 * (E4 * k7 - E1 * k8) + E3 * k4dprime - E2 * k3dprime; 
J_NaCaCa_ss = E2 * k2 - E1 * k1; 

% Currents 
I_NaCa_ss = Gbar_NaCa * 0.2 * allo_ss * (z_Na * J_NaCaNa_ss + z_Ca * J_NaCaCa_ss); 
I_NaCa = I_NaCa_i + I_NaCa_ss; 

%% I_NaK - Na/K ATPase current (page 14)
%{ 
This model is from Smith and Crampin 2004, which has an extensive setup. I
have made indications of what each of these components are with respect to
their original formulation, which was not given in the Supplemental
Materials of the O'Hara Rudy 2011 paper. 

Smith, N.P. and Crampin, E.J.. "Development of models of active ion
transport for whole-cell modelling: cardiac sodium-potassium pump as a case
study." Biophys Mol Biol, 85, 387-405. 2004.
%}

% Voltage-dependent half-saturation vlaues 
K_Nai   = KoNai * exp(Delta * Vm * F / (3 * (R*1000) * T));
K_Nao   = KoNao * exp((1 - Delta) * Vm * F / (3 * (R*1000) * T)); 

% Inorganic phosphate concentration
P = SigmaP / (1 + H / K_HP + Na_i / K_NaP + K_i / K_KP); 

% Lumped parameter pseudo-first-order rate constants for metabolite binding
% and unbinding reactions in the forward (clockwise) direction 
alpha1 = k_1plus * (Na_i / K_Nai)^3 / ( ... 
    (1 + Na_i / K_Nai)^3 + (1 + K_i / K_Ki)^2 - 1); 
alpha2 = k_2plus; 
alpha3 = k_3plus * (Ko / K_Ko)^2 / ( ... 
    (1 + Nao / K_Nao)^3 + (1 + Ko / K_Ko)^2 - 1); 
alpha4 = k_4plus * (MgATP / K_MgATP) / (1 + MgATP / K_MgATP); 

% Lumped parameter pseudo-first-order rate constants for metabolite binding
% and unbinding reactions in the reverse (counterclockwise) direction
beta1 = k_1neg * MgADP; 
beta2 = k_2neg * (Nao / K_Nao)^3 / ( ...
    (1 + Nao / K_Nao)^3 + (1 + Ko / K_Ko)^2 - 1); 
beta3 = k_3neg * P * H / (1 + MgATP / K_MgATP); 
beta4 = k_4neg * (K_i / K_Ki)^2 / ( ... 
    (1 + Na_i / K_Nai)^3 + (1 + K_i / K_Ki)^2 - 1); 

% Used to calculate ATPase states (similar to NCX)
x1 = alpha4 * alpha1 * alpha2 + beta2 * beta4 * beta3 + ...
    alpha2 * beta4 * beta3 + beta3 * alpha1 * alpha2; 
x2 = beta2 * beta1 * beta4 + alpha1 * alpha2 * alpha3 + ...
    alpha3 * beta1 * beta4 + alpha2 * alpha3 * beta4; 
x3 = alpha2 * alpha3 * alpha4 + beta3 * beta2 * beta1 + ...
    beta2 * beta1 * alpha4 + alpha3 * alpha4 * beta1; 
x4 = beta4 + beta3 * beta2 + alpha3 * alpha4 * alpha1 + ...
    beta2  * alpha4 * alpha1 + beta3 * beta2 * alpha1; 

% Four N/K ATPase transporter states 
E1 = x1 / (x1 + x2 + x3 + x4);
E2 = x2 / (x1 + x2 + x3 + x4);
E3 = x3 / (x1 + x2 + x3 + x4);
E4 = x4 / (x1 + x2 + x3 + x4);

% Fluxes of Na and K across membrane 
J_NaKNa = 3 * (E1 * alpha3 - E2 * beta3); 
J_NaKK = 2 * (E4 * beta1 - E3 * alpha1); 

% Current
I_NaK = G_NaK * 30 * (z_Na * J_NaKNa + z_K * J_NaKK); 

%% I_Nab - Background Na current (page 15)

I_Nab = P_Nab * z_Na^2 * (Vm * F^2 / ((R*1000) * T)) * (... 
    (Na_i * exp(z_Na * Vm * F / ((R*1000) * T)) - Nao) / ...
    (exp(z_Na * Vm * F / ((R*1000) * T)) - 1)); 
%{
This forumation doesn't have any gammas in it... is this a typo? 
%}

%% I_Cab - Background Ca current (page 15)

I_Cab = P_Cab * z_Ca^2 * (Vm * F^2 / ((R*1000) * T)) * (... 
    (gamma_Cai * Ca_i * exp(z_Ca * Vm * F / ((R*1000) * T)) - gamma_Cao * Cao) / ...
    (exp(z_Ca * Vm * F / ((R*1000) * T)) - 1));

%% I_Kb - Background K current (page 15) 

x_Kb = 1 / (1 + exp( -(Vm - 14.48) / 18.34)); 
I_Kb = Gbar_Kb * x_Kb * (Vm - E_K); 

%% I_pCa - Sarcolemmal Ca pump (page 16)

I_pCa = Gbar_pCa * Ca_i / (0.0005 + Ca_i); 

%% Itot - Total currents (page 18, separated out from the differential equations)

I_Na_i  = I_Na + 3 * I_NaCa_i + 3 * I_NaK + I_Nab; 
%{
I think there's a typo in the paper because it originally had the
following for the total Na current
I_Na_i = I_Na + I_NaL + 3*I_NaCa_i + 3*I_NaK + I_Nab; 
But there was no I_NaL above. After looking back at the Hund paper, I_Na is
the fast Na channel and I_NaL is the late Na channel. However, this is not
how the ORd model has been formulated. In this model, I_Na = I_Naf + I_Nal,
the sum of the two Na channels. I'm going to leave this formulation. 

Hund, TJ; Decker, KF; Kanter, E; Mohler, PJ; Boyden, PA; Schuessler, RB;
Yamada, KA; and Rudy, Y. "Role of activated CaMKII in abnormal calcium
homeostasis and I_Na remodeling after myocardial infarction: insights from
mathematical modeling. J Mol Cell Cardiol. 2008; 45(3):420-428. 
%}
I_Na_ss = I_CaLNa + 3 * I_NaCa_ss; 
I_K_i   = I_to + I_Kr + I_Ks + I_K1 + I_Kb - 2 * I_NaK; 
%{
I think there's a typo in the paper because it originally had the
following for the total K current
I_K_i = I_to + I_Kr + I_Ks + I_K1 + I_Kur - 2*I_NaK;  
But there was no I_Kur above. After looking back at the Sridhar paper, 
I_Kur is the ultra fast rectifying K channel. However, this is not defined
in the documentation for the ORd model, nor is a formulation defined in the
variations of the Hund-Rudy model or the original paper. Thus, I_Kur was 
omitted. 

Sridhar A, da Cunha DN, Lacombe VA, Zhou Q, Fox JJ, Hamlin RL, Carnes CA. 
The plateau outward current in canine ventricle, sensitive to 
4-aminopyridine, is a constitutive contributor to ventricular 
repolarization. Br J Pharmacol. 2007;152(6):870-879
%}
I_K_ss  = I_CaLK; 
I_Ca_i  = I_pCa + I_Cab - 2 * I_NaCa_i; 
I_Ca_ss = I_CaLCa - 2 * I_NaCa_ss; 

Itot = I_Naf + I_Nal + I_to + I_CaLCa + I_CaLNa + I_CaLK + I_Kr + I_Ks + I_K1 + ... 
    I_NaCa_i + I_NaCa_ss + I_NaK + I_pCa + I_Cab + I_Nab  + I_Kb ; 

%% Metabolite Diffusion Fluxes 

% Fluxes (mM ms^(-1))
J_diffNa = (Na_ss - Na_i) / tau_diffNa; 
J_diffK  = (K_ss - K_i) / tau_diffK; 
J_diffCa = (Ca_ss - Ca_i) / tau_diffCa; 

%% J_rel - SR Ca release Flux via Ryanodine receptor 
%{ 
This model is from Livshitz and Rudy 2007. I
have made indications of what each of these components are with respect to
their original formulation, which was not given in the Supplemental
Materials of the O'Hara Rudy 2011 paper. 

Livshitz, L.M. and Rudy, Y. "Regulation of Ca2+ and electrcal alternans i9n
cardiac myocytes: role of CaMKII and repolarizing currents." Am J Physiol
Heart Circ Physiol, 292, H2854-H2866. 2007. 
%}

% Steady-state values
J_relNPinf   = alpha_rel * (-I_CaLCa) / (1 + (1.5 / Ca_jsr)^8); 
J_relCaMKinf = alpha_relCaMK * (-I_CaLCa) / (1 + (1.5 / Ca_jsr)^8); 

% Time-scales 
tau_relNP   = beta_tau / (1 + 0.0123 / Ca_jsr); 
tau_relCaMK = beta_tauCaMK / (1 + 0.0123 / Ca_jsr); 

phi_relCaMK = 1 / (1 + K_mCaMK_active / CaMK_active); 

% Flux 
J_rel = G_RyR * (1 - phi_relCaMK)*J_relNP + phi_relCaMK * J_relCaMK; 

%% J_up - Ca uptake via SERCA pump

% Fluxes (dimensionless)
J_upNP = G_up_epi * (0.004375 * Ca_i / (0.00092 + Ca_i)); 
J_upCaMK = G_up_epi * ((1 + DeltaJbar_upCaMK) * 0.004375 * Ca_i / (0.00092 - DeltaKbar_mPLB + Ca_i)); 

phi_upCaMK = 1 / (1 + K_mCaMK_active / CaMK_active); 

J_leak = 0.0039375 * Ca_nsr / 15; 
J_up = G_SERCA * ((1 - phi_upCaMK) * J_upNP + phi_upCaMK * J_upCaMK - J_leak); 

%% J_tr - Ca translocation from NSR to JSR 

% Flux (mM ms^(-1))
J_tr = (Ca_nsr - Ca_jsr) / tau_tr; 

%% Buffers 

beta_Cai = 1 / (1 + ...
    cmdnbar * K_mcmdn / (K_mcmdn + Ca_i)^2 + ...
    trpnbar * K_mtrpn / (K_mtrpn + Ca_i)^2); 
beta_Cass = 1 / (1 + ... 
    bsrbar * K_mbsr / (K_mbsr + Ca_ss)^2 + ...
    bslbar * K_mbsl / (K_mbsl + Ca_ss)^2); 
beta_Cajsr = 1 / (1 + ... 
    csqnbar * K_mcsqn / (K_mcsqn + Ca_jsr)^2); 

%% Differential equations 

% Membrane voltage (mV/ms)
dVm = (I_stim - Itot) / Cm; 

% Concentrations (mM)
dNa_i   = -I_Na_i *A_cap / (F * v_myo) + J_diffNa * v_ss / v_myo;
dNa_ss  = -I_Na_ss * A_cap / (F * v_ss) - J_diffNa;
dK_i    = (I_stim - I_K_i) * A_cap / (F * v_myo) + J_diffK * v_ss / v_myo;
dK_ss   = -I_K_ss * A_cap / (F * v_ss) - J_diffK; 
dCa_i   = beta_Cai * (-I_Ca_i * A_cap / (2 * F * v_myo) - ... 
    J_up * v_nsr / v_myo + J_diffCa * v_ss / v_myo); 
dCa_ss  = beta_Cass * (-I_Ca_ss * A_cap / (2 * F * v_ss) + ... 
    J_rel * v_jsr / v_ss - J_diffCa); 
dCa_nsr = J_up - J_tr * v_jsr / v_nsr; 
dCa_jsr = beta_Cajsr * (J_tr - J_rel); 

% I_Na gating variables 
dm       = (m_inf - m) / tau_m; 
dh_f     = (h_inf - h_f) / tau_hf; 
dh_s     = (h_inf - h_s) / tau_hs; 
dj       = (j_inf - j) / tau_j; 
dh_CaMKs = (h_CaMKinf - h_CaMKs) / tau_hCaMKs; 
dj_CaMK  = (j_CaMKinf - j_CaMK) / tau_jCaMK; 
dm_l     = (m_linf - m_l) / tau_ml; 
dh_l     = (h_linf - h_l) / tau_hl; 
dh_lCaMK = (h_lCaMKinf - h_lCaMK) / tau_hLCaMK; 

% I_to gating variables 
da       = (a_inf - a) / tau_a; 
di_f     = (i_inf - i_f) / tau_if; 
di_s     = (i_inf - i_s) / tau_is; 
da_CaMK  = (a_CaMKinf - a_CaMK) / tau_aCaMK; 
di_CaMKf = (i_CaMKinf - i_CaMKf) / tau_iCaMKf; 
di_CaMKs = (i_CaMKinf - i_CaMKs) / tau_iCaMKs; 

% I_CaL gating variables 
dd         = (d_inf - d) / tau_d; 
df_f       = (f_inf - f_f) / tau_ff; 
df_s       = (f_inf - f_s) / tau_fs; 
df_Caf     = (f_Cainf - f_Caf) / tau_fCaf; 
df_Cas     = (f_Cainf - f_Cas) / tau_fCas; 
dj_Ca      = (j_Cainf - j_Ca) / tau_jCa; 
df_CaMKf   = (f_CaMKinf - f_CaMKf) / tau_fCaMKf; 
df_CaCaMKf = (f_CaCaMKinf - f_CaCaMKf) / tau_fCaCaMKf; 
dn         = alpha_n * 1000 - n * j_Ca * 1; 

% I_Kr gating variables 
dx_rf = (x_rinf - x_rf) / tau_xrf; 
dx_rs = (x_rinf - x_rs) / tau_xrs; 

% I_Ks gating variables 
dx_s1 = (x_s1inf - x_s1) / tau_xs1; 
dx_s2 = (x_s2inf - x_s2) / tau_xs2; 

% I_K1 gating variables 
dx_K1 = (x_K1inf - x_K1) / tau_xK1; 

% SR Ca release fluxes 
dJ_relNP   = (J_relNPinf - J_relNP) / tau_relNP; 
dJ_relCaMK = (J_relCaMKinf - J_relCaMK) / tau_relCaMK; 

% Fraction of autonomous CaMK binding sites 
dCaMK_trap = alpha_CaMK * CaMK_bound * CaMK_active - ...
    beta_CaMK * CaMK_trap; 

dI_stim = 0; 

dYdt = [dVm; 
    dNa_i; dNa_ss; dK_i; dK_ss; dCa_i; dCa_ss; dCa_nsr; dCa_jsr; 
    dm; dh_f; dh_s; dj; dh_CaMKs; dj_CaMK; dm_l; dh_l; dh_lCaMK; 
    da; di_f; di_s; da_CaMK; di_CaMKf; di_CaMKs; 
    dd; df_f; df_s; df_Caf; df_Cas; dj_Ca; df_CaMKf; df_CaCaMKf; dn; 
    dx_rf; dx_rs; 
    dx_s1; dx_s2; 
    dx_K1; 
    dJ_relNP; dJ_relCaMK; 
    dCaMK_trap; 
    dI_stim]; 
dYdt = dYdt * 1000; % Corrects units from s^(-1) to ms^(-1) 

data = [I_Na, I_to, I_CaLCa, I_CaLNa, I_CaLK, ...
    I_Kr, I_Ks, I_K1, I_NaCa, I_NaK, ...
    I_Nab, I_Cab, I_Kb, I_pCa, ...
    I_stim]; 
end