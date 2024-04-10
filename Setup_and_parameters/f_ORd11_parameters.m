function pars = f_ORd11_parameters(celltype)


%%  Maximal conductances (mS microF^(-1))

Gbar_Naf  = 75; 
Gbar_Nal  = 0.0075; 
Gbar_to   = 0.02; 
Gbar_Kr   = 0.046; 
Gbar_Ks   = 0.0034; 
Gbar_K1   = 0.1908; 
Gbar_NaCa = 0.0008; % (microA / microF)
Gbar_Kb   = 0.003;
Gbar_pCa  = 0.0005; 
G_RyR     = 1;
G_SERCA   = 1;
G_NaK     = 1;
G_up_epi  = 1;

%% Permeabilities (cm s^(-1))

PR_NaK     = 0.01833; %(dimensionless)  %Permeability ratio for Nernst potential

P_Ca       = 0.0001; 
P_CaNa     = 0.00125 * P_Ca; 
P_CaK      = 3.574e-4 * P_Ca; 
P_CaCaMK   = 1.1 * P_Ca;
P_CaNaCaMK = 0.00125 * P_CaCaMK; 
P_CaKCaMK  = 3.574e-4 * P_CaCaMK; 
P_Nab      = 3.75e-10; 
P_Cab      = 2.5e-8; 

%% Activity coefficients (dimensionless)

gamma_Nai = 0.75; 
gamma_Nao = 0.75; 
gamma_Ki  = 0.75; 
gamma_Ko  = 0.75; 
gamma_Cai = 1; 
gamma_Cao = 0.341; 

%% Average concentrations under physiological conditions (mM)

% Compounds
MgADP  = 0.05; 
MgATP  = 9.8; 
H      = 1e-7; 
SigmaP = 4.2; %Total measurable inorganic phosphate concentration (e.g. by NMR)

% Buffers
cmdnbar = 0.05; 
trpnbar = 0.07;
bsrbar  = 0.047;
bslbar  = 1.124; 
csqnbar = 10; 

%% Half-saturation values (mM)

K_mCaM     = 0.0015; %Half-saturation value for CaMK activation
K_mCaMK_active  = 0.15; %dimensionless 
K_mCass  = 0.002; %called K_mn in original paper 
K_mCaact = 150e-6; 

K_mcmdn = 0.00238; 
K_mtrpn = 0.0005; 
K_mbsr  = 0.00087; 
K_mbsl  = 0.0087; 
K_mcsqn = 0.8; 

%% Time-scales (ms)

tau_diffNa = 2;
tau_diffK  = 2; 
tau_diffCa = 0.2; 
tau_tr     = 100; 

%% Channel percentages (dimensionless)

% I_Na
A_hf = 0.99; %percentage of fast h channels
A_hs = 0.01; %percentage of slow h channels

% I_CaL
% Percentage of fast and slow channels 
A_ff       = 0.6; 
A_fs       = 1 - A_ff; 
A_fCaMKf   = A_ff; 
A_fCaMKs   = A_fs; 

%% CaMK Parameters 

alpha_CaMK = 0.05; %ms^(-1)     %Phosphorylation rate of CaMK 
beta_CaMK  = 0.00068; %ms^(-1)  %Dephosphorylation rate of CaMK
CaMK_o     = 0.05;              %Fraction of CaMK binding sites at equilibrium 

%% I_NaCa parameters

k_Na1      = 15; %mM    
k_Na2      = 5; %mM     
k_Na3      = 88.12; %mM  
k_asymm    = 12.5; %Asymmetry factor for Na binding inside vs outside 
omega_Na   = 6e4; %Hz
omega_Ca   = 6e4; %Hz
omega_NaCa = 5e3; %Hz 
k_Caon     = 1.5e6; %mM/ms
k_Caoff    = 5e3; %Hz
q_Na       = 0.5224; 
q_Ca       = 0.167; 

%% I_NaK parameters

k_1plus = 949.5; %Hz
k_1neg  = 182.4; %mM^(-1)
k_2plus = 687.2; %Hz
k_2neg  = 39.4; %Hz
k_3plus = 1899; %Hz
k_3neg  = 79300; %Hz mM^(-2)
k_4plus = 639; %Hz
k_4neg  = 40; %Hz

% Dissociation constants for ion binding 
Delta   = -0.155; %Percentage of voltage-dependence between Na_i and Na_o dissociation reactions
%{
This doesn't make any sense to me to have a negative percentage. In the
original Smith and Crampin paper, the equation for K_Nao was 
K_Nao = K_0Nao * exp((1 + Delta) * V * F / (3 * R * T));
which, if Delta < 0, would be appropriate due to the 1 + Delta sign. In the
ORd model though, they switch up the sign. So I don't understand what it's
supposed to signify now. 
%} 
K_oNai  = 9.073; %mM %Na_i dissociation when V = 0 
K_oNao  = 27.78; %mM %Na_o dissociation when V = 0 
K_Ki    = 0.5; %mM
K_Ko    = 0.3582; %mM 
K_MgATP = 1.698e-7; %mM
K_HP    = 1.698e-7; %mM %Dissociation for inorganic P bound to H
K_NaP   = 224; %mM %Dissociation for inorganic P bound to Na
K_KP    = 292; %mM %Dissociation for inorganic P bound to K

%% J_rel parameters

% Maximal CaMK-dependent value of tau_rel 
beta_tau     = 4.75; %ms
beta_tauCaMK = 1.25 * beta_tau; 

% Amplitude coefficient 
alpha_rel     = 0.5 * beta_tau; 
alpha_relCaMK = 0.5 * beta_tauCaMK; 

%% J_up parameters

DeltaKbar_mPLB   = 0.00017; %mM %Maximal half-saturation value change due to phospholamban 
DeltaJbar_upCaMK = 1.75; 

if celltype == 1
    Gbar_Nal    = Gbar_Nal * 0.6;
    Gbar_to     = Gbar_to * 4.0;
    P_Ca        = P_Ca * 1.2;
    P_CaNa      = P_CaNa * 1.2;
    P_CaK       = P_CaK * 1.2;
    Gbar_Kr     = Gbar_Kr * 1.3; 
    Gbar_Ks     = Gbar_Ks * 1.4; 
    Gbar_K1     = Gbar_K1 * 1.2; 
    Gbar_NaCa   = Gbar_NaCa * 1.1;
    G_NaK       = G_NaK * 0.9;
    Gbar_Kb     = Gbar_Kb * 0.6;
    G_up_epi    = G_up_epi * 1.3;
    cmdnbar     = cmdnbar * 1.3;
end

%% Output

pars = [Gbar_Naf; Gbar_Nal; Gbar_to; Gbar_Kr; Gbar_Ks; Gbar_K1; Gbar_NaCa; Gbar_Kb; Gbar_pCa; G_RyR; G_SERCA; G_NaK; G_up_epi;  % "conductances" 
    PR_NaK; P_Ca; P_CaNa; P_CaK; P_CaCaMK; P_CaNaCaMK; P_CaKCaMK; P_Nab; P_Cab;                          %Permeabilities
    gamma_Nai; gamma_Nao; gamma_Ki; gamma_Ko; gamma_Cai; gamma_Cao;                                     %Gammas
    MgADP; MgATP; H; SigmaP; cmdnbar; trpnbar; bsrbar; bslbar; csqnbar;                                 %Average conconcentrations
    K_mCaM; K_mCaMK_active; K_mCass; K_mCaact; K_mcmdn; K_mtrpn; K_mbsr; K_mbsl; K_mcsqn;               %Half-saturation values
    tau_diffNa; tau_diffK; tau_diffCa; tau_tr;                                                          %Time-scales
    A_hf; A_hs; A_ff; A_fs; A_fCaMKf; A_fCaMKs;                                    %Percentages
    alpha_CaMK; beta_CaMK; CaMK_o;                                                                      %CaMK parameters
    k_Na1; k_Na2; k_Na3; k_asymm; omega_Na; omega_Ca; omega_NaCa; k_Caon; k_Caoff; q_Na; q_Ca;          %I_NaCa parameters
    k_1plus; k_1neg; k_2plus; k_2neg; k_3plus; k_3neg; k_4plus; k_4neg;                                 %I_NaK parameters
    Delta;                                                              
    K_oNai; K_oNao; K_Ki; K_Ko; K_MgATP; K_HP; K_NaP; K_KP; 
    beta_tau; beta_tauCaMK; alpha_rel; alpha_relCaMK;                                                   %J_rel parameters
    DeltaKbar_mPLB; DeltaJbar_upCaMK;                                                                   %J_up parameters
    ];

end 