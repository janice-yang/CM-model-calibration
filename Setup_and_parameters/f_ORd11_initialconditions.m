function init = f_ORd11_initialconditions

%% Initial conditions 

% Voltage (mV) 
V0 = -87.84; 

% Concentrations (mM)
Na_i0  = 7.23; 
Na_ss0 = 7.23;
K_i0   = 143.79;
K_ss0  = 143.79;
Ca_i0  = 8.54e-5; 
Ca_ss0 = 8.43e-5; 
Ca_nsr = 1.61; 
Ca_jsr = 1.56; 

% I_Na gating variables
m0 = 0;
h_f0 = 0.7;
h_s0 = 0.7;
j0 = 0.7;
h_CaMKs0 = 0.4; 
j_CaMK0 = 0.7; 
m_L0 = 0; 
h_L0 = 0.5; 
h_LCaMK = 0.3; 

% I_to gating variables
a0 = 0; 
i_f0 = 1;
i_s0 = 0.6; 
a_CaMK0 = 0; 
i_CaMKf0 = 1; 
i_CaMKs0 = 0.6;

% I_LCa gating variables
d0 = 0; 
f_f0 = 1;
f_s0 = 0.9; 
f_Caf0 = 1; 
f_Cas0 = 1; 
j_Ca0  = 1;
f_CaMKf0 = 1; 
f_CaCaMKf0 = 1;
n0 = 0; 

% I_Kr gating variables
x_rf0 = 0; 
x_rs0 = 0.5; 

% I_Ks gating variables
x_s10 = 0.3; 
x_s20 = 0; 

% I_K1 gating variables 
x_K10 = 1; 

% Ca Fluxes (mM ms^(-1))
J_relNP0 = 2.53943e-5; 
J_relCaMK0 = 3.17262e-7; 

CaMK_trap = 0.0124065; 

stim = 0;

init = [V0; 
    Na_i0; Na_ss0; K_i0; K_ss0; Ca_i0; Ca_ss0; Ca_nsr; Ca_jsr;  %Voltage and concentrations
    m0; h_f0; h_s0; j0; h_CaMKs0; j_CaMK0; m_L0; h_L0; h_LCaMK;         %I_Na gating variables
    a0; i_f0; i_s0; a_CaMK0; i_CaMKf0; i_CaMKs0;                        %I_to gating variables
    d0; f_f0; f_s0; f_Caf0; f_Cas0; j_Ca0; f_CaMKf0; f_CaCaMKf0; n0;    %I_LCa gating variables
    x_rf0; x_rs0;                                                       %I_Kr gating variables
    x_s10; x_s20;                                                       %I_Ks gating variables  
    x_K10;                                                              %I_K1 gating variables 
    J_relNP0; J_relCaMK0;                                               %Ca fluxes
    CaMK_trap; 
    stim
    ]; 

end 

