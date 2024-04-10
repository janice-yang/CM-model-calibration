%===============================================================================
% Version Frontiers 2018 (New RyR and optimization params)
%==========================================================================

function [dY, data] = f_Paci18(time, Y, geometry, conductances, parameters)
global F R T Nao Ko Cao

%--------------------------------------------------------------------------
% Geometry
%--------------------------------------------------------------------------
V_SR = geometry.V_SR ;
Vc = geometry.Vc ;
Cm  = geometry.Cm ;
%--------------------------------------------------------------------------
% Conductances
%--------------------------------------------------------------------------
c = conductances.baseline .* conductances.scaling ;
%{ 
the following comparison to time happens inside the f_Paci18 function because
the times when drugs are added & removed may lie within intervals set up 
during the pacing protocol. 
%}
for i = 1: length(c)
    if conductances.applicationTimes(1,i) <= time && conductances.applicationTimes(2,i) > time
        c(i) = c(i) .* conductances.drugEffects(i) ;
    end
end
g_Na = c(1) ;
g_f = c(2) ; 
g_CaL = c(3) ; 
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
GNaLmax = c(16) ; 
g_serca = c(17) ; 

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
p = parameters.baseline .* parameters.scaling ;
RyRa1		= p(1) ;	% 然
RyRa2		= p(2) ; 	% 然
RyRahalf	= p(3) ;   	% 然
RyRohalf	= p(4) ;   	% 然
RyRchalf	= p(5) ;  	% 然
Kup			= p(6) ; 	% mM    (in calcium_dynamics)
alpha		= p(7) ;	% -     (in i_NaCa)

%--------------------------------------------------------------------------
% Initial conditions
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
Ki = 150.0;   % mM (in model_parameters)
% Cai  = 0.0002 mM Y(3)
% caSR = 0.3    mM Y(2)

% time (second)

%% Nernst potential
E_Na = R*T./F.*log(Nao./Y(18)) ;
E_Ca = 0.5*R*T./F.*log(Cao./Y(3)) ;
E_K  = R*T/F*log(Ko/Ki);
PkNa = 0.03;                   % - (in electric_potentials)
E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa.*Y(18)));

%% INa
i_Na        = g_Na*Y(14).^3.0.*Y(12).*Y(13).*(Y(1)-E_Na);
% normal before drug is applied, impaired after drug is applied
h_inf       = 1.0/sqrt(1.0+exp((Y(1)*1000.0+72.1)/5.7));
alpha_h     = 0.057*exp(-(Y(1)*1000.0+80.0)/6.8);
beta_h      = 2.7*exp(0.079*Y(1)*1000.0)+3.1*10.0^5.0*exp(0.3485*Y(1)*1000.0);
if (Y(1) < -0.0385)
    tau_h   = 1.5/((alpha_h+beta_h)*1000.0);
else
    tau_h   = 1.5*1.6947/1000.0;
end
dY(12, 1)   = (h_inf-Y(12))./tau_h;

j_inf       = 1.0/sqrt(1.0+exp((Y(1)*1000.0+72.1)/5.7));
if (Y(1) < -0.04)
    alpha_j = (-25428.0*exp(0.2444*Y(1)*1000.0)-6.948*10.0^-6.0*exp(-0.04391*Y(1)*1000.0)).*(Y(1)*1000.0+37.78)./(1.0+exp(0.311*(Y(1)*1000.0+79.23))) ;
else
    alpha_j = 0.0;
end
if (Y(1) < -0.04)
    beta_j  = (0.02424*exp(-0.01052*Y(1)*1000)./(1+exp(-0.1378*(Y(1)*1000+40.14))));
else
    beta_j  = (0.6*exp((0.057)*Y(1)*1000)./(1+exp(-0.1*(Y(1)*1000+32))));
end
tau_j       = 7.0/((alpha_j+beta_j)*1000.0);
dY(13, 1)   = (j_inf-Y(13))./tau_j;

m_inf       = 1.0/(1.0+exp((-Y(1)*1000.0-34.1)/5.9))^(1.0/3.0);
alpha_m     = 1.0/(1.0+exp((-Y(1)*1000.0-60.0)/5.0));
beta_m      = 0.1/(1.0+exp((Y(1)*1000.0+35.0)/5.0))+0.1/(1.0+exp((Y(1)*1000.0-50.0)/200.0));
tau_m       = 1.0*alpha_m.*beta_m/1000.0;
dY(14, 1)   = (m_inf-Y(14))./tau_m;

%% INaL
myCoefTauM  = 1;
tauINaL     = 200;      % ms
Vh_hLate    = 87.61;
i_NaL       = GNaLmax* Y(19).^(3).*Y(20).*(Y(1)-E_Na);

m_inf_L     = 1/(1+exp(-(Y(1)*1000+42.85)/(5.264)));
alpha_m_L   = 1/(1+exp((-60-Y(1)*1000)/5));
beta_m_L    = 0.1/(1+exp((Y(1)*1000+35)/5))+0.1/(1+exp((Y(1)*1000-50)/200));
tau_m_L     = 1/1000 * myCoefTauM*alpha_m_L.*beta_m_L;
dY(19, 1)   = (m_inf_L-Y(19))./tau_m_L;

h_inf_L     = 1/(1+exp((Y(1)*1000+Vh_hLate)/(7.488)));
tau_h_L     = 1/1000 * tauINaL;
dY(20, 1)   = (h_inf_L-Y(20))/tau_h_L;

%% If
E_f         = -0.017;	% V     (in i_f)
i_f         = g_f*Y(15).*(Y(1)-E_f);
i_fNa       = 0.42*g_f*Y(15).*(Y(1)-E_Na);

Xf_infinity = 1.0/(1.0+exp((Y(1)*1000.0+77.85)/5.0));
tau_Xf      = 1900.0/(1.0+exp((Y(1)*1000.0+15.0)/10.0))/1000.0;
dY(15, 1)   = (Xf_infinity-Y(15))./tau_Xf;

%% ICaL
i_CaL       = g_CaL*4.0*Y(1)*F.^2.0./(R*T).*(Y(3).*exp(2.0*Y(1)*F/(R*T))-0.341*Cao)./(exp(2.0*Y(1)*F/(R*T))-1.0).*Y(5).*Y(6).*Y(7).*Y(8);

d_infinity  = 1.0/(1.0+exp(-(Y(1)*1000.0+9.1)/7.0));
alpha_d     = 0.25+1.4/(1.0+exp((-Y(1)*1000.0-35.0)/13.0));
beta_d      = 1.4/(1.0+exp((Y(1)*1000.0+5.0)/5.0));
gamma_d     = 1.0/(1.0+exp((-Y(1)*1000.0+50.0)/20.0));
tau_d       = (alpha_d.*beta_d+gamma_d)*1.0/1000.0;
dY(5, 1)    = (d_infinity-Y(5))/tau_d;

f1_inf      = 1.0/(1.0+exp((Y(1)*1000.0+26.0)/3.0));
if (f1_inf-Y(6) > 0.0)
    constf1 = 1.0+1433.0*(Y(3)-50.0*1.0e-6);
else
    constf1 = 1.0;
end
tau_f1      = (20.0+1102.5*exp(-((Y(1)*1000.0+27.0)^2.0/15.0)^2.0)+200.0/(1.0+exp((13.0-Y(1)*1000.0)/10.0))+180.0/(1.0+exp((30.0+Y(1)*1000.0)/10.0)))*constf1/1000.0;
dY(6, 1)    = (f1_inf-Y(6))/tau_f1;

f2_inf      = 0.33+0.67/(1.0+exp((Y(1)*1000.0+32.0)/4.0));
constf2     = 1.0;
tau_f2      = (600.0*exp(-(Y(1)*1000.0+25.0)^2.0/170.0)+31.0/(1.0+exp((25.0-Y(1)*1000.0)/10.0))+16.0/(1.0+exp((30.0+Y(1)*1000.0)/10.0)))*constf2/1000.0;
dY(7, 1)    = (f2_inf-Y(7))/tau_f2;

alpha_fCa   = 1.0/(1.0+(Y(3)/0.0006)^8.0);
beta_fCa    = 0.1/(1.0+exp((Y(3)-0.0009)/0.0001));
gamma_fCa   = 0.3/(1.0+exp((Y(3)-0.00075)/0.0008));
fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
if ((Y(1) > -0.06) && (fCa_inf > Y(8)))
    constfCa = 0.0;
else
    constfCa = 1.0;
end
tau_fCa     = 0.002;   % s (in i_CaL_fCa_gate)
dY(8, 1)    = constfCa*(fCa_inf-Y(8))/tau_fCa;

%% Ito
i_to        = g_to*(Y(1)-E_K)*Y(16)*Y(17);

q_inf       = 1.0/(1.0+exp((Y(1)*1000.0+53.0)/13.0));
tau_q       = (6.06+39.102/(0.57*exp(-0.08*(Y(1)*1000.0+44.0))+0.065*exp(0.1*(Y(1)*1000.0+45.93))))/1000.0;
dY(16, 1)   = (q_inf-Y(16))/tau_q;

r_inf       = 1.0/(1.0+exp(-(Y(1)*1000.0-22.3)/18.75));
tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(Y(1)*1000.0+30.61))+0.369*exp(-0.12*(Y(1)*1000.0+23.84))))/1000.0;
dY(17, 1)   = (r_inf-Y(17))/tau_r;

%% IKs
i_Ks        = g_Ks*(Y(1)-E_Ks)*Y(11)^2.0*(1.0+0.6/(1.0+(3.8*0.00001/Y(3))^1.4));

Xs_infinity = 1.0/(1.0+exp((-Y(1)*1000.0-20.0)/16.0));
alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-Y(1)*1000.0)/6.0));
beta_Xs     = 1.0/(1.0+exp((-60.0+Y(1)*1000.0)/20.0));
tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
dY(11, 1)   = (Xs_infinity-Y(11))/tau_Xs;

%% IKr
L0           = 0.025;	% -     (in i_Kr_Xr1_gate)
Q            = 2.3;     % -     (in i_Kr_Xr1_gate)
i_Kr         = g_Kr*(Y(1)-E_K)*Y(9)*Y(10)*sqrt(Ko/5.4);

V_half       = 1000.0*(-R*T/(F*Q)*log((1.0+Cao/2.6)^4.0/(L0*(1.0+Cao/0.58)^4.0))-0.019);

Xr1_inf      = 1.0/(1.0+exp((V_half-Y(1)*1000.0)/4.9));
alpha_Xr1    = 450.0/(1.0+exp((-45.0-Y(1)*1000.0)/10.0));
beta_Xr1     = 6.0/(1.0+exp((30.0+Y(1)*1000.0)/11.5));
tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
dY(9, 1)     = (Xr1_inf-Y(9))/tau_Xr1;

Xr2_infinity = 1.0/(1.0+exp((Y(1)*1000.0+88.0)/50.0));
alpha_Xr2    = 3.0/(1.0+exp((-60.0-Y(1)*1000.0)/20.0));
beta_Xr2     = 1.12/(1.0+exp((-60.0+Y(1)*1000.0)/20.0));
tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
dY(10, 1)    = (Xr2_infinity-Y(10))/tau_Xr2;

%% IK1
alpha_K1 = 3.91/(1.0+exp(0.5942*(Y(1)*1000.0-E_K*1000.0-200.0)));
beta_K1  = (-1.509*exp(0.0002*(Y(1)*1000.0-E_K*1000.0+100.0))+exp(0.5886*(Y(1)*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(Y(1)*1000.0-E_K*1000.0)));
XK1_inf  = alpha_K1/(alpha_K1+beta_K1);
i_K1     = g_K1*XK1_inf*(Y(1)-E_K)*sqrt(Ko/5.4);

%% INaCa
KmCa   = 1.38;	% mM  (in i_NaCa)
KmNai  = 87.5;	% mM  (in i_NaCa)
Ksat   = 0.1;   % -   (in i_NaCa)
gamma  = 0.35;	% -   (in i_NaCa)
kNaCa1 = kNaCa;	% A/F (in i_NaCa)
i_NaCa = kNaCa1*(exp(gamma*Y(1)*F/(R*T))*Y(18)^3.0*Cao-exp((gamma-1.0)*Y(1)*F/(R*T))*Nao^3.0*Y(3)*alpha)/((KmNai^3.0+Nao^3.0)*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*Y(1)*F/(R*T))));

%% INaK
Km_K  = 1.0;    % mM  (in i_NaK)
Km_Na = 40.0;   % mM  (in i_NaK)
PNaK1 = PNaK;   % A/F (in i_NaK)
i_NaK = PNaK1*Ko/(Ko+Km_K)*Y(18)/(Y(18)+Km_Na)/(1.0+0.1245*exp(-0.1*Y(1)*F/(R*T))+0.0353*exp(-Y(1)*F/(R*T)));

%% IpCa
KPCa  = 0.0005;   % mM  (in i_PCa)
i_PCa = g_PCa*Y(3)/(Y(3)+KPCa);

%% Background currents
i_b_Na = g_b_Na*(Y(1)-E_Na);

i_b_Ca = g_b_Ca*(Y(1)-E_Ca);

%% Sarcoplasmic reticulum
i_up = g_serca * VmaxUp/(1.0+Kup^2.0/Y(3)^2.0);

i_leak = (Y(2)-Y(3))*V_leak;

dY(4, 1) = 0;                                                              %<-dKi

% RyR
RyRSRCass = (1 - 1/(1 +  exp((Y(2)-0.3)/0.1)));
i_rel = g_irel_max*RyRSRCass*Y(22)*Y(23)*(Y(2)-Y(3));

RyRainfss = RyRa1-RyRa2/(1 + exp((1000*Y(3)-(RyRahalf))/0.0082));
RyRtauadapt = 1; %s
dY(21,1) = (RyRainfss- Y(21))/RyRtauadapt;

RyRoinfss = (1 - 1/(1 +  exp((1000*Y(3)-(Y(21)+ RyRohalf))/0.003)));
if (RyRoinfss>= Y(22))
	RyRtauact = 18.75e-3;       %s
else
	RyRtauact = 0.1*18.75e-3;   %s
end
dY(22,1) = (RyRoinfss- Y(22))/RyRtauact;

RyRcinfss = (1/(1 + exp((1000*Y(3)-(Y(21)+RyRchalf))/0.001)));
if (RyRcinfss>= Y(23))
	RyRtauinact = 2*87.5e-3;    %s
else
    RyRtauinact = 87.5e-3;      %s
end
dY(23,1) = (RyRcinfss- Y(23))/RyRtauinact;

%% Ca2+ buffering
Buf_C       = 0.25;   % mM (in calcium_dynamics)
Buf_SR      = 10.0;   % mM (in calcium_dynamics)
Kbuf_C      = 0.001;  % mM (in calcium_dynamics)
Kbuf_SR     = 0.3;    % mM (in calcium_dynamics)
Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/(Y(3)+Kbuf_C)^2.0);
Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/(Y(2)+Kbuf_SR)^2.0);

%% Ionic concentrations
%Nai
dY(18, 1) = -Cm*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa)/(F*Vc*1.0e-18);
%caSR
dY(3, 1)  = Cai_bufc*(i_leak-i_up+i_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*Cm/(2.0*Vc*F*1.0e-18));
%Cai
dY(2, 1)  = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));

%% Stimulation
i_stim = Y(end) ; 
dY(24, 1) = 0 ;

%% Membrane potential
dY(1, 1) = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);

%% Output variables
IK1     = i_K1;
Ito     = i_to;
IKr     = i_Kr;
IKs     = i_Ks;
ICaL    = i_CaL;
INaK    = i_NaK;
INa     = i_Na;
INaCa   = i_NaCa;
IpCa    = i_PCa;
If      = i_f;
IbNa    = i_b_Na;
IbCa    = i_b_Ca;
Irel    = i_rel;
Iup     = i_up;
Ileak   = i_leak;
Istim   = i_stim;
INaL    = i_NaL;

%% Return
data = [INa, If, ICaL, Ito, IKs, IKr, IK1, INaCa, INaK, IpCa, IbNa, IbCa, Irel, Iup, Ileak, Istim, E_K, E_Na, INaL];
%===============================================================================
% End of file
%===============================================================================
