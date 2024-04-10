function [deriv, currents] = f_Paci13_a(t, statevar, geometry, conductances, parameters)

global F R T Nao Ko Cao

% geometry
Vsr = geometry.V_SR ;
Vc = geometry.Vc ;
Cm  = geometry.Cm ;

% conductances
c = conductances.baseline .* conductances.scaling ;
for i = 1: length(c)
    if conductances.applicationTimes(1,i) <= t && conductances.applicationTimes(2,i) > t
        c(i) = c(i) .* conductances.drugEffects(i) ;
    end
end
GNa =    c(1) ;
Gf =     c(2) ;
GCaL =   c(3) ;
Gto =    c(4) ;
GKs =    c(5) ;
GKr =    c(6) ;
GK1 =    c(7) ;
GpCa =   c(8) ;
GbNa =   c(9) ;
GbCa =   c(10) ;
Vmaxup = c(11) ;
G_RyR =  c(12) ;
KNaCA =  c(13) ;
PNaK =   c(14) ;
Vleak =  c(15) ;

% parameters
p = parameters.baseline .* parameters.scaling ;
arel =      p(1) ;
brel =      p(2) ;
crel =      p(3) ;
Bufc =      p(4) ;
Bufsr =  	p(5) ;
Kbufc =     p(6) ;
Kbufsr =	p(7) ;
Kup =       p(8) ;
KpCa =      p(9) ;
L0 =        p(10) ;
Pkna =      p(11) ;
Ksat =      p(12) ;
KmCa =      p(13) ;
KmNai =     p(14) ;
alpha =     p(15) ;
gamma =     p(16) ;
KmNa =      p(17) ;
Kmk =       p(18) ;

% intracellular ion concentrations -> all change over time except Ki.
Ki = 150.0 ; % mM

%Reversal potentials

ENa = ((R*T)/(F))*(log(Nao/statevar(16))) ; %0.063 V
Ek = ((R*T)/(F))*(log(Ko/Ki)) ; %-0.088 V
Eks = ((R*T)/(F))*log((Ko+(Pkna*Nao))/(Ki+(Pkna*statevar(16)))) ; %-0.072 V
ECa = ((0.5*R*T)/(F))*log(Cao/statevar(17)) ; %0.121 V
Ef = -0.017 ; % V

%% Na Current, INa

INa = GNa*(statevar(3)^3)*statevar(1)*statevar(2)*(statevar(15)-ENa) ;

% INa, h gate equations

hinf = (1+exp((statevar(15)*1000+72.1)/5.7))^(-1/2) ;

if statevar(15) < -0.040
    alphah = 0.057*exp(-(statevar(15)*1000.0+80.0)/6.8) ;
    betah = 2.7*exp(0.079*statevar(15)*1000.0)+3.1*10.0^5.0*exp(0.3485*statevar(15)*1000.0) ;
    tauh = 1.5/((alphah+betah)*1000.0) ;
else
    alphah = 0 ;
    betah = 0.77/(0.13*(1.0+exp((statevar(15)*1000.0+10.66)/-11.1))) ;
    tauh  = 2.542/1000 ;
end

dhdt = (hinf - statevar(1))/tauh ;

% INa, j gate

jinf = 1.0/sqrt(1.0+exp((statevar(15)*1000.0+72.1)/5.7)) ;

if statevar(15) < -40e-3
    %alphaj = (-25428.0*exp(0.2444*statevar(15)*1000.0)-6.948*10.0^-6.0*exp(-0.04391*statevar(15)*1000.0))*(statevar(15)*1000.0+37.78)/(1.0+exp(0.311*(statevar(15)*1000.0+79.23))) ; %ignore supplementals, attempting consistency with Paci code.
    alphaj = (-25428.0*exp(0.2444*statevar(15)*1000.0)-6.948*10.0^-6.0*exp(-0.04391*statevar(15)*1000.0))*(statevar(15)*1000.0+37.78)/(1.0+exp(0.311*(statevar(15)*1000.0+79.23)));
    betaj = 0.02424*exp(-0.01052*statevar(15)*1000.0)/(1.0+exp(-0.1378*(statevar(15)*1000.0+40.14))) ;
    %betaj = (0.02424*exp(-0.01052*statevar(15)*1000))/(1+exp(-0.1378*(statevar(15)*1000+40.14))) ;
else
    alphaj = 0 ;
    %betaj = (0.6*exp(0.057*statevar(15)*1000))/(1+exp(-0.1*(statevar(15)*1000+32))) ;
    betaj = 0.6*exp(0.057*statevar(15)*1000.0)/(1.0+exp(-0.1*(statevar(15)*1000.0+32.0)));
end

tauj = 7.0/((alphaj+betaj)*1000.0) ;

djdt = (jinf - statevar(2))/tauj ;

%INa, m gate

%         minf = (1+exp((-34.1-statevar(15)*1000)/5.9))^(-1/3) ;
%
%         alpham = (1+exp((-60-statevar(15)*1000)/5))^(-1) ;
%
%         betam = (0.1 / (1+exp((statevar(15)*1000+35)/5))) + (0.1 / (1+exp((statevar(15)*1000-50)/200))) ;
%
%         taum = (alpham*betam)/1000 ;
%
%         dmdt = ((minf - statevar(3))/taum) ;

minf = 1.0/(1.0+exp((-statevar(15)*1000.0-34.1)/5.9))^(1.0/3.0);

alpham = 1.0/(1.0+exp((-statevar(15)*1000.0-60.0)/5.0));

betam = 0.1/(1.0+exp((statevar(15)*1000.0+35.0)/5.0))+0.1/(1.0+exp((statevar(15)*1000.0-50.0)/200.0));

taum = 1.0*alpham*betam/1000.0;

dmdt = (minf-statevar(3))/taum;

%% L-type Calcium current

%ICaL = (GCaL*4*statevar(15)*(F^2)/(R*T))*((statevar(17)*exp(2*statevar(15)*F/(R*T)))-(0.341*Cao)) / (exp((2*statevar(15)*F)/(R*T))-1)*statevar(4)*statevar(5)*statevar(6)*statevar(7) ;
ICaL = GCaL*4.0*statevar(15)*F^2.0/(R*T)*(statevar(17)*exp(2.0*statevar(15)*F/(R*T))-0.341*Cao)/(exp(2.0*statevar(15)*F/(R*T))-1.0)*statevar(4)*statevar(5)*statevar(6)*statevar(7) ;

% ICaL, d gate

dinf = (1+exp((-(statevar(15)*1000+5.986))/7))^(-1) ;
alphad = 0.25+1.4/(1.0+exp((-statevar(15)*1000.0-35.0)/13.0));
betad = 1.4/(1.0+exp((statevar(15)*1000.0+5.0)/5.0));
gammad = 1.0/(1.0+exp((-statevar(15)*1000.0+50.0)/20.0));
taud = (alphad*betad+gammad)*1.0/1000.0;

%         alphad = (1.4 / (1+exp((-35-statevar(15)*1000)/13))) + 0.25 ;
%
%         betad = (1.4 / (1+exp((statevar(15)*1000+5)/5))) ;
%
%         gammad = (1 / (1+exp((-statevar(15)*1000+50)/20))) ;
%
%         taud = ((alphad*betad)+gammad)/1000 ;

dddt = ((dinf - statevar(4)) / taud) ;

% ICaL, fCa gate

%alphafCa = (1 + ((statevar(17)/0.0006)^8))^(-1) ;

%betafCa = 0.1/(1 + exp((statevar(17)-0.0009)/0.0001)) ;

%gammafCa = 0.3/(1 + exp((statevar(17)-0.00075)/0.0008)) ;

%fCainf = (alphafCa + betafCa + gammafCa) / (1.3156) ;

alphafCa = 1.0/(1.0+(statevar(17)/0.0006)^8.0);
betafCa = 0.1/(1.0+exp((statevar(17)-0.0009)/0.0001));
gammafCa = 0.3/(1.0+exp((statevar(17)-0.00075)/0.0008));
fCainf = (alphafCa+betafCa+gammafCa)/1.3156;

taufCa = 2/1000 ; % put into seconds

if (fCainf > statevar(5)) && (statevar(15) > -0.060)
    dfCadt = 0 ;
else
    dfCadt = ((fCainf - statevar(5)) / taufCa) ;
end

% ICaL, f1 gateF
f1inf = (1+exp((statevar(15)*1000+25.226)/3))^(-1) ;
if f1inf - statevar(6) > 0.0
    tauf1constant = 1.0 + 1433.0*(statevar(17)-(50e-6)) ;
else
    tauf1constant = 1.0 ;
end

%tauf1 = (20.0+1102.5*exp(-((statevar(15)*1000.0+27.0)^2.0/15.0)^2.0)+200.0/(1.0+exp((13.0-statevar(15)*1000.0)/10.0))+180.0/(1.0+exp((30.0+statevar(15)*1000.0)/10.0)))*tauf1constant/1000.0 ;
tauf1 = (20.0+1102.5*exp(-((statevar(15)*1000.0+27.0)^2.0/15.0)^2.0)+200.0/(1.0+exp((13.0-statevar(15)*1000.0)/10.0))+180.0/(1.0+exp((30.0+statevar(15)*1000.0)/10.0)))*tauf1constant/1000.0;

df1dt = ((f1inf - statevar(6)) / tauf1) ;

% ICaL, f2 gate
f2inf = (0.67/(1+exp((statevar(15)*1000+31.226)/4)))+0.33 ;
tauf2 = (((600*exp(-(((statevar(15)*1000+25)^2)/170))) + (31/(1+exp((-statevar(15)*1000+25)/10))) +  (16/(1+exp((statevar(15)*1000+30)/10))))*2)/1000 ;

df2dt = ((f2inf - statevar(7)) / tauf2) ;

% Transient outward current, Ito

Ito = (Gto*statevar(8)*statevar(9)*(statevar(15) - Ek)) ;

% Ito, r gate

rinf = 1.0/(1.0+exp(-(statevar(15)*1000.0-22.3)/18.75));

taur = (2.75352+14.40516/(1.037*exp(0.09*(statevar(15)*1000.0+30.61))+0.369*exp(-0.12*(statevar(15)*1000.0+23.84))))/1000.0;

drdt = (rinf - statevar(8))/(taur) ;

% Ito, q gate

qinf = (1+exp((statevar(15)*1000+53)/13))^(-1) ;

tauq = ((39.102 / ((0.57*exp(-0.08*(statevar(15)*1000+44))) + (0.065*exp(0.01*(statevar(15)*1000+45.93))))) + 6.06)/1000 ;

dqdt = ((qinf - statevar(9)) / tauq) ;

%% Rapid delayed rectifier K+ current, IKr

IKr = GKr*((Ko/5.4)^.5)*statevar(10)*statevar(11)*(statevar(15)-Ek) ;

% IKr, Xr1 gate

Vonehalf = 1000.0*(-R*T/(F*2.3)*log((1.0+Cao/2.6)^4.0/(L0*(1.0+Cao/0.58)^4.0))-0.019);

xr1inf = 1.0/(1.0+exp((Vonehalf-statevar(15)*1000.0)/4.9));

alphaxr1 = 450.0/(1.0+exp((-45.0-statevar(15)*1000.0)/10.0));

betaxr1 = 6.0/(1.0+exp((30.0+statevar(15)*1000.0)/11.5));

tauxr1 = 1.0*alphaxr1*betaxr1/1000.0;

dXr1dt = ((xr1inf - statevar(10)) / tauxr1) ;

%IKr, Xr2 gate

xr2inf = 1.0/(1.0+exp((statevar(15)*1000.0+88.0)/50.0));

alphaxr2 = 3.0/(1.0+exp((-60.0-statevar(15)*1000.0)/20.0));

betaxr2 = 1.12/(1.0+exp((-60.0+statevar(15)*1000.0)/20.0));

tauxr2 = 1.0*alphaxr2*betaxr2/1000.0;

dXr2dt = ((xr2inf-statevar(11)) / (tauxr2)) ;

%% Slow delayed rectifier K+ current, IKs

IKs = GKs*(statevar(12)^2)*(1+(0.6/(1+(((3.8e-5)/statevar(17))^1.4))))*(statevar(15)-Eks) ;

% IKs, Xs gate

XSinf = (1+exp((-20-statevar(15)*1000)/16))^(-1) ;

alphaxs = 1100/((1+exp((-10-statevar(15)*1000)/6))^.5) ;

betaxs = (1+exp((-60+statevar(15)*1000)/20))^(-1) ;

tauxs = (alphaxs*betaxs)/1000 ;

dXsdt = ((XSinf - statevar(12)) / tauxs) ;

%Inward rectifier K+ current, IK1

alphaK1 = 3.91/(1.0+exp(0.5942*(statevar(15)*1000.0-Ek*1000.0-200.0)));

betaK1 = (-1.509*exp(0.0002*(statevar(15)*1000.0-Ek*1000.0+100.0))+exp(0.5886*(statevar(15)*1000.0-Ek*1000.0-10.0)))/(1.0+exp(0.4547*(statevar(15)*1000.0-Ek*1000.0)));

XK1inf = alphaK1/(alphaK1+betaK1);

IK1 = GK1*XK1inf*(statevar(15)-Ek)*sqrt(Ko/5.4);

% Hyperpolarization activated funny current, If

If = Gf*statevar(13)*(statevar(15)-Ef) ;

%If, Xf gate

xfinf = (1+exp((statevar(15)*1000+77.85)/5))^(-1) ;

tauf = 1900/(1+exp((statevar(15)*1000+15)/10)) ;

dXfdt = ((xfinf - statevar(13)) / tauf) ;

% Na+K+ pump current, INaK

INaK = PNaK*Ko/(Ko+Kmk)*statevar(16)/(statevar(16)+KmNa)/(1.0+0.1245*exp(-0.1*statevar(15)*F/(R*T))+0.0353*exp(-statevar(15)*F/(R*T)));

% Na+/Ca2+ exchanger current, INaCa

INaCa = KNaCA*(exp(gamma*statevar(15)*F/(R*T))*statevar(16)^3.0*Cao-exp((gamma-1.0)*statevar(15)*F/(R*T))*Nao^3.0*statevar(17)*alpha)/((KmNai^3.0+Nao^3.0)*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*statevar(15)*F/(R*T))));
% Ca2+ dynamics
Irel = G_RyR* ( (((arel*(statevar(18)^2))/((brel^2)+(statevar(18)^2)))+crel)*statevar(4)*statevar(14)*0.0556 ) ;

Iup = Vmaxup/(1.0+Kup^2.0/statevar(17)^2.0);

Ileak = (statevar(18)-statevar(17))*Vleak;

if (statevar(17) <= 0.00035)
    ginf = 1.0/(1.0+(statevar(17)/0.00035)^6.0);
else
    ginf = 1.0/(1.0+(statevar(17)/0.00035)^16.0);
end

taug = 0.002 ; %ms-->seconds

if ((ginf > statevar(14)) && (statevar(15) > -0.060))
    dgdt = 0 ;
else
    dgdt = ((ginf - statevar(14))/taug) ;
end

Caibufc = 1.0/(1.0+Bufc*Kbufc/(statevar(17)+Kbufc)^2.0);
Casrbufsr = 1.0/(1.0+Bufsr*Kbufsr/(statevar(18)+Kbufsr)^2.0);

IPCa = (GpCa*statevar(17))/(statevar(17)+KpCa) ;  %Ca2+ pump current, IpCa

IbCa = GbCa*(statevar(15)-ECa) ; %Background current

dCaidt = Caibufc*(Ileak-Iup+Irel-(((ICaL+IbCa+IPCa-(2*INaCa))/(2*Vc*F*1e-18))*Cm)) ; %RANDOM DIVISION BY 18TH ORDERS IN PACI CODE?

dCasrdt = Casrbufsr*Vc/Vsr*(Iup-(Irel+Ileak));

% Background currents

IbNa = GbNa*(statevar(15)-ENa) ;

% Sodium dynamics

dNaidt = -Cm*((INa+IbNa+(3*INaK)+(3*INaCa))/(F*Vc*1e-18)) ; % another random division

% Stimulation
Istim = -statevar(end) ;
dStim = Istim ;

%%%%%% SOLVE DERIVATIVE %%%%%%%

Iion = (IK1 + Ito + IKr + IKs + ICaL + INaK + INa + INaCa + IPCa + If + IbNa + IbCa) ;
dVdt = - (Iion + Istim) ;

currents = [INa, If, ICaL, Ito, IKs, IKr, IK1, INaCa, INaK, IPCa, IbNa, IbCa, Istim, Ek, ENa] ; % things in Paci18 that are missing here: Irel, Iup, Ileak, INaL

deriv = [dhdt; djdt; dmdt; dddt; dfCadt; df1dt; df2dt; drdt; dqdt; dXr1dt; dXr2dt; dXsdt; dXfdt; dgdt; dVdt; dNaidt; dCaidt; dCasrdt; dStim] ;

end


