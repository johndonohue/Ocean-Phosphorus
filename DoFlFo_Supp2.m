%% Matlab supplement to "The development of deep-ocean anoxia in a comprehensive ocean phosphorus model"
% Authors: J. G. Donohue, B. J. Florio, and A. C. Fowler
% Corresponding Author: J. G. Donohue (john.donohue.nui@gmail.com), MACSI, University of Limerick, Limerick, Ireland
% Journal: Computers & Geosciences
%
% The script below computes the value of the critical anoxic parameter (\lambda_6*s_3 - \nu) for a given set of physical parameters.


% Carbon cycle
kCbur_prox  = 0.0905660;     
kCF1        = 38.6597938;      
kCF10       = 0.820594208;       
kCF11       = (496.6/(3600+28.0125));        
kCF12       = 8.8109658*1E-3;
kCF13       = 1.6/496.6;     
kCF4        = 1.0;              
kCF5        = 1.0589392;
kCF6        = 2.19836 ;
kCF7        = 2.7/(560.25+4.66); 
kCF8        = 1.0;              
kCF9        = 0.71893;          
kCrel_prox  = 6.9383815;        

% Phosphorus cycle
kCaP_prox   = 0.056559;   
kFeP_prox   = 0.925;     
korgP_prox  = 1.0;       
Kp          = 0.124;
kPF10       = 2.69703*1e-3;
kPF11       = 1.0;         
kPF12a      = 0.01;        
kPF12b      = 0.5;         
kPF13       = 2.18342;     
kPF14       = 0.00675/(0.044+5.28538);         
kPF15       = 1;
kPF19       = 0.811116;        
kPF20       = 0.01;             
kPF22       = 0.1197461;        
kPF23       = 0.5;              
kPF24       = 8.8267*1e-3;     
kPF25       = 0.5;             

kPF26       = 6.75*1e9;         
kPF27       = 2.8858*1e-3;      
kPF28       = 1.0;              
kPF29       = 0.0;

kPF4        = 1.0;   
kPF5a       = 0.01;   
kPF5b       = 0.5;      
kPF8        = 1.0;       
kPF9        = 0.00135238; 
kPrel_prox  = 7.4338699;

Redfield_CO2= 106/138;
Redfield_CP = 106;
CPoxdeep    = 237.0;
CPandeep    = 1100;

k_evap          = 37e12;    
k_river_water   = 37e12;     
k_river_P       = 0.09e12;   
k_fb_prox       = 0.0;       
k_Pfish_dis_surf= 0.0;       
k_f12c          = 0.0;       
k_O2_surf       = 0.325;     

KmO2            = 0.0001;     
C_NO3_init   = 0.0;  
C_RP_max     = 0.03;  
dcoeff       = 0.0;   
area_deep    = 3.32e14; 
dx           = 0.1;   
sedvel_t0    = 0.000; 
C_bur_deep_t0= 1.6e12;
kredox       = 1e8;   
kprec        = 1e-3;  


%define mixing parameters
vmix_oce = 0.1;
vmix_cstl = 0.1;

%Water constants:
W1=36e12;    
W2=3600e12; 
W3=49830e12;
W4=1297670e12;

%Water fluxes
r = k_river_water;
w12 = k_river_water;
w43 = 3780e12*vmix_oce;
w42 = 378e12*vmix_cstl;
E =  k_evap;
w23 = w12 + w42;
w34 = w43 + w42;

%a constants:
a1 = W4*kCF12/Redfield_CO2;
a2 = W4*kredox;
a3 = (1-kCbur_prox)*kCF1*Redfield_CP*W1;
a4 = kCrel_prox*W1 + kCF4*w12;
a5 = (1 - kCF7)*kCF4*w12;
a6 = (1 - kCF7)*kCF5*Redfield_CP*W2;
a7 = w23*kCF8 + kCF6*W2;
a8 = (1 - kCF11)*kCF9*Redfield_CP*W3;
a9 = (1 - kCF11)*w23*kCF8;
a10 = kCF10*W3;
a11 = kCF11*w23*kCF8;
a12 = kCF11*kCF9*Redfield_CP*W3;
a13 = kCF12*W4;
a14 = kPrel_prox*W1*(1-kCaP_prox);
a15 = kPF5b*W1;
a16 = (kCF1+kFeP_prox)*W1 + kPF4*w12;
a17 = kCF1*W1*(1-(korgP_prox/400)*kCbur_prox*Redfield_CP - kPF5a);
a18 = kPrel_prox*W1 + kPF8*w12;
a19 = kPF5a*kCF1*W1;
a20 = kPF5b*W1;
a21 = kPF4*w12;
a22 = (1-kPF10)*kPF13*W2;
a23 = kPF12b*W2;
a24 = kPF11*w23+W2*(kPF9 + kCF5);
a25 = kCF5*W2*(1-kPF14 - kPF12a);
a26 = kPF8*w12*(1-kPF14);
a27 = kPF15*w23 + kPF13*W2;
a28 = kPF12a*kCF5*W2;
a29 = kPF12b*W2;
a30 = kPF19*W3;
a31 = kPF11*w23;
a32 = kCF9*W3 + w34;
a33 = kCF9*W3*(1 - kPF20 - kCF11);
a34 = kPF15*w23;
a35 = kCF11*kCF8*w23/Redfield_CP;
a36 = kPF19*W3;
a37 = kPF20*kCF9*W3;
a38 = kPF23*W3;
a39 = kPF24*W4*(1 - kPF27);
a40 = kPF25*W4;
a41 = kCF11*kCF8*w23;
a42 = kCF11*kCF9*W3*Redfield_CP;
a43 = kPF24*W4;
a44 = W3*kPF23*(1-kPF29);
a45 = kPF25*W4;
a46 = (kCF12/Redfield_CO2)*W4;
a47 = dcoeff*C_RP_max*kCbur_prox*area_deep/dx;
a49 = kCF8*w23*sedvel_t0/C_bur_deep_t0*area_deep*kCF13*kCF11;
a50 = kCF9*Redfield_CP*W3*sedvel_t0/C_bur_deep_t0*area_deep*kCF13*kCF11;
a51 = kredox*W4;
a52 = kprec*W4;

%miscellanea
h28 = kCF13;
Km = KmO2;
gs = k_O2_surf;
h1 = k_river_P;
h31 = k_fb_prox;
h21 = k_Pfish_dis_surf;
h33 = k_f12c;
h26 = kPF26;

% Further reduced quantities:
b0 = w34/W4;
b1 = a1/W4;
b2 = a2/W4;
b3 = a3/W1;
b4 = a4/W1;
b5 = a5/W2;
b6 = a6/W2;
b7 = a7/W2;
b8 = a8/W3;
b9 = a9/W3;
b10 = a10/W3;
b11 = a11/W4;
b12 = a12/W4;
b13 = a13/W4;
b14 = a14/W1;
b15 = a15/W1;
b16 = a16/W1;
b17 = a17/W1;
b18 = a18/W1;
b19 = a19/W1;
b20 = a20/W1;
b21 = a21/W2;
b22 = a22/W2;
b23 = a23/W2;
b24 = a24/W2;
b25 = a25/W2;
b26 = a26/W2;
b27 = a27/W2;
b28 = a28/W2;
b29 = a29/W2;
b30 = a30/W3;
b31 = a31/W3;
b32 = a32/W3;
b33 = a33/W3;
b34 = a34/W3;
b35 = a35/W3;
b36 = a36/W3;
b37 = a37/W3;
b38 = a38/W3;
b39 = a39/W4;
b40 = a40/W4;
b41 = a41/W4;
b42 = a42/W4;
b43 = a43/W4;
b44 = a44/W4;
b45 = a45/W4;
b46 = a46/W4;
b47 = a47/W4;
b49 = a49/W4;
b50 = a50/W4;
b52 = a52/W4;
b53 = h1/W1; 
b54 = w42/W2;
b55 = h33/W2;
b56 = w43/W3;
b57 = h21/W3;
b58 = w34/W4;
b59 = h26/W4;
b60 = h31/W1;

d16 = b16 - b14*b17/b18;
d24 = b24 - b22*b25/b27;
d32b = b32 - b30*b33/b36 - b50*b39/b58*b42/(b43*Redfield_CP);
f32b = d24*d32b/b31-b54*b42*b39/(b43*Redfield_CP*b58);
e32b = f32b - b23*b28*d32b/(b29*b31);

%scales
S1 = b53/d16;
P1 = b17/b18*S1;
C1 = b3/b4*S1;
F1 = b19/b20*S1;
S3 = b21/e32b*S1;
S2 = d32b/b31*S3;
S4 = b39/b58*b42/(Redfield_CP*b43)*S3;
P4 = b42/(Redfield_CP*b43)*S3;
P2 = b25/b27*S2;
F2 = b28/b29*S2;
F3 = b37/b38*S3;
F4 = b44/b45*F3;
P3 = b33/b36*S3;
C2 = b6/b7*S2;
C3 = b8/b10*S3;
C4 = b12/b13*S3;
R1 = C_RP_max;
G4 = 10^6*b0*gs/(b2*R1);
R2 = 0;

%dimensionless parameters;
%; lambda
lam1 = b53/(b16*S1);
lam2 = b14*P1/(b16*S1);
lam3 = b22*P2/(b24*S2);
lam4 = b30*P3/(b32*S3);
lam5 = b59/(b58*S4);
lam11 = Km/G4;
lam20 = b58*R1/(b46*C4);
%; delta
del1 = b56*S4/(b32*S3);
del2 = b31*S2/(b32*S3);
del3 = b2*R1*G4/(10^6*b52);
del4 = b5*C1/(b7*C2);
del5 = b26*P1/(b27*P2);
%epsilon
eps1 = b15*F1/(b16*S1);
eps2 = b21*S1/(b24*S2);
eps3 = b23*F2/(b24*S2);
eps4 = b54*S4/(b24*S2);
eps6 = b40*F4/(b58*S4);
eps7 = S3/S4;
eps8 = b34*P2/(b36*P3);
eps9 = b35*C2/(b36*P3);
eps10 = h28*b41*C2/(237*b43*P4);
eps11 = h28*b42*S3/(237*b43*P4);
eps13 = b9*C2/(b10*C3);
eps14 = b11*C2/(b13*C4);
eps15 = h28*b11*C2/(b13*C4);
eps16 = h28*b12*S3/(b13*C4);
eps19 = G4/gs;
eps20 = b58/b16;
eps21 = b58/b24;
eps22 = b58/b32;
eps23 = b58/b18;
eps24 = b58/b27;
eps25 = b58/b36;
eps26 = b58/b43;
eps27 = b58/b20;
eps28 = b58/b4;
eps29 = b58/b7;
eps30 = b58/b10;
eps31 = b58/b13;
eps32 = b58*G4/(b0*gs);
eps33 = b41*C2/(Redfield_CP*b43*P4);
eps34 = b58/b29;
eps35 = b58/b38;
eps36 = b58/b45;
eps37 = b58*R1/b52;
eps38 = b46*C4/b52;
eps39 = b1*C4/(b0*gs);


%Rescales:
Sb1 = 1;
Pb1 = 1;
Cb1 = 1;
Fb1 = 1;
Sb3 = 1/((1-lam4-del1)*(1-lam3-eps3));
Sb2 = 1/(1-lam3-eps3);
Sb4 = Sb3;
Pb4 = Sb3;
Pb2 = Sb2;
Fb2 = Sb2;
Fb3 = Sb3;
Fb4 = Sb3;
Pb3 = Sb3;
Cb2 = Sb2;
Cb3 = Sb3;
Cb4 = Sb3;
Rb1 = 1;
Gb4 = 1;
Rb2 = 0;

%rescaled ND coefficients
lam6 = eps39*Cb4;
lam8 = eps38*Sb3;
eps40 = del4/Sb2;
eps41 = eps13*Sb2/Sb3;
eps42 = eps14*Sb2/Sb3;
eps43 = eps15*Sb2/Sb3;
eps44 = eps2/Sb2;
eps45 = eps4*Sb3/Sb2;
eps46 = del5/Sb2;
eps47 = del2*Sb2/Sb3;
eps48 = eps8*Sb2/Sb3;
eps49 = eps9*Sb2/Sb3;
eps50 = lam5/Sb4;
eps51 = eps33*Sb2/Sb3;
eps52 = eps10*Sb2/Sb3;
eps53 = 1- lam3;
eps54 = 1 - del1 - lam4;

del6=  (eps44+eps46)/eps45;
del7 = (eps47+eps48-eps49)/del1;
lam9 = 1 +eps6+eps7;
lam10 = (eps53-eps3)/(eps45);
lam12n = 1+eps54/del1;
eps55=  (eps21+eps24)/(eps45);
eps56=   (eps22+eps25)/(del1);


s_4_approx = (lam9*del6*del7)/(lam10*(lam12n-lam9)-del7*lam9);
s_3_approx = ( (lam10+del7)*s_4_approx+del6*del7)/(lam10*lam12n);

anoxia_parameter = lam6*s_3_approx
