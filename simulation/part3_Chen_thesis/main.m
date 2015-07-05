%%%%%%%%%%%%%%
%
% AlGaAs-GaAs quantum well bandgap energy calculation according to Chen's
% thesis (page 93). Theory bandgap = 1.5331eV.
%
% author: Xin Zhang
% date: May, 2015
%
%%%%%%%%%%%%%%

clear;clc;
global ELEMENTARY_CHARGE;
global ELECTRON_MASS;
global PLANCK_CONSTANT;
global H_BAR;
global delta_x;
global PI;
PI = 3.14159265359;
ELEMENTARY_CHARGE = 1.6e-19;
ELECTRON_MASS = 9.10938291e-31;
PLANCK_CONSTANT = 6.62606957e-34;
H_BAR = PLANCK_CONSTANT/(2*PI);
meV = ELEMENTARY_CHARGE*1e-3;

%%%
offset = 0.7;
well = 0;
barrier = 0.26;
L = 42e-10;
L_barrier = 200e-10;
%%%

delta_x = L/100;
x = (  -(L_barrier+L/2) : delta_x : (L_barrier+L/2)  );
Cx = composotionIII( 0,x,L,well,barrier );

%%% strain
% plot(x, strain(Cx));

Vc  =     offset * Eg(Cx) - (1-offset)*S_biaxial(Cx);
Vhh = (1-offset) * Eg(Cx) - (1-offset)*S_biaxial(Cx) + S_uniaxial(Cx);
Vlh = (1-offset) * Eg(Cx) - (1-offset)*S_biaxial(Cx) + 0.5*( S_uniaxial(Cx)+Delta0(Cx) ) - 0.5*sqrt( 9*S_uniaxial(Cx).^2 + Delta0(Cx).^2 - 2*S_uniaxial(Cx).*Delta0(Cx) );


%%% E
Ec  = TMM( x, Vc-min(Vc),   mc(Cx),  (1:1:200)*meV  ,1) + min(Vc);
Ehh = TMM( x, Vhh-min(Vhh), mHH(Cx), (1:1:200)*meV  ,1) + min(Vhh);
Elh = TMM( x, Vlh-min(Vlh), mLH(Cx), (1:1:200)*meV  ,1) + min(Vlh);


%%% transition E
CHH = ( Ec+Ehh )/meV
CLH = ( Ec+Elh )/meV