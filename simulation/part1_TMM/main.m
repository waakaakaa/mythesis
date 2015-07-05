%%%%%%%%%%%%%%%%%%%
%
% A simulation program to calculate the band structure of AlGaAs quantum
% wells according to the paper of IEEE Journal of Quantum Electronics, vol.
% 26, no. 11, November 1990.
%
% Author: Xin Zhang
% Date: May, 2015
% Thanks to: Yiming Xia
%
%%%%%%%%%%%%%%%%%%%
clear;clc;
global H_BAR;
global delta_x;
global PI;

%%% constants
PI = 3.14159265359;
ELEMENTARY_CHARGE = 1.6e-19;
ELECTRON_MASS = 9.10938291e-31;
PLANCK_CONSTANT = 6.62606957e-34;
H_BAR = PLANCK_CONSTANT/(2*PI);

%%% variable
L = 20e-9; % well width

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% x
delta_x = L/1000;
x = (-2*L:delta_x:2*L);

%%% find -L/2 and L/2
temp = find(x>=-L/2);
well_start = temp(1);
temp = find(x<=L/2);
well_stop = temp(length(temp));

%%% find -L and L
temp = find(x>=-L);
x_start = temp(1);
temp = find(x<=L);
x_stop = temp(length(temp));

%%% potential energy
Vx = x;
Vx(:) = 0.225*ELEMENTARY_CHARGE;
Vx(well_start:well_stop) = 0;

%%% effective mass
Mx = x;
Mx(:) = 0.0919*ELECTRON_MASS;
Mx(well_start:well_stop) = 0.067*ELECTRON_MASS;

%%% get ready for plot wave function
x_new = x(x_start:x_stop)/1e-9;
Vx_new = Vx(x_start:x_stop)/1e-3/ELEMENTARY_CHARGE;
figure(1);
subplot(2,2,1);hold on;
subplot(2,2,2);hold on;
subplot(2,2,3);hold on;
subplot(2,2,4);hold on;

%%% energy level 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EE = (9.9 : 0.001 : 10)*ELEMENTARY_CHARGE*1e-3;
eigenvalue = TMM(x,Vx,Mx,EE);
eigenvalue/1e-3/ELEMENTARY_CHARGE

W = wave_function( x,Vx,Mx,eigenvalue );
W_new = abs(W(x_start:x_stop));
W_new = W_new/max(W_new);
subplot(2,2,1);
[AX,H1,H2] = plotyy(x_new,Vx_new,   x_new,W_new);
set(AX,'fontsize',24,'fontname','Times New Roman')
set(get(AX(1),'Ylabel'),'String','potential energy (meV)','fontsize',24,'fontname','Times New Roman' )
set(get(AX(2),'Ylabel'),'String','normalized |\psi |^2','fontsize',24,'fontname','Times New Roman')
xlabel('x (nm)');

%%% energy level 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EE = (39 : 0.01 : 40)*ELEMENTARY_CHARGE*1e-3;
eigenvalue = TMM(x,Vx,Mx,EE);
eigenvalue/1e-3/ELEMENTARY_CHARGE

W = wave_function( x,Vx,Mx,eigenvalue );
W_new = abs(W(x_start:x_stop));
W_new = W_new/max(W_new);
subplot(2,2,2);
[AX,H1,H2] = plotyy(x_new,Vx_new,   x_new,W_new);
set(AX,'fontsize',24,'fontname','Times New Roman')
set(get(AX(1),'Ylabel'),'String','potential energy (meV)','fontsize',24,'fontname','Times New Roman' )
set(get(AX(2),'Ylabel'),'String','normalized |\psi |^2','fontsize',24,'fontname','Times New Roman')
xlabel('x (nm)');

%%% energy level 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EE = (88 : 0.01 : 89)*ELEMENTARY_CHARGE*1e-3;
eigenvalue = TMM(x,Vx,Mx,EE);
eigenvalue/1e-3/ELEMENTARY_CHARGE

W = wave_function( x,Vx,Mx,eigenvalue );
W_new = abs(W(x_start:x_stop));
W_new = W_new/max(W_new);
subplot(2,2,3);
[AX,H1,H2] = plotyy(x_new,Vx_new,   x_new,W_new);
set(AX,'fontsize',24,'fontname','Times New Roman')
set(get(AX(1),'Ylabel'),'String','potential energy (meV)','fontsize',24,'fontname','Times New Roman' )
set(get(AX(2),'Ylabel'),'String','normalized |\psi |^2','fontsize',24,'fontname','Times New Roman')
xlabel('x (nm)');

%%% energy level 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EE = (155 : 0.01 : 156)*ELEMENTARY_CHARGE*1e-3;
eigenvalue = TMM(x,Vx,Mx,EE);
eigenvalue/1e-3/ELEMENTARY_CHARGE

W = wave_function( x,Vx,Mx,eigenvalue );
W_new = abs(W(x_start:x_stop));
W_new = W_new/max(W_new);
subplot(2,2,4);
[AX,H1,H2] = plotyy(x_new,Vx_new,   x_new,W_new);
set(AX,'fontsize',24,'fontname','Times New Roman')
set(get(AX(1),'Ylabel'),'String','potential energy (meV)','fontsize',24,'fontname','Times New Roman' )
set(get(AX(2),'Ylabel'),'String','normalized |\psi |^2','fontsize',24,'fontname','Times New Roman')
xlabel('x (nm)');