%%%%%%%%%%%%
%
% A simple simulation to repeat some results of IEEE JOURNAL OF QUANTUM
% ELECTRONICS, VOL. 34, NO. 1, JANUARY 1998. The conduction band energies
% of quantum wells changed with Ld are calculated.
%
% author: Xin Zhang
% date: May, 2015
%
%%%%%%%%%%%%

clc;

global ELEMENTARY_CHARGE;
global ELECTRON_MASS;
global PLANCK_CONSTANT;
global H_BAR;
global delta_x;

%%% constants
ELEMENTARY_CHARGE = 1.6e-19;
ELECTRON_MASS = 9.11e-31;
PLANCK_CONSTANT = 6.63e-34;
H_BAR = PLANCK_CONSTANT/(2*pi);

%%% variable
L = 10e-9; % well width
well = 0; % GaAs
barrier = 0.3; % Al_0.3Ga_0.7As

%%% x
delta_x = L/100;
x = (-2*L:delta_x:2*L);

%%% find -L/2 and L/2
temp = find(x>=-L/2);
well_start = temp(1);
temp = find(x<=L/2);
well_stop = temp(length(temp));

Ld = (0:5:50)*1e-10;
levels = 2;
eigenvalues = zeros(length(Ld),levels);
for i=Ld
   %%% Ld
   Cx = composotionIII( i,x,L,well,barrier );

   %%% potential energy
   % Vx = 0.7*(  Eg(Cx) - min(Eg(Cx))  ); 
   Vx = 0.3 * (Eg(Cx) - min(Eg(Cx)));
   % figure(1);plot(x,Vx);hold on;
   
   %%% effective mass
   % Mx = mc(Cx); 
   Mx = mlh(Cx);
   % figure(2);plot(x, Mx);hold on;
   
   %%% energy level
   EE = (1 : 1 : 200)*ELEMENTARY_CHARGE*1e-3;
   eigenvalues(Ld==i,:) = TMM(x,Vx,Mx,EE,levels);
   % eigenvalues(Ld==i,:)/1e-3/ELEMENTARY_CHARGE
end

eigenvalues = eigenvalues/1e-3/ELEMENTARY_CHARGE;
Ld_plot = Ld/1e-10;
figure(3);
plot(Ld_plot,eigenvalues(:,1),  Ld_plot,eigenvalues(:,2))

Elh = eigenvalues;