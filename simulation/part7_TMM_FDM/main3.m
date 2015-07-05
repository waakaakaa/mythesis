clear;clc;
global h_bar;
global h;
global H;
global length_z;
global bx;
global by;
global wx;
global wy;
global offset;
global z;
global speed_of_light;
h_bar = 6.6261e-34 /2 /pi;
speed_of_light = 2.9979e8;
offset = 0.6;
%% QW dimension
length_z = 0.01e-9;
% h = 6e-9;
% H = 20e-9;
% wx = 0.53;
% wy = 1;
% % [bx,by] = q2xy(1.25);
% bx = 0.74;
% by = 0.57;
h = 3.5e-9;
H = 20e-9;
wx = 0.71;
wy = 0.61;
bx = 1;
by = 0;
z = -(H+h)/2 : length_z : (H+h)/2;
%% variables
LdIII = (0:1:6)*1e-9;
k = 0.25;
wavelengthCHH = [];
wavelengthCLH = [];
%% figure
figure(1);
xlabel('L_{III} (nm)');
ylabel('Transition Wavelength(nm)');
xlim([LdIII(1) LdIII(length(LdIII))]*1e9);
box on;hold on;
pause(0.1);
%% for
for Ld=LdIII
    %% Composition
    CIn = compositionIII(Ld);
    CAs = compositionV(Ld*k);
    %% Potential
    Pc = vc(CIn,CAs);
    Phh = vhh(CIn,CAs);
    Plh = vlh(CIn,CAs);
    %% Energy
    EE = (1:1:100)*1e-3*1.6e-19;
    
    Ec  = TMM( z,Pc-min(Pc),  mc(CIn,CAs),      EE ,1)+min(Pc);
    Ehh = TMM( z,Phh-min(Phh),mhh_perp(CIn,CAs),EE ,1)+min(Phh);
    Elh = TMM( z,Plh-min(Plh),mlh_perp(CIn,CAs),EE ,1)+min(Plh);
    %% Transition energy
    chh = Ec + Ehh;
    clh = Ec + Elh;
    LD(LdIII==Ld)=Ld;
    wavelengthCHH(LdIII==Ld) = speed_of_light*h_bar*2*pi/chh*1e9;
    wavelengthCLH(LdIII==Ld) = speed_of_light*h_bar*2*pi/clh*1e9;
    plot(LD*1e9,wavelengthCHH,'m.-',  LD*1e9,wavelengthCLH,'b.-');
    pause(0.1);
end