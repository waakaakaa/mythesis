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
h_bar = 6.63e-34 /2 /pi;
speed_of_light = 3e8;
offset = 0.6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
length_z = 0.05e-9;
h = 5.5e-9;
H = 10e-9;
wx = 0.8;
wy = 0.8;
[bx,by] = q2xy(1.25);
z = -(H+h)/2 : length_z : (H+h)/2;
LdIII = 0;
k = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[chh,clh] = qwi_fdm(LdIII,k);
wavelengthCHH = speed_of_light*h_bar*2*pi/chh*1e9;
wavelengthCLH = speed_of_light*h_bar*2*pi/clh*1e9;