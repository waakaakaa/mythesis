clear;
h = 6e-9; % 阱宽
H = 12e-9; % 垒宽

length = 0.01e-9; %向量间隔

%bx = 0.82;
%by = 0.38;
%wx = 0.53;
%wy = 1;

bx = 1; % InGaAsP  ——> x,1-x,y,1-y
by = 0;
wx = 0.53;
wy = 1;

N = 1; % 量子阱数目
offset = 0.6;
z = (-N*(H+h)/2):length:(N*(H+h)/2); % 定义域
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CIn0 = composition(length,h,H,bx,wx,N,1e-15);
CIn1 = composition(length,h,H,bx,wx,N,0.5e-9);
CIn2 = composition(length,h,H,bx,wx,N,1e-9);
CIn3 = composition(length,h,H,bx,wx,N,2e-9);

CGa0 = composition(length,h,H,1-bx,1-wx,N,1e-14);
CGa1 = composition(length,h,H,1-bx,1-wx,N,0.5e-9);
CGa2 = composition(length,h,H,1-bx,1-wx,N,1e-9);
CGa3 = composition(length,h,H,1-bx,1-wx,N,2e-9);

CAs0 = composition(length,h,H,by,wy,N,k*1e-14);
CAs1 = composition(length,h,H,by,wy,N,k*0.5e-9);
CAs2 = composition(length,h,H,by,wy,N,k*1e-9);
CAs3 = composition(length,h,H,by,wy,N,k*2e-9);

CP0 = composition(length,h,H,1-by,1-wy,N,k*1e-14);
CP1 = composition(length,h,H,1-by,1-wy,N,k*0.5e-9);
CP2 = composition(length,h,H,1-by,1-wy,N,k*1e-9);
CP3 = composition(length,h,H,1-by,1-wy,N,k*2e-9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(z,CIn2,z,CGa2,z,CAs2,z,CP2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E0 = potential(CIn0,CAs0);
E1 = potential(CIn1,CAs1);
E2 = potential(CIn2,CAs2);
E3 = potential(CIn3,CAs3);

E0 = E0 - E0(1,N*(H+h)/2/length); %阱数目为偶数的时候需要修正！！！
E1 = E1 - E1(1,N*(H+h)/2/length);
E2 = E2 - E2(1,N*(H+h)/2/length);
E3 = E3 - E3(1,N*(H+h)/2/length);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(z,E0,z,E1,z,E2,z,E3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str0 = strain(CIn0,CAs0);
str1 = strain(CIn1,CAs1);
str2 = strain(CIn2,CAs2);
str3 = strain(CIn3,CAs3);
%plot(z,str2);
%plot(z,str0,z,str1,z,str2,z,str3);