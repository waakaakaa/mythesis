clear;
h = 6e-9; % 阱宽
H = 34e-9; % 垒宽

length = 0.01e-9; %向量间隔

%wx = 0.53; % 中文小论文
%wy = 1;
%bx = 0.82;
%by = 0.38;

%bx = 1;   %T.~Mei,"k值与红移"
%by = 0;
%wx = 0.71;
%wy = 0.29;

wx = 0.53; % 教科书
wy = 1;
bx = 1;
by = 0;

N = 1; % 量子阱数目
offset = 0.6;
z = (-N*(H+h)/2):length:(N*(H+h)/2); % 定义域
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CIn0 = composition(length,h,H,bx,wx,N,1e-15);
CIn1 = composition(length,h,H,bx,wx,N,0.5e-9);
CIn2 = composition(length,h,H,bx,wx,N,1e-9);
CIn3 = composition(length,h,H,bx,wx,N,2e-9);

CAs0 = composition(length,h,H,by,wy,N,k*1e-14);
CAs1 = composition(length,h,H,by,wy,N,k*0.5e-9);
CAs2 = composition(length,h,H,by,wy,N,k*1e-9);
CAs3 = composition(length,h,H,by,wy,N,k*2e-9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(z,CIn0,z,CAs0,z,CIn1,z,CAs1,z,CIn2,z,CAs2,z,CIn3,z,CAs3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E0 = potential(CIn0,CAs0);
E1 = potential(CIn1,CAs1);
E2 = potential(CIn2,CAs2);
E3 = potential(CIn3,CAs3);

%E0 = E0 - E0(1,N*(H+h)/2/length); %阱数目为偶数的时候需要修正！！！
%E1 = E1 - E1(1,N*(H+h)/2/length);
E2 = (1-offset)*(E2 - E2(1,round(N*(H+h)/2/length)));
%E3 = E3 - E3(1,N*(H+h)/2/length);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(z,E0,z,E1,z,E2,z,E3);
plot(z,-E2/1.6e-19,'r.');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str0 = strain(CIn0,CAs0);
str1 = strain(CIn1,CAs1);
str2 = strain(CIn2,CAs2);
str3 = strain(CIn3,CAs3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(z,str2);
%plot(z,str0,z,str1,z,str2,z,str3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vlh0 = Vlh(CIn0,CAs0,offset);
%Vlh1 = Vlh(CIn1,CAs1,offset);
Vlh2 = Vlh(length,h,H,CIn2,CAs2,N,offset);
%Vlh3 = Vlh(CIn3,CAs3,offset);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(z,-Vlh2/1.6e-19);
%plot(z,Vc0,z,Vc1,z,Vc2,z,Vc3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%