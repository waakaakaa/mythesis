clear;
h = 3.5e-9; % 阱宽
H = 7e-9; % 垒宽

length = 0.1e-9; %向量间隔

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
k = 0.25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CIn0 = composition(length,h,H,bx,wx,N,1e-15);
%CIn1 = composition(length,h,H,bx,wx,N,0.5e-9);
CIn2 = composition(length,h,H,bx,wx,N,1e-9);
%CIn3 = composition(length,h,H,bx,wx,N,2e-9);

%CAs0 = composition(length,h,H,by,wy,N,k*1e-14);
%CAs1 = composition(length,h,H,by,wy,N,k*0.5e-9);
CAs2 = composition(length,h,H,by,wy,N,k*1e-9);
%CAs3 = composition(length,h,H,by,wy,N,k*2e-9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(z,CIn0,z,CAs0,z,CIn1,z,CAs1,z,CIn2,z,CAs2,z,CIn3,z,CAs3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%E0 = offset*potential(CIn0,CAs0);
%E1 = offset*potential(CIn1,CAs1);
E2 = offset*potential(CIn2,CAs2);
%E3 = offset*potential(CIn3,CAs3);

%E0 = E0 - E0(1,N*(H+h)/2/length); %阱数目为偶数的时候需要修正！！！
%E1 = E1 - E1(1,N*(H+h)/2/length);
%E2 = E2 - E2(1,N*(H+h)/2/length);
%E3 = E3 - E3(1,N*(H+h)/2/length);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(z,E0,z,E1,z,E2,z,E3);
%plot(z,E2);
%hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%str0 = strain(CIn0,CAs0);
%str1 = strain(CIn1,CAs1);
%str2 = strain(CIn2,CAs2);
%str3 = strain(CIn3,CAs3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(z,str2);
%plot(z,str0,z,str1,z,str2,z,str3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vc0 = Vc(CIn0,CAs0,offset);
%Vc1 = Vc(CIn1,CAs1,offset);
Vc2 = Vc(CIn2,CAs2,offset);
%Vc3 = Vc(CIn3,CAs3,offset);

%Vhh0 = Vhh(CIn0,CAs0,offset);
%Vhh1 = Vhh(CIn1,CAs1,offset);
%Vhh2 = Vhh(CIn2,CAs2,offset);
%Vhh3 = Vhh(CIn3,CAs3,offset);

%Vlh0 = Vlh(CIn0,CAs0,offset);
%Vlh1 = Vlh(CIn1,CAs1,offset);
%Vlh2 = Vlh(CIn2,CAs2,offset);
%Vlh3 = Vlh(CIn3,CAs3,offset);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(z,Vc2,z,-Vhh2,z,-Vlh2);
%plot(z,Vc0,z,Vc1,z,Vc2,z,Vc3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[vectc,valuec] = clevel(length,N,H,h,Vc0,CIn0,CAs0);
%E = sum(valuec);
%plot(z(1,1:4),E(1,1:4),'r.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
In = CIn2;
As = CAs2;
V = Vc2;

h_bar = 6.626196e-34 / 2 / 3.1415926;
m = 0.0632*(1-In).*As + 0.0213*In.*As + 0.17*(1-In).*(1-As) + 0.077*In.*(1-As);%InGaAsP电子的有效质量
m = 9.109534e-31 * m;

x = N*(H+h)/length+1

A = m;
B = m;
A(1,x) = h_bar^2 / (2*length^2*m(1,x));
B(1,1) = h_bar^2 / (2*length^2*m(1,1));
for i = 1 : (x-1)
    A(1,i) = h_bar^2 / (length^2 * (m(1,i)+m(1,i+1)));
end
for i = 2 : x
    B(1,i) = h_bar^2 / (length^2 * (m(1,i)+m(1,i-1)));
end
C = A+B+V;

W = zeros(x,x);
W(1,1) = C(1,1);
W(1,2) = -A(1,1);
W(x,x) = C(1,x);
W(x,x-1) = -B(1,x);

for i = 2:(x-1)
    W(i,i) = C(1,i);
    W(i,i+1) = -A(1,i);
    W(i,i-1) = -B(1,i);
end

[vect,value] = eig(W);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E(1,:) = sum(value);
plot(z,E(1,:),'r.');