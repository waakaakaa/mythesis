clear;
clc;

h = 5e-9; % 阱宽
H = 10e-9; % 垒宽

length = 0.1e-9;

bx = 0.82;
by = 0.38;
wx = 0.53;
wy = 1;

N = 1; % 量子阱数目

offset = 0.6;

z = (-N*(H+h)/2):length:(N*(H+h)/2); % 定义域

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CIn0 = composition(length,h,H,bx,wx,N,1e-15);
%CIn1 = composition(length,h,H,bx,wx,N,0.5e-9);
%CIn2 = composition(length,h,H,bx,wx,N,2.5e-9);
%CIn3 = composition(length,h,H,bx,wx,N,5e-9);

CAs0 = composition(length,h,H,by,wy,N,1e-15);
%CAs1 = composition(length,h,H,by,wy,N,0.5e-9);
%CAs2 = composition(length,h,H,by,wy,N,2.5e-9);
%CAs3 = composition(length,h,H,by,wy,N,5e-9);

E0 = offset*potential(CIn0,CAs0);
%E1 = offset*potential(CIn1,CAs1);
%E2 = offset*potential(CIn2,CAs2);
%E3 = offset*potential(CIn3,CAs3);

E0 = E0 - E0(1,round(N*(H+h)/2/length)); %阱数目为偶数的时候需要修正！！！
%E1 = E1 - E1(1,N*(H+h)/2/length);
%E2 = E2 - E2(1,N*(H+h)/2/length);
%E3 = E3 - E3(1,N*(H+h)/2/length);

%plot(z,E0,z,E1,z,E2,z,E3);
%title('Potential Curve ');
%legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
%xlabel('z/nm');
%ylabel('Potential');
%axis([-N*(H+h)/2 N*(H+h)/2 0 1]);
x = N*(H+h)/length+1;
In = CIn0;
As = CAs0;
h_bar = 6.626196e-34 / 2 / 3.1415926;
m = 0.0632*(1-In).*As + 0.0213*In.*As + 0.17*(1-In).*(1-As) + 0.077*In.*(1-As);%InGaAsP电子的有效质量
m = 9.109534e-31 * m;
E_infinite = 3.14^2 * h_bar^2 / 2 / h^2 / m(1,round((x+1)/2))   %无限深势阱模型的能级

[F,E] = clevel(length,N,H,h,E0,CIn0,CAs0);

%F = F.*F;
E = sum(E);

%subplot(1,2,1);
%plot(z,F(:,1),z,F(:,2),z,F(:,3),z,F(:,4));
%subplot(1,2,2);
%plot(z(1,1:4),E(1,1:4),'r.');
E(1,1)
E(1,2)
E(1,3)