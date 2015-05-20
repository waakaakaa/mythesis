clear;
h = 3.5e-9; % 阱宽
H = 8e-9; % 垒宽

length = 0.01e-9; %向量间隔

bx = 1; % InGaAsP  ——> x,1-x,y,1-y
by = 0;
wx = 0.71;
wy = 0.61;

N = 1; % 量子阱数目

z = (-N*(H+h)/2):length:(N*(H+h)/2); % 定义域

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CIn0 = composition(length,h,H,bx,wx,N,1e-14);
CIn1 = composition(length,h,H,bx,wx,N,0.5e-9);
CIn2 = composition(length,h,H,bx,wx,N,2.5e-9);
CIn3 = composition(length,h,H,bx,wx,N,5e-9);

CGa0 = composition(length,h,H,1-bx,1-wx,N,1e-14);
CGa1 = composition(length,h,H,1-bx,1-wx,N,0.5e-9);
CGa2 = composition(length,h,H,1-bx,1-wx,N,2.5e-9);
CGa3 = composition(length,h,H,1-bx,1-wx,N,5e-9);

CAs0 = composition(length,h,H,by,wy,N,1e-14);
CAs1 = composition(length,h,H,by,wy,N,0.5e-9);
CAs2 = composition(length,h,H,by,wy,N,2.5e-9);
CAs3 = composition(length,h,H,by,wy,N,5e-9);

CP0 = composition(length,h,H,1-by,1-wy,N,1e-14);
CP1 = composition(length,h,H,1-by,1-wy,N,0.5e-9);
CP2 = composition(length,h,H,1-by,1-wy,N,2.5e-9);
CP3 = composition(length,h,H,1-by,1-wy,N,5e-9);


z = z*1e9;
subplot(2,2,1);plot(z,CIn0,z,CIn1,z,CIn2,z,CIn3);
title('In Composotions Before And After QWI');
legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
xlabel('z/nm');
ylabel('In compositions');

subplot(2,2,2);plot(z,CGa0,z,CGa1,z,CGa2,z,CGa3);
title('Ga Composotions Before And After QWI');
legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
xlabel('z/nm');
ylabel('Ga compositions');

subplot(2,2,3);plot(z,CAs0,z,CAs1,z,CAs2,z,CAs3);
title('As Composotions Before And After QWI');
legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
xlabel('z/nm');
ylabel('As compositions');


subplot(2,2,4);plot(z,CP0,z,CP1,z,CP2,z,CP3);
title('P Composotions Before And After QWI');
legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
xlabel('z/nm');
ylabel('P compositions');

%axis([-N*(H+h)/2 N*(H+h)/2 0 1]);