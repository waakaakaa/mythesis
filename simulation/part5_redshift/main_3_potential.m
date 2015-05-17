h = 5e-9; % ���
H = 10e-9; % �ݿ�

length = 0.01e-9;

bx = 0.82;
by = 0.38;
wx = 0.53;
wy = 1;

N = 1; % ��������Ŀ

offset = 0.6;

z = (-N*(H+h)/2):length:(N*(H+h)/2); % ������

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CIn0 = composition(length,h,H,bx,wx,N,1e-15);
CIn1 = composition(length,h,H,bx,wx,N,0.5e-9);
CIn2 = composition(length,h,H,bx,wx,N,2.5e-9);
CIn3 = composition(length,h,H,bx,wx,N,5e-9);

CAs0 = composition(length,h,H,by,wy,N,1e-15);
CAs1 = composition(length,h,H,by,wy,N,0.5e-9);
CAs2 = composition(length,h,H,by,wy,N,2.5e-9);
CAs3 = composition(length,h,H,by,wy,N,5e-9);

E0 = offset*potential(CIn0,CAs0);
E1 = offset*potential(CIn1,CAs1);
E2 = offset*potential(CIn2,CAs2);
E3 = offset*potential(CIn3,CAs3);

E0 = E0 - E0(1,round(N*(H+h)/2/length)); %����ĿΪż����ʱ����Ҫ����������
E1 = E1 - E1(1,round(N*(H+h)/2/length));
E2 = E2 - E2(1,round(N*(H+h)/2/length));
E3 = E3 - E3(1,round(N*(H+h)/2/length));

plot(z,E0,z,E1,z,E2,z,E3);
title('Potential Curve ');
legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
xlabel('z/nm');
ylabel('Potential/eV');
%axis([-N*(H+h)/2 N*(H+h)/2 0 1]);