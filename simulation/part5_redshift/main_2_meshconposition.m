clear;
h = 5e-9; % ���
H = 10e-9; % �ݿ�

length = 0.01e-9; %�������

bx = 0.82; % InGaAsP  ����> x,1-x,y,1-y
by = 0.38;
wx = 0.53;
wy = 1;

N = 1; % ��������Ŀ

z = (-N*(H+h)/2):length:(N*(H+h)/2); % ������

Ld_min = 0;
Ld_length = 0.1e-9;
Ld_max = 2.5e-9;
Ld = Ld_min:Ld_length:Ld_max;

%[zz,LdLd] = meshgrid(z,Ld);

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
end

plot(z,CIn(:,:));
