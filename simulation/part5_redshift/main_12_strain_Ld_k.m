clear;
h = 5e-9; % ���
H = 10e-9; % �ݿ�

length_z = 0.01e-9; %�������

wx = 0.53; % ����С����
wy = 1;
bx = 0.82;
by = 0.38;

%bx = 1;   %T.~Mei,"kֵ�����"
%by = 0;
%wx = 0.71;
%wy = 0.29;

%wx = 0.53; % �̿���
%wy = 1;
%bx = 1;
%by = 0;

N = 1; % ��������Ŀ
offset = 0.6;
z = (-N*(H+h)/2):length_z:(N*(H+h)/2); % ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0;
Ld_min = 0;
Ld_length = 0.1e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    str(i,:) = strain(CIn(i,:),CAs(i,:));
end
plot(Ld,str(:,round(N*(H+h)/length_z/2+1)),'.');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.25;
Ld_min = 0;
Ld_length = 0.1e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    str(i,:) = strain(CIn(i,:),CAs(i,:));
end
plot(Ld,str(:,round(N*(H+h)/length_z/2+1)),'.');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.5;
Ld_min = 0;
Ld_length = 0.1e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    str(i,:) = strain(CIn(i,:),CAs(i,:));
end
plot(Ld,str(:,round(N*(H+h)/length_z/2+1)),'.');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
Ld_min = 0;
Ld_length = 0.1e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    str(i,:) = strain(CIn(i,:),CAs(i,:));
end
plot(Ld,str(:,round(N*(H+h)/length_z/2+1)),'.');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 2;
Ld_min = 0;
Ld_length = 0.1e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    str(i,:) = strain(CIn(i,:),CAs(i,:));
end
plot(Ld,str(:,round(N*(H+h)/length_z/2+1)),'.');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 5;
Ld_min = 0;
Ld_length = 0.1e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    str(i,:) = strain(CIn(i,:),CAs(i,:));
end
plot(Ld,str(:,round(N*(H+h)/length_z/2+1)),'.');
title('strain');
xlabel('Ld_III/m');
ylabel('strain/%');