clear;
h = 5e-9; % ���
H = 10e-9; % �ݿ�

length_z = 0.05e-9;

bx = 0.82;
by = 0.38;
wx = 0.53;
wy = 1;

N = 1; % ��������Ŀ

offset = 0.6;

z = (-N*(H+h)/2):length_z:(N*(H+h)/2); % ������

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.25;
Ld_min = 0;
Ld_length = 0.02e-9;
Ld_max = 1e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    E(i,:) = E(i,:) - E(i,round(N*(H+h)/2/length_z));
end

%plot(z,E(:,:));
%title('Potential Curve ');
%legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
%xlabel('z/nm');
%ylabel('Potential');
%axis([-N*(H+h)/2 N*(H+h)/2 0 1]);

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectc,valuec] = clevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
end
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectv,valuev] = hhlevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ev(i,:) = sum(valuev);
end
EHH = offset*Ec(:,1)-(offset-1)*Ev(:,1);
wavelength = 3e8 * 6.63e-34 ./ EHH;
plot(Ld,wavelength);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.5;
Ld_min = 0;
Ld_length = 0.02e-9;
Ld_max = 1e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    E(i,:) = E(i,:) - E(i,round(N*(H+h)/2/length_z));
end

%plot(z,E(:,:));
%title('Potential Curve ');
%legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
%xlabel('z/nm');
%ylabel('Potential');
%axis([-N*(H+h)/2 N*(H+h)/2 0 1]);

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectc,valuec] = clevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
end
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectv,valuev] = hhlevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ev(i,:) = sum(valuev);
end
EHH = offset*Ec(:,1)-(offset-1)*Ev(:,1);
wavelength = 3e8 * 6.63e-34 ./ EHH;
plot(Ld,wavelength);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
Ld_min = 0;
Ld_length = 0.02e-9;
Ld_max = 1e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    E(i,:) = E(i,:) - E(i,round(N*(H+h)/2/length_z));
end

%plot(z,E(:,:));
%title('Potential Curve ');
%legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
%xlabel('z/nm');
%ylabel('Potential');
%axis([-N*(H+h)/2 N*(H+h)/2 0 1]);

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectc,valuec] = clevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
end
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectv,valuev] = hhlevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ev(i,:) = sum(valuev);
end
EHH = offset*Ec(:,1)-(offset-1)*Ev(:,1);
wavelength = 3e8 * 6.63e-34 ./ EHH;
plot(Ld,wavelength);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1.5;
Ld_min = 0;
Ld_length = 0.02e-9;
Ld_max = 1e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    E(i,:) = E(i,:) - E(i,round(N*(H+h)/2/length_z));
end

%plot(z,E(:,:));
%title('Potential Curve ');
%legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
%xlabel('z/nm');
%ylabel('Potential');
%axis([-N*(H+h)/2 N*(H+h)/2 0 1]);

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectc,valuec] = clevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
end
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectv,valuev] = hhlevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ev(i,:) = sum(valuev);
end
EHH = offset*Ec(:,1)-(offset-1)*Ev(:,1)
wavelength = 3e8 * 6.63e-34 ./ EHH;
plot(Ld,wavelength);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 2;
Ld_min = 0;
Ld_length = 0.02e-9;
Ld_max = 1e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    E(i,:) = E(i,:) - E(i,round(N*(H+h)/2/length_z));
end

%plot(z,E(:,:));
%title('Potential Curve ');
%legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
%xlabel('z/nm');
%ylabel('Potential');
%axis([-N*(H+h)/2 N*(H+h)/2 0 1]);

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectc,valuec] = clevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
end
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectv,valuev] = hhlevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ev(i,:) = sum(valuev);
end
EHH = offset*Ec(:,1)-(offset-1)*Ev(:,1);
wavelength = 3e8 * 6.63e-34 ./ EHH;
plot(Ld,wavelength);