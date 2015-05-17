clear;
h = 5e-9; % 阱宽
H = 10e-9; % 垒宽

length_z = 0.05e-9;

bx = 0.82;
by = 0.38;
wx = 0.53;
wy = 1;

N = 1; % 量子阱数目

offset = 0.6;

z = (-N*(H+h)/2):length_z:(N*(H+h)/2); % 定义域

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
Ld_min = 0;
Ld_length = 0.02e-9;
Ld_max = 1e-9;
Ld = Ld_min:Ld_length:Ld_max;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = offset*potential(CIn(i,:),CAs(i,:));
    E(i,:) = E(i,:) - E(i,round(N*(H+h)/2/length_z));
end
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectc,valuec] = clevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = (1-offset)*potential(CIn(i,:),CAs(i,:));
    E(i,:) = E(i,:) - E(i,round(N*(H+h)/2/length_z));
end
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vectv,valuev] = hhlevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    Ev(i,:) = sum(valuev);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EHH = Ec(:,1)-Ev(:,1);

plot(Ld,Ec(:,1),Ld,Ev(:,1));
%plot(Ld,Ec(:,1),Ld,Ev(:,1),Ld,EHH);

%wavelength = 3e8 * 6.63e-34 ./ EHH;
%for i = 2:((Ld_max-Ld_min)/Ld_length+1)
%    wavelength(i,1) = wavelength(1,1)-wavelength(i,1);
%end
%wavelength(1,1) = 0;
%plot(Ld,wavelength);