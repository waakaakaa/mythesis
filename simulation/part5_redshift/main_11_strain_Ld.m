clear;
h = 6e-9; % 阱宽
H = 12e-9; % 垒宽

length_z = 0.01e-9; %向量间隔

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
z = (-N*(H+h)/2):length_z:(N*(H+h)/2); % 定义域
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.25;
Ld_min = 0;
Ld_length = 0.1e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    str(i,:) = strain(CIn(i,:),CAs(i,:));
end
plot(Ld,str(:,N*(H+h)/length_z/2+1)*100,'.');
title('strain');
legend('k=0.25');
xlabel('Ld_III/m');
ylabel('strain/%');