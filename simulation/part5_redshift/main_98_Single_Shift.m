clear;clc;
meV = 1e-3 * 1.6e-19;
qcolor = 'brgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcy';

h = 3.5e-9; % 阱宽
H = 7e-9; % 垒宽

length = 0.1e-9; %向量间隔

bx = 1;   % T.~Mei,"k值与红移"
by = 0;
wx = 0.71;
wy = 0.61;

N = 1; % 量子阱数目
offset = 0.6;
z = (-N*(H+h)/2):length:(N*(H+h)/2); % 定义域
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
Ld_min = 0e-9;
Ld_length = 0.1e-9;
Ld_max = 1e-9;
Ld = Ld_min:Ld_length:Ld_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(1,2,1);
xlabel('z (nm)');
ylabel('indium composition');
ylim([0 0.3]);
box on;hold on;
subplot(1,2,2);
xlabel('diffusion length L_{dIII} (nm)');
ylabel('band energy levels (meV)');
box on;hold on;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    
    V_c(i,:) = Vc(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_hh(i,:) = Vhh(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_lh(i,:) = Vlh(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    
    subplot(1,2,1);plot(z*1e9, CIn(i,:)-CIn(i,round(N*(H+h)/2/length)),qcolor(i));
    
    [vectc,valuec] = clevel(length,N,H,h,V_c(i,:),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
    
    [vecthh,valuehh] = hhlevel(length,N,H,h,V_hh(i,:),CIn(i,:),CAs(i,:));
    Ehh(i,:) = sum(valuehh);
    
    [vectlh,valuelh] = lhlevel(length,N,H,h,V_lh(i,:),CIn(i,:),CAs(i,:));
    Elh(i,:) = sum(valuelh);
end
wavelengthhh(:,1) = 1./(Ec(:,1)+Ehh(:,1)) * 6.63e-34 * 3e8/1e-9;
wavelengthhh0 = 1./(Ec(1,1)+Ehh(1,1)) * 6.63e-34 * 3e8/1e-9;

wavelengthlh(:,1) = 1./(Ec(:,1)+Elh(:,1)) * 6.63e-34 * 3e8/1e-9;
wavelengthlh0 = 1./(Ec(1,1)+Elh(1,1)) * 6.63e-34 * 3e8/1e-9;

subplot(1,2,2);
% plot(Ld/1e-9,Ec(:,1)/meV,'.b-',Ld/1e-9,Ec(:,2)/meV,'.g-',Ld/1e-9,Ec(:,3)/meV,'.r-');
plot(Ld/1e-9,Ec(:,1)/meV,'.b-', Ld/1e-9,Ehh(:,1)/meV,'.g-', Ld/1e-9,Elh(:,1)/meV,'.r-');
legend('C1','HH1','LH1');
hold off;

subplot(1,2,1);
legend('0','0.1nm','0.2nm','0.3nm','0.4nm','0.5nm','0.6nm','0.7nm','0.8nm','0.9nm','1nm');
hold off;
%%%%%%%%%%%%%%%%%%%

figure(2);
subplot(1,2,1);
xlabel('z (nm)');
ylabel('indium composition');
% ylim([0 0.3]);
box on;hold on;
subplot(1,2,2);
xlabel('diffusion length L_{dIII} (nm)');
ylabel('band energy levels (meV)');
box on;hold on;

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    
    V_c(i,:) = Vc_2(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_hh(i,:) = Vhh_2(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_lh(i,:) = Vlh_2(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    
    subplot(1,2,1);plot(z*1e9, CIn(i,:), qcolor(i));
    
    [vectc,valuec] = clevel(length,N,H,h,V_c(i,:),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
    
    [vecthh,valuehh] = hhlevel(length,N,H,h,V_hh(i,:),CIn(i,:),CAs(i,:));
    Ehh(i,:) = sum(valuehh);
    
    [vectlh,valuelh] = lhlevel(length,N,H,h,V_lh(i,:),CIn(i,:),CAs(i,:));
    Elh(i,:) = sum(valuelh);
end

subplot(1,2,2);
plot(Ld/1e-9,Ec(:,1)/meV,'.b-', Ld/1e-9,Ehh(:,1)/meV,'.g-', Ld/1e-9,Elh(:,1)/meV,'.r-');
legend('C1','HH1','LH1');

subplot(1,2,1);
legend('0','0.1nm','0.2nm','0.3nm','0.4nm','0.5nm','0.6nm','0.7nm','0.8nm','0.9nm','1nm');

%%%%%%%%%%

figure(3);
xlabel('diffusion length L_{dIII} (nm)');
ylabel('wavelength (nm)');

wavelengthhh(:,1) = 1./(Ec(:,1)+Ehh(:,1)) * 6.63e-34 * 3e8/1e-9;
wavelengthlh(:,1) = 1./(Ec(:,1)+Elh(:,1)) * 6.63e-34 * 3e8/1e-9;
plot(Ld/1e-9,wavelengthhh,'.b-',   Ld/1e-9,wavelengthlh,'.g-');
legend('C1-HH1','C1-LH1');