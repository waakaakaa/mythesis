clc;
clear;
x = 0.53;
y = 1;
d = 6e-9;
%%%%%%%%%%%%%%%%%%%%  华丽丽的常量  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 300;
e = 1.6e-19;
Eg = e * (1.35 - 1.17*y + 0.668*(1-x) - 0.069*y.*(1-x) + 0.18*y.^2 + 0.03*(1-x).*y.^2 + 0.758*(1-x).^2 - 0.322*y.*(1-x).^2); %InGaAsP
m0 = 9.11e-31;
KB = 1.38e-23;
h = 6.626e-34;
h_bar = h /2 /pi;
mc = m0 * (0.0632*(1-x).*y + 0.0213*x.*y + 0.17*(1-x).*(1-y) + 0.077*x.*(1-y));%InGaAsP电子的有效质量
mv = m0 * (0.5*(1-x).*y + 0.41*x.*y + 0.54*(1-x).*(1-y) + 0.12*x.*(1-y));% HH !
mr = 1/(1/mc+1/mv);
n = 3.41;
% tor_r = 2e-9;
c = 2.998e8;
tor_in = 0.1e-12;
epsilon0 = 8.854e-12;
%%%%%%%%%%%%%%%%%%%  注入  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_n = 2.0e18 * 1e6;
%%%%%%%%%%%%%%%%%%%  T=0K时的准费米能级的Eg  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EfEg0 = ( Eg + (3*pi^2)^(2/3)*h_bar^2/2/mr*(delta_n)^(2/3) );
%%%%%%%%%%%%%%%%%%%  T=300K的准费米能级  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 0:Eg/1000:Eg*10;
y = -Eg*10:Eg/1000:0;
rou_c = (2*mc)^1.5 / (2*pi^2*h_bar^3) * (x).^0.5;
rou_v = (2*mv)^1.5 / (2*pi^2*h_bar^3) * (-y).^0.5;

tmp = 1;
rou_c_QW = rou_c-rou_c;
E1 = (h_bar^2*pi^2) / (2*mc*d^2);
En = E1;
while En<=Eg*10
    rou_c_QW = rou_c_QW + mc / (pi*h_bar^2*d) .* ( x>En );
    tmp = tmp + 1;
    En = E1 * tmp^2;
end
tmp = 1;
rou_v_QW = rou_v-rou_v;
E1 = -(h_bar^2*pi^2) / (2*mv*d^2);
En = E1;
while En>=-Eg*10
    rou_v_QW = rou_v_QW + mv / (pi*h_bar^2*d) .* ( y<En );
    tmp = tmp + 1;
    En = E1 * tmp^2;
end
% plot(x,rou_c,  x,rou_c_QW,  y,rou_v,  y,rou_v_QW);

Efc = 0;
while 1
    nc = trapz(x, rou_c_QW./(1+exp((x-Efc)/(KB*T))));
    if(nc>delta_n)
        break;
    end
    Efc = Efc + EfEg0/10000;
end
Efv = 0;
while 1
    nv = trapz(y, rou_v_QW./(1+exp((Efv-y)/(KB*T))));
    if(nv<delta_n)
        break;
    end
    Efv = Efv + EfEg0/10000;
end
EfEg = Efc - Efv + Eg;
%%%%%%%%%%%%%%%%%%%  态密度啦啦啦啦  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = Eg/h:(EfEg-Eg)/h/1000:EfEg/h;% 频率 华丽丽的自变量
rou = (2*mr)^1.5 / (pi*h_bar^2) * (h*nu-Eg).^0.5;
tmp = 1;
rouQW = rou - rou;
E1 = (h_bar^2*pi^2) / (2*mc*d^2);
E2 = (h_bar^2*pi^2) / (2*mv*d^2);
En = E1 + E2;
while En<=EfEg
    rouQW = rouQW + (2*mr) / (h_bar*d) .*((h*nu)>(En+Eg));
    tmp = tmp + 1;
    En = ( E1 + E2 ) * tmp^2;
end
% plot(nu,rou,  nu,rouQW);
%%%%%%%%%%%%%%%%%%%  前置项啦啦啦啦  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
front = (h*e^2) ./ (2*n*epsilon0*m0^2*c*h*nu);
%%%%%%%%%%%%%%%%%%%  光学矩阵元！！  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = (m0^2*Eg) / (12*mc);
%%%%%%%%%%%%%%%%%%%  计算增益啦啦啦  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EE2 = (mr/mc) * (h*nu-Eg) + Eg;
EE1 = EE2 - h*nu;
fc = 1 ./ ( 1 + exp( (EE2 - Efc -Eg) /KB/T ) );
fv = 1 ./ ( 1 + exp( (EE1 - Efv) /KB/T ) );
fg = fc-fv;
% plot(nu,fc,  nu,fv,  nu,fg);
% plot(nu,fg);

% gama = (c/n)^2 / (8*pi*tor_r) .* rouQW .* fg ./nu ./nu;
gama = front .* M.^2 .* rouQW .* fg;
%%%%%%%%%%%%%%%%%%%  线性增宽啦啦啦  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % nunu = -Eg/h/10:Eg/h/10000:Eg/h/10;
% % L = (1/pi) * (h_bar/tor_in) ./ ((h_bar/tor_in)^2+(h*nunu).^2);
% % plot(nunu*h/e,L);
% nunu = 0:EfEg/h/10000:EfEg/h*1.2;
% gamaFinal = nunu-nunu;
% i = 1;
% for x = 0:EfEg/h/10000:EfEg/h*1.2
%     gamaFinal(1,i) = trapz(nu, gama * (1/pi) * (h_bar/tor_in) ./ ( (h_bar/tor_in)^2+(h*x-h*nu).^2) );
%     i = i + 1;
% end
% % plot(nunu*h/e,gamaFinal/1e2);
% plot(c./nunu/1e-9,gamaFinal/1e2);
%%%%%%%%%%%%%%%%%%%%  绘图是一门艺术  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(nu*h/e,gama/1e2);
% plot(c./nu/1e-9,gama/1e2);
% plot(nu*h/e,gama/1e2,nu*h/e,gamaFinal/1e2);
xlabel('h\nu (eV)');
% xlabel('\lambda (nm)');
ylabel('Gain coefficient \gamma_0(\nu) (cm^{-1})');