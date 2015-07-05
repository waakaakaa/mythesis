function [EEE,gamaTE,gamaTM] = gain(delta_n, step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The final version of gain imformation for bulk materials at T = 300k.
%
% @author Frank
% @date 2010-08-10
%%%%%%%%%%%%%%%%%%%%  华丽丽的常量  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global wx;
global wy;
global m0;
global KB;
global h;
global h_bar;
global mc;
global mhh;
global mlh;
global mrhh;
global mrlh;
global n;
global c;
global epsilon;
global d;
e = 1.6e-19;
Eg = Eg(wx,wy);
T = 300;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 态密度公式
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = 0:Eg/1000:Eg*10;

t = 1;
rouc = tmp-tmp;
E1 = (h_bar^2*pi^2) / (2*mc*d^2);
En = E1;
while En<=Eg*10
    rouc = rouc + mc / (pi*h_bar^2*d) .* ( tmp>En );
    t = t + 1;
    En = E1 * t^2;
end
t = 1;
rouhh = tmp-tmp;
E1 = (h_bar^2*pi^2) / (2*mhh*d^2);
En = E1;
while En<=Eg*10
    rouhh = rouhh + mhh / (pi*h_bar^2*d) .* ( tmp>En );
    t = t + 1;
    En = E1 * t^2;
end
t = 1;
roulh = tmp-tmp;
E1 = (h_bar^2*pi^2) / (2*mlh*d^2);
En = E1;
while En<=Eg*10
    roulh = roulh + mlh / (pi*h_bar^2*d) .* ( tmp>En );
    t = t + 1;
    En = E1 * t^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 准费米能级
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EfcminusEc = -Eg/2;
while 1
    nc = trapz(tmp, rouc./(1+exp((tmp-EfcminusEc)/KB/T)));
    if(nc>delta_n)
        break;
    end
    EfcminusEc = EfcminusEc + Eg/1000;
end
EhhminusEfv = Eg/2;
while 1
    nhh = trapz(tmp, rouhh./(1+exp((tmp-EhhminusEfv)/KB/T)));
    if(nhh<delta_n)
        break;
    end
    EhhminusEfv = EhhminusEfv - Eg/1000;
end
ElhminusEfv = Eg/2;
while 1
    nlh = trapz(tmp, roulh./(1+exp((tmp-ElhminusEfv)/KB/T)));
    if(nlh<delta_n)
        break;
    end
    ElhminusEfv = ElhminusEfv - Eg/1000;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 自变量 E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EndPoint = 1.1*max([(EfcminusEc+EhhminusEfv+Eg) (EfcminusEc+ElhminusEfv+Eg)]);
E = Eg:step:EndPoint;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 费米反转因子
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E2minusEc = (mrhh/mc) * (E-Eg);
fc = 1 ./ ( 1 + exp( (E2minusEc - EfcminusEc) /KB/T ) );
EhhminusE1 = (mrhh/mhh) * (E-Eg);
fv = 1 ./ ( 1 + exp( (EhhminusEfv - EhhminusE1) /KB/T ) );
fghh = fc-fv;

E2minusEc = (mrlh/mc) * (E-Eg);
fc = 1 ./ ( 1 + exp( (E2minusEc - EfcminusEc) /KB/T ) );
ElhminusE1 = (mrlh/mlh) * (E-Eg);
fv = 1 ./ ( 1 + exp( (ElhminusEfv - ElhminusE1) /KB/T ) );
fglh = fc-fv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 折合态密度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1;
rouhh = E-E;
E1 = (h_bar^2*pi^2) / (2*mrhh*d^2);
En = E1;
while En<=Eg*10
    rouhh = rouhh + mrhh / (pi*h_bar^2*d) .* ( (E-Eg)>En );
    t = t + 1;
    En = E1 * t^2;
end
t = 1;
roulh = E-E;
E1 = (h_bar^2*pi^2) / (2*mrlh*d^2);
En = E1;
while En<=Eg*10
    roulh = roulh + mrlh / (pi*h_bar^2*d) .* ( (E-Eg)>En );
    t = t + 1;
    En = E1 * t^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 光学矩阵元
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OME = (m0/mc-1)*m0*Eg/2;
% OME = (19.7+5.6*y)*m0/2*e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 静增益
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamahh = (pi*e^2*h_bar)./(m0^2*epsilon*n*E*c) .* OME .* rouhh .*fghh;
gamalh = (pi*e^2*h_bar)./(m0^2*epsilon*n*E*c) .* OME .* roulh .*fglh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 洛伦兹增宽
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
torr_in = 0.1e-12;
EE = -0.1*e : step : 0.1*e;
L = 1/pi * (h_bar/torr_in) ./ ( (h_bar/torr_in)^2 + EE.^2);

EEE = (Eg - 0.1*e) : step : (EndPoint + 0.1*e);
gamahhb = conv(gamahh,L) * step;
gamalhb = conv(gamalh,L) * step;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 量子阱增益
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamaTE = gamahhb/2 + gamalhb/6;
gamaTM = gamalhb*2/3;