function V = Vhh(length,h,H,x,y,N,offset)

Eg = 1.35 - 1.17*y + 0.668*(1-x) - 0.069*y.*(1-x) + 0.18*y.^2 + 0.03*(1-x).*y.^2 + 0.758*(1-x).^2 - 0.322*y.*(1-x).^2; %InGaAsP
Eg = 1.6e-19 * Eg;

V = (1-offset) * ( Eg - Eg(1,round(N*(H+h)/2/length)) - s_1(x,y) ) + s_2(x,y);