function V = Vlh(length,h,H,x,y,N,offset)

Eg = 1.35 - 1.17*y + 0.668*(1-x) - 0.069*y.*(1-x) + 0.18*y.^2 + 0.03*(1-x).*y.^2 + 0.758*(1-x).^2 - 0.322*y.*(1-x).^2; %InGaAsP
Eg = 1.6e-19 * Eg;

delta = 0.34*(1-x).*y + 0.43*x.*y + 0.10*(1-x).*(1-y) + 0.10*x.*(1-y);
delta = delta * 1.6e-19;

V = (1-offset) * ( Eg - Eg(1,round(N*(H+h)/2/length)) - s_1(x,y) ) + 0.5 * ( s_2(x,y) + delta ) - 0.5*sqrt( delta.^2 - 2*s_2(x,y).*delta + 9*s_2(x,y).^2 );