function V = vlh(x,y)

global offset;
global z;

E = Eg(x,y);
mysize = max(size(z));
center = round(mysize/2);

V = (1-offset) * ( E - E(1,center) - s_1(x,y) )...
    + 0.5 * ( s_2(x,y) + delta(x,y) )...
    - 0.5 * sqrt( delta(x,y).^2 - 2*s_2(x,y).*delta(x,y) + 9*s_2(x,y).^2 );