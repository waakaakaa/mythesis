function V = vhh(x,y)

global offset;
global z;

E = Eg(x,y);
mysize = max(size(z));
center = round(mysize/2);
V = (1-offset) * ( E - E(1,center) - s_1(x,y) ) + s_2(x,y);