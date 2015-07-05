function V = vc(x,y)

global offset;
global z;

E = Eg(x,y);
mysize = max(size(z));
center = round(mysize/2);
V = offset * ( E - E(1,center) ) - (1-offset) * s_1(x,y) + E(1,center);