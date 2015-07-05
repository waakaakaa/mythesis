function C = compositionV(Ld)

global h;
global z;
global by;
global wy;

temp = 2 - erf((h/2-z)/Ld/2) - erf((h/2+z)/Ld/2);
C = wy + (by-wy)*temp/2;