function C = compositionIII(Ld)

global h;
global z;
global bx;
global wx;

temp = 2 - erf((h/2-z)/Ld/2) - erf((h/2+z)/Ld/2);
C = wx + (bx-wx)*temp/2;