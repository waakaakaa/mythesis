function [ C ] = composotionIII( Ld,x,L,well,barrier )

temp = 2 - erf((L/2-x)/Ld/2) - erf((L/2+x)/Ld/2);
C = well + (barrier-well)*temp/2;

end