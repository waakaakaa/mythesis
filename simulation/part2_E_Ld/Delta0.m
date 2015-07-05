function [ out ] = Delta0( x )

global ELEMENTARY_CHARGE;
out = 0.34 - 0.065*x;
out = out * ELEMENTARY_CHARGE;

end

