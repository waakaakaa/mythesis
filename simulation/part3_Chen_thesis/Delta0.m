function [ out ] = Delta0( w )

global ELEMENTARY_CHARGE;
out = 0.34 - 0.065*w;
out = out*ELEMENTARY_CHARGE;

end

