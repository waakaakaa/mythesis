function [ out ] = b( w )

global ELEMENTARY_CHARGE;
out = -1.7 + 0.2*w;
out = out*ELEMENTARY_CHARGE;

end

