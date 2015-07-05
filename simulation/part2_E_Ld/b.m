function [ out ] = b( x )

global ELEMENTARY_CHARGE;

out = -1.7 + 0.2*x;
out = out*ELEMENTARY_CHARGE;

end

