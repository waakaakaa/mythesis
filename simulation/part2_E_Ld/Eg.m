function [ out ] = Eg( x )

global ELEMENTARY_CHARGE;

out = 1.424 + 1.594*x + x.*(1-x).*(0.127-1.31*x);
out = out*ELEMENTARY_CHARGE;

end

