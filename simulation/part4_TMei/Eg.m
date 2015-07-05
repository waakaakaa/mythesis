function [ out ] = Eg( w )

global ELEMENTARY_CHARGE;
out = 1.423 + 1.594*w + w.*(1-w).*(0.127-1.31*w);
out = out*ELEMENTARY_CHARGE;

end

