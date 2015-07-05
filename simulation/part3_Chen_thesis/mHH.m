function [ out ] = mHH( w )

global ELECTRON_MASS;
out = 0.5 + 0.2*w;
out = out*ELECTRON_MASS;


end

