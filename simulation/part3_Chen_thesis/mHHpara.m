function [ out ] = mHHpara( w )

global ELECTRON_MASS;
out = 0.11 + 0.0668*w;
out = out*ELECTRON_MASS;

end

