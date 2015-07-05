function [ out ] = mLHpara( w )

global ELECTRON_MASS;
out = 0.23 + 0.1223*w;
out = out*ELECTRON_MASS;

end

