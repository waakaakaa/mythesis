function [ out ] = mLH( w )

global ELECTRON_MASS;
out = 0.088 + 0.0372*w + 0.0162*w.^2;
out = out*ELECTRON_MASS;

end

