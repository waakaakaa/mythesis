function [ out ] = mhh( x )

global ELECTRON_MASS;
out = 0.50 + 0.2*x;
out = out*ELECTRON_MASS;

end