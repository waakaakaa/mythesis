function [ out ] = mc( x )

global ELECTRON_MASS;
out = 0.0632 + 0.0856*x + 0.0231*x.^2;
out = out*ELECTRON_MASS;

end