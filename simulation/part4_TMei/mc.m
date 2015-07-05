function [ out ] = mc( w )

global ELECTRON_MASS;
out = 0.0632 + 0.0856*w + 0.0231*w.^2;
out = out*ELECTRON_MASS;

end

