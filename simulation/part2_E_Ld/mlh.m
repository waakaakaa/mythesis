function [ out ] = mlh( x )

global ELECTRON_MASS;
out = 0.088 + 0.0372*x + 0.0162*x.^2;
out = out*ELECTRON_MASS;

end

