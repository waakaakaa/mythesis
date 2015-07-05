function [ out ] = C11( x )

out = 11.9 + 0.12*x;
out = out * 1e11 * 0.1; % 1dyn/cm2 = 0.1pa

end

