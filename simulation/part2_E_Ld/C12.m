function [ out ] = C12( x )

out = 5.38 - 0.08*x;
out = out * 1e11 * 0.1; % 1dyn/cm2 = 0.1pa

end

