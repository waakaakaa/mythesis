function [ out ] = dEgdP( x )

out = 11.5 - 1.3*x;
out = out * 1e-6 * 1.60e-24; % 1eV/bar = 1.60217657 ¡Á 10-24m^3

end

