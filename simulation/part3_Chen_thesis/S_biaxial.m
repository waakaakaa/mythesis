function [ out ] = S_biaxial( w )

a = -1/3 * dEgdP(w) .* (  c11(w) + 2*c12(w)  );
out = -2 * a .* strain(w) .* (  1 - c12(w)/c11(w)  );

end