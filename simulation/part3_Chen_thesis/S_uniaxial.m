function [ out ] = S_uniaxial( w )

out = -b(w) .* (  1 + 2 * c12(w)/c11(w)  ) .* strain(w);

end

