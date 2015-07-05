function m = mlh_perp(In,As)

m = 0.88*(1-In).*As + 0.024*In.*As + 0.16*(1-In).*(1-As) + 0.12*In.*(1-As);
m = 9.109534e-31 * m;