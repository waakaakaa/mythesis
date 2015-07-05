function m = mhh_perp(In,As)

m = 0.5*(1-In).*As + 0.41*In.*As + 0.54*(1-In).*(1-As) + 0.12*In.*(1-As);
m = 9.109534e-31 * m;