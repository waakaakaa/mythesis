function m = mc(In,As)

m = 0.0632*(1-In).*As + 0.0213*In.*As + 0.17*(1-In).*(1-As) + 0.077*In.*(1-As);%InGaAsP电子的有效质量
%m = 0.0632 + 0.0856*x + 0.0231*x.^2;%AlGaAs电子的有效质量
m = 9.109534e-31 * m;