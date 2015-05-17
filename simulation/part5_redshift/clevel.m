function [vect,value] =  clevel(length,N,H,h,V,In,As)

h_bar = 6.626196e-34 / 2 / 3.1415926;
m = 0.0632*(1-In).*As + 0.0213*In.*As + 0.17*(1-In).*(1-As) + 0.077*In.*(1-As);%InGaAsP电子的有效质量
%m = 0.0632 + 0.0856*x + 0.0231*x.^2;%AlGaAs电子的有效质量
m = 9.109534e-31 * m;

x = round(N*(H+h)/length+1);

A(1,x) = h_bar^2 / (2*length^2*m(1,x));
B(1,1) = h_bar^2 / (2*length^2*m(1,1));
for i = 1:(x-1)
    A(1,i) = h_bar^2 / (length^2 * (m(1,i)+m(1,i+1)));
end
for i = 2:x
    B(1,i) = h_bar^2 / (length^2 * (m(1,i)+m(1,i-1)));
end
C = A+B+V;

W = zeros(x,x);
W(1,1) = C(1,1);
W(1,2) = -A(1,1);
W(x,x) = C(1,x);
W(x,x-1) = -B(1,x);

for i = 2:(x-1)
    W(i,i) = C(1,i);
    W(i,i+1) = -A(1,i);
    W(i,i-1) = -B(1,i);
end

[vect,value] = eig(W);