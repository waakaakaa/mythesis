function [ magnitude_square ] = wave_function( x,Vx,Mx,eigenvalue )

global H_BAR;
global delta_x;

A = x;
B = x;
A(1) = 0;
B(1) = 1;
k = 1/H_BAR * sqrt( 2*Mx .* (eigenvalue - Vx) );
beta = k ./ Mx;
for j=1:length(x)-1
    T1 = (1/(2*beta(j+1))) * [beta(j+1)+beta(j)  beta(j+1)-beta(j);beta(j+1)-beta(j)  beta(j+1)+beta(j)];
    T2 = [exp(1i*k(j)*delta_x) 0;0 exp(-1i*k(j)*delta_x)];
    M = T1*T2*[A(j);B(j)];
    A(j+1) = M(1,1);
    B(j+1) = M(2,1);
end

magnitude_square = (abs(A+B)).^2;

end