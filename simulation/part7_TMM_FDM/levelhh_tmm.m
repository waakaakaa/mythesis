function E = levelhh_tmm(V,CIn,CAs)

global h_bar;
global z;
global length_z;

x = numel(z);
m = mhh_perp(CIn,CAs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = 0;
Told = 1;
while 1
    k = 1/h_bar * sqrt(2*m.*(E-V));
    beta = k./m;
    T = [1,0;0,1];
    for j=1:(x-1)
        A = zeros(2,2);
        B = zeros(2,2);
        A(1,1) = beta(1,j+1) + beta(1,j);
        A(2,2) = beta(1,j+1) + beta(1,j);
        A(1,2) = beta(1,j+1) - beta(1,j);
        A(2,1) = beta(1,j+1) - beta(1,j);
        B(1,1) = exp(1i*k(1,j)*length_z);
        B(2,2) = exp(-1i*k(1,j)*length_z);
        T = T * 0.5/beta(1,j+1) * A * B;
    end
    Tnew = real(T(2,2));
    if Told*Tnew<=0
        break;
    else
        Told = Tnew;
    end
    E = E + 1.6e-19/1000;
end