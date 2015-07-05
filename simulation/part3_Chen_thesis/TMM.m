function [ eigenvalues ] = TMM( x,Vx,Mx,EE ,how_many_solutions)

global H_BAR;
global delta_x;

TT = EE;
for E=EE
    %%% wavenumber
    k = 1/H_BAR * sqrt( 2*Mx .* (E - Vx) );
    beta = k ./ Mx;
    
    %%% transfer matrix
    T = [1 0;0 1];
    for j=1:length(x)-1
        T1 = (1/(2*beta(j+1))) * [beta(j+1)+beta(j)  beta(j+1)-beta(j);beta(j+1)-beta(j)  beta(j+1)+beta(j)];
        T2 = [exp(1i*k(j)*delta_x) 0;0 exp(-1i*k(j)*delta_x)];
        T = T1*T2 * T;
    end
    TT(EE==E) = abs(T(2,2));
end

eigenvalues=zeros(1,how_many_solutions);
i=1;
for index=2:length(EE)-1
    if TT(index)<TT(index-1) && TT(index)<TT(index+1)
        eigenvalues(1,i) = EE(index);
        i =i+1;
        if i>how_many_solutions
            break;
        end
    end
end


end