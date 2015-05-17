function C = composition(length,h,H,Nb,Nw,N,Ld)

z = (-N*(H+h)/2):length:(N*(H+h)/2); % 定义域

C = 0;  % 组分初始值

for i = 1:N
    j = ((N-i)*(H+h) + H*(1+(-1)^i)/2 + h*(1-(-1)^i)/2 )  /  2;
    C = C + (-1)^(i+1)*   (  2  -  erf((j-z)/Ld/2)  -  erf((j+z)/Ld/2)  );
end

C = Nb*(1+(-1)^N)/2 + Nw*(1-(-1)^N)/2  +  (Nb-Nw)*C/2;