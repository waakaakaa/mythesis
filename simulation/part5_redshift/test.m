clear;
h = 5e-9; % 阱宽
H = 10e-9; % 垒宽

length = 0.01e-9;

bx = 0.82;
by = 0.38;
wx = 0.53;
wy = 1;

N = 1; % 量子阱数目

%offset = 0.6;

z = (-N*(H+h)/2):length:(N*(H+h)/2); % 定义域

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CIn0 = composition(length,h,H,bx,wx,N,0);
CIn1 = composition(length,h,H,bx,wx,N,0.5e-9);
CIn2 = composition(length,h,H,bx,wx,N,2.5e-9);
CIn3 = composition(length,h,H,bx,wx,N,5e-9);

CAs0 = composition(length,h,H,by,wy,N,0);
CAs1 = composition(length,h,H,by,wy,N,0.5e-9);
CAs2 = composition(length,h,H,by,wy,N,2.5e-9);
CAs3 = composition(length,h,H,by,wy,N,5e-9);

E0 = potential(CIn0,CAs0); % *offset?
E1 = potential(CIn1,CAs1);
E2 = potential(CIn2,CAs2);
E3 = potential(CIn3,CAs3);

E0 = E0 - E0(1,N*(H+h)/2/length); %阱数目为偶数的时候需要修正！！！
E1 = E1 - E1(1,N*(H+h)/2/length);
E2 = E2 - E2(1,N*(H+h)/2/length);
E3 = E3 - E3(1,N*(H+h)/2/length);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
In = CIn0;
As = CAs0;

h_bar = 6.63e-34 / 2 / 3.14;

m = 0.0632*(1-In).*As + 0.0213*In.*As + 0.17*(1-In).*(1-As) + 0.077*In.*(1-As);
m = 9.109534e-31 * m;

x = N*(H+h)/length+1;


A(1,x) = h_bar^2 / (2*length^2*m(1,x));
B(1,1) = h_bar^2 / (2*length^2*m(1,1));
for i = 1:(x-1)
    A(1,i) = h_bar^2 / ( length^2 * ( m(1,i) + m(1,i+1) ) );
end
for i = 2:x
    B(1,i) = h_bar^2 / ( length^2 * ( m(1,i) + m(1,i-1) ) );
end
C = A+B+E0;

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

[F,E] = eig(W);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%F = F.*F;
E = sum(E);


subplot(2,1,1);
plot(z,F(:,1),z,F(:,2));
subplot(2,1,2);
plot(z(1,1:4),E(1,1:4),'r.');hold on;
plot(z,E0,'b.');

%t = 1:x;
%plot(t,A,t,B,t,C);

E(1,1)
E(1,2)
E(1,3)
E_infinite = 3.14^2 * h_bar^2 / 2 / h^2 / m(1,(x+1)/2)