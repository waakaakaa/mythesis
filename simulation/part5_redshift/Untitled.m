clear;
length = 0.01; % 间隔
N = 2000; % 点子总数
x = -(N-1)*length/2:length:(N-1)*length/2;
V = x.^2 * 1e-22; % meV量级

h_bar = 6.63e-34 / 2 / 3.14;
%h_bar = 1;
m = 9.109534e-34; %电子量级

R = h_bar^2  /  ( (length)^2 * m );
A = 2*R + V;
W = zeros(N,N);
W(1,1) = A(1,1);
W(1,2) = -R;
W(N,N) = A(1,N);
W(N,N-1) = -R;
for i = 2:(N-1)
W(i,i) = A(1,i);
W(i,i+1) = -R;
W(i,i-1) = -R;
end

[F,E] = eig(W); % 对角化
%F = F.*F;
E = sum(E);
plot(x,F(:,1));
E(1,1:5)