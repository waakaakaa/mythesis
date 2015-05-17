clear;
h = 5e-9; % ���
H = 10e-9; % �ݿ�

length_z = 0.01e-9;

bx = 0.82;
by = 0.38;
wx = 0.53;
wy = 1;

N = 1; % ��������Ŀ

offset = 0.6;

z = (-N*(H+h)/2):length_z:(N*(H+h)/2); % ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
Ld_min = 0;
Ld_length = 0.01e-9;
Ld_max = 2.5e-9;
Ld = Ld_min:Ld_length:Ld_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length_z,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length_z,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = offset*potential(CIn(i,:),CAs(i,:));
    E(i,:) = E(i,:) - E(i,round(N*(H+h)/2/length_z));
end

%plot(z,E(:,:));
%title('Potential Curve ');
%legend('Ld=0','Ld=0.5','Ld=2.5','Ld=5');
%xlabel('z/nm');
%ylabel('Potential');
%axis([-N*(H+h)/2 N*(H+h)/2 0 1]);

for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    [vect,value] = clevel(length_z,N,H,h,E(i,:),CIn(i,:),CAs(i,:));
    E(i,:) = sum(value);
end
plot(Ld,E(:,1),'r.');
%plot(z(1,1:4),E(1,1:4),'r.');