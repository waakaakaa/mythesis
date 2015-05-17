clear;
h = 6e-9; % ���
H = 34e-9; % �ݿ�

length = 0.1e-9; %�������

%wx = 0.53; % ����С����
%wy = 1;
%bx = 0.82;
%by = 0.38;

%bx = 1;   %T.~Mei,"kֵ�����"
%by = 0;
%wx = 0.71;
%wy = 0.29;

wx = 0.53; % �̿���
wy = 1;
bx = 1;
by = 0;

N = 1; % ��������Ŀ
offset = 0.6;
z = (-N*(H+h)/2):length:(N*(H+h)/2); % ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0;
Ld_min = 0.5e-9;
Ld_length = 0.5e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    V_c(i,:) = Vc(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_hh(i,:) = Vhh(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_lh(i,:) = Vlh(length,h,H,CIn(i,:),CAs(i,:),N,offset);

    [vectc,valuec] = clevel(length,N,H,h,V_c(i,:)+E(i,N*(H+h)/2/length),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
    
    [vecthh,valuehh] = hhlevel(length,N,H,h,V_hh(i,:),CIn(i,:),CAs(i,:));
    Ehh(i,:) = sum(valuehh);
    
    [vectlh,valuelh] = lhlevel(length,N,H,h,V_lh(i,:),CIn(i,:),CAs(i,:));
    Elh(i,:) = sum(valuelh);
end
plot(Ld,(Ec(:,1)+Elh(:,1))/1.6e-19);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.25;
Ld_min = 0.5e-9;
Ld_length = 0.5e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    V_c(i,:) = Vc(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_hh(i,:) = Vhh(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_lh(i,:) = Vlh(length,h,H,CIn(i,:),CAs(i,:),N,offset);

    [vectc,valuec] = clevel(length,N,H,h,V_c(i,:)+E(i,N*(H+h)/2/length),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
    
    [vecthh,valuehh] = hhlevel(length,N,H,h,V_hh(i,:),CIn(i,:),CAs(i,:));
    Ehh(i,:) = sum(valuehh);
    
    [vectlh,valuelh] = lhlevel(length,N,H,h,V_lh(i,:),CIn(i,:),CAs(i,:));
    Elh(i,:) = sum(valuelh);
end
plot(Ld,(Ec(:,1)+Elh(:,1))/1.6e-19);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.5;
Ld_min = 0.5e-9;
Ld_length = 0.5e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    V_c(i,:) = Vc(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_hh(i,:) = Vhh(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_lh(i,:) = Vlh(length,h,H,CIn(i,:),CAs(i,:),N,offset);

    [vectc,valuec] = clevel(length,N,H,h,V_c(i,:)+E(i,N*(H+h)/2/length),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
    
    [vecthh,valuehh] = hhlevel(length,N,H,h,V_hh(i,:),CIn(i,:),CAs(i,:));
    Ehh(i,:) = sum(valuehh);
    
    [vectlh,valuelh] = lhlevel(length,N,H,h,V_lh(i,:),CIn(i,:),CAs(i,:));
    Elh(i,:) = sum(valuelh);
end
plot(Ld,(Ec(:,1)+Elh(:,1))/1.6e-19);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
Ld_min = 0.5e-9;
Ld_length = 0.5e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    V_c(i,:) = Vc(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_hh(i,:) = Vhh(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_lh(i,:) = Vlh(length,h,H,CIn(i,:),CAs(i,:),N,offset);

    [vectc,valuec] = clevel(length,N,H,h,V_c(i,:)+E(i,N*(H+h)/2/length),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
    
    [vecthh,valuehh] = hhlevel(length,N,H,h,V_hh(i,:),CIn(i,:),CAs(i,:));
    Ehh(i,:) = sum(valuehh);
    
    [vectlh,valuelh] = lhlevel(length,N,H,h,V_lh(i,:),CIn(i,:),CAs(i,:));
    Elh(i,:) = sum(valuelh);
end
plot(Ld,(Ec(:,1)+Elh(:,1))/1.6e-19);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 2;
Ld_min = 0.5e-9;
Ld_length = 0.5e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    V_c(i,:) = Vc(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_hh(i,:) = Vhh(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_lh(i,:) = Vlh(length,h,H,CIn(i,:),CAs(i,:),N,offset);

    [vectc,valuec] = clevel(length,N,H,h,V_c(i,:)+E(i,N*(H+h)/2/length),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
    
    [vecthh,valuehh] = hhlevel(length,N,H,h,V_hh(i,:),CIn(i,:),CAs(i,:));
    Ehh(i,:) = sum(valuehh);
    
    [vectlh,valuelh] = lhlevel(length,N,H,h,V_lh(i,:),CIn(i,:),CAs(i,:));
    Elh(i,:) = sum(valuelh);
end
plot(Ld,(Ec(:,1)+Elh(:,1))/1.6e-19);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 6;
Ld_min = 0.5e-9;
Ld_length = 0.5e-9;
Ld_max = 10e-9;
Ld = Ld_min:Ld_length:Ld_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:((Ld_max-Ld_min)/Ld_length+1)
    CIn(i,:) = composition(length,h,H,bx,wx,N,Ld_min+(i-1)*Ld_length);
    CAs(i,:) = composition(length,h,H,by,wy,N,(Ld_min+(i-1)*Ld_length)*k);
    E(i,:) = potential(CIn(i,:),CAs(i,:));
    V_c(i,:) = Vc(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_hh(i,:) = Vhh(length,h,H,CIn(i,:),CAs(i,:),N,offset);
    V_lh(i,:) = Vlh(length,h,H,CIn(i,:),CAs(i,:),N,offset);

    [vectc,valuec] = clevel(length,N,H,h,V_c(i,:)+E(i,N*(H+h)/2/length),CIn(i,:),CAs(i,:));
    Ec(i,:) = sum(valuec);
    
    [vecthh,valuehh] = hhlevel(length,N,H,h,V_hh(i,:),CIn(i,:),CAs(i,:));
    Ehh(i,:) = sum(valuehh);
    
    [vectlh,valuelh] = lhlevel(length,N,H,h,V_lh(i,:),CIn(i,:),CAs(i,:));
    Elh(i,:) = sum(valuelh);
end
plot(Ld,(Ec(:,1)+Elh(:,1))/1.6e-19);