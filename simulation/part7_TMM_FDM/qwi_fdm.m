function [chh,clh] = qwi_fdm(LdIII,k)
%  qwi function (single InP-QW version)
%
%  @author   X.Zhang
%  @date     July 3rd, 2010
%  @Email    zju.zhangxin@gmail.com
%  注意      LdIII避免用0，否则会引发0/0运算
if LdIII==0
    LdIII=1e-12;
end
% Composition
CIn = compositionIII(LdIII);
CAs = compositionV(LdIII*k);
% Potential
Pc = vc(CIn,CAs);
Phh = vhh(CIn,CAs);
Plh = vlh(CIn,CAs);
% Energy
[Fc,Ec] = levelc_fdm(Pc,CIn,CAs);
[Fhh,Ehh] = levelhh_fdm(Phh,CIn,CAs);
[Flh,Elh] = levellh_fdm(Plh,CIn,CAs);
Ec = sum(Ec);
Ehh = sum(Ehh);
Elh = sum(Elh);
% Transition energy
chh = Ec(1,1)+Ehh(1,1);
clh = Ec(1,1)+Elh(1,1);