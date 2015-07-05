function [chh clh] = qwi_tmm(LdIII,k)
%  qwi function
%
%  @author   X.Zhang
%  @date     July 3rd, 2013
%  @Email    zju.zhangxin@gmail.com
%  注意      入参LdIII避免用0，否则会引发0/0运算
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy
Ec = levelc_tmm(Pc,CIn,CAs);
Ehh = levelhh_tmm(Phh,CIn,CAs);
Elh = levellh_tmm(Plh,CIn,CAs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition energy
chh = Ec + Ehh;
clh = Ec + Elh;