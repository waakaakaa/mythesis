%% constants
const = 1.239841905763867;
global ELEMENTARY_CHARGE;
global ELECTRON_MASS;
global PLANCK_CONSTANT;
global H_BAR;
global delta_x;
global PI;
PI = 3.14159265359;
ELEMENTARY_CHARGE = 1.6e-19;
ELECTRON_MASS = 9.10938291e-31;
PLANCK_CONSTANT = 6.62606957e-34;
H_BAR = PLANCK_CONSTANT/(2*PI);
meV = ELEMENTARY_CHARGE*1e-3;
%% AlGaAs-GaAs QW params
barrier = 0.26;
L = 42;
L_barrier = 200;
%% calculate material parameters
[QW, strain, parameterTable] = HeteroIIIV('(AlGaAs,200,0.26)(GaAs,42)(AlGaAs,200,0.26)',300,'GaAs',100);
ee = strain(1,1,:)+strain(2,2,:)+strain(3,3,:);ee=reshape(ee,length(ee),1);

%% x
delta_x = L/100;
x = (  -(L_barrier+L/2) : delta_x : (L_barrier+L/2)  );

temp = find(x>=-L/2);
well_start = temp(1);
temp = find(x<=L/2);
well_stop = temp(length(temp));

%% C
Vc = x;
Vc(:) = QW(1,30);
Vc(well_start:well_stop) = QW(2,30);
Vc = Vc*ELEMENTARY_CHARGE;
Mc = x;
Mc(:) = QW(1,2);
Mc(well_start:well_stop) = QW(2,2);
Mc = Mc*ELECTRON_MASS;
Ec = TMM( x,Vc-min(Vc),Mc,(1:1:100)*meV ,1)+min(Vc);


%% HH
Vhh = x;
Vhh(:) = QW(1,31);
Vhh(well_start:well_stop) = QW(2,31);
Vhh = abs(Vhh)*ELEMENTARY_CHARGE;
Mhh = x;
Mhh(:) = QW(1,21);
Mhh(well_start:well_stop) = QW(2,21);
Mhh = Mhh*ELECTRON_MASS;
Ehh = TMM( x,Vhh-min(Vhh),Mhh,(1:1:100)*meV ,1)+min(Vhh);

Ec/meV+Ehh/meV
