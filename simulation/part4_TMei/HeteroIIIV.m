function [QW, strain, parameterTable] = HeteroIIIV(layers,T,substrate,orientation,material,rework)
% function [QW, strain, parameterTable] = HeteroIIIV(layers,T,substrate,orientation,material)
% Band parameter calculation for heterostructures of zinc-blende III-V compounds contaning Ga,Al,In,As,P,Sb,N 
%
% This software is written by Dr. Ting Mei, Nanyang Techonogical University, Singapore, March 2005.
% Please get the permission before using this program. The author preserves all rights.
% The author is grateful if being acknowledged when you use the program for publication  
% If you have any question, feedback, or comment, please contact Dr. Mei at:
%     Tel: +65 67904387; Email: etmei@ntu.edu.sg
%
% Revisions
%     Ver 2 (8/2005): add global variables for control
%     Ver 3 (9/2005): incorporate changes form ZincBlendeIIIV to make it faster for large number of layers
%
% References: 
% 1. I. Vurgaftman, J.R. Meyer and L.R. Ram-Mohan, J. Appl. Phys. 89, 5815 (2001)
% 2. K. Boujadaria, etc., "Luttinger-like parameter calculations", Phys Rev B 63, 235302(2001)
% 3. T.B. Bahder, "Eight-band k.p model of strained zinc-blende crystals", Phys Rev B 41, 11992(1990)
% 4. Jasprit Singh, "Electronic and optoelectronic properties of semiconductor structures", Cambridge University Press, 2003. 
%
%
% Input variables
%   layers: 
%        1. a cell array. each cell element contains {'Symbol',double(thickness),double(x),double(y) or 'latticematched'; ...} 
%        2. a string of formula : See InterpFormula.m; e.g. '2((AlGaAs,40,0.45)(InGaAsP,30,0.25,0.3)(GaAs,48)(InAlAsP,40,0.45,latticematched))'
%        3. a numeric array: [d,x;...], or [d,x,y; ...] where d is thickness; x,y is molar fraction; y = 9 for a lattice-matched layer
%   T: substrate temperature (K)
%   substrate: symbol of substrate, e.g. 'GaAs', 'InP', etc.
%   orientation: orientation of the substrate, e.g. 100, 111
%   material: 'InGaAsP', etc. or cell strings;  used only when layers is a numeric array (case 3).
%
% Output variables
%   QW: an [nx34] array containg band parameters of all layers. 
%       Note that the hole mass calculation does not alway make sense as it is based on Luttinger parameters. See Ref[1]
%   strain: a [3x3xn] array storing strain tensors (3x3) of all layers. Each tensor contains
%           [exx,exy,exy;
%            exy,exx,exy;
%            exy,exy,exx]; 
%   parameterTable: parameter table with parameter names corresponding to the array elements in QW'
%   rework: {nx2} cell array, containing bandparameter names to be recalcluated by user-defined formulas 
%           {'parameter name1','0.2*X+0.3*Y-0.1*.X*.Y';...}
%
%
% Examples:
%   [QW, strain, parameterTable] = HeteroIIIV('2((AlGaAs,40,0.45)(InGaAsP,30,0.25,0.1))(GaAs,48)(InGaAsP,40,0.45,latticematched)',0,'GaAs',100)
%   [QW, strain, parameterTable] = HeteroIIIV([40,0.1,0.3;50,0.15,0.35;60,0.2,9],0,'InP',100,'InGaAsP')
%   [QW, strain, parameterTable] = HeteroIIIV([40,0.1;50,0;60,0.2],0,'GaAs',100,'AlGaAs')
%   [QW, strain, parameterTable] = HeteroIIIV({'AlGaAs',40,0.1,[];'GaAs',50,[],[];'InGaAsP',30,0.1,0.2;'InGaAsP',30,0.1,'latticematched'},0,'InP',111)
%   [QW, strain, parameterTable] = HeteroIIIV({'AlGaAs',40,0.1;'InP',50,[];'GaAs',60,[]},0,'InP',100)
%   [QW, strain, parameterTable] = HeteroIIIV({'GaAs',40;'InP',50;'InAs',60},0,'InP',111)
%   ParameterNames = HeteroIIIV(); % return parametername list if no input variable
%   [QW, strain, parameterTable] = HeteroIIIV([40,0.1,0.3;50,0.15,0.35;60,0.2,9],0,'InP',100,'InGaAsP',{'mhh[110]','0.5*X+0.3*Y';'mlh[110]','0.5*X+0.3*Y.^2'})
%
%
% Additional control parameters (optional) to be announced and assigned in the calling program
% global VBOratio: use to determine the way to calculate valence band offset.
%   if VBORatio is assigned, QW(:,9) = - Eg*VBOratio;
%   if VBORatio is empty [], QW(:,9) = VBO, as calculated in Ref[1].
%
% HeteroIIIV calls ZincBlendeIIIV.m, InterpFormula.m 
%
%
%
% Try the following codes:
% 
% clear all; [QW, strain, parameterTable] = HeteroIIIV('(AlGaAs,100,0.3)(InGaAs,30,0.30)(GaAs,50)(AlGaAs,100,0.3)',0,'GaAs',100);
% ee = strain(1,1,:)+strain(2,2,:)+strain(3,3,:);ee=reshape(ee,length(ee),1);
% figure(1);clf;
% subplot(2,1,1);ladderplot(QW(:,34),QW(:,30),'k');xlabel('z (A)'); ylabel('E_c(eV)');legend('E_c');% Ec
% subplot(2,1,2);ladderplot(QW(:,34),QW(:,31),'b');hold on;ladderplot(QW(:,34),QW(:,32),'r');ladderplot(QW(:,34),QW(:,33),'g');% Evh,Evl,Evso
% xlabel('z (A)'); ylabel('E_v(eV)'); legend('E_v_h_h','E_v_l_h','E_v_s_o'); hold off;
% figure(2);clf;subplot(2,1,1);
% ladderplot(QW(:,34),QW(:,30)-QW(:,9)-QW(:,18),'k');hold on; % DEc
% ladderplot(QW(:,34),QW(:,31)-QW(:,9),'b'); %DEvh
% ladderplot(QW(:,34),QW(:,32)-QW(:,9),'r'); %DEvl
% ladderplot(QW(:,34),QW(:,33)-QW(:,9)+QW(:,1),'g'); %DEvso
% xlabel('z (A)'); ylabel('Strain-induced \DeltaE(eV)');legend('\DeltaE_c','\DeltaE_v_h_h','\DeltaE_v_l_h','\DeltaE_v_s_o'); hold off;
% subplot(2,1,2);ladderplot(QW(:,34),ee,'b');xlabel('z (A)');ylabel('Strain');


global binaryCompoundTable bowingCTable T_effect_on_VBO EqABCD BF4D VBOratio

if nargout > 2 ||nargin == 0
    parameterTable ={ ...
        '01. Spin-orbit splitting Dltaso';        '02. Effective electron mass (Gamma valley) meG';'03. Luttinger parameter gama1'; ...
        '04. Luttinger parameter gama2';          '05. Luttinger parameter gama3';                 '06. SO band effective mass mso'; ...
        '07. Interaction energy of CB&VB Ep(eV)'; '08. Kane parameter F';                          '09. Valence bandoffset VBO(eV)'; ...
        '10. Hydrostatic deformation potential on CB ac(eV)'; '11. Hydrostatic deformation potential on VB av(eV)';'12. Sheer deformation potential b(eV)'; ...
        '13. Sheer deformation potential d(eV)';  '14. Elastic stiffness coefficient c11(GPa)';    '15. Elastic stiffness coefficient c12(GPa)'; ...
        '16. Elastic stiffness coefficient c44(GPa)'; '17. Lattice constant alc(T)';               '18. Bandgap at Gamma Eg_Gama'; ...
        '19. Bandgap energy at X Eg_X';           '20. Bandgap energy at L Eg_L';                  '21. Effective hole mass mhh[100]'; ...
        '22. Effective heavy hole mass mhh[110]'; '23. Effective heavy hole mass mhh[111]';        '24. Effective light hole mass mlh[100]'; ...
        '25. Effective light hole mass mlh[110]'; '26. Effective light hole mass mlh[111]';        '27. Modified Luttinger parameter mGama1';...
        '28. Modified Luttinger parameter mGama2';'29. Modified Luttinger parameter mGama3';       '30. Conduction bandedge Ec(eV)'; ...
        '31. Heavy hole valence bandedge Evhh(eV)';'32. Light hole valence bandedge Evlh(eV)';     '33. SO bandedge Evso(eV)'; ...
        '34. Layer thickness'};
    if nargin == 0
        QW = parameterTable;
        return
    end
end

switch nargin
    case 3
        orientation = 100;
        rework = [];
    case 4
        if iscell(orientation)
            rework = orientation;
        else
            rework = [];
        end
    case 5
        if iscell(material)
            rework = material;
        else
            rework = [];
        end
    case 6
        ;
    otherwise
        error('incorrect number of input variables');
end

if ~isnumeric(layers)
    if ischar(layers)
        [layers, expandedformul] = InterpFormula(layers);
    end
    if ~iscell(layers)
        error('Input argument layers error');
    end
    [rlayers,clayers] = size(layers);
    material = layers(:,1);
    d = [layers{:,2}];
    X = NaN*ones(rlayers,1); Y = X;
    xnull = logical(ones(rlayers,1)); ynull = xnull; islm = ~xnull;
    if clayers > 2
        xnull=cellfun('isempty',layers(:,3));
        X(~xnull) = [layers{:,3}]';
    end
    if clayers > 3
        nlength = cellfun('length',layers(:,4));
        ynull= nlength == 0;
        islm = nlength == 14;
        Y(~ynull&~islm) = [layers{~ynull&~islm,4}]';
        Y(islm) = 9;
    end
else
    [rlayers,clayers] = size(layers);
    X = NaN*ones(rlayers,1); Y = X;
    xnull = logical(ones(rlayers,1)); ynull = xnull; islm = ~xnull;
    d = layers(:,1);
    if clayers > 1
        X = layers(:,2); 
        xnull = isnan(X);
    end
    if clayers > 2
        Y = layers(:,3); 
        ynull = isnan(Y);
        islm = Y==9;
   end
    if ischar(material)
        material = cellstr(repmat(material,rlayers,1));
    end
end

isbinary = xnull;
isternary = ynull & ~isbinary;
isquaternary = ~ynull & ~islm;
subs_para = ZincBlendeIIIV(substrate,T);
altc_subs = subs_para(34);
none = logical(zeros(rlayers,1));
done = none;
QW = NaN*ones(43,rlayers);
strain(1:3,1:3,1:rlayers)=0;
ilayer = 1;
while 1
    calc_strain = 1;
    samematerial = [strmatch(material{ilayer},material,'exact')];
    todo = none; todo(samematerial) = 1;
    if isbinary(ilayer)
        todo = todo & isbinary;
        [materialParameter, materialFormula, parameterTablex] = ZincBlendeIIIV(material{ilayer},T);
        materialParameter = repmat(materialParameter,1,sum(todo));
    elseif islm(ilayer)
        todo = todo & islm;
        [materialParameter, materialFormula, parameterTablex] = ZincBlendeIIIV(material{ilayer},T,X(todo),substrate);
        calc_strain = 0;
    elseif isternary(ilayer)
        todo = todo & isternary;
        [materialParameter, materialFormula, parameterTablex] = ZincBlendeIIIV(material{ilayer},T,X(todo));
    elseif isquaternary(ilayer)
        todo = todo & isquaternary;
        [materialParameter, materialFormula, parameterTablex] = ZincBlendeIIIV(material{ilayer},T,[X(todo),Y(todo)]);
    else
        error('The input layers is wrong');
    end
    QW(:,todo) = materialParameter;
    if calc_strain
        ntodo = sum(todo);
        e_para = (altc_subs - QW(34,todo))./QW(34,todo);
        switch orientation
            case 100
                strain(1,1,todo) = reshape(e_para,1,1,ntodo);
                strain(2,2,todo) = reshape(e_para,1,1,ntodo);
                strain(3,3,todo) = reshape(-2*e_para.*QW(32,todo)./QW(31,todo),1,1,ntodo); 
            case 111
                exx = e_para.*(2/3-1/3*(2*QW(31,todo)+4*QW(32,todo)-4*QW(34,todo))./(QW(31,todo)+2*QW(32,todo)+4*QW(34,todo)));
                exy = e_para.*(-1/3-1/3*(2*QW(31,todo)+4*QW(32,todo)-4*QW(34,todo))/(QW(31,todo)+2*QW(32,todo)+4*QW(34,todo)));
                strain(1,1,todo) = reshape(exx, 1,1,ntodo);
                strain(2,2,todo) = reshape(exx, 1,1,ntodo);
                strain(3,3,todo) = reshape(exx, 1,1,ntodo);
                strain(1,2,todo) = reshape(exy, 1,1,ntodo);
                strain(1,2,todo) = reshape(exy, 1,1,ntodo);
                strain(2,1,todo) = reshape(exy, 1,1,ntodo);
                strain(2,3,todo) = reshape(exy, 1,1,ntodo);
                strain(3,1,todo) = reshape(exy, 1,1,ntodo);
                strain(3,2,todo) = reshape(exy, 1,1,ntodo);
                %             [exx,exy,exy;
                %              exy,exx,exy;
                %              exy,exy,exx]; 
            otherwise
                error(['substrate orientation ',nm2str(orientation),' is not considered']);
                return
        end
    end
    done(todo) = 1;
    [vidx,idx] = sort(done);
    if vidx(1)
        break
    else
        ilayer = idx(1);
    end
end    
 
QW=QW';
QW(:,14:19) = [];
QW(:,1:11) = [];
QW(:,27) = QW(:,3)-QW(:,7)./QW(:,18)/3; % modified Luttinger parameter mGama1 = Gama1-Ep/3Eg
QW(:,28) = QW(:,3)-QW(:,7)./QW(:,18)/6; % modified Luttinger parameter mGama2 = Gama2-Ep/6Eg
QW(:,29) = QW(:,3)-QW(:,7)./QW(:,18)/6; % modified Luttinger parameter mGama3 = Gama3-Ep/6Eg
ee = reshape(strain(1,1,:)+strain(2,2,:)+strain(3,3,:),rlayers,1); % exx+eyy+ezz
p = QW(:,11).*ee;
q = QW(:,12).*reshape(strain(3,3,:)-(strain(1,1,:)+strain(2,2,:))/2,rlayers,1);

if ~isempty(VBOratio)
    QW(:,9) = -QW(:,18)*VBOratio;
end
QW(:,30) = QW(:,9)+QW(:,18)+QW(:,10).*ee; % Ec band profile = Eg+VBO+ac*(exx+eyy+ezz) 
QW(:,31) = QW(:,9)-p-q;                   % HH band profile = VBO-p-q;  p=av(exx+eyy+ezz)
QW(:,32) = QW(:,9)-p+q;                   % LH band profile = VBO-p+q;  q=b(ezz-(exx+eyy)/2)
QW(:,33) = QW(:,9)-QW(:,1)-p;             % SO band profile = VBO-DSO-p
QW(:,34) = d';

if ~isempty(rework)
    npara = size(rework,1);
    for ipara = 1:npara
        itable = strmatch(rework{ipara,1},parameterTable);
        QW(:,itable) = eval(rework{ipara,2});
    end
end

if nargout > 2
    parameterTable(:,2:rlayers+1) = num2cell(QW');
end
return
