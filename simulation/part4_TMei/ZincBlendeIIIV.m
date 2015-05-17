function [materialParameter, materialFormula, parameterTable, BF3] = ZincBlendeIIIV(material, T, composition, substrate)
% function [materialParameter, materialFormula, parameterTable, BF3] = ZincBlendeIIIV(material, T, composition, substrate)
% Band parameter calculation for the zinc-blende III-V compounds contaning Ga,Al,In,As,P,Sb,N 
%
% This software is written by Dr. Ting Mei, Nanyang Techonogical University, Singapore, March 2005.
% Please get the permission before using this program. The author preserves all rights.
% The author is grateful if being acknowledged when you use the program for publication  
% If you have any question, feedback, or comment, please contact Dr. Mei at:
%     Tel: +65 67904387; Email: etmei@ntu.edu.sg
%
% Revisions  
%   Ver 2(8/2005): this version includes the used parameters in calculatin in parameterTable for alloy materials
%   Ver 3(9/2005): this version accepts a series of compositions [x], or [x,y] as input for calculation
%   Ver 4(7/2006): add calculation for dilute nitride (based on Ref[7])
%
% References: 
% 1. I. Vurgaftman, J.R. Meyer and L.R. Ram-Mohan, J. Appl. Phys. 89, 5815 (2001) (Major reference for the program)
% 2. M.A. Littlejohn, JR. Hauser, and T.H. Glisson, Appl Phys Lett 30, 242(1977)
% 3. T. Mei, unpublished
% 4. G.P. Donati, R. Kaspi, K.J.Malloy, J. Appl. Phys. 94, 5814(2003)
% 5. Adachi, Sadao, 1950-, Physical properties of III-V semiconductor compounds : InP, InAs, GaAs, GaP, InGaAs, and InGaAsP, Wiley, c1992
% 6. I.Shikawa, J.E.Bowers, IEEE J. Quantum Elect. 30, 562(1994)
% 7. I. Vurgaftman and J.R. Meyer, J. Appl. Phys. 94, 3675 (2003)
%
%
% Input variables:
%
% material: material formula for binary, ternary, and quaternary compounds. Dilute nitrides are included after Version 4.
% - Standard atomic symbols like Ga,As,P,etc. It is case sensitive. Forms like GA or ga are not recognized.
% - Group V elements must be placed after group III elements.
% - Example: 'GaAs', 'InGaAs', 'InAsP', 'InGaAsP', 'AlInGaSb', etc.
%
% T: 0, 77, 300, etc. (temperature K).
%
% composition:
% - For binary compounds AB: the variable should be void.
% - For ternary compounds A(x)B(1-x)C, or AB(x)C(1-x): 
%   A scalar x or a [nx1] column vector X, assigned to 1st element in the same group. e.g. x-> Al for AlGaAs but ->Ga for GaAlAs 
% - For quaternary compounds A(x)B(1-x)C(y)D(1-y), A(x)B(y)C(1-x-y)D, or AB(x)C(y)D(1-x-y): 
%   1D vector [x,y], e.g. [0.1,0.2]; or [nx2] array: [X,Y] 
%   x,y are assigned to the element shown in the material formula by sequence. e.g.  x->In,y->As for InGaAsP; x->In,y->Al for InAlGaAs
% - For quaternary compounds latticed-matched to substrate:
%   A scalar z, or a [nx1] column vector Z specific to the following
%   - Lattice-matched to GaAs:
%       - AlGaInP: [Al(0.52)In(0.48)P](z)/[Ga(0.51)In(0.49)P](1-z)
%       - GaInAsP: [Ga(0.51)In(0.49)P](z)/[GaAs](1-z)
%       - AlGaInAs: not programmed
%   - Lattice-matched to InP:
%       - GaInAsP: [Ga(0.47)In(0.53)As](z)/[InP](1-z)
%       - AlGaInAs: [Al(0.48)In(0.52)As](z)/[Ga(0.47)In(0.53)As](1-z)
%       - GaInAsSb: [Ga(0.47)In(0.53)As](z)/[GaAs(0.5)Sb(0.5)](1-z)
%   - Lattice-matched to InAs:
%       - GaInAsSb: [GaAs(0.08)Sb(0.92)](z)/[InAs](1-z)
%       - AlGaAsSb: [AlAs(0.16)Sb(0.84](z)/[GaAs(0.08)Sb(0.92)](1-z)
%       - InAsPSb: [InSb(0.31)P(0.69)](z)/[InAs](1-z)
%   - Lattice-matched to GaSb:
%       - GaInAsSb: [InAs(0.91)Sb(0.09)](z)/[GaSb](1-z)
%       - AlGaAsSb: [AlAs(0.08)Sb(0.92](z)/[GaSb](1-z)
%   - Lattice-matched to GaAs(0.61)P(0.39):
%       - AlGaAsP: [AlAs(0.61)P(0.39)](z)/[GaAs(0.61)P(0.39)](1-z)
%
% substrate: 'GaAs', 'InP', 'InAs', 'GaSb', or 'GaAsP'.
%
%
% Output variables:
%
% materialParameter: 1D array  containing the following parameters
%   01. a1c1(A)  02. a2c1(1e-5 A/K) 03. Eg0G(eV) 04. alfaG(meV/K) 05. betaG(K) 06. Eg0X(eV) 07. alfaX(meV/K) 08. betaX(K)     
%   09. Eg0L(eV) 10. alfaL(meV/K)   11. betaL(K) 12. Dltaso       13. meG      14. mlL      15. mtL          16. mDOSL
%   17. mlX      18. mtX            19. mDOSX    20. gama1        21. gama2    22. gama3    23. mso          24. Ep(eV)  
%   25. F        26. VBO(eV)        27. ac(eV)   28. av(eV)       29. b(eV)    30. d(eV)    31. c11(GPa)     32. c12(GPa)
%   33. c44(GPa) 34. alc = alc1+alc2*1e-5*(T-300)                 
%   35. EgG=Eg0-alfa*T^2/(T+beta)   36. EgX=Eg0-alfa*T^2/(T+beta)              37. EgL=Eg0-alfa*T^2/(T+beta) 
%   38. mhh100 = 1/(gama1-2*gama2)  39. mhh110 = 2/(2*gama1-gama2-3*gama3)     40. mhh111 = 1/(gama1-2*gama3)
%   41. mlh100 = 1/(gama1+2*gama2)  42. mlh110 = 2/(2*gama1+gama2+3*gama3)     43. mlh111 = 1/(gama1+2*gama3)
%   VBO: valence band offset w.r.t. InSb valence band maxima
% - Elements #01 ~ #11 make sense only for binary compounds 
% - Elements are trival if their values are NaN or 1e-100
%
% materialFormula: some information about material formula and composition expressed in a string format
%
% parameterTable: a table of parameters and all parameters used in calculation; help user read data in materialParameter
%
% BF3: cases other than the two below will make error if BF3 is included
%  (1) bowing factors of a ternary compound. The author used it for recursive call in quaternary compound calculation
%  (2) BAC model parameters for dilute nitride [EN, V]. See Table V in Ref[7] 
%      You should do post-calculation dilute nitride using BAC model. The Eg(Gama) value of the calculated material in materialParameter is invalid
%                                              __________________  
%         E(+/-)(k) = (1/2){ [Ec(k)+EN] (+/-) /[Ec(k)-EN]^2+4xV^2 }, x: nitrogen molar fraction
%      where Ec(k) is the conduction-band dispersion of the unperturbed no-nitride compound, obtainable from parameterTable 
%      Validity of other band parameters need to be checked in future
%
%
% Examples:
%   [materialParameter, materialFormula, parameterTable] = ZincBlendeIIIV('GaAs',300)
%   [materialParameter, materialFormula] = ZincBlendeIIIV('GaInP',77,0.3)
%   [materialParameter, materialFormula] = ZincBlendeIIIV('AlGaInAs',0,[0.2,0.3])
%   materialParameter = ZincBlendeIIIV('GaInAsP',300,0.3,'InP')
%   [materialParameter, materialFormula, parameterTable, BF3] = ZincBlendeIIIV('GaAsN',300,0.95)
%
%
% Additional control parameters (optional) to be announced and assigned in the calling program 
%
% global T_effect_on_VBO: set as 1 if you want to consider temperature varying effect of VBO.
%
% global EqABCD: select an equation for calculating a band parameter QA(x)B(1-x)C(y)D(1-y) of quaternary alloys 
%   - EqABCD is empty or 0 (default): use weighted sum of ternary values with no surface bowing consideration. See  Eq.5.2 in Ref[1] or Eq.1. in Ref[2]
%   - EqABCD = 1: use weighted sum of ternary values with surface bowing, modified from Eq.5.2 in Ref[1] or Eq.1 in Ref[2] as proposed by Mei(Ref[3])
%   - EqABCD = 2: use biquardratic formula with D factor in Eq.7&8 in Ref[4](Note the C terms are opposite in sign against that in Ref[1]). See also Eq.A.8 in Ref[5], and Eq.2.1 in Ref[6]
%   - EqABCD should be assigned as an [43x1] array in case you use different equations to calculate different parameters. 
%
% global BF4D: Bowing factor(s) D when EqABCD = 1,2 for quaternary alloys
%   - BF4D is empty or 0 (default): DIII = 0, DV = 0; or D = 0. Used for most cases
%   - BF4D is a [43x1] array: DIII=BF4D,DV=BF4D for the EqABCD = 1 case (Ref[3]); 
%                             D=BF4D for the EqABCD = 2 case (Ref[4])                           
%   - BF4D is a [43x2] array: DIII=BF4D(:,1),DV=BF4D(:,2) for the EqABCD = 1 case (Ref[3]);
%                             D=BF4D(:,1) for the EqABCD = 2 case (Ref[4])
%
% Example: include the following lines in the main program for GaInAsP (Ref[3]) after declaring those globals 
%   EqABCD= zeros(43,1); BF4D = zeros(43,2); 
%   EqABCD(12)=1; BF4D(12,:) = [1.2,0.3]; % SO, weighted sum of ternary values with surface bowing 
%   EqABCD(35)=1;                         % Eg, weighted sum of ternary values with surface bowing
%
% global binaryCompoundTable: use your own set of parameters for binary compounds GaAs, InAs, etc. if not empty
%
% global bowingCTable: use your own set of bowing parameters for tenery alloys InGaAs, AlInAs, etc. if not empty
%                      Write an element in complex form c1 + i*c2 if it is composition dependent, C = c1+c2*x1  
% 
% ZincBlendeIIIV.m calls SetBandParameter.m

global binaryCompoundTable bowingCTable T_effect_on_VBO EqABCD BF4D
global recursion_level

alloy = '';
elementIIIsym = {'Al','Ga','In'};
elementVsym = {'As','P','Sb','N'};
elementIIIcount = 0;            %get the information of how many elements in the input material
elementVcount = 0;
for ielement = 1:3
    symPosition = strfind(material,elementIIIsym{ielement});
    if ~isempty(symPosition) 
        alloy = [alloy, elementIIIsym{ielement}]; 
        elementIIIcount = elementIIIcount + 1;
        elementIIIposition(elementIIIcount) = symPosition;
    end
end
for ielement = 1:4
    symPosition = strfind(material,elementVsym{ielement});
    if ~isempty(symPosition) 
        alloy = [alloy, elementVsym{ielement}]; 
        elementVcount = elementVcount + 1;
        elementVposition(elementVcount) = symPosition;
    end
end                            %**************

if elementIIIcount == 0 || elementVcount == 0
    error([alloy,' - material formula error - missing either group III or V element!']);
    materialParameter = []; materialFormula = '';
    return
elseif elementIIIcount + elementVcount >=5
    error([alloy,' - pentanary and hexanary alloys are not included in this program.']);
    materialParameter = []; materialFormula = '';
    return
else
    alloyType =[num2str(elementIIIcount), num2str(elementVcount)];
    switch alloyType
        case '11'
            if nargin ~= 2         %nargin: Number of function input arguments.
                error([alloy,' - incorrect input! input parameters include material and T only'])
                materialParameter = []; materialFormula = '';
                return
            end                
        case {'12','21'}
            if nargin ~= 3 || size(composition,2) ~= 1
                error([alloy, ' - incorrect input! input parameters include material, T, composition; a scalar or column vector  for composition x in A(x)B(1-x)C or AB(x)C(1-x)']);
                materialParameter = []; materialFormula = '';
                return
            else
                totalpoints = size(composition,1);
            end
        case '22'
            switch nargin
                case 3
                    if size(composition,2) ~= 2
                        error([alloy, ' - incorrect input! input parameters include material, T, composition; an nx2 array or 1x2 row vector for composition [X,Y] in A(x)B(1-x)C(y)D(1-y)'])
                        materialParameter = []; materialFormula = '';
                        return
                    else
                        totalpoints = size(composition,1);
                    end
                case 4
                    if size(composition,2) ~= 1 || ~ischar(substrate)
                        error([alloy, ' - incorrect input! Only a lattice-matched case has 4 input parameters; a scalar or column vector for composition z in ABCD and a string for substrate'])
                        materialParameter = []; materialFormula = '';
                        return
                    else
                        alloyType = 'lattice matched';
                        totalpoints = size(composition,1);
                    end
                otherwise
                    error([alloy,' - wrong number of input parameters (more than 4)!']);
                    materialParameter = []; materialFormula = '';
                    return
            end
        case {'31','13'}
            switch nargin
                case 3
                    if size(composition,2) == 1
                        error([alloy,' - incorrect input! composition [X,Y] should be an nx2 array  for A(x)B(y)C(1-x-y)D or AB(x)C(y)D(1-x-y)'])
                        materialParameter = []; materialFormula = '';
                        return
                    else
                        totalpoints = size(composition,1);
                    end
                case 4
                    if size(composition,2) ~= 1 || ~ischar(substrate)
                        error([alloy,' - incorrect input! Only a lattice-matched case has 4 input parameters; a scalar or column vector for composition z in ABCD and a string for substrate'])
                        materialParameter = []; materialFormula = '';
                        return
                    else
                        totalpoints = size(composition,1);
                        alloyType = 'lattice matched';
                    end
                otherwise
                    error([alloy,' - wrong number of input parameters (more than 4)!']);
                    materialParameter = []; materialFormula = '';
                    return
            end            
    end
end

if ~isempty(strmatch('lattice matched', alloyType))
    z = composition;
else
    [positionInStrIII, indexElementIII] = sort(elementIIIposition);
    switch elementIIIcount
        case 2 % two group III elements, A(x)B(1-x)C or A(x)B(1-x)C(y)D(1-y)
            x(:,indexElementIII(1)) = composition(:,1);
            x(:,indexElementIII(2)) = 1-composition(:,1);    
        case 3 % three group III elements, A(x)B(y)C(1-x-y)D
            x(:,indexElementIII(1)) = composition(:,1);
            x(:,indexElementIII(2)) = composition(:,2);
            x(:,indexElementIII(3)) = 1-composition(:,1)-composition(:,2);
    end 
    [positionInStrV, indexElementV] = sort(elementVposition);
    switch elementVcount
        case 2 % two group V elements
            if elementIIIcount == 1 % AB(x)C(1-x)
                x(:,indexElementV(1)) = composition(:,1);
                x(:,indexElementV(2)) = 1-composition(:,1);  
            else  % A(x)B(1-x)C(y)D(1-y)        
                y(:,indexElementV(1)) = composition(:,2);
                y(:,indexElementV(2)) = 1-composition(:,2);  
            end
        case 3 % three group V elements DA(x)B(y)C(1-x-y)
            x(:,indexElementV(1)) = composition(:,1);
            x(:,indexElementV(2)) = composition(:,2);
            x(:,indexElementV(3)) = 1-composition(:,1)-composition(:,2);
    end
end

if isempty(recursion_level) 
    recursion_level = 0;
else
    recursion_level = recursion_level + 1;
end

binaryTable = {'GaAs', 'AlAs', 'InAs', 'GaP', 'AlP', 'InP', 'GaSb', 'AlSb', 'InSb', 'GaN', 'AlN', 'InN'};
ternaryTable = {'AlGaAs', 'GaInAs', 'AlInAs', 'AlGaP', 'GaInP', 'AlInP','AlGaSb', 'GaInSb', 'AlInSb', 'AlAsP', 'GaAsP', 'InAsP',... 
        'AlPSb', 'GaPSb', 'InPSb', 'AlAsSb', 'GaAsSb', 'InAsSb', 'GaInN', 'AlGaN', 'AlInN', 'GaAsN','GaPN','InPN','InAsN','GaSbN','InSbN'};
switch alloyType
    case '11' % binary compounds
        if ~isempty(binaryCompoundTable)
            binaryCompound = binaryCompoundTable;
        else
%           1.GaAs  2.AlAs  3.InAs  4.GaP   5.AlP   6.InP   7.GaSb  8.AlSb  9.InSb  10.GaN 11.AlN 12.InN    Band parameters
		    binaryCompound = [...
            5.65325,5.6611, 6.0583, 5.4505, 5.4672, 5.8697, 6.0959, 6.1355, 6.4794, 4.50,  4.38,  4.98;...  % 01. a1c1(A) 
            3.88,   2.90,   2.74    2.92,   2.92,   2.79,   4.72,   2.60,   3.48,   0,     0,     0;...     % 02. a2c1(1e-5 A/K)
            1.519,  3.099,  0.417,  2.886,  3.63,   1.4236, 0.812,  2.386,  0.235,  3.299, 5.4,   0.78;...  % 03. Eg0G(eV)
            0.5405, 0.885,  0.276,  0.1081, 0.5771, 0.363,  0.417,  0.42,   0.32,   0.593, 0.593, 0.245;... % 04. alfaG(meV/K)
            204,    530,    83,     164,    372,    162,    140,    140,    170,    600,   600,   624;...   % 05. betaG(K)
%                           93 in Vurgaftman, but a mistake 
            1.981,  2.24,   1.433,  2.35,   2.52,   2.384,  1.141,  1.696,  0.63,   4.52,  4.9,   2.51;...  % 06. Eg0X(eV)
            0.460,  0.70,   0.276,  0.5771, 0.318,  3.7e-4, 0.475,  0.39,   0,      0.593, 0.593, 0.245;... % 07. alfaX(meV/K)
            204,    530,    93,     372,    588,    NaN,    94,     140,    1e-100, 600,   600,   624;...   % 08. betaX(K)
            1.815,  2.46,   1.133,  2.72,   3.57,   2.014,  0.875,  2.329,  0.93,   5.59,  9.3,   5.82;...  % 09. Eg0L(eV)
            0.605,  0.605,  0.276,  0.5771, 0.318,  0.363,  0.597,  0.58,   0,      0.593, 0.593, 0.245;... % 10. alfaL(meV/K)
            204,    204,    93,     372,    588,    162,    140,    140,    1e-100, 600,   600,   624;...   % 11. betaL(K)
            0.341,  0.28,   0.39,   0.08,   0.07,   0.120,  0.76,   0.676,  0.81,   0.017, 0.019, 0.005;... % 12. Dltaso
            0.067,  0.15,   0.026,  0.13,   0.22,   0.0795, 0.039,  0.14,   0.0135, 0.15,  0.25,  0.07;...  % 13. meG
            1.9,    1.32,   0.64,   2.0,    NaN,    NaN,    1.3,    1.64,   NaN,    NaN,   NaN,   NaN;...   % 14. mlL
            0.0754, 0.15,   0.05,   0.253,  NaN,    NaN,    0.10,   0.23,   NaN,    NaN,   NaN,   NaN;...   % 15. mtL
            0.56,   NaN,    0.29,   NaN,    NaN,    0.47,   NaN,    NaN,    0.25,   NaN,   NaN,   NaN;...   % 16. mDOSL
            1.3,    0.97,   1.13,   1.2,    2.68,   NaN,    1.51,   1.357,  NaN,    0.5,   0.53,  0.48;...  % 17. mlX
            0.23,   0.22,   0.16,   0.15,   0.155,  NaN,    0.22,   0.123,  NaN,    0.3,   0.31,  0.27;...  % 18. mtX
            0.85,   NaN,    0.64,   NaN,    NaN,    0.88,   NaN,    NaN,    NaN,    NaN,   NaN,   NaN;...   % 19. mDOSX
            6.98,   3.76,   20.0,   4.05,   3.35,   5.08,   13.4,   5.18,   34.8,   2.70,  1.92,  3.72;...  % 20. gama1
            2.06,   0.82,   8.5,    0.49,   0.71,   1.60,   4.7,    1.19,   15.5,   0.76,  0.47,  1.26;...  % 21. gama2
            2.93,   1.42,   9.2,    2.93,   1.23,   2.10,   6.0,    1.97,   16.5,   1.11,  0.85,  1.63;...  % 22. gama3
            0.172,  0.28,   0.14,   0.25,   0.30,   0.21,   0.12,   0.22,   0.11,   0.29,  0.47,  0.3;...   % 23. mso
%                                                           0.109 in Vurgaftmen, but change to be in accord with Madeluong, Perea
            28.8,   21.1,   21.5,   31.4,   17.7,   20.7,   27.0,   18.7,   23.3,   25.0,  27.1,  17.2;...  % 24. Ep(eV)
            -1.94,  -0.48,  -2.90,  -2.04,  -0.65,  -1.31,  -1.63,  -0.56,  -0.23,  -0.95, -1.01, -4.36;... % 25. F
            -0.80,  -1.33,  -0.59,  -1.27,  -1.74,  -0.94,  -0.03,  -0.41,  0,      -2.64, -3.44, -2.34;... % 26. VBO(eV)(valence band offset w.r.t. InSb valence band maxima)
            -7.17,  -5.64,  -5.08,  -8.2,   -5.7,   -6.0,   -7.5,   -4.5,   -6.94,  -6.71, -4.5,  -2.65;... % 27. ac(eV)
            -1.16,  -2.47,  -1.00,  -1.7,   -3.0,   -0.6,   -0.8,   -1.4,   -0.36,  -0.69, -4.9,  -0.7;...  % 28. av(eV)
            -2.0,   -2.3,   -1.8,   -1.6,   -1.5,   -2.0,   -2.0,   -1.35,  -2.0,   -2.0,  -1.7,  -1.2;...  % 29. b(eV)
            -4.8,   -3.4,    -3.6,  -4.6,   -4.6,   -5.0,   -4.7,   -4.3,   -4.7,   -3.7,  -5.5,  -9.3;...  % 30. d(eV)
            1221,   1250,   832.9,  1405,   1330,   1011,   884.2,  876.9,  684.7,  293,   304,   187;...   % 31. c11(GPa)
            566,    534,    452.6,  620.3,  630,    561,    402.6,  434.1,  373.5,  159,   160,   125;...   % 32. c12(GPa)
            600,    542,    395.9,  703.3,  615,    456,    432.2,  407.6,  311.1,  155,   193,   86;...    % 33. c44(GPa)
            ];
%           1.GaAs  2.AlAs  3.InAs  4.GaP   5.AlP   6.InP   7.GaSb  8.AlSb  9.InSb  10.GaN 11.AlN 12.InN    Band structure parameters
        end
		binaryCompound(34,:) = binaryCompound(1,:) + binaryCompound(2,:)*1e-5*(T-300);                      % 34. alc = alc1+alc2*1e-5*(T-300) 
		binaryCompound(35,:) = binaryCompound(3,:) - 1e-3*binaryCompound(4,:)./(T+binaryCompound(5,:))*T^2; % 35. EgG=Eg0-alfa*T^2/(T+beta) 
		binaryCompound(36,:) = binaryCompound(6,:) - 1e-3*binaryCompound(7,:)./(T+binaryCompound(8,:))*T^2; % 36. EgX=Eg0-alfa*T^2/(T+beta) 
		binaryCompound(37,:) = binaryCompound(9,:) - 1e-3*binaryCompound(10,:)./(T+binaryCompound(11,:))*T^2;% 37. EgL=Eg0-alfa*T^2/(T+beta) 
        binaryCompound(38,:) = 1./(binaryCompound(20,:)-2*binaryCompound(21,:));                             % 38. mhh001 = 1/(gama1-2*gama2) 
        binaryCompound(39,:) = 2./(2*binaryCompound(20,:)-binaryCompound(21,:)-3*binaryCompound(22,:));      % 39. mhh110 = 2/(2*gama1-gama2-3*gama3) 
        binaryCompound(40,:) = 1./(binaryCompound(20,:)-2*binaryCompound(22,:));                             % 40. mhh111 = 1/(gama1-2*gama3)
        binaryCompound(41,:) = 1./(binaryCompound(20,:)+2*binaryCompound(21,:));                             % 41. mlh001 = 1/(gama1+2*gama2) 
        binaryCompound(42,:) = 2./(2*binaryCompound(20,:)+binaryCompound(21,:)+3*binaryCompound(22,:));      % 42. mlh110 = 2/(2*gama1+gama2+3*gama3) 
        binaryCompound(43,:) = 1./(binaryCompound(20,:)+2*binaryCompound(22,:));                             % 43. mlh111 = 1/(gama1+2*gama3)
		
		binaryCompound(35,4) = binaryCompound(3,4) + binaryCompound(4,4)*(1-coth(binaryCompound(5,4)/(T+1e-200)));% 35. GaP: EgG=2.886+0.1081(1-coth(164/T)) 
		binaryCompound(36,6) = binaryCompound(6,6) - binaryCompound(7,6)*T;                                  % 36. InP: EgX=2.384-3.7e-4*T 
        
        if ~isempty(T_effect_on_VBO) && T_effect_on_VBO == 1 % Consider temperature varying effect of VBO, not being reported.
            EgG300K=[1.4225, 3.003, 0.35379, 2.777, 3.5527, 1.3529, 0.7267, 2.3001, 0.17372, 3.2397, 4.8407, 1.9155]; 
            binaryCompound(26,:) = binaryCompound(26,:).*binaryCompound(35,:)./EgG300K; 
        end

        nBinary = strmatch(alloy, binaryTable);
        materialParameter = binaryCompound(:,nBinary); % set output
        materialFormula = binaryTable{nBinary};        % set output
        parameterHistory ={};
    case {'21','12'} % ternary alloys
        nTernary = strmatch(alloy, ternaryTable);        
%                      1      2      3      4     5     6     7      8      9      10    11    12    13    14    15    16     17     18     19    20    21    22    23   24   25
%       select binary components for 
%                      AlGaAs GaInAs AlInAs AlGaP GaInP AlInP AlGaSb GaInSb AlInSb AlAsP GaAsP InAsP AlPSb GaPSb InPSb AlAsSb GaAsSb InAsSb GaInN AlGaN AlInN GaAsN GaPN InPN InAsN      
        binaryIndex = [2,1;   1,3;   2,3;   5,4;  4,6;  5,6;  8,7;   7,9;   8,9;   2,5;  1,4;  3,6;  5,8;  4,7;  6,9;  2,8;   1,7;   3,9;   10,12;11,10;11,12;1,10; 4,10;6,12;3,12]; 
%       where the number means 1.GaAs  2.AlAs  3.InAs  4.GaP   5.AlP   6.InP   7.GaSb  8.AlSb  9.InSb  10.GaN 11.AlN 12.InN    Band structure parameters
        if ~isempty(bowingCTable)
            bowingC = bowingCTable;
        else
%           1      2      3      4     5     6     7      8      9      10    11    12    13    14    15    16     17     18     19    20    21    22    23   24   25    26    27
%           AlGaAs GaInAs AlInAs AlGaP GaInP AlInP AlGaSb GaInSb AlInSb AlAsP GaAsP InAsP AlPSb GaPSb InPSb AlAsSb GaAsSb InAsSb GaInN AlGaN AlInN GaAsN GaPN InPN InAsN GaSbN InSbN     
		    bowingC = [...
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 01. a1c1(A) 
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 02. a2c1(1e-5 A/K)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 03. Eg0G(eV)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 04. alfaG(meV/K)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 05. betaG(K)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 06. Eg0X(eV)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 07. alfaX(meV/K)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 08. betaX(K)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 09. Eg0L(eV)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 10. alfaL(meV/K)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 11. betaL(K)
            0,     0.15,  0.15,  0,    0,    -0.19,0.3,   0.1,   0.25,  0,    0,    0.16, 0,    0,    0.75, 0.15   0.6,   1.2,   0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 12. Dltaso
            0,     0.0091,0.049, 0,    0.051,0.22, 0,     0.0092,0,     0,    0,    0,    0,    0,    0,    0,     0,     0.35,  0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 13. meG
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 14. mlL
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 15. mtL
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 16. mDOSL
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 17. mlX
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 18. mtX
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 19. mDOSX
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 20. gama1
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 21. gama2
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 22. gama3
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 23. mso
            0,     -1.48, -4.81, 0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 24. Ep(eV)
            0,     1.77,  -4.44, 0,    0.78, 0,    0,     -6.84, 0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 25. F
            0,     -0.38, -0.64, NaN,  0,    0,    0,     NaN,   NaN,   NaN,  NaN,  NaN,  NaN,  0,    0,    -1.71, -1.2   0,     NaN,  NaN,  NaN,  NaN,  NaN, NaN, NaN,  NaN,  NaN; ... % 26. VBO(eV)
            0,     2.61,  -1.4,  0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 27. ac(eV)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 28. av(eV)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 29. b(eV)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 30. d(eV)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 31. c11(GPa)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 32. c12(GPa)
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 33. c44(GPa)
		    0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 34. alc = alc1+alc2*1e-5*(T-300) 
		    -0.127+1.310i,...
                   0.477, 0.70,  0,    0.65, -0.48,-0.044+1.22i,...
                                                          0.415, 0.43,  0.22, 0.19, 0.10, 2.7,  2.7,  1.9,  0.8,   1.43,  0.67,  1.4,  0.7,  2.5,  -79.6+100i,...
                                                                                                                                                         3.9, 15   4.22, 0,    0;   ... % 35. EgG=Eg0-alfa*T^2/(T+beta) 
		    0.055, 1.4,   0,     0.13, 0.20, 0.38, 0,     0.33,  0,     0.22, 0.24, 0.27, 2.7,  2.7,  1.9,  0.28,  1.2,   0.6,   0.69, 0.61, 0.61, 0,    10,  0,   0,    0,    0;   ... % 36. EgX=Eg0-alfa*T^2/(T+beta) 
		    0,     0.33,  0,     0,    1.03, 0,    0,     0.4,   0,     0.22, 0.16, 0.27, 2.7,  2.7,  1.9,  0.28,  1.2,   0.6,   1.84, 0.80, 0.80, 0,    0,   0,   0,    0,    0;   ... % 37. EgL=Eg0-alfa*T^2/(T+beta) 
            0,     -0.145,0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 38. mhh001 = 1/(gama1-2*gama2) 
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 39. mhh110 = 2/(2*gama1-gama2-3*gama3) 
            0,     0.0202,0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 40. mhh111 = 1/(gama1-2*gama3)
            0,     0,     0,     0,    0,    0,    0,     0.011, 0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 41. mlh001 = 1/(gama1+2*gama2) 
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 42. mlh110 = 2/(2*gama1+gama2+3*gama3) 
            0,     0,     0,     0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,    0,    0,     0,     0,     0,    0,    0,    0,    0,   0,   0,    0,    0;   ... % 43. mlh111 = 1/(gama1+2*gama3)
            ];
%           1      2      3      4     5     6     7      8      9      10    11    12    13    14    15    16     17     18     19    20    21    22    23   24   25    26    27
%           AlGaAs GaInAs AlInAs AlGaP GaInP AlInP AlGaSb GaInSb AlInSb AlAsP GaAsP InAsP AlPSb GaPSb InPSb AlAsSb GaAsSb InAsSb GaInN AlGaN AlInN GaAsN GaPN InPN InAsN GaSbN InSbN      
        end
        BF3 = repmat(real(bowingC(:,nTernary)),1,totalpoints) + imag(bowingC(:,nTernary))*x(:,1)';
        unknownVBOC = isnan(BF3(26,:));  
        BF3(26, unknownVBOC) = 0; % Bowing parameters denoted by NaN stand for unknown and are assumed 0

        [X1, fX1, tX1] = ZincBlendeIIIV(binaryTable{binaryIndex(nTernary,1)}, T); % recursive call
        [X2, fX2, tX2] = ZincBlendeIIIV(binaryTable{binaryIndex(nTernary,2)}, T); % recursive call
        materialParameter = X1*x(:,1)' + X2*x(:,2)' - BF3.*repmat((x(:,1).*x(:,2))',43,1);
        if nTernary == 2  % GaInAs
            materialParameter(20,:) = (1./materialParameter(41,:) + 1./materialParameter(38,:))/2;  % Gama1
            materialParameter(21,:) = (1./materialParameter(41,:) - 1./materialParameter(38,:))/4;  % Gama2
            materialParameter(22,:) = materialParameter(21,:) + (X1(22)-X1(21))*x(:,1)' + (X2(22)-X2(21))*x(:,2)' - 0.481*(x(:,1).*x(:,2))';
            materialParameter(39,:) = 2./(2*materialParameter(20,:)-materialParameter(21,:)-3*materialParameter(22,:)); % 39. mhh110 = 2/(2*gama1-gama2-3*gama3) 
            materialParameter(40,:) = 1./(materialParameter(20,:)-2*materialParameter(22,:));                           % 40. mhh111 = 1/(gama1-2*gama3)
            materialParameter(42,:) = 2./(2*materialParameter(20,:)+materialParameter(21,:)+3*materialParameter(22,:)); % 42. mlh110 = 2/(2*gama1+gama2+3*gama3) 
            materialParameter(43,:) = 1./(materialParameter(20,:)+2*materialParameter(22,:));                           % 43. mlh111 = 1/(gama1+2*gama3)
        end
        materialFormula = [alloy,', x = ', num2str(x(1,1))]; 
        parameterHistory = [[{[alloy,' BF']}; num2cell(BF3(:,1))],tX1(:,2),tX2(:,2)]; 
        if any(unknownVBOC) 
            parameterHistory(27,1) = {'Assumed 0'}; 
        end
    case '22' % quaternary alloys A(x)B(1-x)C(y)D(1-y)
        A = alloy(1:2);
        B = alloy(3:4);
        P_position = strfind(alloy,'P');
        if isempty(P_position)
            C = alloy(5:6); D = alloy(7:end);
        else
            C = alloy(5:min(6,P_position)); D = alloy(max(6,P_position):end);
        end
        
        if isempty(EqABCD)
            QABCD(1:43,1) = 0; 
        elseif length(EqABCD)==1 
            QABCD(1:43,1) = EqABCD; 
        elseif size(EqABCD,1)==43 && size(EqABCD,2)==1 
            QABCD = EqABCD; 
        else
            error('EqABCD for selecting quarternary parameter equations only allows [], a scalar, or a [43x1] array');
        end
        if isempty(BF4D)
            DIII(1:43,1) = 0; DV(1:43,1) = 0; DQ(1:43,1) = 0;
        elseif length(BF4D)==1
            DIII(1:43,1) = BF4D; DV(1:43,1) = BF4D; DQ(1:43,1) = BF4D; 
        elseif size(BF4D,1) == 1 && size(BF4D,2) == 2
            DIII(1:43,1) = BF4D(1); DV(1:43,1) = BF4D(2); DQ(1:43,1) = BF4D(1); 
        elseif size(BF4D,1) == 43 && size(BF4D,2) == 1
            DIII(1:43,1) = BF4D; DV(1:43,1) = BF4D; DQ(1:43,1) = BF4D; 
        elseif size(BF4D,1) == 43 && size(BF4D,2) == 2
            DIII(1:43,1) = BF4D(:,1); DV(1:43,1) = BF4D(1:43,2); DQ(1:43,1) = BF4D(1:43,1); 
        else 
            error('Surface bowing parameter BF4 for quarternary calculation only allows [], a scalar, a [1x2], [43x1] or [43x2] array');
        end
        idx0=(QABCD==0); idx1=(QABCD==1); idx2=(QABCD==2);
        if sum(idx0)+sum(idx1)+sum(idx2)~=43
            error('wrong assignment for BF4D'); 
        end
        
        [GABC,fABC,tABC,BF3ABC] = ZincBlendeIIIV([A B C],T,x(:,1));
        [GABD,fABD,tABD,BF3ABD] = ZincBlendeIIIV([A B D],T,x(:,1));
        [GACD,fACD,tACD,BF3ACD] = ZincBlendeIIIV([A C D],T,y(:,1));
        [GBCD,fBCD,tBCD,BF3BCD] = ZincBlendeIIIV([B C D],T,y(:,1));
        parameterHistory = [tABC(:,2),tABD(:,2),tACD(:,2),tBCD(:,2),tABC(:,3),tABD(:,3),tACD(:,3),tBCD(:,3),tABC(:,4:5),tABD(:,4:5)];
        BFrec = cell(43,1);
        BFrec(idx0) = cellstr([num2str(QABCD(idx0)),repmat(',-',sum(idx0),1)]); 
        BFrec(idx1) = cellstr([num2str(QABCD(idx1)),repmat(',',sum(idx1),1),num2str(DIII(idx1)),repmat(',',sum(idx1),1),num2str(DV(idx1))]);
        BFrec(idx2) = cellstr([num2str(QABCD(idx2)),repmat(',',sum(idx2),1),num2str(DQ(idx2))]);
        parameterHistory =  [[{'Eq, BF D'};cellstr(BFrec)],parameterHistory];
        materialFormula = [alloy,',', num2str(x(1,1)),',', num2str(y(1,1))];
        
        if any(idx0) || any(idx1) % Weighted sum of ternary values, Eq.1.(Littlejohn), Eq.5.2(Vurgaftman)
                                  % Wighted sum of ternary values with surface bowing, Eq.2,3(Mei)
            idx = idx0|idx1; 
            idy00 = y(:,2)==0; materialParameter(idx,idy00) = GABC(idx,idy00);
            idy01 = y(:,1)==0; materialParameter(idx,idy01) = GABD(idx,idy01);
            idy10 = x(:,2)==0; materialParameter(idx,idy10) = GACD(idx,idy10);
            idy11 = x(:,1)==0; materialParameter(idx,idy11) = GBCD(idx,idy11);
            BF3CD = zeros(43,totalpoints); BF3AB = BF3CD;
            BF = BF3ACD.*repmat(x(:,1)',43,1) + BF3BCD.*repmat(x(:,2)',43,1) - DV*(x(:,1).*x(:,2))';
            BF3CD(idx1,:) = BF(idx1,:);
            BF = BF3ABC.*repmat(y(:,1)',43,1) + BF3ABD.*repmat(y(:,2)',43,1) - DIII*(y(:,1).*y(:,2))'; 
            BF3AB(idx1,:) = BF(idx1,:);
            idy = ~(idy00|idy01|idy10|idy11); 
            x1 = repmat(x(idy,1)',sum(idx),1); x2 = 1-x1;
            y1 = repmat(y(idy,1)',sum(idx),1); y2 = 1-y1;
            materialParameter(idx,idy) = (x1.*x2.*(y1.*GABC(idx,idy)+y2.*GABD(idx,idy)-y1.*y2.*BF3CD(idx,idy)) + y1.*y2.*(x1.*GACD(idx,idy)+x2.*GBCD(idx,idy)-x1.*x2.*BF3AB(idx,idy)))./(x1.*x2 + y1.*y2); 
        end
        if any(idx2) % Biquadratic formula, Eq.7,8(Donati), Eq.A.8(Adachi), Eq.2.1.(Shikawa)
            nrep = sum(idx2);
            x1 = repmat(x(:,1)',nrep,1); x2 = 1-x1;
            y1 = repmat(y(:,1)',nrep,1); y2 = 1-y1;
            DQQ = repmat(DQ(idx2),1,totalpoints);
            materialParameter(idx2,:) = y1.*GABC(idx2,:) + y2.*GABD(idx2,:) ...
                - y1.*y2.*(BF3ACD(idx2,:).*repmat(x(:,1)',nrep,1)+BF3BCD(idx2,:).*repmat(x(:,2)',nrep,1)-x1.*x2.*DQQ);
        end        

    case {'13','31'} 
        switch alloyType
            case '13' % quaternary alloys DA(x)B(y)C(1-x-y)
                D = alloy(1:2); A = alloy(3:4); B = alloy(5); C = alloy(6:7);     
            case '31' % quaternary alloys A(x)B(y)C(1-x-y)D
                A = alloy(1:2); B = alloy(3:4); C = alloy(5:6); D = alloy(7:length(alloy));   
        end
        
        if isempty(EqABCD)
            QABCD(1:43,1) = 0; 
        elseif length(EqABCD)==1 
            QABCD(1:43,1) = EqABCD; 
        elseif size(EqABCD,1)==43 && size(EqABCD,2)==1 
            QABCD = EqABCD; 
        else
            error('EqABCD for selecting quarternary parameter equations only allows [], a scalar, or a [43x1] array');
        end
        if isempty(BF4D)
            DIII(1:43,1) = 0; DV(1:43,1) = 0; DQ(1:43,1) = 0;
        elseif length(BF4D)==1
            DIII(1:43,1) = BF4D; DV(1:43,1) = BF4D; DQ(1:43,1) = BF4D; 
        elseif size(BF4D,1) == 1 && size(BF4D,2) == 2
            DIII(1:43,1) = BF4D(1); DV(1:43,1) = BF4D(2); DQ(1:43,1) = BF4D(1); 
        elseif size(BF4D,1) == 43 && size(BF4D,2) == 1
            DIII(1:43,1) = BF4D; DV(1:43,1) = BF4D; DQ(1:43,1) = BF4D; 
        elseif size(BF4D,1) == 43 && size(BF4D,2) == 2
            DIII(1:43,1) = BF4D(:,1); DV(1:43,1) = BF4D(1:43,2); DQ(1:43,1) = BF4D(1:43,1); 
        else 
            error('Surface bowing parameter BF4 for quarternary calculation only allows [], a scalar, a [1x2], [43x1] or [43x2] array');
        end
        
        idx0=(QABCD==0); idx1=(QABCD==1); idx2=(QABCD==2);
        if sum(idx0)+sum(idx1)+sum(idx2)~=43
            error('wrong assignment for BF4D');  
        end

        [GABD,fABD,tABD,BF3ABD] = ZincBlendeIIIV([A B D],T,1-(1-x(:,1)+x(:,2))/2);
        [GBCD,fBCD,tBCD,BF3BCD] = ZincBlendeIIIV([B C D],T,1-(1-x(:,2)+x(:,3))/2);
        [GACD,fACD,tACD,BF3ACD] = ZincBlendeIIIV([A C D],T,1-(1-x(:,1)+x(:,3))/2);
        parameterHistory = [tABD(:,2),tBCD(:,2),tACD(:,2),tABD(:,3),tBCD(:,3),tACD(:,3),tABD(:,4:5),tBCD(:,5)];
        BFrec = cell(43,1);
        BFrec(idx0) = cellstr([num2str(QABCD(idx0)),repmat(',-',sum(idx0),1)]); 
        BFrec(idx1) = cellstr([num2str(QABCD(idx1)),repmat(',',sum(idx1),1),num2str(DIII(idx1)),repmat(',',sum(idx1),1),num2str(DV(idx1))]);
        BFrec(idx2) = cellstr([num2str(QABCD(idx2)),repmat(',',sum(idx2),1),num2str(DQ(idx2))]);
        parameterHistory =  [[{'Eq, BF D'};cellstr(BFrec)],parameterHistory];
        
        
       
        if any(idx0) || any(idx1) % Weighted sum of ternary values, Eq.5. C.K. Williams, J. Electron. Mat. 7, 639(1978); Eq.5.2.(Vurgaftman)
                                 % Wighted sum of ternary values with surface bowing, modified for the last (Mei), note D =0 used
            idx = idx0|idx1; 
            idy01 = x(:,1)==0; materialParameter(idx,idy01) = GBCD(idx,idy01);
            idy02 = x(:,2)==0; materialParameter(idx,idy02) = GACD(idx,idy02);
            idy03 = x(:,3)==0; materialParameter(idx,idy03) = GABD(idx,idy03);
            idy = ~(idy01|idy02|idy03); 
            BF4ABCD = zeros(43,totalpoints);
            nrep = sum(idx1);
            BF4ABCD(idx1,idy) = BF3ABD(idx1,idy).*repmat((1./(x(idy,1)+x(idy,2)))',nrep,1) ...
                              + BF3BCD(idx1,idy).*repmat((1./(x(idy,2)+x(idy,3)))',nrep,1) ...
                              + BF3ACD(idx1,idy).*repmat((1./(x(idy,1)+x(idy,3)))',nrep,1);
            nrep = sum(idx);
            x123 = repmat((x(idy,1).*x(idy,2).*x(idy,3))',nrep,1);
            x1 = repmat(x(idy,1)',nrep,1); x2 = repmat(x(idy,2)',nrep,1); x3 = repmat(x(idy,3)',nrep,1);
            materialParameter(idx,idy) = (x1.*x2.*GABD(idx,idy) + x2.*x3.*GBCD(idx,idy) + x1.*x3.*GACD(idx,idy))./(x1.*x2 + x2.*x3 + x1.*x3); ...
                                         - x123.*BF4ABCD(idx,idy);
        end
        if any(idx2) % Biquadratic formula, Eq.7,8(Donati)
%           A(x1)B(x2)C(x3)D ==> B(X){A(Y)C(1-Y)}(1-X)D ==> B2(X){B1(Y)B4(1-Y)}(1-X)D
%           1 - A, 2 - B, 3 - B, 4 - C
            X = x(:,2); Y = zeros(size(X));
            Y(X~=0) = x(X~=0,1)./(1-x(X~=0,2));
            [GB1,fB1,tB1] = ZincBlendeIIIV([A D],T);
            [GB2,fB2,tB2] = ZincBlendeIIIV([B D],T);
            [GB4,fB4,tB4] = ZincBlendeIIIV([C D],T);
            GB3 = GB2; C23 = zeros(43,totalpoints); C12 = -BF3ABD; C34 = -BF3BCD; C14 = -BF3ACD;
            nrep = sum(idx2);
            materialParameter(idx2,:) = GB1(idx2)*(Y.*(1-X))' + GB2(idx2)*(X.*Y)' + GB3(idx2)*((1-Y).*X)' + GB4(idx2)*((1-X).*(1-Y))' ...
                + C34(idx2,:).*repmat((X.*(1-X).*(1-Y))',nrep,1) + C12(idx2,:).*repmat((X.*(1-X).*Y)',nrep,1) ...
                + C14(idx2,:).*repmat(((1-X).*Y.*(1-Y))',nrep,1) + C23(idx2,:).*repmat((X.*Y.*(1-Y))',nrep,1) ...
                + DQ(idx2)*(X.*(1-X).*Y.*(1-Y))';
        end
        materialFormula = [alloy,',', num2str(x(1,1)),',', num2str(x(1,2))];    
     case 'lattice matched' % lattice matched cases
        bowingC = zeros(43,1);
        switch substrate
            case 'GaAs'
                switch alloy
                    case 'AlGaInP' % [Al(0.52)In(0.48)P](z)/[Ga(0.51)In(0.49)P](1-z)
                        [endA,fA,tA] = ZincBlendeIIIV('AlInP',T,0.52);
                        [endB,fB,tB] = ZincBlendeIIIV('GaInP',T,0.51);
                        bowingC(35) = 0.18; % eV, for EgG; linear interpolation for the rest
                        materialFormula = ['[Al(0.52)In(0.48)P](z)/[Ga(0.51)In(0.49)P](1-z), z=', num2str(z(1)), ', latticed matched to GaAs'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';
                        parameterHistory = [[{'AlGaInP BF'};num2cell(bowingC)],tA(:,2),tB(:,2),tA(:,3),tB(:,3),tA(:,4:5),tB(:,4)];
                    case 'GaInAsP' % [Ga(0.51)In(0.49)P](z)/[GaAs](1-z). Note: the bowing factor is not reliable!!
                        [endA,fA,tA] = ZincBlendeIIIV('GaInP',T,0.51);
                        [endB,fB,tB] = ZincBlendeIIIV('GaAs',T);
                        bowingC(35) = -0.62; % eV, for EgG; 
                        bowingC(36) = 0.53; % eV, for EgX; linear interpolation for the rest
                        materialFormula = ['[Ga(0.51)In(0.49)P](z)/[GaAs](1-z), z=', num2str(z(1)), ', latticed matched to GaAs'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';
                        parameterHistory = [[{'GaInAsP BF'};num2cell(bowingC)],tA(:,2:5),tB(:,2)];
                    case 'AlGaInAs'% not programmed yet
                        error([alloy, ' lattice matched to GaAs is not programmed yet']); 
                        materialParameter = [];
                        materialFormula = '';
                        return
                    otherwise
                        error([alloy, ' lattice matched to GaAs is not programmed']); 
                        materialParameter = [];
                        materialFormula = '';
                        return
                end
            case 'InP'
                switch alloy
                    case 'GaInAsP' % [Ga(0.47)In(0.53)As](z)/[InP](1-z)
%                         [endA,fA,tA] = ZincBlendeIIIV('InGaAs',T,0.5343+(0.5343-0.5378)*(T-300)/300); % 0K 0.5378 300K 0.5343
                        [endA,fA,tA] = ZincBlendeIIIV('InGaAs',T,0.53); 
                        endA = SetBandParameter(endA,'gama1',11.01,'gama2',4.18,'gama3',4.84); % To be consistent with Alavi, Ref437 in Ref[1]
                        [endB,fB,tB] = ZincBlendeIIIV('InP',T);
                        bowingC(35) = 0.13; % eV, for EgG; 
                        bowingC(12) = -0.06; % eV, for Dltaso; 
                        bowingC(27) = -6.7; % ev, for hydrostatic deformation potential ac; 
                        bowingC(41) = 0.032; % for mlh001; 
                        materialFormula = ['[Ga(0.47)In(0.53)As](z)/[InP](1-z), z=', num2str(z(1)), ', latticed matched to InP'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';
                        materialParameter(23,:) = 1./(materialParameter(20,:)-materialParameter(24,:).*materialParameter(12,:)...
                                                /3./materialParameter(35,:)./(materialParameter(35,:)+materialParameter(12,:))); % mso
                        parameterHistory = [[{'GaInAsP BF'};num2cell(bowingC)],tA(:,2:5),tB(:,2)];
                    case 'AlGaInAs' % [Al(0.48)In(0.52)As](z)/[Ga(0.47)In(0.53)As](1-z)
                        [endA,fA,tA] = ZincBlendeIIIV('AlInAs',T,0.48);
                        [endB,fB,tB] = ZincBlendeIIIV('GaInAs',T,0.4657);
                        bowingC(35) = 0.22; % eV, for EgG; 
                        bowingC(13) = -0.016; % m0, for me; 
                        bowingC(24) = -5.68; % ev, for Ep; 
                        materialFormula = ['[Al(0.48)In(0.52)As](z)/[Ga(0.47)In(0.53)As](1-z), z=', num2str(z(1)), ', latticed matched to InP'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';
                        parameterHistory = [[{'AlGaInAs BF'};num2cell(bowingC)],tA(:,2),tB(:,2),tA(:,3),tB(:,3),tA(:,4:5),tB(:,4)];
                    case 'GaInAsSb' % [Ga(0.47)In(0.53)As](z)/[GaAs(0.5)Sb(0.5)](1-z)
                        [endA,fA,tA] = ZincBlendeIIIV('GaInAs',T,0.47);
                        [endB,fB,tB] = ZincBlendeIIIV('GaAsSb',T,0.5);
                        bowingC(35) = 0.22; % eV, for EgG; 
                        materialFormula = ['[Ga(0.47)In(0.53)As](z)/[GaAs(0.5)Sb(0.5)](1-z), z=', num2str(z(1)), ', latticed matched to InP'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';
                        parameterHistory = [[{'GaInAsSb BF'};num2cell(bowingC)],tA(:,2),tB(:,2),tA(:,3),tB(:,3),tA(:,4:5),tB(:,4)];
                    otherwise
                        error([alloy, ' lattice matched to ', substrate, ' is not programmed']); 
                        materialParameter = [];
                        materialFormula = '';
                        return
                end
            case 'InAs'
                switch alloy
                    case 'GaInAsSb' % [GaAs(0.08)Sb(0.92)](z)/[InAs](1-z)
                        [endA,fA,tA] = ZincBlendeIIIV('GaAsSb',T,0.08);
                        [endB,fB,tB] = ZincBlendeIIIV('InAs',T);
                        bowingC(35) = 0.6; % eV, for EgG; 
                        bowingC(36) = 0.15; % eV, for EgX; 
                        bowingC(37) = 0.6; % eV, for EgL; 
                        materialFormula = ['[GaAs(0.08)Sb(0.92)](z)/[InAs](1-z), z=', num2str(z(1)), ', latticed matched to InAs'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';
                        parameterHistory = [[{'GaInAsSb BF'};num2cell(bowingC)],tA(:,2:5),tB(:,2)];
                    case 'AlGaAsSb' % [AlAs(0.16)Sb(0.84](z)/[GaAs(0.08)Sb(0.92)](1-z)
                        [endA,fA,tA] = ZincBlendeIIIV('AlAsSb',T,0.16);
                        [endB,fB,tB] = ZincBlendeIIIV('GaAsSb',T,0.08);
                        bowingC(35) = 0.48; % eV, for EgG; use parameters of AlGaAsSb on GaSb
                        bowingC(36) = 1.454; % eV, for EgX; need verification
                        bowingC(37) = 0.807; % eV, for EgL; need verification
                        materialFormula = ['[AlAs(0.16)Sb(0.84](z)/[GaAs(0.08)Sb(0.92)](1-z), z=', num2str(z(1)), ', latticed matched to InAs. Bowing factors for indirect Eg need verification'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';                        
                        parameterHistory = [[{'AlGaAsSb BF'};num2cell(bowingC)],tA(:,2),tB(:,2),tA(:,3),tB(:,3),tA(:,4:5),tB(:,4:5)];
                    case 'InAsPSb' % [InSb(0.31)P(0.69)](z)/[InAs](1-z)
                        [endA,fA,tA] = ZincBlendeIIIV('InPSb',T,0.69);
                        [endB,fB,tB] = ZincBlendeIIIV('InAs',T);
                        bowingC(12) = -0.75; % eV, for Dltaso; 
                        materialFormula = ['[InSb(0.31)P(0.69)](z)/[InAs](1-z), z=', num2str(z(1)), ', latticed matched to InAs'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';
                        parameterHistory = [[{'InAsPSb BF'};num2cell(bowingC)],tA(:,2:5),tB(:,2)];
                    otherwise
                        error([alloy, ' lattice matched to ', substrate, ' is not programmed']); 
                        materialParameter = [];
                        materialFormula = '';
                        return
                end
            case 'GaSb'
                switch alloy
                    case 'GaInAsSb' % [InAs(0.91)Sb(0.09)](z)/[GaSb](1-z)
                        [endA,fA,tA] = ZincBlendeIIIV('InAsSb',T,0.91);
                        [endB,fB,tB] = ZincBlendeIIIV('GaSb',T);
                        bowingC(35) = 0.75; % eV, for EgG; 
                        bowingC(36) = 0.43; % eV, for EgX; 
                        bowingC(37) = 0.85; % eV, for EgL; 
                        bowingC(12) = -0.26; % eV, for Dltaso; 
                        materialFormula = ['[InAs(0.91)Sb(0.09)](z)/[GaSb](1-z), z=', num2str(z(1)), ', latticed matched to GaSb'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';
                        parameterHistory = [[{'GaInAsSb BF'};num2cell(bowingC)],tA(:,2:5),tB(:,2)];
                    case 'AlGaAsSb' % [AlAs(0.08)Sb(0.92](z)/[GaSb](1-z)
                        [endA,fA,tA] = ZincBlendeIIIV('AlAsSb',T,0.08);
                        [endB,fB,tB] = ZincBlendeIIIV('GaSb',T);
                        bowingC(35) = 0.48; % eV, for EgG; 
                        bowingC(36) = 1.454; % eV, for EgX; need verification
                        bowingC(37) = 0.807; % eV, for EgL; need verification
                        materialFormula = ['[AlAs(0.08)Sb(0.92](z)/[GaSb](1-z), z=', num2str(z(1)), ', latticed matched to GaSb. Bowing factors for indirect Eg need verification'];
                        materialParameter = endA*z' + endB*(1-z)' - bowingC*((1-z).*z)';                        
                        parameterHistory = [[{'AlGaAsSb BF'};num2cell(bowingC)],tA(:,2:5),tB(:,2)];
                    otherwise
                        error([alloy, ' lattice matched to ', substrate, ' is not programmed']); 
                        materialParameter = [];
                        materialFormula = '';
                        return
                end
            case 'GaAsP' % GaAs(0.61)P(0.39)
                switch alloy 
                    case 'AlGaAsP' % [AlAs(0.61)P(0.39)](z)/[GaAs(0.61)P(0.39)](1-z)
                        [endA,fA,tA] = ZincBlendeIIIV('AlAsP',T,0.61);
                        [endB,fB,tB] = ZincBlendeIIIV('GaAsP',T,0.61);
                        bowingC(35) = 1.3; % eV, for EgG; need verification 
                        materialFormula = ['[AlAs(0.61)P(0.39)](z)/[GaAs(0.61)P(0.39)](1-z), z=', num2str(z), ', latticed matched to GaSb. Bowing factor for direct Eg needs verification'];
                        materialParameter = endA*z + endB*(1-z) - bowingC*(1-z)*z;
                        parameterHistory = [[{'AlGaAsP BF'};num2cell(bowingC)],tA(:,2),tB(:,2),tA(:,3),tB(:,3),tA(:,4:5),tB(:,4:5)];
                    otherwise
                        error([alloy, ' lattice matched to ', substrate, ' is not programmed']); 
                        materialParameter = [];
                        materialFormula = '';
                        return
                end
            otherwise
                error([alloy,' - the substrate ',substrate,' is not recognized or is not programmed!']);
                materialParameter = [];
                materialFormula = '';
                return
        end
        tVBObf = '';
        if length(tA) == 44, tVBObf = [tVBObf tA{44,2}, ' ']; end
        if length(tB) == 44, tVBObf = [tVBObf tB{44,2}]; end
        if ~isempty(tVBObf), parameterTable(44,1:2) = {'VBO B.F. assumed 0', tVBObf}; end 
end
parameterTable(1:44,1) ={['  T=',num2str(T),' K  '],'01. a1c1(A)','02. a2c1(1e-5 A/K)','03. Eg0G(eV)','04. alfaG(meV/K)','05. betaG(K)','06. Eg0X(eV)','07. alfaX(meV/K)',...
        '08. betaX(K)','09. Eg0L(eV)','10. alfaL(meV/K)','11. betaL(K)','12. Dltaso','13. meG','14. mlL','15. mtL','16. mDOSL',...
        '17. mlX','18. mtX','19. mDOSX','20. gama1','21. gama2','22. gama3','23. mso','24. Ep(eV)','25. F','26. VBO(eV)',...
        '27. ac(eV)','28. av(eV)','29. b(eV)','30. d(eV)','31. c11(GPa)','32. c12(GPa)','33. c44(GPa)','34. alc(A)','35. Eg_Gama',...
        '36. Eg_X','37. Eg_L','38. mhh[100]','39. mhh[110]','40. mhh[111]','41. mlh[100]','42. mlh[110]','43. mlh[111]'};
parameterTable(1,2) = {materialFormula};
parameterTable(2:44,2) = num2cell(materialParameter(:,1));
parameterTable=[parameterTable,parameterHistory];

% Considering dilute nitrides, version 4, T. Mei, July 2006
recursion_level = recursion_level - 1;
if recursion_level == 0
    dilutenitrideTable = {'GaAsN','InAsN','GaPN','InPN','InSbN','GaInAsN','GaInPN'};
    ndnitride = strmatch(alloy, dilutenitrideTable, 'exact');
    if elementVcount ~= 1 && ~isempty(ndnitride)
        switch ndnitride
            case 6 % GaInAsN
                BF3(:,1) = 1.65*x(:,1)+1.44*x(:,2)-0.38*x(:,1).*x(:,2);
                BF3(:,2) = 2.7*x(:,1) +2.0*x(:,2) -3.5*x(:,1).*x(:,2);
            case 7 % GaInPN
                BF3(:,1) = 2.18*x(:,1)+1.79*x(:,2);
                BF3(:,2) = 3.05*x(:,1)+3.0*x(:,2) -3.3*x(:,1).*x(:,2);
            otherwise
                dilutenitrideParameter =[...
                    1.65,   1.44,   2.18,  1.79,  0.65;...  % 01. EN w.r.t. to VBM (eV)
                    2.7,    2.0,    3.05,  3.0,   3.0;]';   % 02. V (eV)
                %   1 GaAsN 2 InAsN 3 GaPN 4 InPN 5 InSbN
                BF3 = dilutenitrideParameter(ndnitride,:);
        end
        materialFormula =[materialFormula,', dilute nitrides: BAC model parameters [EN,V] stored in BF3'];
    end
end
return
