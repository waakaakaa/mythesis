function [layers, expandedformul] = InterpFormula(formul)
% function [layers, expandedformula] = InterpFormula(formula)
% This function translates a heterostructure formula expression into a parameter array to help easy programming
%
% This software is written by Dr. Ting Mei, Nanyang Techonogical University, Singapore, March 2005.
% Please get the permission before using this program. The author preserves all rights.
% The author is grateful if being acknowledged when you use the program for publication  
% If you have any question, feedback, or comment, please contact Dr. Mei at:
%     Tel: +65 67904387; Email: etmei@ntu.edu.sg
%
% 
% Syntax of formula (string):
%   1. Express a layer
%       Binary: (material_symbol,layer_thickness) 
%       Ternary: (material_symbol,layer_thickness,x) 
%       Qaurternary: (material_symbol,layer_thickness,x,y) 
%       Lattice Matched: (material_symbol,layer_thickness,x,latticematched)
%       Refer to the corresponding parts in ZincBlendeIIIV.m for definitions of x,y 
%   2. Nest segments
%       n((layer1)(layer2))
%   3. material_symbol: GaAs, AlGaAs, InGaAs, etc. It is case sensitive
%
%
%  Examples: [layers, expandedformul] = InterpFormula('2((AlGaAs,40,0.45)(InGaAsP,30,0.25,0.3)(GaAs,48)(InAlAsP,40,0.45,latticematched))')


% remove blanks
formul = strrep(formul, ' ', '');

% expand the formula
while 1
    bra = findstr('(',formul);
    ket = findstr(')',formul);
    b = size(bra,2);
    k = size(ket,2);
    if b ~= k
        error('Symbolc formula error (3): Brackets mismatch');
    end
    if isempty(bra)
        break;
    end
    kx = 1;
    while ket(kx) < bra(b) 
        kx = kx+1;
        if kx > k
            error('Symbolc formula error (3): Brackets mismatch');
        end
    end
    braketed_str = formul(bra(b)+1:ket(kx)-1);
    repeat_num_str = '';
    for r = bra(b)-1:-1:1 
        if formul(r)<'0' || formul(r)>'9'
            break;
        end
        repeat_num_str = [formul(r),repeat_num_str];
    end
    if isempty(repeat_num_str)
        repeat_num = 1;
    else
        repeat_num = str2double(repeat_num_str);
    end
    replacement = '';
    for r = 1:repeat_num
        replacement = ['{', braketed_str, '}', replacement];
    end
%    replacement
    formul = strrep(formul, [repeat_num_str,'(',braketed_str,')'], replacement);
%     formul = strrep(formul, [repeat_num_str,'(',braketed_str,')'], ['{',replacement,'}']);
    formul = strrep(formul,'{}','');
    formul = strrep(formul, '{{','{');
    formul = strrep(formul, '}}','}');
end
% formul

% form layer structure
bra = findstr('{',formul);
ket = findstr('}',formul);
for kx = 1:length(bra)
    braketed_str = formul(bra(kx)+1:ket(kx)-1);
    comma = findstr(',',braketed_str);
    layers{kx,1} = braketed_str(1:comma(1)-1);
    cn = length(comma);
    for cx = 2:cn
        layers{kx,cx} = braketed_str(comma(cx-1)+1:comma(cx)-1);
        if isempty(strmatch('latticematched',layers{kx,cx}))
            layers{kx,cx} = str2double(layers{kx,cx});
        end
    end
    layers{kx,cn+1} = braketed_str(comma(cn)+1:length(braketed_str));
    if isempty(strmatch('latticematched',layers{kx,cn+1}))
        layers{kx,cn+1} = str2double(layers{kx,cn+1});
    end
end
expandedformul = strrep(formul,'{','(');
expandedformul = strrep(expandedformul,'}',')');
return