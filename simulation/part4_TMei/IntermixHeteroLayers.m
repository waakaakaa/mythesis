function [XY,d_step] = IntermixHeteroLayers(layerD,layerXY,Ldiff,Nnodes,Pmode)
% function [XY,d_step] = IntermixHeteroLayers(layerD,layerXY,Ldiff,Nnodes,Pmode)
% This program calculate composition distribution after interdiffusion for a heterostructure
%
% This software is written by Dr. Ting Mei, Nanyang Techonogical University, Singapore, March 2005.
% This software is free to use, but the author preserves all rights.
% The author is grateful if being acknowledged when you use the program for publication  
% If you have any question, feedback, or comment, please contact Dr. Mei at:
%     Tel: +65 67904387; Email: etmei@ntu.edu.sg
%
%
% Input variables
% layerD: layer thicknesses of a heterostructure, 1D array 
% layerXY: layer compositions X([X;Y]) of the ternary(quaternary) heterostructure, 1D(2D) array 
% Ldiff: diffusion length(s) in sublattice(s), a scalar or a two-element array - Ldiff or [Ldiff1, Lidiff2]
% Nnodes: total number of meshpoints for the calculation
% Pmode: 
%   default: treat as infite on both sides of the heterostructure
%   1: treat as a periodic structure
%
% Output variables
% XY: composition distribution X or [X,Y] on meshpoints after interdiffusion
% d_step: the interval of the meshpoints
%
% IntermixHeteroLayers.m and MeshHeteroLayers.m both accept input of a heterostructure 


[nlayers, m] = size(layerXY);
if nlayers ~= length(layerD) || m ~= length(Ldiff)
    error('Input arrays error: the length or dimension does not match');
end

if size(layerD,2)<length(layerD)
    layerD=layerD';
end
z_int = cumsum(layerD);
Lperiod = z_int(end);
d_step = Lperiod/Nnodes;
z = d_step*(0:Nnodes-1)';
if nargin == 5 && Pmode ==1 % periodic structure
    z_int_1=[0,z_int(1:end-1)];
    k=0;
    Xlast=zeros(Nnodes,1);
    Ylast=zeros(Nnodes,1);
    while 1
        for ilayer=1:nlayers
            X(:,ilayer) = layerXY(ilayer,1)/2*(ferf(z-z_int_1(ilayer)-k*Lperiod,Ldiff(1)) ...
                -ferf(z-z_int(ilayer)-k*Lperiod,Ldiff(1)));
            if k ~=0
                X(:,ilayer) = X(:,ilayer)+ layerXY(ilayer,1)/2*(ferf(z-z_int_1(ilayer)+k*Lperiod,Ldiff(1)) ...
                    -ferf(z-z_int(ilayer)+k*Lperiod,Ldiff(1)));
            end
            if m > 1 
                Y(:,ilayer) = layerXY(ilayer,2)/2*(ferf(z-z_int_1(ilayer)-k*Lperiod,Ldiff(2)) ...
                    -ferf(z-z_int(ilayer)-k*Lperiod,Ldiff(2)));
                if k ~=0
                    Y(:,ilayer) = Y(:,ilayer)+ layerXY(ilayer,2)/2*(ferf(z-z_int_1(ilayer)+k*Lperiod,Ldiff(2)) ...
                        -ferf(z-z_int(ilayer)+k*Lperiod,Ldiff(2)));
                end
            end
        end
        
        if m > 1
            Xs = sum(X,2); Ys = sum(Y,2);
            Xlast=Xlast+Xs;Ylast=Ylast+Ys; 
            if max(abs(Xs))<1e-4 && max(abs(Ys))<1e-4
                XY=[Xlast,Ylast];
                break
            end
        else   
            Xs = sum(X,2); 
            Xlast=Xlast+Xs;
            if max(abs(Xs))<1e-4 
                XY = Xlast;
                break
            end
        end
        k=k+1;
    end
else % infinite struture
    for ilayer=1:nlayers-1
        XRightSideErf(:,ilayer) = layerXY(ilayer,1)*(1-ferf(z-z_int(ilayer),Ldiff(1)))/2;
        XLeftSideErf(:,ilayer) = layerXY(ilayer+1,1)*(1+ferf(z-z_int(ilayer),Ldiff(1)))/2;
        if m > 1
            YRightSideErf(:,ilayer) = layerXY(ilayer,2)*(1-ferf(z-z_int(ilayer),Ldiff(2)))/2;
            YLeftSideErf(:,ilayer) = layerXY(ilayer+1,2)*(1+ferf(z-z_int(ilayer),Ldiff(2)))/2;
        end
    end
    XY = sum(XRightSideErf,2)+sum(XLeftSideErf,2)-sum(layerXY(2:nlayers-1,1));
    if m > 1
        Y = sum(YRightSideErf,2)+sum(YLeftSideErf,2)-sum(layerXY(2:nlayers-1,2));
        XY = [XY,Y];
    end
end
return

% 
% x(i)
% i =      1        2         3          . . .       nlayer
%               +------+               +------+
%               |      |               |      |
%               |      |               |      |
%        -------+      +---------------+      +-----------
% 
% z_int(i)  
% i=            1      2              . . .  nlayer-1
%                     (i)
%                     
%        nlayer-1   x(i)            z-z_int(i)       x(i+1)          z-z_int(i)
% y(z) =   sum    -------(1+ erf(- -----------)) + --------(1 + erf(-----------))
%          i=1       2               2*Ldiff           2              2*Ldiff
%                           ---.                                .---      
%                               \                              /
%                                ---                        ---
%                          right edge                      left edge
%     - sum( x(2:nlayer-1))
%
% This is dervived from
%
%         nlayer   x(i)           z-z_int(i-1)        x(i)             z-z_int(i)
% y(z) =   sum   -------(1+ erf( --------------)) + -------(1 + erf(- -----------)) - x(i)
%          i=1      2               2*Ldiff            2                2*Ldiff
%                           .---                                ---.      
%                          /                                        \
%                      ___/                                          \____
%                      left edge                                right edge
%                _____________________________       _____________________________   ______
%                            ^                                   ^                    ^
%                            |                                   |                    |
%                 absent for the left end            absent for the right end   absent for both ends        
%                 
%                
% but if handing a peroidic structure, i.e.  x(i) = x(i+nlayer) and y(z) = y(z+Lp)   
%        k=+inf  nlayer    x(i)      z-z_int(i-1)-kLp          z-z_int(i)-kLp
% y(z) =  sum     sum   -------(erf(------------------) - erf(----------------)) 
%        k=-inf   i=1       2             2*Ldiff                 2*Ldiff
%                           .---                                ---.      
%                          /                                        \
%                      ___/                                          \____
%                      left edge                                right edge
%

function Y=ferf(Z,Ld)
% step function => error function with difussion length Ldiff; allows Ld=0
if Ld == 0
    Y = double(Z>0)-double(Z<0);
else
    Y = erf(Z/Ld/2);
end
return
