function ladderplot(varargin)

qcolor = 'brgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcybrgkmcy';

k = ishold;

fault = 0;
icurve = 1;
iargin = 1;
holdplot = 0;
while iargin <= nargin
    if isnumeric(varargin{iargin}) && isnumeric(varargin{iargin+1})
        d = varargin{iargin};
        Y = varargin{iargin+1};
        npoints = length(Y);
        if length(d) ~= npoints
            fault = 1; break
        end
        X = []; X(2:npoints+1) = cumsum(d);
        Z = Y;  Z(end+1) = Y(end);
        if iargin+2<=nargin && ischar(varargin{iargin+2})
            linestr = varargin{iargin+2};
            stairs(X,Z,linestr);
            if ~holdplot, hold on, holdplot = 1; end
            iargin = iargin + 3;
        else
            stairs(X,Z,qcolor(icurve));
            if ~holdplot, hold on, holdplot = 1; end
            iargin = iargin + 2;
            icurve = icurve + 1;
        end
    else
        fault = 1;
        break
    end
end

if ~k, hold off, end
if fault == 1, error('wrong input in ladderplot'), end

return
