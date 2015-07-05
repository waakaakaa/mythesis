function del = delta(x,y)

del = 0.34*(1-x).*y + 0.43*x.*y + 0.10*(1-x).*(1-y) + 0.10*x.*(1-y);
del = del * 1.6e-19;