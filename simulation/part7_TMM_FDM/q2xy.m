function [x,y]=q2xy(Q)

global h_bar;
global speed_of_light;

for y = (0:0.01:1)
    x = 1-0.47*y;
    if(Eg(x,y)<=(h_bar*2*pi*speed_of_light/(Q*1e-6)))
        break;
    end
end