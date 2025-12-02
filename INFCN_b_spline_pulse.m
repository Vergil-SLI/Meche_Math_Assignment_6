%% B-spline pulse function

% INPUTS:
%   t: current time
%   w: width of pulse (starts at t=0, ends at t=h)
%   h: height of pulse
% OUTPUTS:
%   res: pulse evaluated at t

function res = INFCN_b_spline_pulse(t, w, h)
    
    t = 4*t/w;
    b3 = (0<=t).*(t<1).*(t.^3)/4;
    t = t-1;
    b2 = (0<=t).*(t<1).*(-3*t.^3 + 3*t.^2 + 3*t + 1)/4;
    t = t-1;
    b1 = (0<=t).*(t<1).*(3*t.^3 - 6*t.^2 + 4)/4;
    t = t-1;
    b0 = (0<=t).*(t<1).*(-t.^3 + 3*t.^2 - 3*t + 1)/4;
    
    res = h*(b0 + b1 + b2 + b3);
end