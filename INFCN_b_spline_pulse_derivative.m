%% B-spline pulse function (derivative)

% INPUTS:
%   t: current time
%   w: width of pulse (starts at t=0, ends at t=h)
%   h: height of pulse
% OUTPUTS:
%   res: derivative of pulse evaluated at t

function res = INFCN_b_spline_pulse_derivative(t, w, h)
    
    t = 4*t/w;
    b3 = (0<=t).*(t<1).*(3*t.^2)/4;
    t = t-1;
    b2 = (0<=t).*(t<1).*(-9*t.^2+6*t+3)/4;
    t = t-1;
    b1 = (0<=t).*(t<1).*(9*t.^2-12*t)/4;
    t = t-1;
    b0 = (0<=t).*(t<1).*(-3*t.^2+6*t-3)/4;
    
    res = (4*h/w)*(b0 + b1 + b2 + b3);
end