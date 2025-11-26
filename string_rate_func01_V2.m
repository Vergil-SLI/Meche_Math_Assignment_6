

% INPUTS
%   t: current time
%   V: system state. V = [U; dUdt] where U and dUdt are nx1 column vectors
%   string_params: a struct containing the system parameters describing the string
%       string_params.n: number of masses
%       string_params.M: total mass attached to the string
%       string_params.Uf_func: function describing motion of end point
%       string_params.dUfdt_func: time derivative of Uf
%       string_params.Tf: tension in string
%       string_params.L: length of string
%       string_params.c: damping coefficient
%       string_params.dx: horizontal spacing between masses
% OUTPUTS
%   dVdt: time derivative of the state vector V at time t, V = [dUdt; d2Udt2]
%           has size nx1 column vector

function dVdt = string_rate_func01_V2(t, V, string_params)
    
    % unpack string_params
    n = string_params.n; % number of masses
    M_total = string_params.M; % total mass attached to the string
    Uf_func = string_params.Uf_func; % function describing motion of end point
    dUfdt_func = string_params.dUfdt_func; % time derivative of Uf
    Tf = string_params.Tf; % tension in string
    c = string_params.c; % damping coefficient
    dx = string_params.dx; % horizontal spacing between masses
    
    u0 = 0; % by sys definition
    
    % unpack state variable
    U = V(1:n);
    dUdt = V((n+1):(2*n));
    
    % compute acceleration
    % M_total/n * u'' = Tf/dx * (u_i-1 - 2*u_i + u_i+1)     % normal sys
    %                  + c/dx * (u_i-1' - 2*u_i' + u_i+1')  % + damping
    % [M]*[u''] = Tf/dx*[K]*[u] + Tf/dx*[B_u]               % matrix vers
    %           + c/dx*[K]*[u'] + c/dx*[B_u']
    
    % construct M
    mi = M_total/n; I = eye(n); 
    M = mi*I;
    
    % construct K
    K = -2*I + circshift(I, [1, 0]) + circshift(I, [0, 1]);
    K(end, 1) = 0; K(1, end) = 0;

    % construct B (end conditions)
    B_u = zeros(n, 1); B_dudt = zeros(n, 1);
    B_u(1) = u0; B_u(end) = Uf_func(t);
    B_dudt(1) = u0; B_dudt(end) = dUfdt_func(t);

    % calculate 
    % d2Udt2 = (M\K)*U;
    partial_calc = (Tf/dx)*K*U + (Tf/dx)*B_u + (c/dx)*K*dUdt + (c/dx)*B_dudt;
    d2Udt2 = M\partial_calc;
    
    % assemble state derivative
    dVdt = [dUdt; d2Udt2];
end
