%INPUTS

%t: current time
%V: system state. V = [U;dUdt] where
% U and dUdt are n x 1 column vectors
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
function dVdt = string_rate_func01(t,V,string_params)
n = string_params.n; %number of masses
M = string_params.M; %total mass attached to the string
Uf_func = string_params.Uf_func; %function describing motion of end point
dUfdt_func = string_params.dUfdt_func; %time derivative of Uf
Tf = string_params.Tf; %tension in string
L = string_params.L; %length of string
c = string_params.c; %damping coefficient
dx = string_params.dx; %horizontal spacing between masses
%unpack state variable
U = V(1:n);
dUdt = V((n+1):(2*n));
Uf = Uf_func(t);
dUfdt = dUfdt_func(t);
%compute acceleration
Twos = -2.*eye(n,n);
Ones1 = eye(n-1,n-1);
Ones2 = eye(n-1,n-1);
col1 = zeros(n,1);
col2 = zeros(n,1);
row1 = zeros(1,n-1);
row2 = zeros(1,n-1);
Ones1 = [Ones1;row1];
Ones1 = [col1, Ones1];
Ones2 = [row2;Ones2];
Ones2 = [Ones2,col2];
K = Twos+Ones1+Ones2
K
d2Udt2 = (M\K)*U;
%assemble state derivative
dVdt = [dUdt;d2Udt2];
end