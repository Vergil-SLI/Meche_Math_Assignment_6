function day_2025_12_9_coding
 %exp1()
 exp2()
end
function exp1
n = 5;
mode_index = 3;

string_params.n = 5;
string_params.M = 10;
string_params.Tf = 2; %tension in string
string_params.L = 7;%length of string
string_params.c = 0.0001;  %damping coefficient
string_params.dx = string_params.L/(n+1);%horizontal spacing between masses

[M_mat, K_mat] = construct_2nd_order_matrices(string_params);

[Ur_mat,lambda_mat] = eig(K_mat,M_mat);
mode_shape_LA = Ur_mat(:,mode_index);
omega_n = sqrt(-lambda_mat(mode_index,mode_index));

omega = omega_n;
A = 3;

string_params.Uf_func = @(t_in) A*sin(omega*t_in);
string_params.dUfdt_func = @(t_in) -omega*A*cos(omega*t_in);
string_params.Uf_amplitude = A;

U0 = zeros(n,1);
dUdt0 = zeros(n,1);

V0 = [U0;dUdt0];

t_list = linspace(0,20*(2*pi)/omega,1000+1);

my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);

[t_list,vresult] = ode45(my_rate_func,t_list,V0);

animate_string_with_mode(t_list,vresult,string_params)
end

function exp2
n = 200;


string_params.n = 5;
string_params.M = 10;
string_params.Tf = 2; %tension in string
string_params.L = 7;%length of string
string_params.c = 0.0001;  %damping coefficient
string_params.dx = string_params.L/(n+1);%horizontal spacing between masses


string_params.Uf_func = @(t_in) 0;
string_params.dUfdt_func = @(t_in) 0;
string_params.Uf_amplitude = 0;

rho = string_params.M/string_params.L;
c = sqrt(string_params.Tf/rho);

x_list = linspace(0,string_params.L,n+2)';
x_list = (x_list(2:end-1));

w_pulse = string_params.L/20;
h_pulse = 7;

U0 = INFCN_triangle_pulse(x_list,w_pulse,h_pulse);
dUdt0 = -c*INFCN_triangle_pulse_derivative(x_list,w_pulse,h_pulse);

hold on
plot(x_list,U0,'o-','Color','r');
plot(x_list,dUdt0,'o-','Color','b');

V0 = [U0;dUdt0];

t_list = linspace(0,100,5000+1);

my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);

[t_list,vresult] = ode45(my_rate_func,t_list,V0);

animate_string_travelling_wave(t_list,vresult,string_params,x_centroid0,c)
end

function animate_string_with_mode(t_list,vresult,string_params,mode_shape_LA)
n = string_params.n;
L = string_params.L;

X_list = linspace(0,L,n+2);

string_plot = plot(0,0,'o-','Color','k','LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4);

mode_shape_plot = plot(0,0,'o-','Color','b','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4);

mode_shape_padded = [0;mode_shape_LA;0];
scalefactor = max(abs(maxU))
set(mode_shape_plot,'xdata',X_list,'ydata',mode_shape_padded)

maxU = max(max(vresult(:,1:n)));
minU = min(min(vresult(:,1:n)));

h = max(abs(maxU),abs(minU));

axis([0,L,-1.1*h,1.1*h])
hold on
xlabel('x')
ylabel('U(x,t)')
title("Plot of Vibrating String")

for k = 1:length(t_list)
    t = t_list(k);
    U = vresult(k,1:n);
    Uf = string_params.Uf_func(t);

    U_padded = [0,U,Uf];
    X_list = linspace(0,L,n+2);

    set(string_plot,'xdata',X_list,'ydata',U_padded)
    drawnow;
end
end

function animate_string_travelling_wave(t_list,vresult,string_params,x_centroid0,c)
n = string_params.n;
L = string_params.L;

X_list = linspace(0,L,n+2);

string_plot = plot(0,0,'o-','Color','k','LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4);

centroid_plot = plot(0,0,'b','LineWidth',2);

maxU = max(max(vresult(:,1:n)));
minU = min(min(vresult(:,1:n)));

h = max(abs(maxU),abs(minU));

axis([0,L,-1.1*h,1.1*h])
hold on
xlabel('x')
ylabel('U(x,t)')
title("Plot of Vibrating String")

for k = 1:length(t_list)

    %short blurb showing how to find x-coord of tracking line
    %x = x-coord of tracking line, t = current time
    %c = wave speed, w = pulse width (in time), L = string length
    x_centroid = c*t+x_centroid0;
    x_centroid = mod(x_centroid,2*L);
    if x_centroid > L
        x_centroid = 2*L - x_centroid;
    end

    t = t_list(k);
    U = vresult(k,1:n);
    Uf = string_params.Uf_func(t);

    U_padded = [0,U,Uf];
    X_list = linspace(0,L,n+2);

    set(string_plot,'xdata',X_list,'ydata',U_padded)
    set(centroid_plot,'xdata',[x_centroid,x_centroid],'ydata',[0,L,-1.1*h,1.1*h])
    drawnow;
end

end

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
    dUdt = dUfdt_func(t);
    
    u_left = [0;(U(1:end-1))];
    u_right = [U(2:end);Uf];
    
    term1 = (Tf/dx)*(u_left -2*U+u_right);
    term2 = (c/dx)*(dUdt_left -2*dUdt+dUdt_right);
    %compute acceleration
    d2Udt2 = term1+term2;
    %assemble state derivative
    dVdt = [dUdt;d2Udt2];
end

%build the mass and stiffness matrices that describe the 2nd order system.
%INPUTS
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
%OUTPUTS
%M_mat: the n x n mass (inertia) matrix
%K_mat: the n x n stiffness matrix
function [M_mat,K_mat] = construct_2nd_order_matrices(string_params)

n = string_params.n;
Tf = string_params.Tf;
M = string_params.M;
dx = string_params.dx;

I_n = eye(n);
M_left = circshift(I_n,[0,-1]);
M_right = circshift(I_n,[0,1]);
my_Laplacian = M_left-2*I_n+M_right;
my_Laplacian(1,end) = my_Laplacian(1,end)-1; %delete unwanted 1 in top right corner
my_Laplacian(end,1) = my_Laplacian(end,1)-1; %delete unwanted 1 in bottom right corner

K_mat = (Tf/dx)*my_Laplacian;
M_mat = (M/n)*I_n;
end
