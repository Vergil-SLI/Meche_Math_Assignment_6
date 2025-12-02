%% Modal system analysis
% explanations for each step included in the collapsed comments

%% Establish system 

string_params.n = 3; % num masses on string
string_params.M = 5; % total mass on string;
string_params.Tf = 5; % tension force
string_params.L = 5; % string length 
string_params.dx = string_params.L/string_params.n;

%% Solve the system 

% Get M and K
% If we assume that the damping is zero, and that there is no forcing
% function input, then: 
% MU'' + CU' + KU = B_1*uf(t) + B_2*uf'(t) => MU'' + KU = 0 
[M_mat, K_mat] = construct_2nd_order_matrices(string_params);

% Use MATLAB to solve the generalized eigenvalue problem
% find solutions of the form U(t) = U_r * cos(ω_r*t + φ), where 
% the mode shape vector, U_r, is a vector of constants. => K*U_r = λ*M*U_r 
[Ur_mat, lambda_mat] = eig(K_mat, M_mat);

% Compute the mode shapes and resonant frequencies of the string.
% note: the ith col of Ur_mat and ith diagonal element of lambda_mat  
% correspond to a vector U_r and scalar λ = ω_r^2
mode_shapes = Ur_mat; % by col
resonant_freqs = sqrt(sum(lambda_mat, 1));

%% Choose a mode shape and a resonant frequency. Run the vibrating string
% simulation (with the initial displacements and velocities of zero), using
% an input frequency for uf(t) that is near or at the corresponding
% resonant frequency.

i = 2; omega = resonant_freqs(i);
damping = 0.01; 
tspan = [0, 25];
save_vid = false;

% simulate
simulate_mode(omega, damping, tspan, save_vid)

%% Animate a comparison of the predicted mode shape to the actual result

% two subplots, one with the mode shape and the other with the animation of the vibrating string.
subplot(2, 1, 1);
mode_x = linspace(0, string_params.L, string_params.n+2);
mode_y = 
plot()