%% Simulate the string system via differential equations for a specified
% environment / input. Here, the simulate code is refactored so that
% the user can simply input 
% each mass moves acoording to a sine fucntion**

function simulate_mode_soln(U_r, omega, damping, tspan, string_params, save_vid)
    
    % unpack string_params
    num_masses = string_params.n;
    string_params.c = damping; % set damping coeff
    
    % solution has form
    % U(t) = U_r*cos(ω_r*t + φ)  

    string_params.Uf_func = @(t_in) 0*t_in;
    string_params.dUfdt_func = @(t_in) 0*t_in;

    % generate a tlist
    tlist = linspace(tspan(1), tspan(2), 100*diff(tspan));
    
    % construct Vlist (according to the ideal cosine movement)
    Vlist = zeros(length(tlist), num_masses);
    for i = 1:num_masses
        Vlist(:,i) = U_r(i)*cos(omega*tlist - pi/2);
    end

    % animate result
    animate_mode_comparison(Vlist, tlist, string_params, save_vid)

end