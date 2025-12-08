%% simulate string and weights system

function string_simulation()
    %Normal or Wave Animation 0 = Normal # of Masses, 1 = Wave Test
    Selection = 1;
    %Wave Test Pulse Shape: 0 = Triangle 1 = BSpline
    pulse = 1;

    if Selection == 0 
        num_masses = 20; % per later analysis
        damping_coeff = 0.01;% selected per Orion's suggestion 0.1-0.01
    elseif Selection == 1
        num_masses = 300; % per later analysis
        damping_coeff = 0; % 0 for wave
    else
        disp('Bad Input Selection')
    end

    % assign system params
    total_mass = 4.5; % arbitrary, may need to change
    tension_force = 5; % arbitrary, may need to change
    string_length = 5; % selected arbitrarily
    
    dx = string_length/(num_masses+1);
    amplitude_Uf = 0.5; % CHOSEN
    omega_Uf = 2; % CHOSEN

    if Selection == 0 
        % construct the forcing function and its derivative:
        % uf = A*cos(w*t), uf' = -w*A*sin(w*t)
        Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
        dUfdt_func = @(t_in) - omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);
    elseif Selection == 1
        %Pulse Function as Boundaries
        w = .75; %Chosen
        if pulse == 0 
            Uf_func = @(t_in) INFCN_triangle_pulse(t_in, w, amplitude_Uf);
            dUfdt_func = @(t_in) INFCN_triangle_pulse_derivative(t_in, w, amplitude_Uf);
        elseif pulse == 1
            Uf_func = @(t_in) INFCN_b_spline_pulse(t_in, w, amplitude_Uf);
            dUfdt_func = @(t_in) INFCN_b_spline_pulse_derivative(t_in, w, amplitude_Uf);
        else
            disp('Bad Input Pulse')
        end
    else
        disp('Bad Input Selection')
    end
    
    % generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    
    % create anonymized function
    my_rate_func = @(t_in, V_in) string_rate_func01(t_in, V_in, string_params);
    
    % initial conditions
    U0 = zeros(num_masses, 1); % CHOSEN
    dUdt0 = zeros(num_masses, 1); % CHOSEN
    V0 = [U0; dUdt0];
    tspan = [0, 15]; % CHOSEN
    
    % run the integration
    [tlist, Vlist] = ode45(my_rate_func, tspan, V0); % may want to specify abstol?
    
    % animate the system
    save_vid = false;
    animate_string(Vlist, tlist, string_params, save_vid)
end

