
function simulate_mode(U_r, omega, damping, tspan, string_params, save_vid)
    

    % unpack string_params
    num_masses = string_params.n;
    % solution has form
    % U(t) = U_r*cos(ω_r*t + φ)  

    tlist = linspace(tspan(1), tspan(2), 100*diff(tspan));
    
    Vlist = zeros(length(tlist), num_masses);
    for i = 1:num_masses
        Vlist(:,i) = U_r(i)*cos(omega*tlist - pi/2);
    end


    animate_string(Vlist, tlist, string_params, save_vid)


end