
function sim_mode_comparison(U_r, omega, damping, tlist, string_params, save_vid)
    
    % extract system params
    num_masses = string_params.n;
    omega_Uf = omega; % use the resonant frequency
    A = 1; % ignore amplitude coefficient
    
    % construct the forcing function and its derivative:
    % uf = A*cos(w*t), uf' = -w*A*sin(w*t)
    string_params.Uf_func = @(t_in) A*cos(omega_Uf*t_in);
    string_params.dUfdt_func =  @(t_in) -omega_Uf*A*sin(omega_Uf*t_in);
    string_params.c = damping;

    % create anonymized function
    my_rate_func = @(t_in, V_in) string_rate_func01(t_in, V_in, string_params);
    
    % initial conditions
    V0 = zeros(2*num_masses, 1);
    
    % run the integration
    [t_list, Vlist] = ode45(my_rate_func, tlist, V0);

    % construct Mlist (according to the ideal cosine movement)
    Mlist = zeros(length(t_list), num_masses+2);
    for i = 1:num_masses
        Mlist(:,i+1) = U_r(i)*cos(omega*t_list - pi/2);
    end

    animate_mode_comparison(Vlist, Mlist, t_list, string_params, save_vid)

end

%% functions to animate the string-weight and the ideal mode system 
% (with color coordination)
% mlist is ideal mode movement w/ fixed ends [0, m1, m2, ..., mn, 0]

function animate_mode_comparison(Vlist, Mlist, tlist, string_params, save_vid)

    % File set up:
    if save_vid == true
        % define location and filename where video will be stored
        mypath1 = 'C:\Users\lodio\OneDrive - Olin College of Engineering\Desktop\';
        mypath2 = 'Classes\Junior Fall\Orion Math\Assignment-06\';
        fname = 'comparison.avi';
        input_fname = [mypath1, mypath2, fname];

        % create a videowriter to write frames to the animation file
        writerObj = VideoWriter(input_fname);
        open(writerObj);
    end
    
    % initialize system vars
    L = string_params.L; % total horizontal length
    n = string_params.n; % number of masses on the string
    Uf_func = string_params.Uf_func; % func describing the motion of uf
    
    % precompute lists
    uf_list = Uf_func(tlist);
    x_data = linspace(0, L, n+2);
    y_data = [0*tlist, Vlist(:, 1:n), uf_list]; % y pos data [u0, mass pos, uf]

    % scale Mlist 
    ymax_sim = max(abs(y_data), [], "all");
    Mlist = (ymax_sim/max(abs(Mlist), [], "all"))*Mlist;

    % initialize figure
    fig1 = figure(1); hold on; 

    % determine axis limits
    axis_lims = [-0.1*L, 1.1*L, -1.1*ymax_sim, 1.1*ymax_sim];
    
    % initialize string & mode plots 
    string_plot_struct = initialize_string(x_data, y_data(1,:), axis_lims);
    mode_plot_struct = initialize_mode(x_data, Mlist(1,:), axis_lims);

    % Animation: loop and plot each timestep
    for i = 2:length(tlist)

        % update ydata in string plot to be U vals
        update_string(string_plot_struct, y_data(i,:)); 

        % update mode plot
        update_mode(mode_plot_struct, Mlist(i,:));

        % redraw
        drawnow;

        if save_vid == true
            % capture a frame (what is currently plotted)
            current_frame = getframe(fig1);
            % write the frame to the video
            writeVideo(writerObj, current_frame);
        end
    end

    % Clean up: close writer object if recording
    if save_vid == true
        close(writerObj)
    end
end

%% initialize string plots

function string_plot_struct = initialize_string(xdata, ydata, axis_lims)

    C = orderedcolors("gem");
    string_plot_struct = struct();

    % plot just the lines between points (blue continuous)
    subplot(2, 1, 1); 
    string_plot_struct.line_plot = plot(xdata, ydata, "-", Color=C(1,:), LineWidth=1.5);
    hold on; axis(axis_lims);

    % plot markers for the masses (filled red circles)
    string_plot_struct.mass_plot = plot(xdata(2:end-1), ydata(2:end-1), ...
        "o", MarkerSize=5, MarkerFaceColor=C(2,:), MarkerEdgeColor=C(2,:));

    % plot the marker for the fixed end (black x)
    string_plot_struct.fixed_end = plot(xdata(1), ydata(1), "kx", MarkerSize=7, LineWidth=2);

    % plot the marker for the input end (blue diamond)
    string_plot_struct.input_end = plot(xdata(end), ydata(end), ...
        Marker="diamond", MarkerSize=5, MarkerFaceColor=C(3,:), MarkerEdgeColor=C(3,:));

    title("String simulation")

end

%% update plots

function update_string(string_plot_struct, ydata)

    set(string_plot_struct.line_plot, 'ydata', ydata);
    set(string_plot_struct.mass_plot, 'ydata', ydata(2:end-1));
    set(string_plot_struct.input_end, 'ydata', ydata(end));

end

%% initialize mode plot

function mode_plot_struct = initialize_mode(xdata, ydata, axis_lims)

    C = orderedcolors("gem");
    mode_plot_struct = struct();

    % plot just the lines between points (blue continuous)
    subplot(2, 1, 2); 
    mode_plot_struct.line_plot = plot(xdata, ydata, "--", Color=C(1,:), LineWidth=1.5);
    hold on; axis(axis_lims);

    % plot markers for the masses (filled red circles)
    mode_plot_struct.mass_plot = plot(xdata(2:end-1), ydata(2:end-1), ...
        "o", MarkerSize=5, MarkerFaceColor=C(2,:), MarkerEdgeColor=C(2,:));

    % plot the marker for the fixed ends (black x)
    mode_plot_struct.ends = plot([xdata(1), xdata(end)], [ydata(1), ydata(end)], ...
        "kx", MarkerSize=7, LineWidth=2);

    title("Mode simulation")
    
end

%% update mode plot

function update_mode(mode_plot_struct, ydata)

    set(mode_plot_struct.line_plot, 'ydata', ydata);
    set(mode_plot_struct.mass_plot, 'ydata', ydata(2:end-1));
    set(mode_plot_struct.ends, 'ydata', [ydata(1), ydata(end)]);

end
