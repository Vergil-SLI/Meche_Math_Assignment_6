%% functions to animate the string-weight system (with color coordination)

function animate_string_sys_V2(Vlist, tlist, string_params, save_vid)

    % File set up:
    if save_vid == true
        % define location and filename where video will be stored
        mypath1 = 'C:\Users\lodio\OneDrive - Olin College of Engineering\Desktop\';
        mypath2 = 'Classes\Junior Fall\Orion Math\Assignment-06\';
        fname = 'video.avi';
        input_fname = [mypath1, mypath2, fname];

        % create a videowriter to write frames to the animation file
        writerObj = VideoWriter(input_fname);
        open(writerObj);
    end
    
    % initialize system vars
    string_length = string_params.L; % total horizontal length
    n = string_params.n; % number of masses on the string
    Uf_func = string_params.Uf_func; % func describing the motion of uf
    
    % precompute lists
    uf_list = Uf_func(tlist);
    tdiff = diff(tlist);

    % initialize figure
    fig1 = figure(1); clf(fig1); hold on; 
    axis([-0.1*string_length, 1.1*string_length, -1, 1]);
    
    % initialize plot 
    x_data = linspace(0, string_length, n+2);
    V_data = Vlist(1, 1:n); u0 = 0; uf = uf_list(1);
    y_data = [u0, V_data, uf];
    string_plot_struct = initialize_plot(x_data, y_data);

    % Animation: loop and plot each timestep
    for i = 2:length(tlist)
        % delay according to timestep
        pause(0.3*tdiff(i-1));

        % update ydata in string plot to be U vals
        ydata_new = [u0, Vlist(i, 1:n), uf_list(i)];
        update_plots(string_plot_struct, ydata_new); 
        
        % redraw
        drawnow;
        axis equal; ylim([-1, 1]);

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

function string_plot_struct = initialize_plot(xdata, ydata)

    C = orderedcolors("gem");
    string_plot_struct = struct();

    % plot just the lines between points (blue continuous)
    string_plot_struct.line_plot = plot(xdata, ydata, "-", Color=C(1,:), LineWidth=1.5);
    hold on;

    % plot markers for the masses (filled red circles)
    string_plot_struct.mass_plot = plot(xdata(2:end-1), ydata(2:end-1), ...
        "o", MarkerSize=5, MarkerFaceColor=C(2,:), MarkerEdgeColor=C(2,:));

    % plot the marker for the fixed end (black x)
    string_plot_struct.fixed_end = plot(xdata(1), ydata(1), "kx", MarkerSize=7, LineWidth=2);

    % plot the marker for the input end (blue diamond)
    string_plot_struct.input_end = plot(xdata(end), ydata(end), ...
        Marker="diamond", MarkerSize=5, MarkerFaceColor=C(3,:), MarkerEdgeColor=C(3,:));

end

%% update plots

function update_plots(string_plot_struct, ydata)

    set(string_plot_struct.line_plot, 'ydata', ydata);
    set(string_plot_struct.mass_plot, 'ydata', ydata(2:end-1));
    set(string_plot_struct.input_end, 'ydata', ydata(end));

end