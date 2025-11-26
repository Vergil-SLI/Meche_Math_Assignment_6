%% functions to animate the string-weight system

function animate_string_sys(Vlist, tlist, string_params, save_vid)

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
    y_data = Vlist(1, 1:n); u0 = 0; uf = uf_list(1);
    string_plot = plot(x_data, [u0, y_data, uf], "o-", ...
        MarkerSize=5, MarkerFaceColor="r", MarkerEdgeColor="r", LineWidth=1); 
    hold on;

    % Animation: loop and plot each timestep
    for i = 2:length(tlist)
        % delay according to timestep
        pause(0.1*tdiff(i-1));

        % update ydata in string plot to be U vals
        set(string_plot, 'ydata', [u0, Vlist(i, 1:n), uf_list(i)]);
        
        % redraw
        drawnow;
        axis equal;

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
