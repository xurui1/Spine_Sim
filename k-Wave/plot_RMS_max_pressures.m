function plot_RMS_max_pressures(kgrid,sensor_data_fluid,sensor_data_elastic,freq_data)
% plot and save the rms recorded pressure and HM line
    % open new figure
    figure;
    subplot(2, 2, 1), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data_fluid.p_rms,[-1 1]);
    title('RMS Pressure fluid simulation');
    colormap(getColorMap);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
    
    subplot(2, 2, 2), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data_fluid.p_max,[-1 1]);
    title('Max Pressure fluid simulation');
    colormap(getColorMap);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
    
    subplot(2, 2, 3), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data_elastic.p_rms, [-1 1]);
    title('RMS Pressure elastic simulation');
    colormap(getColorMap);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
        
    subplot(2, 2, 4), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data_elastic.p_max, [-1 1]);
    title('Max Pressure elastic simulation');
    colormap(getColorMap);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
    fname = sprintf('./Spine_Sim/HM_freq%d.fig', freq_data.frequency(freq_tracker));
    savefig(fname);