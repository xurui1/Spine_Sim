% Simulating Transducer Field Patterns Example
%
% In this code I create two curved transducers that are driven at different
% frequencies

clear all;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 216;           % number of grid points in the x (row) direction
Ny = 216;           % number of grid points in the y (column) direction
x_length = 100e-3; 
dx = x_length/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
% medium.sound_speed = 1500;  % [m/s] ok for homogeneous medium
% define the properties of the propagation medium    
vsound = 1500; 
medium.sound_speed = makemediaDisk(Nx,Ny,round(Nx/2),round(Nx/2),round(Nx/10),vsound,vsound);         % [m/s]
medium.density = makedensityDisc(Nx,Ny,round(Nx/2),round(Nx/2),round(Nx/10),1000,1000);         % [kg/m^3]
medium.alpha_power = 1.5;   % [dB/(MHz^y cm)]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]

% create the time array
cfl   = 0.1;
t_end = 12e-5;
kgrid.t_array= makeTime(kgrid, 1500, cfl, t_end);

% define two curved transducer elements
source.p_mask = makedoubleSemiCircle(Nx, Ny,round(Nx/2),round(Ny/2),round(Nx/2),...
    pi/12, pi/3+pi/12,round(Nx/2),round(Ny/2),round(Nx/2), pi/12+pi/2, pi/3+pi/12+pi/2);

% define constant varying sinusoidal sources
source_freq1 = 0.25e6;       % [Hz]
source_mag1 = 0.5;           % [Pa]

% define changing source (source 2) 
initial = 0.2e6;            % [Hz]
final = 0.3e6;               % [Hz]
step_size = 0.01e6;          % [Hz]
Num_steps = (final - initial)/step_size;

% define struct of frequencies, max rms pressure locations
freq_data = struct('frequency',[],'x_max',[],'y_max',[],'HM_x_max',[],'HM_y_max',[],'FWHM',[]);

step_tracker = 1;

for source_freq2 = initial:step_size:final 
    source_mag2 = 0.5;           % [Pa]
    freq_data.frequency(step_tracker) = source_freq2;

    %Here I define the frequencies of my source points
    p_number = 1;
    for j = 1:Ny/2
      for i = 1:Nx
          if source.p_mask(i,j) == 1
               source.p(p_number,:) = source_mag1*sin(2*pi*source_freq1*kgrid.t_array);
               p_number = p_number + 1;
          end
      end
    end
    for j = (Ny/2+1):Ny
        for i = 1:Nx
          if source.p_mask(i,j) == 1
               source.p(p_number,:) = source_mag2*sin(2*pi*source_freq2*kgrid.t_array);
               p_number = p_number + 1;
          end
        end
    end

% filter the source to remove high frequencies not supported by the grid
    for i = 1:p_number-1
        source.p(i,:) = filterTimeSeries(kgrid, medium, source.p(i,:));
    end
% create a display mask to display the transducer
    display_mask = source.p_mask;

% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
    sensor.mask = [1, 1, Nx, Ny].';

% set the record mode capture the final wave-field and the statistics at
% each sensor point 
    sensor.record = {'p_final', 'p_max', 'p_rms'};
% define time to start recording data by calculating the time it takes for
% the US wave to propagate to the furthest point in the US field, in this
% case from a transducer to the opposite corner (5cm from focal point plus
% distance between focal point and far corner.
    max_distance = 0.05 + dx*sqrt((Nx/2)^2+(Ny/2)^2)
    travel_time = max_distance /vsound;
    i = 1;
    while kgrid.t_array(i) < travel_time
        i = i+1;
    end
    sensor.record_start_index = i;



% assign the input options
    input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

    
    
% run the simulation
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

%Find location of peak intensity
    i_max = 0;
    j_max = 0;
    max_rms_p = 0;
    for i = 1:Nx
        for j = 1:Ny
            if max_rms_p < sensor_data.p_rms(i,j)
                i_max = i;
                j_max = j;
                max_rms_p = sensor_data.p_rms(i,j);
            end
        end
    end
    freq_data.x_max(step_tracker) = i_max;
    freq_data.y_max(step_tracker) = j_max;

%Find FWHM (making assumption that the profile is approximately circular)
    half_max_rms_p = max_rms_p/2;
    half_max_coords = struct('x',[],'y',[]);
    num_half_max_pts = 1;
    for i = 2:Nx-1
        for j = 2:Ny-2
            %compute max,min at local 3-by-3 grid
            local_max = sensor_data.p_rms(i,j);
            local_min = sensor_data.p_rms(i,j);
            for local_i = (i-1):(i+1)
                for local_j = (j-1):(j+1)
                    if sensor_data.p_rms(local_i,local_j) < local_min
                        local_min = sensor_data.p_rms(local_i,local_j);
                    end
                    if sensor_data.p_rms(local_i,local_j) > local_max
                        local_max = sensor_data.p_rms(local_i,local_j);
                    end
                end
            end
            %record location of HM point
            if local_max >= half_max_rms_p  && local_min <=half_max_rms_p
                half_max_coords.x(num_half_max_pts) = i;
                half_max_coords.y(num_half_max_pts) = j;
                num_half_max_pts = num_half_max_pts + 1;
            end
        end
    end
    %find centre of HM circle 
    HM_circle_x = 0;
    HM_circle_y = 0;
     for i = 1:num_half_max_pts-1
          HM_circle_x = HM_circle_x+half_max_coords.x(i);
          HM_circle_y = HM_circle_y+half_max_coords.y(i);
     end
     HM_circle_x = HM_circle_x/(num_half_max_pts-1);
     HM_circle_y = HM_circle_y/(num_half_max_pts-1);
     freq_data.HM_x_max(step_tracker)= HM_circle_x;
     freq_data.HM_y_max(step_tracker)= HM_circle_y;


     %find mean distance from HM point to centre of HM circle
     mean_norm=0; 
     for i = 1:num_half_max_pts-1
         mean_norm = mean_norm + sqrt((half_max_coords.x(i)-HM_circle_x)^2 + (half_max_coords.y(i)-HM_circle_y)^2);
     end
     mean_norm = mean_norm/(num_half_max_pts-1);
     freq_data.FWHM(step_tracker) = 2.0*dx*mean_norm;
     
     %plot HM points to see where the points actually are
     FWHM_matrix = zeros(Nx,Ny);
     for i = 1:num_half_max_pts-1
         FWHM_matrix(half_max_coords.x(i),half_max_coords.y(i)) = 1;
     end
     sensor_data.p_rms(source.p_mask ~= 0) = 1;
     sensor_data.p_rms(FWHM_matrix == 1) = 1;
    % plot and save the rms recorded pressure and HM line
    figure;
    subplot(1, 1, 1), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data.p_rms, [-1 1]);
    colormap(getColorMap);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
    title('RMS Pressure');
    scaleFig(2, 1);   
    fname = sprintf('./Spine_Sim/HM_freq%d.fig', freq_data.frequency(step_tracker));
    %savefig(fname);
    
    
    %Update step tracker
    step_tracker = step_tracker+1;
end

% =========================================================================
% VISUALISATION
% =========================================================================

% add the source mask onto the recorded wave-field
sensor_data.p_final(source.p_mask ~= 0) = 1;
sensor_data.p_max(source.p_mask ~= 0) = 1;
sensor_data.p_rms(source.p_mask ~= 0) = 1;

% plot the final wave-field
figure;
subplot(1, 4, 1), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data.p_final, [-1 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Final Wave Field');

% plot the maximum recorded pressure
%subplot(1, 4, 2), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data.p_max, [-1 1]);
%colormap(getColorMap);
%ylabel('x-position [mm]');
%xlabel('y-position [mm]');
%axis image;
%title('Maximum Pressure');

%plot the change in location of the centre of the FWHM 'circle'
subplot(1, 4, 2), scatter( freq_data.HM_x_max,freq_data.HM_y_max );
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Centre of FWHM circle for RMS Pressure');
scaleFig(2, 1);

% plot the rms recorded pressure
%subplot(1, 4, 3), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data.p_rms, [-1 1]);
%colormap(getColorMap);
%ylabel('x-position [mm]');
%xlabel('y-position [mm]');
%axis image;
%title('RMS Pressure');
%scaleFig(2, 1);

%plot the change in FWHM as a function of frequency
subplot(1, 4, 3), plot( freq_data.frequency, freq_data.FWHM );
ylabel('FWHM');
xlabel('Frequency of source');
axis image;
title('FWHM as a function of frequency');
scaleFig(2, 1);

% plot the change in location of max rms recorded pressure for different
% frequencies
subplot(1, 4, 4), scatter( freq_data.x_max,freq_data.y_max );
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Max RMS Pressure location for different frequencies');
scaleFig(2, 1);