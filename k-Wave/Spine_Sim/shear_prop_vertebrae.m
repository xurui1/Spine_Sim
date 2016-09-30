%Simulating two transducers, an object, and accounting for shear waves
% Shear Waves And Critical Angle Reflection Example


clear all;

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% change scale to 2 to reproduce the higher resolution figures used in the
% help file
scale = 1;

% create the computational grid
Nx = 472;           % number of grid points in the x (row) direction
Ny = 472;           % number of grid points in the y (column) direction
PML_size = 20;
x_length = 100e-3; 
dx = x_length/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
vsound = 1580; %[m/s]

% create the time array
cfl   = 0.1;
t_end = 12e-5;
kgrid.t_array= makeTime(kgrid, vsound, cfl, t_end);

% define two curved transducer elements
source_mask = makedoubleSemiCircle(Nx, Ny,round((Nx)/2),round(Ny/2),...
    round(Nx/2), pi/12, pi/3+pi/12,round((Nx)/2),round(Ny/2),round(Nx/2), pi/12+pi/2, pi/3+pi/12+pi/2);
source.p_mask = source_mask;


% define the source properties
% define constant varying sinusoidal sources
source_freq1 = 0.25e6;       % [Hz]
source_mag1 = 0.5;           % [Pa]
source_mag2 = 0.5;           % [Pa]

% define changing source (source 2) 
initial = 0.25e6;            % [Hz]
final = 0.25e6;               % [Hz]
step_size = 0.01e6;          % [Hz]
Num_steps = (final - initial)/step_size;

% define the sensor to record the maximum particle velocity everywhere
sensor.record = {'u_max_all'};

% define the medium properties structure
    I = loadImage('vertebrae.png');
    scale = Nx/745;
    medium_img = imresize(I, scale);
    medium_img = im2bw(medium_img,0.01);

step_tracker = 1;

for source_freq2 = initial:step_size:final 
    
    % define struct of frequencies, max rms pressure locations
    freq_data = struct('frequency',[],'x_max_elastic',[],'y_max_elastic',[],...
    'HM_x_max_elastic',[],'HM_y_max_elastic',[],'FWHM_elastic',[],...
    'x_max_fluid',[],'y_max_fluid',[],...
    'HM_x_max_fluid',[],'HM_y_max_fluid',[],'FWHM_fluid',[]);
    freq_data.frequency(step_tracker) = source_freq2;
    
    
    % =========================================================================
    % FLUID SIMULATION
    % =========================================================================
    
    for i = 1:Nx
        for j = 1:Ny
            if medium_img(i,j) == 0 
                medium.sound_speed(i,j) = 1580; %[m/s]
                medium.density(i,j) = 1000; %[kg/m^3]
                medium.alpha_coeff(i,j) = 0.57; %[dB/(MHz^y cm)]
            elseif medium_img(i,j) == 1
                medium.sound_speed(i,j) = 2820; %[m/s]
                medium.density(i,j) = 1800;    %[kg/m^3]
                medium.alpha_coeff(i,j) = 9; %[dB/(MHz^y cm)]
            end
        end
    end
    medium.alpha_power = 2;
   
    %Here I define the frequencies of my source points for fluid simulation
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
    % create a display mask to display the transducer for fluid simulation
    display_mask = source.p_mask;

    % create a sensor mask covering the entire computational domain using the
    % opposing corners of a rectangle
    sensor.mask = [1, 1, Nx, Ny].';

    % set the record mode capture the final wave-field and the statistics at
    % each sensor point 
    sensor.record = {'p_final', 'p_max', 'p_rms','u_max_all'};
    max_distance = 0.05 + dx*sqrt((Nx/2)^2+(Ny/2)^2)
    travel_time = max_distance /vsound;
    i = 1;
    while kgrid.t_array(i) < travel_time
        i = i+1;
    end
    sensor.record_start_index = i;

    % assign the input options for fluid simulation
    input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

    % run the fluid simulation
    sensor_data_fluid = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    %Find location of peak intensity for fluid simulation
    i_max = 0;
    j_max = 0;
    max_rms_p = 0;
    for i = 1:Nx
        for j = 1:Ny
            if max_rms_p < sensor_data_fluid.p_rms(i,j)
                i_max = i;
                j_max = j;
                max_rms_p = sensor_data_fluid.p_rms(i,j);
            end
        end
    end
    freq_data.x_max_fluid(step_tracker) = i_max;
    freq_data.y_max_fluid(step_tracker) = j_max;
    
    %Find FWHM (making assumption that the profile is approximately circular)
    half_max_rms_p = max_rms_p/2;
    half_max_coords_fluid = struct('x',[],'y',[]);
    num_half_max_pts = 1;
    for i = 2:Nx-1
        for j = 2:Ny-2
            %compute max,min at local 3-by-3 grid
            local_max = sensor_data_fluid.p_rms(i,j);
            local_min = sensor_data_fluid.p_rms(i,j);
            for local_i = (i-1):(i+1)
                for local_j = (j-1):(j+1)
                    if sensor_data_fluid.p_rms(local_i,local_j) < local_min
                        local_min = sensor_data_fluid.p_rms(local_i,local_j);
                    end
                    if sensor_data_fluid.p_rms(local_i,local_j) > local_max
                        local_max = sensor_data_fluid.p_rms(local_i,local_j);
                    end
                end
            end
            %record location of HM point
            if local_max >= half_max_rms_p  && local_min <=half_max_rms_p
                half_max_coords_fluid.x(num_half_max_pts) = i;
                half_max_coords_fluid.y(num_half_max_pts) = j;
                num_half_max_pts = num_half_max_pts + 1;
            end
        end
    end
    %find centre of HM circle for fluid simulation
    HM_circle_x = 0;
    HM_circle_y = 0;
     for i = 1:num_half_max_pts-1
          HM_circle_x = HM_circle_x+half_max_coords_fluid.x(i);
          HM_circle_y = HM_circle_y+half_max_coords_fluid.y(i);
     end
     HM_circle_x = HM_circle_x/(num_half_max_pts-1);
     HM_circle_y = HM_circle_y/(num_half_max_pts-1);
     freq_data.HM_x_max_fluid(step_tracker)= HM_circle_x;
     freq_data.HM_y_max_fluid(step_tracker)= HM_circle_y;


     %find mean distance from HM point to centre of HM circle
     mean_norm=0; 
     for i = 1:num_half_max_pts-1
         mean_norm = mean_norm + sqrt((half_max_coords_fluid.x(i)-HM_circle_x)^2 + (half_max_coords_fluid.y(i)-HM_circle_y)^2);
     end
     mean_norm = mean_norm/(num_half_max_pts-1);
     freq_data.FWHM_fluid(step_tracker) = 2.0*dx*mean_norm;
     
     %plot HM points to see where the points actually are in fluid
     %simulation
     FWHM_matrix = zeros(Nx,Ny);
     for i = 1:num_half_max_pts-1
         FWHM_matrix(half_max_coords_fluid.x(i),half_max_coords_fluid.y(i)) = 1;
     end
     sensor_data_fluid.p_rms(source.p_mask ~= 0) = 1;
     sensor_data_fluid.p_rms(FWHM_matrix == 1) = 1;

    % =========================================================================
    % ELASTIC SIMULATION
    % =========================================================================

    % clear medium properties for the elastic simulation
    clear medium
    % define the medium properties for the elastic simulation
    
    for i = 1:Nx
        for j = 1:Ny
            if medium_img(i,j) == 0
                medium.sound_speed_compression(i,j) = 1580; %[m/s]
                medium.sound_speed_shear(i,j) = 0; %[m/s]
                medium.density(i,j) = 1000; %[kg/m^3]
                medium.alpha_coeff_compression(i,j) = 0.57; %[dB/(MHz^y cm)]
                medium.alpha_coeff_shear(i,j) = 0; %[dB/(MHz^y cm)]
            elseif medium_img(i,j) == 1
                medium.sound_speed_compression(i,j) = 2820; %[m/s]
                medium.sound_speed_shear(i,j) = 1500; %[m/s]
                medium.density(i,j) = 1800;    %[kg/m^3]
                medium.alpha_coeff_compression(i,j) = 9; %[dB/(MHz^y cm)]
                medium.alpha_coeff_shear(i,j) = 20; %[dB/(MHz^y cm)]
            end
        end
    end
   % medium.alpha_power = 2;
    
    
    % assign the source for the elastic simulation
    clear source
    source.s_mask = source_mask;

     % create a sensor mask covering the entire computational domain using the
    % opposing corners of a rectangle
    clear sensor;
    sensor.mask = [1, 1, Nx, Ny].';
    max_distance = 0.05 + dx*sqrt((Nx/2)^2+(Ny/2)^2)
    travel_time = max_distance /vsound;
    i = 1;
    while kgrid.t_array(i) < travel_time
        i = i+1;
    end
    sensor.record_start_index = i;

    % set the record mode capture the final wave-field and the statistics at
    % each sensor point 
    sensor.record = {'p_final', 'p_max', 'p_rms','u_max_all'};


    %Here I define the frequencies of my source points for the elastic
    %simulation
    p_number = 1;
    for j = 1:Ny/2
      for i = 1:Nx
          if source.s_mask(i,j) == 1
               source.sxx(p_number,:) = source_mag1*sin(2*pi*source_freq1*kgrid.t_array);
               %source.sxy(p_number,:) = source_mag1*sin(2*pi*source_freq1*kgrid.t_array);
               source.syy(p_number,:) = source_mag1*sin(2*pi*source_freq1*kgrid.t_array);
               p_number = p_number + 1;
          end
      end
    end
    for j = (Ny/2+1):Ny
        for i = 1:Nx
          if source.s_mask(i,j) == 1
               source.sxx(p_number,:) = source_mag2*sin(2*pi*source_freq2*kgrid.t_array);
               %source.sxy(p_number,:) = source_mag2*sin(2*pi*source_freq2*kgrid.t_array);
               source.syy(p_number,:) = source_mag2*sin(2*pi*source_freq2*kgrid.t_array);
               p_number = p_number + 1;
          end
        end
    end
    
    % filter the source to remove high frequencies not supported by the grid
    for i = 1:p_number-1
        source.sxx(i,:) = filterTimeSeries(kgrid, medium, source.sxx(i,:));
        %source.sxy(i,:) = filterTimeSeries(kgrid, medium, source.sxy(i,:));
        source.syy(i,:) = filterTimeSeries(kgrid, medium, source.syy(i,:));
    end
    % create a display mask to display the transducer for the elastic
    % simulation
    display_mask = source.s_mask;

    % assign the input options
    input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

    % run the fluid simulation
    sensor_data_elastic = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});
    
    %Find location of peak intensity
    i_max = 0;
    j_max = 0;
    max_rms_p = 0;
    for i = 1:Nx
        for j = 1:Ny
            if max_rms_p < sensor_data_elastic.p_rms(i,j)
                i_max = i;
                j_max = j;
                max_rms_p = sensor_data_elastic.p_rms(i,j);
            end
        end
    end
    freq_data.x_max_elastic(step_tracker) = i_max;
    freq_data.y_max_elastic(step_tracker) = j_max;
    
    %Find FWHM (making assumption that the profile is approximately circular)
    half_max_rms_p = max_rms_p/2;
    half_max_coords_elastic = struct('x',[],'y',[]);
    num_half_max_pts = 1;
    for i = 2:Nx-1
        for j = 2:Ny-2
            %compute max,min at local 3-by-3 grid
            local_max = sensor_data_elastic.p_rms(i,j);
            local_min = sensor_data_elastic.p_rms(i,j);
            for local_i = (i-1):(i+1)
                for local_j = (j-1):(j+1)
                    if sensor_data_elastic.p_rms(local_i,local_j) < local_min
                        local_min = sensor_data_elastic.p_rms(local_i,local_j);
                    end
                    if sensor_data_elastic.p_rms(local_i,local_j) > local_max
                        local_max = sensor_data_elastic.p_rms(local_i,local_j);
                    end
                end
            end
            %record location of HM point
            if local_max >= half_max_rms_p  && local_min <=half_max_rms_p
                half_max_coords_elastic.x(num_half_max_pts) = i;
                half_max_coords_elastic.y(num_half_max_pts) = j;
                num_half_max_pts = num_half_max_pts + 1;
            end
        end
    end
    %find centre of HM circle 
    HM_circle_x = 0;
    HM_circle_y = 0;
     for i = 1:num_half_max_pts-1
          HM_circle_x = HM_circle_x+half_max_coords_elastic.x(i);
          HM_circle_y = HM_circle_y+half_max_coords_elastic.y(i);
     end
     HM_circle_x = HM_circle_x/(num_half_max_pts-1);
     HM_circle_y = HM_circle_y/(num_half_max_pts-1);
     freq_data.HM_x_max_elastic(step_tracker)= HM_circle_x;
     freq_data.HM_y_max_elastic(step_tracker)= HM_circle_y;


     %find mean distance from HM point to centre of HM circle
     mean_norm=0; 
     for i = 1:num_half_max_pts-1
         mean_norm = mean_norm + sqrt((half_max_coords_elastic.x(i)-HM_circle_x)^2 + (half_max_coords_elastic.y(i)-HM_circle_y)^2);
     end
     mean_norm = mean_norm/(num_half_max_pts-1);
     freq_data.FWHM_elastic(step_tracker) = 2.0*dx*mean_norm;
     
     %plot HM points to see where the points actually are
     FWHM_matrix = zeros(Nx,Ny);
     for i = 1:num_half_max_pts-1
         FWHM_matrix(half_max_coords_elastic.x(i),half_max_coords_elastic.y(i)) = 1;
     end
     sensor_data_elastic.p_rms(source.s_mask ~= 0) = 1;
     sensor_data_elastic.p_rms(FWHM_matrix == 1) = 1;
     %sensor_data_elastic.p_rms(medium_img == 1) = 1;
     % plot and save the rms recorded pressure and HM line
    %open new figure
    figure;
    subplot(1, 2, 1), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data_fluid.p_rms);
    colormap(getColorMap);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
    title('RMS Pressure fluid simulation');
    scaleFig(2, 1);   
    
    subplot(1, 2, 2), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data_elastic.p_rms);
    colormap(getColorMap);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
    title('RMS Pressure elastic simulation');
    scaleFig(2, 1);   
    fname = sprintf('./Spine_Sim/HM_freq%d.fig', freq_data.frequency(step_tracker));
    %savefig(fname);

end


% =========================================================================
% VISUALISATION
% =========================================================================

% define plot vector
x_vec = kgrid.x_vec(1 + PML_size:end - PML_size)*1e3;
y_vec = kgrid.y_vec(1 + PML_size:end - PML_size)*1e3;

% calculate square of velocity magnitude
u_e = sensor_data_elastic.ux_max_all.^2 + sensor_data_elastic.uy_max_all.^2;
u_f = sensor_data_fluid.ux_max_all.^2 + sensor_data_fluid.uy_max_all.^2;

% plot layout
figure;
imagesc(y_vec, x_vec, double(source_mask ));
xlabel('y [mm]');
ylabel('x [mm]');
axis image;
colormap(flipud(gray));

% plot beam patterns
figure;
subplot(2, 1, 1);
imagesc(y_vec, x_vec, 20*log10(u_f./max(u_f(:))));
xlabel('y [mm]');
ylabel('x [mm]');
axis image;
colorbar;
caxis([-50, 0]);
title('Fluid Model');

subplot(2, 1, 2);
imagesc(y_vec, x_vec, 20*log10(u_e./max(u_e(:))));
xlabel('y [mm]');
ylabel('x [mm]');
axis image;
colorbar;
caxis([-50, 0]);
title('Elastic Model');
colormap(jet(256));