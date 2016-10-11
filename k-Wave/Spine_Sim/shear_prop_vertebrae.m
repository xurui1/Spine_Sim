%Simulating two transducers, an object from image fild, and accounting for shear waves
%This code is somewhat based on the Shear Waves And Critical Angle Reflection Example

clear all;

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================
scale = 1;

% create the computational grid
Nx = 216;           % number of grid points in the x (row) direction
Ny = 216;           % number of grid points in the y (column) direction
PML_size = 20;
x_length = 110e-3;  %length in x-direction [mm]
dx = x_length/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
vsound = 1580; %[m/s]
plot_type = 2;      % 0-RMS data, 1-max pressure data, 2-Both

% create the time array
cfl   = 0.1;
t_end = 12e-5;
kgrid.t_array= makeTime(kgrid, vsound, cfl, t_end);

% define the source properties
source_freq1 = 0.25e6;       % [Hz]
source_mag1 = 0.5;           % [Pa]
source_mag2 = 0.5;           % [Pa]
Theta_zero = pi/48;          % [rad]
dTheta = pi/48;              % [rad]
Theta_final = pi/4-pi/48;          % [rad]

% define changing source (source 2) 
initial = 0.25e6;            % [Hz]
final = 0.25e6;              % [Hz]
step_size = 0.01e6;          % [Hz]
Num_steps = (final - initial)/step_size;

% define the medium properties structure
I = loadImage('vertebrae.png');
scale = Nx/745;
medium_img = imresize(I, scale);
medium_img = im2bw(medium_img,0.01);
% load contour of vertebra for plotting
C = loadImage('vert_contour_img.png');
contour_img = imresize(C, scale);
contour_img = im2bw(contour_img,0.01);
    
% define struct of frequencies, max rms pressure locations, etc
freq_data = struct('frequency',[],'x_max_elastic',[],'y_max_elastic',[],...
'HM_x_max_elastic',[],'HM_y_max_elastic',[],'FWHM_elastic',[],...
'x_max_fluid',[],'y_max_fluid',[],...
'HM_x_max_fluid',[],'HM_y_max_fluid',[],'FWHM_fluid',[]);

% define struct of angular data 
theta_data = struct('angle',[],'spinalcord_fluid_RMS_p',[],'spinalcord_fluid_Max_p',[],...
    'spinalcord_fluid_RMS_pnorm',[],'spinalcord_fluid_Max_pnorm',[],...
    'spinalcord_elastic_RMS_p',[],'spinalcord_elastic_Max_p',[],...
    'spinalcord_elastic_RMS_pnorm',[],'spinalcord_elastic_Max_pnorm',[]);

freq_tracker = 1;
theta_tracker = 1;

for theta = Theta_zero:dTheta:Theta_final
    % define two curved transducer elements
    clear source; 
    
    source_mask = makedoubleSemiCircle(Nx, Ny,round((Nx)/2),round(Ny/2),...
        round(0.05*Nx/x_length), pi/4 - theta,pi/4+theta,round((Nx)/2),round(Ny/2),round(0.05*Nx/x_length),...
        3*pi/4 - theta, 3*pi/4+theta); %want the radius of the spherical transducer to be 5cm.
    source.p_mask = source_mask;
    
    theta_data.angle(theta_tracker) = theta;

for source_freq2 = initial:step_size:final 
    
    %Save frequency of source 2 in freq_data
    freq_data.frequency(freq_tracker) = source_freq2;
    
    % =========================================================================
    % FLUID SIMULATION
    % =========================================================================    
    
    %Set medium properties for fluid simulation
    clear medium;
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
    [source.p] = two_transducers_fluid(Nx,Ny,source.p_mask,...
        source_freq1,source_freq2,source_mag1,source_mag2,kgrid.t_array,kgrid,medium );
    
    % create a display mask to display the transducer for fluid simulation
    display_mask = source.p_mask;

    % create a sensor mask covering the entire computational domain using the
    % opposing corners of a rectangle
    clear sensor;
    sensor.mask = [1, 1, Nx, Ny].';

    % set the record mode capture the final wave-field and the statistics at
    % each sensor point 
    sensor.record = {'p_final', 'p_max', 'p_rms','u_max_all'};
    max_distance = 0.05 + dx*sqrt((Nx/2)^2+(Ny/2)^2);
    travel_time = max_distance /vsound;
    i = 1;
    while kgrid.t_array(i) < travel_time
        i = i+1;
    end
    sensor.record_start_index = i;

    % assign the input options for fluid simulation
    input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false,'Smooth', [true, true, true]};

    % run the fluid simulation
    sensor_data_fluid = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
     if plot_type == 0 || plot_type == 2%rms plot
        %Find location of peak RMS intensity for fluid simulation
        [freq_data.x_max_fluid(freq_tracker),freq_data.y_max_fluid(freq_tracker), max_rms_p] =...
            max_coords(sensor_data_fluid.p_rms,Nx,Ny);
        
        %Find FWHM (making assumption that the profile is approximately circular)
        [FWHM_matrix,freq_data.HM_x_max_fluid(freq_tracker), freq_data.HM_y_max_fluid(freq_tracker),...
            freq_data.FWHM_fluid(freq_tracker)] = FWHM_calc(sensor_data_fluid.p_rms,Nx,Ny,max_rms_p,dx);
        
        %Add masks to fluid data for plotting purposes
        sensor_data_fluid.p_rms(source.p_mask ~= 0) = 1;
        sensor_data_fluid.p_rms(FWHM_matrix == 1) = 1;
        sensor_data_fluid.p_rms(contour_img == 1) = 1;
         
     elseif plot_type == 1 || plot_type == 2 %max pressure plot
         %Find location of peak intensity for fluid simulation
         [freq_data.x_max_fluid(freq_tracker),freq_data.y_max_fluid(freq_tracker), max_p] =...
            max_coords(sensor_data_fluid.p_max,Nx,Ny);
         
        %Find FWHM (making assumption that the profile is approximately circular)
        [FWHM_matrix,freq_data.HM_x_max_fluid(freq_tracker), freq_data.HM_y_max_fluid(freq_tracker),...
            freq_data.FWHM_fluid(freq_tracker)] = FWHM_calc(sensor_data_fluid.p_max,Nx,Ny,max_p,dx);
        
        %Add masks to fluid data for plotting purposes
        sensor_data_fluid.p_max(source.p_mask ~= 0) = 1;
        sensor_data_fluid.p_max(FWHM_matrix == 1) = 1;
        sensor_data_fluid.p_max(contour_img == 1) = 1;
     end
    

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
    
    % assign the source for the elastic simulation
    clear source;
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


    % Here I define the frequencies of my source points for the elastic
    % simulation
    source.s_mode = 'additive';
    %[source.sxx, source.sxy, source.syy] = two_transducers_elastic(Nx,Ny,source.s_mask,...
        %source_freq1,source_freq2,source_mag1,source_mag2,kgrid.t_array,kgrid,medium );
   [source.sxx, source.syy] = two_transducers_elastic(Nx,Ny,source.s_mask,...
        source_freq1,source_freq2,source_mag1,source_mag2,kgrid.t_array,kgrid,medium );
    
    % create a display mask to display the transducer for the elastic
    % simulation
    display_mask = source.s_mask;

    % assign the input options
    input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false, 'Smooth', [true, true, true]};

    % run the fluid simulation
    sensor_data_elastic = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});
    
    if plot_type == 0 || plot_type == 2 % plot RMS pressure data
        % Find location of peak intensity
        [freq_data.x_max_elastic(freq_tracker),freq_data.y_max_elastic(freq_tracker), max_rms_p] =...
            max_coords(sensor_data_elastic.p_rms,Nx,Ny);
    
        % Find FWHM (making assumption that the profile is approximately circular)
        [FWHM_matrix,freq_data.HM_x_max_elastic(freq_tracker), freq_data.HM_y_max_elastic(freq_tracker),...
            freq_data.FWHM_elastic(freq_tracker)] = FWHM_calc(sensor_data_elastic.p_rms,Nx,Ny,max_rms_p,dx);
    
        % Add masks to data for plots
        sensor_data_elastic.p_rms(source.s_mask ~= 0) = 1;
        sensor_data_elastic.p_rms(FWHM_matrix == 1) = 1;
        sensor_data_elastic.p_rms(contour_img == 1) = 1;
    elseif plot_type == 1 || plot_type == 2% plot max pressure data
        % Find location of peak intensity
        [freq_data.x_max_elastic(freq_tracker),freq_data.y_max_elastic(freq_tracker), max_p] =...
            max_coords(sensor_data_elastic.p_max,Nx,Ny);
    
        %Find FWHM (making assumption that the profile is approximately circular)
        [FWHM_matrix,freq_data.HM_x_max_elastic(freq_tracker), freq_data.HM_y_max_elastic(freq_tracker),...
            freq_data.FWHM_elastic(freq_tracker)] = FWHM_calc(sensor_data_elastic.p_max,Nx,Ny,max_p,dx);
    
        %Add masks to data for plots
        sensor_data_elastic.p_max(source.s_mask ~= 0) = 1;
        sensor_data_elastic.p_max(FWHM_matrix == 1) = 1;
        sensor_data_elastic.p_max(contour_img == 1) = 1;
    end

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

end

% calculate RMS and max pressures within the spinal canal
[theta_data.spinalcord_elastic_RMS_p(theta_tracker),theta_data.spinalcord_elastic_Max_p(theta_tracker),...
   theta_data.spinalcord_elastic_RMS_pnorm(theta_tracker),theta_data.spinalcord_elastic_Max_pnorm(theta_tracker) ] ...
= spinaldata(Nx,Ny,sensor_data_elastic.p_max,sensor_data_elastic.p_rms,medium_img);

[theta_data.spinalcord_fluid_RMS_p(theta_tracker),theta_data.spinalcord_fluid_Max_p(theta_tracker),...
   theta_data.spinalcord_fluid_RMS_pnorm(theta_tracker),theta_data.spinalcord_fluid_Max_pnorm(theta_tracker) ] ...
= spinaldata(Nx,Ny,sensor_data_fluid.p_max,sensor_data_fluid.p_rms,medium_img);

theta_tracker = theta_tracker+1;

end

% =========================================================================
% VISUALISATION
% =========================================================================

figure;
subplot(2,2,1);
plot(2*theta_data.angle,theta_data.spinalcord_fluid_RMS_p);
title('Fluid RMS pressure in spinal cord');
xlabel('Angle subtended [rad]');
ylabel('Pressure [Pa]');
subplot(2,2,2);
plot(2*theta_data.angle,theta_data.spinalcord_fluid_Max_p);
title('Fluid Max pressure in spinal cord');
xlabel('Angle subtended [rad]');
ylabel('Pressure [Pa]');
subplot(2,2,3);
plot(2*theta_data.angle,theta_data.spinalcord_elastic_RMS_p);
title('Elastic RMS pressure in spinal cord');
xlabel('Angle subtended [rad]');
ylabel('Pressure [Pa]');
subplot(2,2,4);
plot(2*theta_data.angle,theta_data.spinalcord_elastic_Max_p);
title('Elastic Max pressure in spinal cord');
xlabel('Angle subtended [rad]');
ylabel('Pressure [Pa]');

figure;
subplot(2,2,1);
plot(2*theta_data.angle,theta_data.spinalcord_fluid_RMS_pnorm);
title('Fluid RMS pressure in spinal cord');
xlabel('Angle subtended [rad]');
ylabel('Normalized Pressure');
subplot(2,2,2);
plot(2*theta_data.angle,theta_data.spinalcord_fluid_Max_pnorm);
title('Fluid Max pressure in spinal cord');
xlabel('Angle subtended [rad]');
ylabel('Normalized Pressure');
subplot(2,2,3);
plot(2*theta_data.angle,theta_data.spinalcord_elastic_RMS_pnorm);
title('Elastic RMS pressure in spinal cord');
xlabel('Angle subtended [rad]');
ylabel('Normalized Pressure');
subplot(2,2,4);
plot(2*theta_data.angle,theta_data.spinalcord_elastic_Max_pnorm);
title('Elastic Max pressure in spinal cord');
xlabel('Angle subtended [rad]');
ylabel('Normalized Pressure');

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