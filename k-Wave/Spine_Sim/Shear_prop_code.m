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
% create the computational grid
Nx = 216;           % number of grid points in the x (row) direction
Ny = 216;           % number of grid points in the y (column) direction
PML_size = 20;
x_length = 100e-3; 
dx = x_length/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% create the time array
cfl   = 0.1;
t_end = 12e-5;
kgrid.t_array= makeTime(kgrid, 1580, cfl, t_end);

% define two curved transducer elements
source.p_mask = makedoubleSemiCircle(Nx, Ny,round((Nx+40)/2),round(Ny/2),...
    round(Nx/2), pi/12, pi/3+pi/12,round((Nx+40)/2),round(Ny/2),round(Nx/2), pi/12+pi/2, pi/3+pi/12+pi/2);


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

% set the input arguments
%input_args = {'PMLSize', PML_size, 'PMLAlpha', 2, 'PlotPML', false, ...
%    'PMLInside', false, 'PlotScale', [-1, 1]*source_strength, ...
%    'DisplayMask', 'off', 'DataCast', 'single'};

% =========================================================================
% FLUID SIMULATION
% =========================================================================

% define the medium properties
medium.sound_speed = makemediaDisk(Nx,Ny,round(Nx/2),round(Ny/2),round(Nx/10),1580,2820);         % [m/s]
medium.density = makedensityDisc(Nx,Ny,round(Nx/2),round(Ny/2),round(Nx/10),1000,1800);         % [kg/m^3]
medium.alpha_coeff = makemediaDisk(Nx,Ny,round(Nx/2),round(Ny/2),round(Nx/10),0.57,9);   % [dB/(MHz^y cm)]
medium.alpha_power = 2;
   
% define struct of frequencies, max rms pressure locations
freq_data = struct('frequency',[],'x_max',[],'y_max',[],'HM_x_max',[],'HM_y_max',[],'FWHM',[]);

step_tracker = 1;

for source_freq2 = initial:step_size:final 
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
    sensor.record = {'p_final', 'p_max', 'p_rms','u_max_all'};

% assign the input options
    input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

    % run the fluid simulation
    sensor_data_fluid = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

end


% =========================================================================
% ELASTIC SIMULATION
% =========================================================================

% define the medium properties
clear medium
% define the medium properties
medium.sound_speed_compression = makemediaDisk(Nx,Ny,round(Nx/2),round(Ny/2),round(Nx/10),1580,2820);         % [m/s]
medium.sound_speed_shear = makemediaDisk(Nx,Ny,round(Nx/2),round(Ny/2),round(Nx/10),0,1500);         % [m/s]
medium.density = makedensityDisc(Nx,Ny,round(Nx/2),round(Ny/2),round(Nx/10),1000,1800);         % [kg/m^3]
medium.alpha_coeff_compression = makemediaDisk(Nx,Ny,round(Nx/2),round(Ny/2),round(Nx/10),0.57,9);   % [dB/(MHz^y cm)]
medium.alpha_coeff_shear = makemediaDisk(Nx,Ny,round(Nx/2),round(Ny/2),round(Nx/10),0,20);   % [dB/(MHz^y cm)]

% assign the source
clear source
source_mask = makedoubleSemiCircle(Nx, Ny,round((Nx+40)/2),round(Ny/2),...
    round(Nx/2), pi/12, pi/3+pi/12,round((Nx+40)/2),round(Ny/2),round(Nx/2), pi/12+pi/2, pi/3+pi/12+pi/2);
source.s_mask = makedoubleSemiCircle(Nx, Ny,round((Nx+40)/2),round(Ny/2),...
    round(Nx/2), pi/12, pi/3+pi/12,round((Nx+40)/2),round(Ny/2),round(Nx/2), pi/12+pi/2, pi/3+pi/12+pi/2);
%source.s_mask = source_mask;

clear sensor;
step_tracker = 1;

for source_freq2 = initial:step_size:final 
    freq_data.frequency(step_tracker) = source_freq2;

    %Here I define the frequencies of my source points
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
% create a display mask to display the transducer
    display_mask = source.s_mask;

% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
    sensor.mask = [1, 1, Nx, Ny].';

% set the record mode capture the final wave-field and the statistics at
% each sensor point 
    sensor.record = {'p_final', 'p_max', 'p_rms','u_max_all'};

% assign the input options
    input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

    % run the fluid simulation
    sensor_data_elastic = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

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