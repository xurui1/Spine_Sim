function sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, varargin)
%KSPACEFIRSTORDER1D     1D time-domain simulation of wave propagation.
%
% DESCRIPTION:
%       kspaceFirstOrder1D simulates the time-domain propagation of
%       compressional waves through a one-dimensional homogeneous or
%       heterogeneous acoustic medium given four input structures: kgrid,
%       medium, source, and sensor. The computation is based on a
%       first-order k-space model which accounts for power law absorption
%       and a heterogeneous sound speed and density. If medium.BonA is
%       specified, cumulative nonlinear effects are also modelled. At each
%       time-step (defined by kgrid.t_array), the acoustic field parameters
%       at the positions defined by sensor.mask are recorded and stored. If
%       kgrid.t_array is set to 'auto', this array is automatically
%       generated using makeTime. An anisotropic absorbing boundary layer
%       called a perfectly matched layer (PML) is implemented to prevent
%       waves that leave one side of the domain being reintroduced from the
%       opposite side (a consequence of using the FFT to compute the
%       spatial derivatives in the wave equation). This allows infinite
%       domain simulations to be computed using small computational grids.
%        
%       For a homogeneous medium the formulation is exact and the
%       time-steps are only limited by the effectiveness of the perfectly
%       matched layer. For a heterogeneous medium, the solution represents
%       a leap-frog pseudospectral method with a k-space correction that
%       improves the accuracy of computing the temporal derivatives. This
%       allows larger time-steps to be taken for the same level of accuracy
%       compared to conventional pseudospectral time-domain methods. The
%       computational grids are staggered both spatially and temporally.
%        
%       An initial pressure distribution can be specified by assigning a
%       matrix (the same size as the computational grid) of arbitrary
%       numeric values to source.p0. A time varying pressure source can
%       similarly be specified by assigning a binary matrix (i.e., a matrix
%       of 1's and 0's with the same dimensions as the computational grid)
%       to source.p_mask where the 1's represent the grid points that form
%       part of the source. The time varying input signals are then
%       assigned to source.p. This can be a single time series (in which
%       case it is applied to all source elements), or a matrix of time
%       series following the source elements using MATLAB's standard
%       column-wise linear matrix index ordering. A time varying velocity
%       source can be specified in an analogous fashion, where the source
%       location is specified by source.u_mask, and the time varying input
%       velocity is assigned to source.ux.             
%
%       The field values are returned as arrays of time series at the
%       sensor locations defined by sensor.mask. This can be defined in
%       three different ways. (1) As a binary matrix (i.e., a matrix of 1's
%       and 0's with the same dimensions as the computational grid)
%       representing the grid points within the computational grid that
%       will collect the data. (2) As the grid coordinates of two opposing
%       ends of a line in the form [x1; x2]. This is equivalent to using a
%       binary sensor mask covering the same region, however, the output is
%       indexed differently as discussed below. (3) As a series of
%       Cartesian coordinates within the grid which specify the location of
%       the pressure values stored at each time step. If the Cartesian
%       coordinates don't exactly match the coordinates of a grid point,
%       the output values are calculated via interpolation. The Cartesian
%       points must be given as a 1 by N matrix corresponding to the x
%       positions, where the Cartesian origin is assumed to be in the
%       center of the grid. If no output is required, the sensor input can
%       be replaced with an empty array [].                 
%
%       If sensor.mask is given as a set of Cartesian coordinates, the
%       computed sensor_data is returned in the same order. If sensor.mask
%       is given as a binary matrix, sensor_data is returned using MATLAB's
%       standard column-wise linear matrix index ordering. In both cases,
%       the recorded data is indexed as sensor_data(sensor_point_index,
%       time_index). For a binary sensor mask, the field values at a
%       particular time can be restored to the sensor positions within the
%       computation grid using unmaskSensorData. If sensor.mask is given as
%       a list of opposing ends of a line, the recorded data is indexed as
%       sensor_data(line_index).p(x_index, time_index), where x_index
%       corresponds to the grid index within the line, and line_index
%       corresponds to the number of lines if more than one is specified.
%        
%       By default, the recorded acoustic pressure field is passed directly
%       to the output sensor_data. However, other acoustic parameters can
%       also be recorded by setting sensor.record to a cell array of the
%       form {'p', 'u', 'p_max', ...}. For example, both the particle
%       velocity and the acoustic pressure can be returned by setting
%       sensor.record = {'p', 'u'}. If sensor.record is given, the output
%       sensor_data is returned as a structure with the different outputs
%       appended as structure fields. For example, if sensor.record = {'p',
%       'p_final', 'p_max', 'u'}, the output would contain fields
%       sensor_data.p, sensor_data.p_final, sensor_data.p_max, and
%       sensor_data.ux. Most of the output parameters are recorded at the
%       given sensor positions and are indexed as
%       sensor_data.field(sensor_point_index, time_index) or
%       sensor_data(line_index).field(x_index, time_index) if using a
%       sensor mask defined as opposing ends of a line. The exceptions are
%       the averaged quantities ('p_max', 'p_rms', 'u_max', 'p_rms',
%       'I_avg'), the 'all' quantities ('p_max_all', 'p_min_all',
%       'u_max_all', 'u_min_all'), and the final quantities ('p_final',
%       'u_final'). The averaged quantities are indexed as
%       sensor_data.p_max(sensor_point_index) or
%       sensor_data(line_index).p_max(x_index) if using line ends, while
%       the final and 'all' quantities are returned over the entire grid
%       and are always indexed as sensor_data.p_final(n), regardless of the
%       type of sensor mask.                        
%
%       kspaceFirstOrder1D may also be used for time reversal image
%       reconstruction by assigning the time varying pressure recorded over
%       an arbitrary sensor surface to the input field
%       sensor.time_reversal_boundary_data. This data is then enforced in
%       time reversed order as a time varying Dirichlet boundary condition
%       over the sensor surface given by sensor.mask. The boundary data
%       must be indexed as
%       sensor.time_reversal_boundary_data(sensor_point_index, time_index).
%       If sensor.mask is given as a set of Cartesian coordinates, the
%       boundary data must be given in the same order. An equivalent binary
%       sensor mask (computed using nearest neighbour interpolation) is
%       then used to place the pressure values into the computational grid
%       at each time step. If sensor.mask is given as a binary matrix of
%       sensor points, the boundary data must be ordered using MATLAB's
%       standard column-wise linear matrix indexing. If no additional
%       inputs are required, the source input can be replaced with an empty
%       array [].                
%
%       Acoustic attenuation compensation can also be included during time
%       reversal image reconstruction by assigning the absorption
%       parameters medium.alpha_coeff and medium.alpha_power and reversing
%       the sign of the absorption term by setting medium.alpha_sign = [-1,
%       1]. This forces the propagating waves to grow according to the
%       absorption parameters instead of decay. The reconstruction should
%       then be regularised by assigning a filter to medium.alpha_filter
%       (this can be created using getAlphaFilter).         
%
%       Note: To run a simple photoacoustic image reconstruction example
%       using time reversal (that commits the 'inverse crime' of using the
%       same numerical parameters and model for data simulation and image
%       reconstruction), the sensor_data returned from a k-Wave simulation
%       can be passed directly to sensor.time_reversal_boundary_data with
%       the input fields source.p0 and source.p removed or set to zero.
%        
% USAGE:
%       sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor)
%       sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, ...) 
%
% INPUTS:
% The minimum fields that must be assigned to run an initial value problem
% (for example, a photoacoustic forward simulation) are marked with a *. 
%
%       kgrid*              - k-space grid structure returned by makeGrid
%                             containing Cartesian and k-space grid fields  
%       kgrid.t_array*      - evenly spaced array of time values [s] (set
%                             to 'auto' by makeGrid) 
%
%
%       medium.sound_speed* - sound speed distribution within the acoustic
%                             medium [m/s] 
%       medium.sound_speed_ref - reference sound speed used within the
%                             k-space operator (phase correction term)
%                             [m/s]
%       medium.density*     - density distribution within the acoustic
%                             medium [kg/m^3] 
%       medium.BonA         - parameter of nonlinearity
%       medium.alpha_power  - power law absorption exponent
%       medium.alpha_coeff  - power law absorption coefficient 
%                             [dB/(MHz^y cm)] 
%       medium.alpha_mode   - optional input to force either the absorption
%                             or dispersion terms in the equation of state
%                             to be excluded; valid inputs are
%                             'no_absorption' or 'no_dispersion' 
%       medium.alpha_filter - frequency domain filter applied to the
%                             absorption and dispersion terms in the
%                             equation of state 
%       medium.alpha_sign   - two element array used to control the sign of
%                             absorption and dispersion terms in the
%                             equation of state  
%
%
%       source.p0*          - initial pressure within the acoustic medium
%       source.p            - time varying pressure at each of the source
%                             positions given by source.p_mask 
%       source.p_mask       - binary matrix specifying the positions of the
%                             time varying pressure source distribution
%       source.p_mode       - optional input to control whether the input
%                             pressure is injected as a mass source or
%                             enforced as a dirichlet boundary condition;
%                             valid inputs are 'additive' (the default) or
%                             'dirichlet'    
%       source.ux           - time varying particle velocity in the
%                             x-direction at each of the source positions
%                             given by source.u_mask 
%       source.u_mask       - binary matrix specifying the positions of the
%                             time varying particle velocity distribution 
%       source.u_mode       - optional input to control whether the input
%                             velocity is applied as a force source or
%                             enforced as a dirichlet boundary condition;
%                             valid inputs are 'additive' (the default) or
%                             'dirichlet'
%
%
%       sensor.mask*        - binary matrix or a set of Cartesian points
%                             where the pressure is recorded at each
%                             time-step  
%       sensor.record       - cell array of the acoustic parameters to
%                             record in the form sensor.record = {'p', 'u',
%                             ...}; valid inputs are:  
%                               'p' (acoustic pressure)
%                               'p_max' (maximum pressure)
%                               'p_min' (minimum pressure)
%                               'p_rms' (RMS pressure)
%                               'p_final' (final pressure field at all grid points)
%                               'p_max_all' (maximum pressure at all grid points)
%                               'p_min_all' (minimum pressure at all grid points)
%                               'u' (particle velocity)
%                               'u_max' (maximum particle velocity)
%                               'u_min' (minimum particle velocity)
%                               'u_rms' (RMS particle velocity)
%                               'u_final' (final particle velocity field at all grid points)
%                               'u_max_all' (maximum particle velocity at all grid points)
%                               'u_min_all' (minimum particle velocity at all grid points)
%                               'u_non_staggered' (particle velocity on non-staggered grid)
%                               'I' (time varying acoustic intensity)
%                               'I_avg' (average acoustic intensity) 
%       sensor.record_start_index - time index at which the sensor should
%                             start recording the data specified by
%                             sensor.record (default = 1) 
%       sensor.time_reversal_boundary_data - time varying pressure
%                             enforced as a Dirichlet boundary condition
%                             over sensor.mask  
%       sensor.frequency_response - two element array specifying the center
%                             frequency and percentage bandwidth of a
%                             frequency domain Gaussian filter applied to
%                             the sensor_data
%
% Note: For heterogeneous medium parameters, medium.sound_speed and
% medium.density must be given in matrix form with the same dimensions as
% kgrid. For homogeneous medium parameters, these can be given as single
% numeric values. If the medium is homogeneous and velocity inputs or
% outputs are not required, it is not necessary to specify medium.density.
%
% OPTIONAL INPUTS:
%       Optional 'string', value pairs that may be used to modify the
%       default computational settings.
%
%       'CartInterp'- Interpolation mode used to extract the pressure when
%                     a Cartesian sensor mask is given. If set to 'nearest'
%                     and more than one Cartesian point maps to the same
%                     grid point, duplicated data points are discarded and
%                     sensor_data will be returned with less points than
%                     that specified by sensor.mask.
%       'CreateLog' - Boolean controlling whether the command line output
%                     is saved using the diary function with a date and
%                     time stamped filename (default = false). 
%       'DataCast'  - String input of the data type that variables are cast
%                     to before computation. For example, setting to
%                     'single' will speed up the computation time (due to
%                     the improved efficiency of fftn and ifftn for this
%                     data type) at the expense of a loss in precision.
%                     This variable is also useful for utilising GPU
%                     parallelisation through libraries such as the
%                     Parallel Computing Toolbox by setting 'DataCast' to
%                     'gpuArray-single' (default = 'off'). 
%       'DataRecast'- Boolean controlling whether the output data is cast
%                     back to double precision. If set to false,
%                     sensor_data will be returned in the data format set
%                     using the 'DataCast' option.
%       'DisplayMask' - Binary matrix overlayed onto the animated
%                     simulation display. Elements set to 1 within the
%                     display mask are set to black within the display
%                     (default = sensor.mask).
%       'LogScale'  - Boolean controlling whether the pressure field is log
%                     compressed before display (default = false). The data
%                     is compressed by scaling both the positive and
%                     negative values between 0 and 1 (truncating the data
%                     to the given plot scale), adding a scalar value
%                     (compression factor) and then using the corresponding
%                     portion of a log10 plot for the compression (the
%                     negative parts are remapped to be negative thus the
%                     default color scale will appear unchanged). The
%                     amount of compression can be controlled by adjusting
%                     the compression factor which can be given in place of
%                     the Boolean input. The closer the compression factor
%                     is to zero, the steeper the corresponding part of the
%                     log10 plot used, and the greater the compression (the
%                     default compression factor is 0.02).
%       'MovieArgs' - Settings for movie2avi. Parameters must be given as
%                     {param, value, ...} pairs within a cell array
%                     (default = {}).
%       'MovieName' - Name of the movie produced when 'RecordMovie' is set
%                     to true (default = 'date-time-kspaceFirstOrder2D').
%       'PlotFreq'  - The number of iterations which must pass before the
%                     simulation plot is updated (default = 10).
%       'PlotLayout'- Boolean controlling whether a four panel plot of the
%                     initial simulation layout is produced (initial
%                     pressure, sensor mask, sound speed, density)
%                     (default = false).
%       'PlotPML'   - Boolean controlling whether the perfectly matched
%                     layer is shown in the simulation plots. If set to
%                     false, the PML is not displayed (default = true).
%       'PlotScale' - [min, max] values used to control the scaling for
%                     imagesc (visualisation). If set to 'auto', a
%                     symmetric plot scale is chosen automatically for each
%                     plot frame.
%       'PlotSim'   - Boolean controlling whether the simulation iterations
%                     are progressively plotted (default = true).
%       'PMLAlpha'  - Absorption within the perfectly matched layer in
%                     Nepers per grid point (default = 2).
%       'PMLInside' - Boolean controlling whether the perfectly matched
%                     layer is inside or outside the grid. If set to false,
%                     the input grids are enlarged by PMLSize before
%                     running the simulation (default = true). 
%       'PMLSize'   - Size of the perfectly matched layer in grid points.
%                     To remove the PML, set the appropriate PMLAlpha to
%                     zero rather than forcing the PML to be of zero size
%                     (default = 20). 
%       'RecordMovie' - Boolean controlling whether the displayed image
%                     frames are captured and stored as a movie using
%                     movie2avi (default = false).  
%       'Smooth'    - Boolean controlling whether source.p0,
%                     medium.sound_speed, and medium.density are smoothed
%                     using smooth before computation. 'Smooth' can either
%                     be given as a single Boolean value or as a 3 element
%                     array to control the smoothing of source.p0,
%                     medium.sound_speed, and medium.density,
%                     independently.  
%
% OUTPUTS:
% If sensor.record is not defined by the user:
%       sensor_data - time varying pressure recorded at the sensor
%                     positions given by sensor.mask
%
% If sensor.record is defined by the user:
%       sensor_data.p         - time varying pressure recorded at the
%                               sensor positions given by sensor.mask
%                               (returned if 'p' is set)  
%       sensor_data.p_max     - maximum pressure recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'p_max' is set)  
%       sensor_data.p_min     - minimum pressure recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'p_min' is set)  
%       sensor_data.p_rms     - rms of the time varying pressure recorded
%                               at the sensor positions given by
%                               sensor.mask (returned if 'p_rms' is set)  
%       sensor_data.p_final   - final pressure field at all grid points
%                               within the domain (returned if 'p_final' is
%                               set)
%       sensor_data.p_max_all - maximum pressure recorded at all grid
%                               points within the domain (returned if
%                               'p_max_all' is set) 
%       sensor_data.p_min_all - minimum pressure recorded at all grid
%                               points within the domain (returned if
%                               'p_min_all' is set)  
%       sensor_data.ux        - time varying particle velocity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u' is set)     
%       sensor_data.ux_max    - maximum particle velocity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_max' is set)   
%       sensor_data.ux_min    - minimum particle velocity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_min' is set)   
%       sensor_data.ux_rms    - rms of the time varying particle velocity
%                               in the x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'u_rms' is set)     
%       sensor_data.ux_final  - final particle velocity field in the
%                               x-direction at all grid points within the
%                               domain (returned if 'u_final' is set) 
%       sensor_data.ux_max_all- maximum particle velocity in the
%                               x-direction recorded at all grid points
%                               within the domain (returned if 'u_max_all'  
%       sensor_data.ux_min_all- minimum particle velocity in the
%                               x-direction recorded at all grid points
%                               within the domain (returned if 'u_min_all'
%                               is set)   
%       sensor_data.ux_non_staggered - time varying particle velocity in
%                               the x-direction recorded at the sensor
%                               positions given by sensor.mask after
%                               shifting to the non-staggered grid
%                               (returned if 'u_non_staggered' is set) 
%       sensor_data.Ix        - time varying acoustic intensity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I' is set) 
%       sensor_data.Ix_avg    - average acoustic intensity in the
%                               x-direction recorded at the sensor
%                               positions given by sensor.mask (returned if
%                               'I_avg' is set)   
%
% ABOUT:
%       author      - Bradley Treeby and Ben Cox
%       date        - 22nd April 2009
%       last update - 27th August 2014
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also fft, ifft, getframe, kspaceFirstOrder2D, kspaceFirstOrder3D,
% makeGrid, makeTime, movie2avi, smooth, unmaskSensorData 

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% suppress mlint warnings that arise from using sub-scripts
%#ok<*NASGU>
%#ok<*COLND>
%#ok<*NODEF>
%#ok<*INUSL>

% =========================================================================
% CHECK INPUT STRUCTURES AND OPTIONAL INPUTS
% =========================================================================

% start the timer and store the start time
start_time = clock;
tic;

% set the name of the simulation code
MFILE = mfilename;

% define literals
kspaceFirstOrder_setDefaults;

% run subscript to check inputs
kspaceFirstOrder_inputChecking;

% =========================================================================
% CALCULATE MEDIUM PROPERTIES ON STAGGERED GRID
% =========================================================================

% interpolate the values of the density at the staggered grid locations
% where sgx = (x + dx/2)
if numel(rho0) > 1 && use_sg
    
    % rho0 is heterogeneous and staggered grids are used
    rho0_sgx = interp1(kgrid.x, rho0, kgrid.x + kgrid.dx/2, '*linear');
    
    % set values outside of the interpolation range to original values 
    rho0_sgx(isnan(rho0_sgx)) = rho0(isnan(rho0_sgx)); 
    
else
    % rho0 is homogeneous or staggered grids are not used
    rho0_sgx = rho0;
end

% =========================================================================
% PREPARE DERIVATIVE AND PML OPERATORS
% =========================================================================

% get the PML operators based on the reference sound speed and PML settings
pml_x     = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, PML_x_size, PML_x_alpha, false, 1);
pml_x_sgx = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, PML_x_size, PML_x_alpha, true && use_sg, 1);

% define the k-space derivative operator
ddx_k = ifftshift(1i*kgrid.kx_vec);

% define the staggered grid shift operators (the option use_sg exists for
% debugging)
if use_sg
    shift_pos = ifftshift( exp(1i*kgrid.kx_vec*kgrid.dx/2) );
    shift_neg = ifftshift( exp(-1i*kgrid.kx_vec*kgrid.dx/2) );
else
    shift_pos = 1;
    shift_neg = 1;
end

% create k-space operator (the option use_kspace exists for debugging)
if use_kspace
    kappa = ifftshift( sinc(c_ref*dt*kgrid.k/2) );
else
    kappa = 1;
end

% =========================================================================
% DATA CASTING
% =========================================================================

% preallocate the loop variables using the castZeros anonymous function
% (this creates a matrix of zeros in the data type specified by data_cast)
p      = castZeros([kgrid.Nx, 1]);
rhox   = castZeros([kgrid.Nx, 1]);
ux_sgx = castZeros([kgrid.Nx, 1]);
p_k    = castZeros([kgrid.Nx, 1]);

% run subscript to cast the remaining loop variables to the data type
% specified by data_cast 
if ~strcmp(data_cast, 'off')
    kspaceFirstOrder_dataCast;
end

% =========================================================================
% CREATE INDEX VARIABLES
% =========================================================================

% setup the time index variable
if ~record.time_rev
    index_start = 1;
    index_step = 1;
    index_end = length(t_array);  
else
    % reverse the order of the input data
    sensor.time_reversal_boundary_data = fliplr(sensor.time_reversal_boundary_data);    
    index_start = 1;
    index_step = 1;

    % stop one time point before the end so the last points are not
    % propagated
    index_end = length(t_array) - 1;      
end

% =========================================================================
% PREPARE VISUALISATIONS
% =========================================================================

% pre-compute suitable axes scaling factor
if plot_layout || plot_sim
    [x_sc, scale, prefix] = scaleSI(max(kgrid.x));  %#ok<ASGLU>
end

% run subscript to plot the simulation layout if 'PlotLayout' is set to true
if plot_layout
    kspaceFirstOrder_plotLayout;
end

% initialise the figure used for animation if 'PlotSim' is set to 'true'
if plot_sim
    kspaceFirstOrder_initialiseFigureWindow;
end  

% initialise movie parameters if 'RecordMovie' is set to 'true'
if record_movie
    kspaceFirstOrder_initialiseMovieParameters;
end

% =========================================================================
% LOOP THROUGH TIME STEPS
% =========================================================================

% update command line status
disp(['  precomputation completed in ' scaleTime(toc)]);
disp('  starting time loop...');

% restart timing variables
loop_start_time = clock;
tic;

% start time loop
for t_index = index_start:index_step:index_end
    
    % enforce time reversal bounday condition
    if record.time_rev   
        
        % load pressure value and enforce as a Dirichlet boundary condition
        p(sensor_mask_index) = sensor.time_reversal_boundary_data(:, t_index);

        % update p_k
        p_k = fft(p);

        % compute rhox using an adiabatic equation of state
        rhox_mod = p./(c.^2);
        rhox(sensor_mask_index) = rhox_mod(sensor_mask_index);
        
    end   
    
    % calculate ux at the next time step using dp/dx at the current time
    % step
    if ~nonuniform_grid && ~use_finite_difference
        
        % calculate gradient using the k-space method on a regular grid
        ux_sgx = pml_x_sgx .* (  pml_x_sgx.*ux_sgx - dt./rho0_sgx .* real(ifft(ddx_k .* shift_pos .* kappa .* p_k))  );
        
    elseif use_finite_difference
        switch use_finite_difference
            case 2
                % calculate gradient using second-order accurate finite
                % difference scheme (including half step forward)
                dpdx = ([p(2:end); 0] - p) /kgrid.dx;    
                ux_sgx = pml_x_sgx .* (  pml_x_sgx.*ux_sgx - dt./rho0_sgx .* dpdx );
            case 4
                % calculate gradient using fourth-order accurate finite
                % difference scheme (including half step forward)
                dpdx = ([0; p(1:(end-1))] - 27*p + 27*[p(2:end); 0] - [p(3:end); 0; 0])/(24*kgrid.dx);
                ux_sgx = pml_x_sgx .* (  pml_x_sgx.*ux_sgx - dt./rho0_sgx .* dpdx );
        end 
    else
        
        % calculate gradient using the k-space method on a non-uniform grid
        % via the mapped pseudospectral method         
        ux_sgx = pml_x_sgx .* (  pml_x_sgx.*ux_sgx - dt./rho0_sgx .* kgrid.dxudxn_sgx .* real(ifft(ddx_k .* shift_pos .* kappa .* p_k))  );
        
    end    
    
    % add in the velocity source term
    if ux_source >= t_index
        if strcmp(source.u_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            ux_sgx(u_source_index) = source.ux(:, t_index);
        else
            % add the source values to the existing field values        
            ux_sgx(u_source_index) = ux_sgx(u_source_index) + source.ux(:, t_index);
        end
    end    

    % calculate du/dx at the next time step
    if ~nonuniform_grid && ~use_finite_difference
        
        % calculate gradient using the k-space method on a regular grid  
        duxdx = real(ifft(ddx_k .* shift_neg .* kappa .* fft(ux_sgx)));   

    elseif use_finite_difference
        switch use_finite_difference
            case 2
                % calculate gradient using second-order accurate finite
                % difference scheme (including half step backward)
                duxdx = (ux_sgx - [0; ux_sgx(1:end-1)])/kgrid.dx;
            case 4
                % calculate gradient using fourth-order accurate finite
                % difference scheme (including half step backward)
                duxdx = ([0; 0; ux_sgx(1:(end-2))] - 27*[0; ux_sgx(1:(end-1))] + 27*ux_sgx -[ux_sgx(2:end); 0])/(24*kgrid.dx);

        end 
    else        
        
        % calculate gradients using a non-uniform grid via the mapped
        % pseudospectral method 
        duxdx = kgrid.dxudxn .* real(ifft(ddx_k .* shift_neg .* kappa .* fft(ux_sgx)));
        
    end    
    
    % calculate rhox at the next time step
    if ~nonlinear
        % use linearised mass conservation equation
        rhox = pml_x .* ( pml_x .* rhox - dt.*rho0 .* duxdx );
    else
        % use nonlinear mass conservation equation (implicit calculation)
        rhox = pml_x .* ( ( pml_x .* rhox - dt.*rho0 .* duxdx) ./ (1 + 2*dt.*duxdx) );
    end      
    
    % add in the pre-scaled pressure source term as a mass source  
    if p_source >= t_index
        if strcmp(source.p_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            rhox(p_source_index) = source.p(:, t_index);
        else
            % add the source values to the existing field values
            rhox(p_source_index) = rhox(p_source_index) + source.p(:, t_index);
        end
    end    
    
    % equation of state
    if ~nonlinear
        switch equation_of_state
            case 'lossless'
                % compute p using an adiabatic equation of state
                p = c.^2.*rhox;
            case 'absorbing'
                % compute p using an absorbing equation of state 
                p = c.^2.*(rhox ...
                     + absorb_tau.*real(ifft( absorb_nabla1.*fft(rho0.*duxdx) )) ...
                     - absorb_eta.*real(ifft( absorb_nabla2.*fft(rhox) )) ...
                     );
        end
    else
        switch equation_of_state
            case 'lossless'
                % compute p using a nonlinear adiabatic equation of state
                p = c.^2.*( rhox + medium.BonA.*rhox.^2./(2*rho0) );
            case 'absorbing'
                % compute p using a nonlinear absorbing equation of state 
                p = c.^2.*( rhox ...
                    + absorb_tau.*real(ifft( absorb_nabla1.*fft(rho0.*duxdx) )) ...
                    - absorb_eta.*real(ifft( absorb_nabla2.*fft(rhox) )) ...                   
                    + medium.BonA.*rhox.^2./(2*rho0) ...
                    );
        end         
    end     
    
    % enforce initial conditions if source.p0 is defined instead of time
    % varying sources
    if t_index == 1 && isfield(source, 'p0')
    
        % add the initial pressure to rho as a mass source
        p = source.p0;
        rhox = source.p0./c.^2;
        
        % compute u(t = t1 - dt/2) based on u(dt/2) = -u(-dt/2) which
        % forces u(t = t1) = 0
        if ~use_finite_difference
            % calculate gradient using the k-space method on a regular grid
            ux_sgx = dt./rho0_sgx .* real(ifft(ddx_k .* shift_pos .* kappa .* fft(p))) / 2;
        else
            switch use_finite_difference
                case 2
                    % calculate gradient using second-order accurate finite
                    % difference scheme (including half step forward)
                    dpdx = ([p(2:end); 0] - p) /kgrid.dx;
                    ux_sgx = dt./rho0_sgx .* dpdx / 2;
                case 4
                    % calculate gradient using fourth-order accurate finite
                    % difference scheme (including half step backward)
                    dpdx = ([p(3:end); 0; 0] - 27*[p(2:end); 0] + 27*p - [0; p(1:(end-1))])/(24*kgrid.dx);
                    ux_sgx = dt./rho0_sgx .* dpdx / 2;
            end
        end
    end    
    
    % precompute fft of p here so p can be modified for visualisation
    p_k = fft(p);    
    
    % extract required sensor data from the pressure and particle velocity
    % fields if the number of time steps elapsed is greater than
    % sensor.record_start_index (defaults to 1) 
    if record.use_sensor && ~record.time_rev && (t_index >= sensor.record_start_index)
    
        % update index for data storage
        file_index = t_index - sensor.record_start_index + 1;
        
        % run sub-function to extract the required data from the acoustic
        % variables
        sensor_data = kspaceFirstOrder_extractSensorData(1, sensor_data, file_index, sensor_mask_index, record, p, ux_sgx, [], []);      
        
    end

    % estimate the time to run the simulation
    if t_index == ESTIMATE_SIM_TIME_STEPS
        disp(['  estimated simulation time ' scaleTime(etime(clock, loop_start_time)*index_end/t_index) '...']);
    end      
    
    % plot data if required
    if plot_sim && (rem(t_index, plot_freq) == 0 || t_index == 1 || t_index == index_end)  

        % update progress bar
        waitbar(t_index/length(t_array), pbar);
        drawnow;

        % ensure p is cast as a CPU variable and remove the PML from the
        % plot if required
        if strcmp(data_cast, 'gpuArray')
            p_plot = double(gather(p(x1:x2)));
        else
            p_plot = double(p(x1:x2));  
        end
                       
        % update plot
        if plot_scale_auto || plot_scale_log || t_index == 1
            
            % update plot scale if set to automatic or log
            if plot_scale_auto || plot_scale_log
                kspaceFirstOrder_adjustPlotScale;
            end
                        
            % replace entire plot
            if ~nonuniform_grid
                img_data = plot(kgrid.x(x1:x2)*scale, p_plot);
            else
                img_data = plot(kgrid.x_size*kgrid.xn(x1:x2)*scale, p_plot);
            end            
            
            % add display mask onto plot
            if ~(strcmp(display_mask, 'default') || strcmp(display_mask, 'off'))
                hold on;
                if ~nonuniform_grid
                    stairs(kgrid.x(x1:x2)*scale, display_mask(x1:x2).*(plot_scale(2) - plot_scale(1)) + plot_scale(1), 'k-');
                else
                    stairs(kgrid.xn(x1:x2)*scale, display_mask(x1:x2).*(plot_scale(2) - plot_scale(1)) + plot_scale(1), 'k-');
                end
                hold off
            end
            
            % set plot options
            xlabel(['x-position [' prefix 'm]']);
            if ~nonuniform_grid
                set(gca, 'YLim', plot_scale, 'XLim', kgrid.x([x1, x2])*scale);      
            else
                set(gca, 'YLim', plot_scale, 'XLim', kgrid.x_size*kgrid.xn([x1, x2])*scale);      
            end            
            
        else
            % just replace the y-data
            set(img_data, 'YData', p_plot);
        end
        
        % force plot update
        drawnow;
        
        % save movie frames if required
        if record_movie

            % set background color to white
            set(gcf, 'Color', [1 1 1]);

            % save the movie frame
            movie_frames(frame_index) = getframe(gcf); %#ok<AGROW>

            % update frame index
            frame_index  = frame_index  + 1;

        end
        
        % update variable used for timing variable to exclude the first
        % time step if plotting is enabled
        if t_index == 1
            loop_start_time = clock;
        end        
    end    
end

% assign the final time reversal values
if record.time_rev
    p(sensor_mask_index) = sensor.time_reversal_boundary_data(:, index_end + 1);
end

% update command line status
disp(['  simulation completed in ' scaleTime(toc)]);

% =========================================================================
% CLEAN UP
% =========================================================================

% clean up used figures
if plot_sim
    close(img);
    close(pbar);
    drawnow;
end

% save the movie frames to disk
if record_movie
    kspaceFirstOrder_saveMovieFile;   
end

% save the final pressure field if in time reversal mode
if record.time_rev
    record.p_final = true;
end

% save the final acoustic pressure if required
if record.p_final
    sensor_data.p_final = p(record.x1_inside:record.x2_inside);
end

% save the final particle velocity if required
if record.u_final
    sensor_data.ux_final = ux_sgx(record.x1_inside:record.x2_inside);    
end

% run subscript to cast variables back to double precision if required
if data_recast
    kspaceFirstOrder_dataRecast;
end

% run subscript to compute and save intensity values
if record.use_sensor && ~record.time_rev && (record.I || record.I_avg)
    save_intensity_matlab_code = true;    
    kspaceFirstOrder_saveIntensity;
end

% reorder the sensor points if a binary sensor mask was used for Cartesian
% sensor mask nearest neighbour interpolation (this is performed after
% recasting as the GPU toolboxes do not all support this subscript)
if record.use_sensor && record.reorder_data
    kspaceFirstOrder_reorderCartData;
end

% filter the recorded time domain pressure signals if transducer filter
% parameters are given 
if record.use_sensor && ~record.time_rev && isfield(sensor, 'frequency_response')
    sensor_data.p = gaussianFilter(sensor_data.p, 1/kgrid.dt, sensor.frequency_response(1), sensor.frequency_response(2));
end

% reorder the sensor points if cuboid corners is used
if record.cuboid_corners
    kspaceFirstOrder_reorderCuboidCorners;
end

if ~record.use_sensor
    % if sensor is not used, return empty sensor data
    sensor_data = [];
elseif record.time_rev
    % if computing time reversal, reassign sensor_data.p_final to
    % sensor_data
    sensor_data = sensor_data.p_final;
elseif ~isfield(sensor, 'record') && ~record.cuboid_corners
    % if sensor.record is not given by the user, and not using a cuboid
    % sensor mask, reassign sensor_data.p to sensor_data
    sensor_data = sensor_data.p;     
end

% update command line status
disp(['  total computation time ' scaleTime(etime(clock, start_time))]);

% switch off log
if create_log
    diary off;
end