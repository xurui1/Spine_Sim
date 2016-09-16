function circle = makedoubleSemiCircle(Nx, Ny, cx1, cy1, radius1, arc_angle1a, arc_angle2a, ...
                    cx2, cy2, radius2, arc_angle1b, arc_angle2b, plot_circle)
%MAKECIRCLE     Create a binary map of a circle within a 2D grid.
%
% DESCRIPTION:
%       makeCircle creates a binary map of a circle or arc (using the
%       midpoint circle algorithm) within a two-dimensional grid (the
%       circle position is denoted by 1's in the matrix with 0's
%       elsewhere). A single grid point is taken as the circle centre thus
%       the total diameter will always be an odd number of grid points.
%
% USAGE:
%       circle = makeCircle(Nx, Ny, cx, cy, radius)
%       circle = makeCircle(Nx, Ny, cx, cy, radius, arc_angle)
%       circle = makeCircle(Nx, Ny, cx, cy, radius, arc_angle, plot_circle)
%
% INPUTS:
%       Nx, Ny          - size of the 2D grid [grid points]
%       cx, cy          - centre of the circle [grid points], if set
%                         to 0, the centre of the grid is used
%       radius          - circle radius [grid points]
%
% OPTIONAL INPUTS:
%       arc_angle       - arc angle for incomplete circle [radians]
%                         (default = 2*pi)
%       plot_circle     - Boolean controlling whether the circle is plotted
%                         using imagesc (default = false)
%
% OUTPUTS:
%       circle          - 2D binary map of a circle
%
% ABOUT:
%       author          - Bradley Treeby
%       date            - 1st May 2009
%       last update     - 20th December 2011
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also makeCartCircle, makeDisc

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

% check for plot_circle input
if nargin < 13
    plot_circle = false;
end

% check for arc_angle input
%if nargin < 6
%    arc_angle2 = 2*pi;
%    arc_angle1 = 0;
%elseif arc_angle2 > 2*pi
%    arc_angle2 = 2*pi;
%elseif arc_angle2 < 0
%    arc_angle2 = 0;
%elseif arc_angle1 < 0
%    arc_angle1 = 0;
%elseif arc_angle1 > 2*pi
%    arc_angle1 = 2*pi;
%end

% force integer values
Nx = round(Nx);
Ny = round(Ny);
cx1 = round(cx1);
cy1 = round(cy1);
radius1 = round(radius1);
cx2 = round(cx2);
cy2 = round(cy2);
radius2 = round(radius2);

% check for zero values
if cx1 == 0
    cx1 = floor(Nx/2) + 1;
end
if cy1 == 0
    cy1 = floor(Ny/2) + 1;
end
if cx2 == 0
    cx2 = floor(Nx/2) + 1;
end
if cy2 == 0
    cy2 = floor(Ny/2) + 1;
end

% check the inputs
if cx1 < 1 || cx1 > Nx || cy1 < 1 || cy1 > Ny
    error('The center of circle 1 must be within the grid');
end
if cx1 < 2 || cx2 > Nx || cy2 < 1 || cy2 > Ny
    error('The center of circle 2 must be within the grid');
end

% define literals
MAGNITUDE = 1;

% create empty matrix
circle = zeros(Nx, Ny);

% initialise loop variables for semicircle1
x1 = 0;
y1 = radius1;
d1 = 1 - radius1;

% draw the first cardinal point
try 
    %circle(cx, cy - y) = MAGNITUDE;
catch
    error('The circle must fit within the grid');
end

% draw the remaining cardinal points
py1 = [cx1, cx1+y1, cx1-y1];
px1 = [cy1+y1, cy1, cy1];
for point_index = 1:length(py1)
    
    % check whether the point is within the arc made by arc_angle
    if (atan2(py1(point_index) - cx1, px1(point_index) - cy1) + pi) <= arc_angle2a
        %circle(py(point_index), px(point_index)) = MAGNITUDE;
    end
end

% loop through the remaining points using the midpoint circle algorithm
%start here tomorrow   
while ( x1 < y1 - 1)
    
    x1 = x1 + 1;
    if ( d1 < 0 ) 
        d1 = d1 + x1 + x1 + 1;
    else 
        y1 = y1 - 1;
        a1 = x1 - y1 + 1;
        d1 = d1 + a1 + a1;
    end
    
    % setup point indices
    py1 = [x1+cx1, y1+cx1, y1+cx1, x1+cx1, -x1+cx1, -y1+cx1, -y1+cx1, -x1+cx1];
    px1 = [y1+cy1, x1+cy1, -x1+cy1, -y1+cy1, -y1+cy1, -x1+cy1, x1+cy1, y1+cy1];
    
    % loop through each point
    for point_index = 1:length(py1)
        
        % check whether the point is within the arc made by arc_angle
        
        angle_val = atan2(py1(point_index) - cx1, px1(point_index) - cy1) + pi;
        if (arc_angle1a<= angle_val && angle_val <= arc_angle2a)
            circle(py1(point_index), px1(point_index)) = MAGNITUDE;
        end
    end
end

% initialise loop variables for semicircle2
x2 = 0;
y2 = radius2;
d2 = 1 - radius2;

% draw the first cardinal point
try 
    %circle(cx, cy - y) = MAGNITUDE;
catch
    error('The circle must fit within the grid');
end

% draw the remaining cardinal points
py2 = [cx2, cx2+y2, cx2-y2];
px2 = [cy2+y2, cy2, cy2];
for point_index = 1:length(py2)
    
    % check whether the point is within the arc made by arc_angle
    if (atan2(py2(point_index) - cx2, px2(point_index) - cy2) + pi) <= arc_angle2b
        %circle(py(point_index), px(point_index)) = MAGNITUDE;
    end
end

% loop through the remaining points using the midpoint circle algorithm
%start here tomorrow   
while ( x2 < y2 - 1)
    
    x2 = x2 + 1;
    if ( d2 < 0 ) 
        d2 = d2 + x2 + x2 + 1;
    else 
        y2 = y2 - 1;
        a2 = x2 - y2 + 1;
        d2 = d2 + a2 + a2;
    end
    
    % setup point indices
    py2 = [x2+cx2, y2+cx2, y2+cx2, x2+cx2, -x2+cx2, -y2+cx2, -y2+cx2, -x2+cx2];
    px2 = [y2+cy2, x2+cy2, -x2+cy2, -y2+cy2, -y2+cy2, -x2+cy2, x2+cy2, y2+cy2];
    
    % loop through each point
    for point_index = 1:length(py2)
        
        % check whether the point is within the arc made by arc_angle
        
        angle_val = atan2(py2(point_index) - cx2, px2(point_index) - cy2) + pi;
        if (arc_angle1b<= angle_val && angle_val <= arc_angle2b)
            circle(py2(point_index), px2(point_index)) = MAGNITUDE;
        end
    end
end


% create the figure
if plot_circle
    figure;
    imagesc(circle, [-1 1]);
    colormap(getColorMap);
    axis image;
    xlabel('y-position [grid points]');
    ylabel('x-position [grid points]');
end