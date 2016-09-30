function [FWHM_matrix, HM_x_max,HM_y_max,FWHM  ] = FWHM_calc(sensor_data,Nx,Ny,max,dx)

    half_max = max/2;
    half_max_coords = struct('x',[],'y',[]);
    num_half_max_pts = 1;
    
    %Loop over all non-edge points
    for i = 2:Nx-1
        for j = 2:Ny-2
            %compute max,min at local 3-by-3 grid
            local_max = sensor_data(i,j);
            local_min = sensor_data(i,j);
            for local_i = (i-1):(i+1)
                for local_j = (j-1):(j+1)
                    if sensor_data(local_i,local_j) < local_min
                        local_min = sensor_data(local_i,local_j);
                    end
                    if sensor_data(local_i,local_j) > local_max
                        local_max = sensor_data(local_i,local_j);
                    end
                end
            end
            %record location of HM point
            if local_max >= half_max  && local_min <=half_max
                half_max_coords.x(num_half_max_pts) = i;
                half_max_coords.y(num_half_max_pts) = j;
                num_half_max_pts = num_half_max_pts + 1;
            end
        end
    end
    %find centre of HM circle for fluid simulation
    HM_circle_x = 0;
    HM_circle_y = 0;
     for i = 1:num_half_max_pts-1
          HM_circle_x = HM_circle_x+half_max_coords.x(i);
          HM_circle_y = HM_circle_y+half_max_coords.y(i);
     end
     HM_circle_x = HM_circle_x/(num_half_max_pts-1);
     HM_circle_y = HM_circle_y/(num_half_max_pts-1);
     HM_x_max= HM_circle_x;
     HM_y_max= HM_circle_y;


     %find mean distance from HM point to centre of HM circle
     mean_norm=0; 
     for i = 1:num_half_max_pts-1
         mean_norm = mean_norm + sqrt((half_max_coords.x(i)-HM_circle_x)^2 + (half_max_coords.y(i)-HM_circle_y)^2);
     end
     mean_norm = mean_norm/(num_half_max_pts-1);
     FWHM = 2.0*dx*mean_norm;
     
     %plot HM points to see where the points actually are in fluid
     %simulation
     FWHM_matrix = zeros(Nx,Ny);
     for i = 1:num_half_max_pts-1
         FWHM_matrix(half_max_coords.x(i),half_max_coords.y(i)) = 1;
     end


end
