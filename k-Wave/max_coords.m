function [i_max,j_max,max_val] = max_coords(sensor_data,Nx,Ny)
i_max = 0;
j_max = 0;
max_val = 0;

for i = 1:Nx
     for j = 1:Ny
         if max_val < sensor_data(i,j)
             i_max = i;
             j_max = j;
             max_val = sensor_data(i,j);
         end
     end
end

end
