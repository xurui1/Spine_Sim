function p = two_transducers_fluid(Nx,Ny,mask,source_freq1,source_freq2,t_array,kgrid,medium)

        p_number = 1;
        for j = 1:Ny/2
            for i = 1:Nx
                if mask(i,j) == 1
                     p(p_number,:) = source_mag1*sin(2*pi*source_freq1*t_array);
                     p_number = p_number + 1;
                end
            end
        end
        for j = (Ny/2+1):Ny
            for i = 1:Nx
                if mask(i,j) == 1
                    p(p_number,:) = source_mag2*sin(2*pi*source_freq2*t_array);
                    p_number = p_number + 1;
                end
            end
        end
    
        % filter the source to remove high frequencies not supported by the grid
        for i = 1:p_number-1
            p(i,:) = filterTimeSeries(kgrid, medium, p(i,:));
        end
end
    
function [sxx,syy] = two_transducers_elastic(Nx,Ny,mask,source_freq1,source_freq2,t_array,kgrid,medium)
        p_number = 1;
        for j = 1:Ny/2
            for i = 1:Nx
                if mask(i,j) == 1
                      sxx(p_number,:) = source_mag1*sin(2*pi*source_freq1*t_array);
                      %sxy(p_number,:) = source_mag1*sin(2*pi*source_freq1*t_array);
                      syy(p_number,:) = source_mag1*sin(2*pi*source_freq1*t_array);
                      p_number = p_number + 1;
                end
            end
        end
        for j = (Ny/2+1):Ny
            for i = 1:Nx
                if mask(i,j) == 1
                    sxx(p_number,:) = source_mag2*sin(2*pi*source_freq2*t_array);
                    %sxy(p_number,:) = source_mag2*sin(2*pi*source_freq2*t_array);
                    syy(p_number,:) = source_mag2*sin(2*pi*source_freq2*t_array);
                    p_number = p_number + 1;
                end
            end
        end
    
    % filter the source to remove high frequencies not supported by the grid
    for i = 1:p_number-1
        sxx(i,:) = filterTimeSeries(kgrid, medium, sxx(i,:));
        %source.sxy(i,:) = filterTimeSeries(kgrid, medium, source.sxy(i,:));
        syy(i,:) = filterTimeSeries(kgrid, medium, syy(i,:));
    end
        
    end
    
    
  
    