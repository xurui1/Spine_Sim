function [RMS_p,Max_p,RMS_pnorm,Max_pnorm] = spinaldata(Nx,Ny,pmax_data,prms_data,medium_img)

Num_spine_pts = 0;
RMS_p = 0;
Max_p = 0;
RMS_pnorm = 0;
Max_pnorm = 0;

%absolute_max_p = max(pmax_data(:));
%absolute_max_prms = max(prms_data(:));

absolute_max_p = 0;
absolute_max_prms = 0;
for i = 2:Nx-1
    for j = 2:Ny-1
        sum = 0;
        for i_local = i-1:i+1
            for j_local = j-1:j+1
                sum = sum+medium_img(i_local,j_local);
            end
        end
        if sum > 0 && sum < 9
           if pmax_data(i,j) > absolute_max_p
               absolute_max_p = pmax_data(i,j);
           end
           if prms_data(i,j) > absolute_max_prms
               absolute_max_prms = prms_data(i,j);
           end
        end
    end
end


for i = 1:Nx
    for j = 1:Ny
        if medium_img(i,j) ==0
           Iplus = 0;
           i_local = i;
           while i_local <= Nx &&  Iplus ==0
               if medium_img(i_local,j) ==1;
                   Iplus = 1;
               end
               i_local = i_local + 1;
           end
           
           Iminus = 0;
           i_local = i;
           while i_local > 0 &&Iminus ==0
                if medium_img(i_local,j) ==1;
                    Iminus = 1;
                end
                i_local = i_local - 1;
           end
           
           Jplus = 0;
           j_local = j;
           while j_local <= Ny &&  Jplus ==0
               if medium_img(i,j_local) ==1;
                   Jplus = 1;
               end
               j_local = j_local + 1;
           end
           
           Jminus = 0;
           j_local = j;
           while j_local > 0 &&Jminus ==0
                if medium_img(i,j_local) ==1;
                    Jminus = 1;
                end
                j_local = j_local - 1;
           end
           
           if Iplus+Iminus+Jplus+Jminus ==4
               RMS_p = RMS_p + prms_data(i,j);
               Max_p = Max_p + pmax_data(i,j);
               Num_spine_pts = Num_spine_pts+1;
           end
           
        end
        
    end
end
RMS_p = RMS_p / Num_spine_pts;
Max_p = Max_p / Num_spine_pts;

RMS_pnorm = RMS_p / (absolute_max_prms);
Max_pnorm = Max_p / (absolute_max_p);
