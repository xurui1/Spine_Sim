figure
scatter(1e3*dx*(freq_data.HM_y_max-(Ny/2)),1e3*dx*(freq_data.HM_x_max-(Nx/2)) );
ylabel('x-position [mm]');
xlabel('y-position [mm]');
set(gca,'YDir','reverse');
title('FWHM centre as a function of frequency');
for i = 1:Num_steps+1
    txt = sprintf(' %i' ,freq_data.frequency(i)/1e3');
    text(1e3*dx*(freq_data.HM_y_max(i)-(Ny/2)),1e3*dx*(freq_data.HM_x_max(i)-(Nx/2)+0.05),txt);
end