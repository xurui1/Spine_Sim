function freq_theta_data(theta_data)

% generate vectors and arrays for plotting
y_vec = 2*theta_data.angle(1:end);
x_vec = (10^(-6))*theta_data.frequency(1,1:end);
fluid_RMS_p = theta_data.spinalcord_fluid_RMS_p;
fluid_Max_p = theta_data.spinalcord_fluid_Max_p;
elastic_RMS_p = theta_data.spinalcord_elastic_RMS_p;
elastic_Max_p = theta_data.spinalcord_elastic_Max_p;
fluid_RMS_pnorm = theta_data.spinalcord_fluid_RMS_pnorm;
fluid_Max_pnorm = theta_data.spinalcord_fluid_Max_pnorm;
elastic_RMS_pnorm = theta_data.spinalcord_elastic_RMS_pnorm;
elastic_Max_pnorm = theta_data.spinalcord_elastic_Max_pnorm;

% make figure for RMS and max pressure from fluid and elastic simulations
figure;
subplot(2,2,1);
imagesc(x_vec,y_vec,fluid_RMS_p);
title('Fluid RMS pressure');
ylabel('Angle subtended [rad]');
xlabel('Frequency [MHz]');
%axis image;
c = colorbar;
minimum = min(fluid_RMS_p(:));
maximum = max(fluid_RMS_p(:));
caxis([minimum, maximum]);
ylabel(c,'Pressure [Pa]');

subplot(2,2,2);
imagesc(x_vec,y_vec,fluid_Max_p);
title('Fluid Max pressure');
ylabel('Angle subtended [rad]');
xlabel('Frequency [MHz]');
%axis image;
c = colorbar;
minimum = min(fluid_Max_p(:));
maximum = max(fluid_Max_p(:));
caxis([minimum, maximum]);
ylabel(c,'Pressure [Pa]');

subplot(2,2,3);
imagesc(x_vec,y_vec,elastic_RMS_p);
title('Elastic RMS pressure');
ylabel('Angle subtended [rad]');
xlabel('Frequency [MHz]');
%axis image;
c = colorbar;
minimum = min(elastic_RMS_p(:));
maximum = max(elastic_RMS_p(:));
caxis([minimum, maximum]);
ylabel(c,'Pressure [Pa]');

subplot(2,2,4);
imagesc(x_vec,y_vec,elastic_Max_p);
title('Elastic Max pressure');
ylabel('Angle subtended [rad]');
xlabel('Frequency [MHz]');
%axis image;
c = colorbar;
minimum = min(elastic_Max_p(:));
maximum = max(elastic_Max_p(:));
caxis([minimum, maximum]);
ylabel(c,'Pressure [Pa]');

% make figure for Max and RMS normalized pressure from fluid and elastic
% simulations
figure;
subplot(2,2,1);
imagesc(x_vec,y_vec,fluid_RMS_pnorm);
title('Fluid RMS normalized pressure');
ylabel('Angle subtended [rad]');
xlabel('Frequency [MHz]');
%axis image;
c = colorbar;
minimum = min(fluid_RMS_pnorm(:));
maximum = max(fluid_RMS_pnorm(:));
caxis([minimum, maximum]);
ylabel(c,'Normalized Pressure');

subplot(2,2,2);
imagesc(x_vec,y_vec,fluid_Max_pnorm);
title('Fluid Max normalized pressure');
ylabel('Angle subtended [rad]');
xlabel('Frequency [MHz]');
%axis image;
c = colorbar;
minimum = min(fluid_Max_pnorm(:));
maximum = max(fluid_Max_pnorm(:));
caxis([minimum, maximum]);
ylabel(c,'Normalized Pressure');

subplot(2,2,3);
imagesc(x_vec,y_vec,elastic_RMS_pnorm);
title('Elastic RMS normalized pressure');
ylabel('Angle subtended [rad]');
xlabel('Frequency [MHz]');
%axis image;
c = colorbar;
minimum = min(elastic_RMS_pnorm(:));
maximum = max(elastic_RMS_pnorm(:));
caxis([minimum, maximum]);
ylabel(c,'Normalized Pressure');

subplot(2,2,4);
imagesc(x_vec,y_vec,elastic_Max_pnorm);
title('Elastic Max normalized pressure');
ylabel('Angle subtended [rad]');
xlabel('Frequency [MHz]');
%axis image;
c = colorbar;
minimum = min(elastic_Max_pnorm(:));
maximum = max(elastic_Max_pnorm(:));
caxis([minimum, maximum]);
ylabel(c,'Normalized Pressure');
