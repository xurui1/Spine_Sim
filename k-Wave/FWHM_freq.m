figure
plot( freq_data.frequency, 1e3*freq_data.FWHM,'b--o' );
ylabel('FWHM (mm)');
xlabel('Frequency of source f_1');
axis([ 19e4 31e4 0 35]);
title('FWHM as a function of frequency');
