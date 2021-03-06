%gainOpPWFSExperiment results 

%% 3PWFSSM

gain=0.1:0.1:1.2;

mag0=[27.55, 41.69, 52.92, 59.3, 63.89, 66.81, 68.56, 69.74, 70, 71.06, 71.74, 71.89];
mag2=[27.55, 41.69, 52.92, 59.3, 63.89, 66.81, 68.56, 69.74, 70, 71.05, 71.73, 71.89];
mag4=[27.55, 41.69, 52.92, 59.3, 63.89, 66.8, 68.55, 69.73, 69.99, 71.04, 71.72, 71.86];
mag6=[27.55, 41.69, 52.91, 59.27, 63.86, 66.77, 68.5, 69.67, 69.91, 70.94, 71.6, 71.73];
mag8=[27.53, 41.64, 52.84, 59.16, 63.86, 66.52, 68.19, 69.29, 69.44, 70.36, 70.9, 70.87];
mag10=[27.53, 41.64, 52.84, 59.16, 63.86, 66.52, 68.19, 69.29, 69.44, 70.36, 70.9, 70.87]; %CHECK THIS
mag12=[26.18, 38.26, 47.75, 52.33, 54.77, 55.34, 54.66, 52.98, 50.69, 48.49, 45.53, 42.24];
mag14=[20.44, 24.43, 28.15, 29.13, 29.33, 28.27, 26.37, 24.51, 22.35, 20.74, 19.04, 17.54];
mag16=[14.96, 13.85, 12.99, 11.8, 11.13, 10.21, 9.55, 8.86, 8.26, 7.63, 7.17, 6.67];

figure; plot(gain, mag0); hold on
plot(gain,mag2); plot(gain,mag4); plot(gain,mag6); plot(gain,mag8); plot(gain,mag10); plot(gain,mag12); plot(gain,mag14); plot(gain,mag16);
ylabel('Strehl', 'FontSize', 24)
xlabel('Gain', 'FontSize', 24)
title('Strehl vs Gain 3PWFS SM', 'FontSize', 24);
lgd=legend('mag0','mag2','mag4','mag6', 'mag8', 'mag10', 'mag12', 'mag14', 'mag16');
lgd.FontSize = 20;