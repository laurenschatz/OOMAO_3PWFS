% Gain Optimization Results

%% 3PWFS Slopes 
SM3mag=[0,2,4,6,8,10,12,16];

%flux is without read noise
SM3fluxNoNoise=[1.0992e+08,1.7421e+07, 2.7610e+06, 4.3759e+05, 6.9354e+04,1.0992e+04,1.7421e+03, 276.1020,6.9354];

SM3fluxNoise=[109920120,17415378,2762307,436958,69579,11011,1738,266,43];


%gain
SM3gain=[0.7, 0.5, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5];



%Mean Strehl
SM3strehl=[68.36, 66.42, 68.35, 63.5, 63.32, 62.13, 54.19, 29.22, 10.25];

%% 4PWFS Slopes

SM4fluxNoise=[109139317,17298805,2740476,433847,68622,10913,1794, 298, 43];

SM4gain=[0.7,0.9,0.6,0.5,0.5,0.5,0.5,0.5];

SM4strehl=[68.59, 70.96, 66.61, 63.65, 63.53, 62.49, 54.89, 28.62, 11.54];

%% 4PWFS RI

RI4gain=[1.2, 0.8, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];

RI4strehl=[72.83, 70.14, 63.87, 63.84, 63.75, 62.78, 55.11, 30.33, 12.29];

%% #3PWFS RI

RI3gain=[0.8, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];

RI3strehl=[69.79, 71.46, 63.62, 63.61, 63.44, 62.35, 54.74, 29.71, 11.5];

%% Plotting

figure; plot(log10(SM3fluxNoise), SM3strehl, 'color', 'red'); hold on
plot(log10(SM4fluxNoise), SM4strehl, 'color', 'blue');
plot(log10(SM4fluxNoise), RI4strehl, 'color', 'black')
plot(log10(SM3fluxNoise), RI3strehl, 'color', 'green')
scatter(log10(SM3fluxNoise), SM3strehl,[], [1 0 0]); 
scatter(log10(SM4fluxNoise), SM4strehl,[], [0 0 1]);
scatter(log10(SM4fluxNoise), RI4strehl,[], [0 0 0]);
scatter(log10(SM3fluxNoise), RI3strehl, [], [0 1 0]);

title('Strehl Vs Flux in Pupil', 'FontSize',24)
ylabel('Strehl Value','FontSize',24)
xlabel('Log of Flux in Pupils','FontSize',24)
lgd =legend('3PWFSM', '4PWFSSM','4PWFSRI', '3PWFSRI');
set (gca, 'xdir', 'reverse' )
lgd.FontSize = 20;

ax=gca;
ax.XAxis.FontSize=18;
ax.YAxis.FontSize=18;

%% 

figure; plot(SM3mag,SM3gain); hold on
plot(SM3mag,SM4gain)
plot(SM3mag, RI4gain)
plot(SM3mag, RI3gain)
title('Loop Gain Vs Guide Star Magnitude', 'FontSize',24);
xlabel('Guide Star Magnitude', 'FontSize',24);
ylabel('Loop Gain', 'FontSize', 24);
axis tight
lgd =legend('3PWFSM', '4PWFSSM','4PWFSRI', '3PWFSRI');

