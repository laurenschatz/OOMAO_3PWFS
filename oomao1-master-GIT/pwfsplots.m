%%Figure generator

%3PWFS Full Frame

PWFS3FFstrND00=[40.81, 43.02, 42.44, 43.48, 43.75, 43.99, 43.55, 43.81, 42.8,44.2];
PWFS3FFstrND10=[42.74,42.33, 42.74,42.37,40.66,41.89,42.14];
%pierre retake
%PWFS3FFstrND10=[37.8, 38.83, 38.28, 38.27, 39.37, 38.4, 38.86, 38, 39.17, 38.83];
PWFS3FFstrND15=[32.04, 32.25, 32.12, 32.2, 32.08, 32.07, 32.29, 32.34, 32.04, 32.14];

%4PWFS Full Frame
PWFS4FFstrND00=[43.49, 43.68, 42.61, 43.3, 42.11, 41.9, 43.03, 43.21, 42.25, 43.37];
PWFS4FFstrND10=[41.46, 41.27, 40.24, 40.67, 42.11, 41.38, 41.02, 40.95, 40.57, 41.97];
PWFS4FFstrND15=[33.05,32.76,31.97,32.03,31.75,31.72,31.21,31.56,31.55,31.79];

%3PWFS SM
PWFS3SMstrND00d1=[37.95, 40.25, 40.82, 38.08, 37.99, 42.81, 43.91, 39.25, 38.01, 37.19];
PWFS3SMstrND10d1=[38.77, 38.66, 37.71, 37.92, 39.17, 40.36, 36.1, 35.96, 38.46];
PWFS3SMstrND15d2=[35.41, 35.96, 35.96, 35.55, 35.32, 35.46, 34.62, 34.75, 35.06, 35.65];
PWFS3SMstrND00d2=[37.74, 38.67, 38.92, 41.12, 38.6, 38.41, 38.8, 38.89, 39, 38.82];
PWFS3SMstrND10d2=[37.33, 37.26, 37.97, 36.38, 37.2, 37.5, 36.86, 36.85, 36.86];
PWFS3SMstrND20d2=[27.15, 27.1, 26.89, 26.67, 26.78, 27.06, 31.2, 26.84, 28.45];

%4PWFS SM
PWFS4SMstrND00=[39.83, 39.99, 39.47, 38.37, 40.07, 39.9, 40.05,40.17, 39.89, 40.07];
PWFS4SMstrND10=[39.37, 40.59, 38.09, 40.31, 39.67, 38.29, 38.94, 38.41, 37.98, 38.67];
PWFS4SMstrND15=[33.56, 34.56, 33.69, 34.02, 31.91, 32.33, 33.41, 34.17, 32.51, 32.68];
PWFS4SMstrND20=[27.41, 27.15, 27.89, 27.04, 27, 27.55, 28.25, 28.34, 27.79, 27.26];

%% Scatter plot
% x=ones(size(PWFS4SMstrND00))*0;
% figure; scatter(x,PWFS4SMstrND00); hold on
% scatter(x, PWFS3SMstrND00d1)
% scatter(x, PWFS3SMstrND00d2)
% x=ones(size(PWFS4SMstrND00))*0;
%% ND 0
figure; plot(PWFS4SMstrND00); hold on
plot(PWFS3SMstrND00d1)
plot(PWFS3SMstrND00d2)
plot(PWFS3FFstrND00);
plot(PWFS4FFstrND00)

x=[0];
figure; scatter(x,mean(PWFS4SMstrND00)); hold on
scatter(x,mean(PWFS3SMstrND00d1))
scatter(x,mean(PWFS3SMstrND00d2))
scatter(x,mean(PWFS3FFstrND00));
scatter(x,mean(PWFS4FFstrND00))

%% ND 1
figure; plot(PWFS4SMstrND10); hold on
plot(PWFS3SMstrND10d1)
plot(PWFS3SMstrND10d2)
plot(PWFS3FFstrND10);
plot(PWFS4FFstrND10)

x=[1];
figure; scatter(x,mean(PWFS4SMstrND10)); hold on
scatter(x,mean(PWFS3SMstrND10d1))
scatter(x,mean(PWFS3SMstrND10d2))
scatter(x,mean(PWFS3FFstrND10));
scatter(x,mean(PWFS4FFstrND10))

%% ND 1.5

figure; plot(PWFS4SMstrND15); hold on
plot(PWFS3SMstrND15d2)
plot(PWFS3FFstrND15);
plot(PWFS4FFstrND15)

x=[1.5];
figure; scatter(x,mean(PWFS4SMstrND10)); hold on
scatter(x,mean(PWFS3SMstrND15d2))
scatter(x,mean(PWFS3FFstrND15));
scatter(x,mean(PWFS4FFstrND15))

%% ND 2

figure; plot(PWFS4SMstrND20); hold on
plot(PWFS3SMstrND20d2)


x=[2];
figure; scatter(x,mean(PWFS4SMstrND20)); hold on
scatter(x,mean(PWFS3SMstrND20d2))

%% Plot of them all

PWFS3SM=[mean(PWFS3SMstrND00d2) ,mean(PWFS3SMstrND10d2),mean(PWFS3SMstrND15d2),mean(PWFS3SMstrND20d2)];
PWFS3SMflux=[23929,2236,641,250];
PWFS3SMerror=[std(PWFS3SMstrND00d2) ,std(PWFS3SMstrND10d2),std(PWFS3SMstrND15d2),std(PWFS3SMstrND20d2)];

PWFS4SM=[mean(PWFS4SMstrND00) ,mean(PWFS4SMstrND10),mean(PWFS4SMstrND10) ,mean(PWFS4SMstrND20)];
PWFS4SMflux=[26643,2310,593,240];
PWFS4SMerror=[std(PWFS4SMstrND00) ,std(PWFS4SMstrND10),std(PWFS4SMstrND10) ,std(PWFS4SMstrND20)];

PWFS3FF=[mean(PWFS3FFstrND00) ,mean(PWFS3FFstrND10),mean(PWFS3FFstrND15)];
PWFS3FFflux=[23970,2211,661];
PWFS3FFerror=[std(PWFS3FFstrND00) ,std(PWFS3FFstrND10),std(PWFS3FFstrND15)];

PWFS4FF=[mean(PWFS4FFstrND00), mean(PWFS4FFstrND10), mean(PWFS4FFstrND15)];
PWFS4FFflux=[25922, 2239, 604];
PWFS4FFerror=[std(PWFS4FFstrND00), std(PWFS4FFstrND10), std(PWFS4FFstrND15)];

figure; errorbar(log10(PWFS3FFflux), PWFS3FF,PWFS3FFerror); hold on
errorbar(log10(PWFS3SMflux),PWFS3SM,PWFS3SMerror)
errorbar(log10(PWFS4SMflux),PWFS4SM,PWFS4SMerror)
errorbar(log10(PWFS4FFflux), PWFS4FF,PWFS4FFerror)
title('Strehl Vs Flux in Pupil', 'FontSize',24)
ylabel('Strehl Value','FontSize',24)
xlabel('Log of Flux in Pupils','FontSize',24)
lgd =legend('3PWFSFF','3PWFSM','PWFS4SM','4PWFSFF');
set (gca, 'xdir', 'reverse' )
lgd.FontSize = 20;

ax=gca;
ax.XAxis.FontSize=18;
ax.YAxis.FontSize=18;





