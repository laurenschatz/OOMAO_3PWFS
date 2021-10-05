gain=0.1:0.1:1.2;

PWFSSM3=[27.55,41.69, 52.92, 59.3, 63.89, 66.81, 68.56, 69.74, 70, 71.06,71.74, 71.89];
PWFSSM4=[27.54,41.62, 53.03, 59.41, 64.05, 67.04, 68.8, 70, 70.27, 71.31, 72.1, 72.19];

PWFSRI3=[27.57, 41.47, 53, 59.39, 63.97,66.89, 68.64,69.81, 70.05, 71.11,71.8, 71.94];
PWFSRI4=[27.58, 41.72, 53.18, 59.61, 64.26, 67.22, 68.98 70.18, 70.43, 71.47, 72.26, 72.34];

%%

figure; plot(gain,PWFSSM3); hold on
plot(gain,PWFSSM4)
plot(gain, PWFSRI3)
plot(gain, PWFSRI4)
title('Strehl vs Loop Gain Mag 0 star', 'FontSize',24);
xlabel('Gain', 'FontSize',24);
ylabel('Strehl', 'FontSize', 24);
axis tight
lgd =legend('3PWFSM', '4PWFSSM','3PWFSRI', '4PWFSRI');