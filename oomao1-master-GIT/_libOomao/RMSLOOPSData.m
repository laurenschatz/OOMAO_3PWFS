%% Data from testbench

%% Slopes Map
PWFS4_SM_flux=[6477809.76, 582800.8752,224616.7368,177301.0944,65802.5424,47281.374,24272.8488];%30311.3016];
PWFS4_SM_RMS=[107,110,120,125,188,220,254];%346];

PWFS3_SM_flux=[6300756.21,551862.3717,223697.3517,171013.7448,58543.5033,45168.95259];
PWFS3_SM_RMS=[114,115.3,128.8,139.9,182,216];

figure; plot(log(PWFS4_SM_flux),PWFS4_SM_RMS); hold on
plot(log(PWFS3_SM_flux),PWFS3_SM_RMS)
xlabel('Log of Flux in PWFS Pupils');
ylabel('RMS nanometers')
legend('4PWFS Slopes Map','3PWFS Slopes Map');
%% Full Frame
PWFS4_FF_flux=[5835187.08, 526669.236, 201462.12, 154521.6264, 57291.354];
PWFS4_FF_RMS=[100,107.2,125.1,138.6,247];

PWFS3_FF_flux=[5843393.46, 518719.869, 201272.2911, 153208.1547,56920.2633];
PWFS3_FF_RMS=[103.9, 110.8, 122.5, 138.1, 243.5];

figure; plot(log(PWFS4_FF_flux),PWFS4_FF_RMS); hold on
plot(log(PWFS3_FF_flux),PWFS3_FF_RMS)
xlabel('Log of Flux in PWFS Pupils');
ylabel('RMS nanometers')
legend('4PWFS Full Frame','3PWFS Full Frame');

%% All
figure; plot(log(PWFS4_SM_flux),PWFS4_SM_RMS); hold on
plot(log(PWFS3_SM_flux),PWFS3_SM_RMS)
plot(log(PWFS4_FF_flux),PWFS4_FF_RMS);
plot(log(PWFS3_FF_flux),PWFS3_FF_RMS);
xlabel('Log of Flux in PWFS Pupils');
ylabel('RMS nanometers')
legend('4PWFS Slopes Map','3PWFS Slopes Map', '4PWFS Full Frame','3PWFS Full Frame');




