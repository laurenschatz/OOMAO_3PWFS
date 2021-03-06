% 
% 
% 
% 
% figure; plot(fmPhase(120,:)/max(fmPhase(120,:)))
% title('Phase Error','FontSize',28)
% xlabel('Position in Pixels', 'FontSize',28)
% ylabel('Amplitude', 'FontSize',28)
% 
% 
% % set(gca,'FontSize', 28)

%%
fmPhase=fmPhase(:,:,14);
figure; imagesc(fmPhase/max(max(fmPhase)));axis equal
title('Phase Error','FontSize',28)

%%

figure; imagesc(wfs.camera.frame); axis equal
title('Pyramid Detector Signal', 'FontSize',28)

%% pyramid dector signal slices
%4PWFS
figure; plot(wfs.camera.frame(20,:));
title('Slice of Pyramid Detector Signal', 'FontSize',26)
xlabel('Position in Pixels', 'FontSize',28)
ylabel('Intensity', 'FontSize',28)

% 3PWFS
figure; plot(wfs.camera.frame(40,:));
title('Slice of Pyramid Detector Signal', 'FontSize',26)
xlabel('Position in Pixels', 'FontSize',28)
ylabel('Intensity', 'FontSize',28)


%% 4pwfs detector signal -ref Img

subSlice=wfs.camera.frame-refImg;
figure; plot(subSlice(20,1:41)); 
title('Slice of Pyramid Detector Signal-I_0', 'FontSize',28)
%%
%3pwfs detector signal -ref Img
subSlice=wfs.camera.frame-refImg;
figure; plot(subSlice(20,40:61)); 
title('Slice of Pyramid Detector Signal-I_0', 'FontSize',28)
%% 4PWFS slice of wfs.camera frame -I_0


c=1:22;
b=0:0.2992:2*pi+0.2992;
x=0:0.0251:2*pi;
a=0:1/250:1;
y=1/4+pi^2-pi*sin(3*b);


z=wfs.camera.frame(21,50:71)/max(max(wfs.camera.frame(21,50:71)));
figure; plot(c,z, 'LineWidth',3); hold on
plot(c,y/max(max(y)), 'LineWidth',3);
xlim([0 23])
title('Slice of Pyramid Detector Signal', 'FontSize',26)
xlabel('Position in Pixels', 'FontSize',28)
ylabel('Scaled Intensity', 'FontSize',28)
legend('Cross section of 4PWFS Pupil', 'Mathematical Prediction','Location','southeast')
set(gca,'FontSize', 32, 'LineWidth',3)

z=subSlice(21,50:71)/max(max(subSlice(21,50:71)));
figure; plot(c,z); hold on
% y=y-mean(y);
% y=y/max(max(y));
y=-sin(3*b);
plot(c,y);
xlim([0 23])
title('Slice of Pyramid Detector Signal - I_{0}', 'FontSize',26)
xlabel('Position in Pixels', 'FontSize',28)
ylabel('Intensity', 'FontSize',28)
legend('Cross section of 4PWFS Pupil', 'Mathematical Prediction','Location','southeast')
set(gca,'FontSize', 32)


%% 3pwfs
c=1:22;
subSlice=wfs.camera.frame-refImg;
z=subSlice(20,40:61)/max(max(subSlice(20,40:61)));
figure; plot(c,z); hold on
% y=y-mean(y);
% y=y/max(max(y));
y=-sin(3*b);
plot(c,y);
xlim([0 23])
title('Slice of Pyramid Detector Signal - I_{0}', 'FontSize',26)
xlabel('Position in Pixels', 'FontSize',28)
ylabel('Intensity', 'FontSize',28)
legend('Cross section of 3PWFS Pupil', 'Mathematical Prediction','Location','southeast')
set(gca,'FontSize', 32)

ylim([-1.5 1])


%%
c=1:22;
b=0:0.2992:2*pi+0.2992;
x=0:0.0251:2*pi;
a=0:1/250:1;
y=1/4+pi^2-pi*sin(3*b);

z=wfs.camera.frame(24,39:60)/max(max(wfs.camera.frame(24,39:60)));
figure; plot(c,z, 'LineWidth',3); hold on
plot(c,y/max(max(y)),'LineWidth',3);
xlim([0 23])
title('Slice of Pyramid Detector Signal', 'FontSize',26)
xlabel('Position in Pixels', 'FontSize',28)
ylabel('Scaled Intensity', 'FontSize',28)
legend('Cross section of 3PWFS Pupil', 'Mathematical Prediction','Location','southeast')
set(gca,'FontSize', 32, 'LineWidth',3)

%%
% 
% slopesMap4=wfs.slopesMap.*wfs.validSlopes;
% 
%  figure; imagesc(slopesMap4); axis equal
%  title('Slopes Map', 'FontSize',26)
 

c=1:22;
slopesMap3=(wfs.slopesMap-wfs.referenceSlopesMap).*wfs.validSlopes;
%slopesMap3=(wfs.slopesMap).*wfs.validSlopes;
figure; imagesc(slopesMap3); axis equal
title('Slopes Map', 'FontSize',26)
figure; plot(c,slopesMap3(20,50:71)/max(slopesMap3(20,50:71))); hold on
figure; plot(slopesMap3(20,:)); hold on
 
%% 

 
%slopesMap4=(wfs.slopesMap-wfs.referenceSlopesMap).*wfs.validSlopes;

slopesMap4=(wfs.slopesMap).*wfs.validSlopes;
figure; imagesc(slopesMap4); axis equal
title('Slopes Map', 'FontSize',26)

figure; plot(slopesMap4(20,1:40))
title('X-Slopes Map Slice', 'FontSize',26)
xlabel('Position in Pixels', 'FontSize',28)
ylabel('Slope', 'FontSize',28)
set(gca,'FontSize', 26)


figure; plot(slopesMap4(20,1:40))
%%
c=1:22;

figure; plot(c,slopesMap4(20,10:31)/max(slopesMap4(20,10:31)),'r', 'LineWidth',3); hold on
title(' 4PWFS X-Slopes Map Cross Section', 'FontSize',16)
xlabel('Position in Pixels', 'FontSize',28)
ylabel('Amplitude', 'FontSize',28)
set(gca,'FontSize', 32, 'LineWidth',3)
ylim([-1.5 1])

b=0:0.2992:2*pi+0.2992;

plot(c,-sin(3*b),'k*-','LineWidth',3)
legend('4PWFS Sx Measurement', 'Mathematical Prediction','Location','southeast')

%% 3pwfs

c=1:22;
figure; plot(c,slopesMap3(20,50:71)/max(slopesMap3(20,50:71)),'r','LineWidth',3); hold on
title('3PWFS X-Slopes Map Cross Section', 'FontSize',16)
xlabel('Position in Pixels', 'FontSize',28)
ylabel('Amplitude', 'FontSize',28) 
set(gca,'FontSize', 32, 'LineWidth',3)
ylim([-1.5 1])

b=0:0.2992:2*pi+0.2992;

plot(c,-sin(3*b),'k*-', 'LineWidth',3)
legend('Scaled 3PWFS Sx Measurement', 'Mathematical Prediction','Location','southeast')

%% errors
b=0:0.2992:2*pi+0.2992;
%direct signal
subSlice=wfs.camera.frame-refImg;
z=subSlice(20,40:61)/max(max(subSlice(20,40:61)));
MSerror=mean((y-z).^2)


% subtracted
slopesMap3=(wfs.slopesMap-wfs.referenceSlopesMap).*wfs.validSlopes;

slopesMap3refSub=slopesMap3(20,50:71)/max(slopesMap3(20,50:71));
b=0:0.2992:2*pi+0.2992;
y=-sin(3*b);
MSerrorSMrefSub=mean((y-slopesMap3refSub).^2)



%not subtracted
slopesMap3=(wfs.slopesMap).*wfs.validSlopes;
slopesMap3SM=slopesMap3(20,50:71)/max(slopesMap3(20,50:71));
MSerrorSM=mean((y-slopesMap3SM).^2)

