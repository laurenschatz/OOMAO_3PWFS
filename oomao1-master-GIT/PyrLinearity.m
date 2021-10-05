%% This Script creates a Shack Hartmann wavefront sensor, and builds an interaction matrix based on Fourier Modes



%% Step  0 parameters
% npix=60; %Number of pixels across a pupil
nLenslet=20;
resolution=200;
lambda=700*10^-9;%wavelength
%stroke=lambda/100; %Deformable mirror stroke. (OPD error 1/100 wave)
%amp=(2*pi/lambda)*stroke; %Amplitude of error of the fourier modes.
Amp=(2*pi/lambda)*1*10^-9; %wavefront error in radians.
%% Step One define the telecsope
nPx = 200;
tel = telescope(3.6,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/250);
cam = imager('diameter',tel.D);


%% Step Two define calibration source
ngs = source('wavelength',photometry.V);


%%
ngs*tel*cam
figure; imagesc(cam.frame/max(max(cam.frame))); axis equal
title('Airy Pattern','FontSize',28)
%%
x=-39:40
figure; plot(x,log(cam.frame(40,:)/max(cam.frame(40,:))), 'LineWidth',4)
title('Airy Pattern Log Scale','FontSize',28)
ylabel('Log of Intensity')
xlabel('Distance from peak in pixels')
ax = gca;
ax.LineWidth = 4; % Make tick marks thicker.
ax.TickLength = [.03, .03];
set(gca,'FontSize',20)
xlim([-39 40])
%% Step Three Define the wavefront sensor

%pi*.3315
wfs = pyramid(10,200,'modulation',0, 'alpha',pi/3, 'nFaces',3, 'src',ngs,'tel',tel, 'c', 2, 'rotation',3*pi/2, 'altSlopes',0,'alternative', 0);
%wfs = pyramid(nLenslet,nPx,'modulation',0, 'alpha', pi/2, 'rotation',0, 'nFaces',3, 'src',ngs,'tel',tel, 'c', 2, 'rotation',pi/2, 'alternative',0);

%wfs = pyramid(nLenslet,nPx,'modulation',0, 'alpha',pi/2, 'rotation',0, 'nFaces',4, 'src',ngs,'tel',tel, 'c', 10, 'rotation',0, 'alternative',1);
%% Step Four Propogate
ngs=ngs.*tel*wfs;


%% YOU NEED THIS GURL DONT DELETE
ngs=ngs.*tel*wfs;
wfs.INIT
+wfs
%%
ngs=ngs.*tel*wfs;
refImg=wfs.camera.frame;
%Display the WFS
figure; imagesc(wfs.camera.frame)
% % WFS slopes display
% subplot(1,2,2)
% slopesDisplay(wfs)


%% Step Six Generate Fourier Modes
% % x spatial frequency
%         l;
% % y spatial frequency
zmode = zernike(1:15, tel.D,'resolution', nPx );
%fmode=fourierModes(15, 200);
%fmode=fourierModes(9, 200);
%
sz=size(zmode.modes);
fmode=zmode.modes(:,9);
fmPhase=reshape(fmode,sqrt(sz(1)),sqrt(sz(1)));
figure; imagesc(fmPhase)

%%
sz=size(fmode.modes);
for i=1:81
    fmPhase(:,:,i)=reshape(fmode.modes(:,i),sqrt(sz(1)),sqrt(sz(1)));
    %figure; imagesc(fmPhase(:,:,i)); 
   
end
%%
theFourierMode={tel.pupil, -fmPhase};
%theFourierMode={tel.pupil, fmPhase(:,:,14).*1};
ngs=ngs.*tel*theFourierMode*wfs;
%figure; imagesc(wfs.slopesMap);
figure; imagesc(wfs.camera.frame); axis equal 
%figure; imagesc(wfs.slopesMap-wfs.referenceSlopesMap); 
%%
a=wfs.slopesMap
figure; plot(a(20,1:40)); axis equal

b=wfs.camera.frame-refImg;
figure; plot(b(20,:)); axis equal
title('Slice of Pyramid Detector Signal Reference Subtracted', 'FontSize',28)
c=wfs.camera.frame
figure; plot(c(20,:)); axis equal
title('Slice of Pyramid Detector Signal', 'FontSize',28)


figure; imagesc(wfs.camera.frame); axis equal
title('Pyramid Detector Signal', 'FontSize',28)


%%

a=wfs.slopesMap;
figure; imagesc(a(:,1:40)); axis equal
title('Pupil Plane Signal', 'FontSize',28)
%%
ngs.*tel*theFourierMode*cam;
figure; imagesc(cam.frame); axis equal
title('PSF with Cosine Phase Error', 'FontSize',28)

%% Propogate
counter=1;
for i=-1000:50:1000
    
theFourierMode={tel.pupil, fmPhase.*Amp.*i};
ngs=ngs.*tel*theFourierMode*wfs;

savecam(:,:,counter)=wfs.camera.frame;
saveslope(:,:,counter)=wfs.slopesMap;
counter=counter+1;
end

%%
sCord=saveslope(100,100,:); %dead center
sCord2=saveslope(97,92,:); % side
sCord3=saveslope(97,90,:); %diffraction

sCord=squeeze(sCord);
sCord2=squeeze(sCord2);
sCord3=squeeze(sCord3);

%% 
figure; plot(sCord);
figure; plot(sCord2);
figure; plot(sCord3);


