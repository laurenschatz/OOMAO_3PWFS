%--------------------------------------------------------------------------
% Script to simulate current bench
% 9x9 DM and 62x62 pixels across the pupil.
%
% Charlotte Bond    05/06/17
%--------------------------------------------------------------------------
%

clear all;
%load('sss.mat')
%% Optic parameters
% Equivalent telescope aperture (for atm params)
D = 2;

% Full resolution
nPx = 240;

% Frequency
freq = 300;

% Set source wave-length (closest to 660nm laser, R2 = 650nm)
ngs = source('wavelength',photometry.R);

% Initiate telescope
tel = telescope(D,'resolution',nPx,'samplingTime',1/freq);

%% Atmosphere
r0 = 6e-2; %7 .8e-2
L0 = 6.25;
v = 7.4;

atm = atmosphere(photometry.V0,r0,L0,'windSpeed',v,'windDirection',0);

%% Deformable mirror
% Zernikes normalised to provide 1 m rms optical path difference with a
% coefficent of 1 (although this will be 2 with the double pass on the DM)

nAct = 9;
d = D/(nAct-1);
fmax = 1/(2*d);
nMax = ceil(D*fmax/0.37 - 1);

nModes = sum(1:nMax+1)-1;
zern = zernike(2:nModes+1,D,'resolution',nPx);

% Modes as produced on the bench
% load('zernikeModes')
% load('influenceFunctions')
% clear modeHR
% for n=1:69
%     modeHR(:,:,n) = interp2(modes(:,:,n),2);
% end
% 
% modeHR = modeHR(1:end-1,1:end-1,:);
% modeHR(isnan(modeHR)) = 0;
% 
% modeRes = size(modeHR,1);
% infFs = reshape(modeHR,modeRes*modeRes,69);
% zmodes = infFs*z2c;

dm = deformableMirror(nModes,'modes',zern.modes,'resolution',nPx,...
    'validActuator',logical(1:nModes));

%% Pyramid parameters
nL = 80;
mod = 5;

pyr = pyramid(nL,nPx,'modulation',mod,'alpha',3*pi/2, 'src',ngs, 'tel',tel);

% % NCPA for reference measurement
% load('NCPAresults')
% zern.c(:) = 0;
% zern.c = 1e-6*Ab_micron(2:end);

ngs = ngs.*tel*pyr;

pyr.INIT;

figure
    subplot(2,1,1)
    imagesc(pyr.camera)
    
    subplot(2,1,2)
    slopesDisplay(pyr)

%% Shack-Hartman parameters
nLsh = 30;
minLight = .01;
sh = shackHartmann(nLsh,nPx,minLight);

ngs = ngs.*tel*sh;

sh.INIT

%% Interaction matrix
ngs = ngs.*tel;
calib = calibration(dm,pyr,ngs,ngs.wavelength/100,10);
% load('zernikeIM_mod1 _rms10nm')
% calib = dmCalib;
%% Science camera
sci = source('wavelength',photometry.R);
cam = imager();

tel = tel - atm;
sci = sci.*tel*cam;

cam.referenceFrame = cam.frame;

figure
    imagesc(cam)

flush(cam)
cam.frame = cam.frame*0;
cam.clockRate    = 1;
exposureTime     = 200;
cam.exposureTime = exposureTime;
startDelay       = 30;
cam.startDelay   = startDelay;

%% Noise
% pyr.camera.photonNoise = false;
% pyr.camera.readOutNoise = 0.1;
% ngs.magnitude = 2;
% 
% ngs = ngs.*tel*pyr;
% 
% figure
%     imagesc(wfs.camera)

%% Closed loop
% CL parameters
tel = tel+atm;
gain = .5;
%nIter = cam.exposureTime+cam.startDelay;
nIter = 1000;
dm.coefs = 0;
calib.nThresholded = 0;
% Recorded results
total = zeros(1,nIter);
residual = zeros(1,nIter);

%%
clear a b c e sr
slps = zeros(length(pyr.slopes),1);

zerzer = zernike(2:2*nModes+1,D,'resolution',nPx);
SumX2_turb=zeros(2*nModes,1);
SumX_turb=zeros(2*nModes,1);
SumX2=zeros(2*nModes,1);
SumX=zeros(2*nModes,1);

STD_turb = [];
STD = [];

figure
for n=1:nIter
    
    +tel;
    ngs = ngs.*tel;
    
    zerzer.\ngs;
    zz=zerzer.c;
    SumX2_turb = SumX2_turb + zz.^2;
    SumX_turb = SumX_turb + zz;
    STD_turb = sqrt(SumX2_turb/n-(SumX_turb/n).^2);
    
    % Record initial phase
    phaseTurb = ngs.meanRmPhase;
    total(n) = var(ngs);
    
    % Propagate to dm/wfs
    ngs = ngs.*tel*dm*pyr;
    
    % Record residual phase
    phaseRes = ngs.meanRmPhase;
    residual(n) = var(ngs);
    
    % Propagate through science camera
    sci = sci.*tel*dm*cam;
    saveStrehl(:,nIter)=cam.strehl;
    % Compute dm commands
    dm.coefs = dm.coefs - gain*(calib.M*pyr.slopes);
    
    slps = pyr.slopes;
    
    % Plotting
    figure(1)
    subplot(2,2,[1 2])
    imagesc(cam)
    drawnow
    subplot(2,2,[3 4])
    imagesc([phaseTurb phaseRes])
    axis equal tight
    colorbar()
    drawnow

    cam.flush()
    a(n) = std2(phaseTurb);
    c(n) = std2(phaseRes);
    sr(n) = exposureTime*cam.strehl;
    
    if n>10
    zerzer.\ngs;
    zz=zerzer.c;
    SumX2 = SumX2 + zz.^2;
    SumX = SumX + zz;
    STD = sqrt(SumX2/n-(SumX/n).^2);
    end
    
    figure(2)
    subplot(1,2,1)
    plot(a*1e9*sci.wavelength/(2*pi))
    hold on
    plot(c*1e9*ngs.wavelength/(2*pi));
    set(gca, 'YScale', 'log')
    hold off
    subplot(1,2,2)
    plot(sr)
    
%     figure(3)
%     plot(STD_turb*1e6,'b-')
%     hold on
%     plot(STD*1e6,'b--')
%     plot(sssOL,'r-')
%     plot(sssCL,'r--')
%     hold off
%     set(gca, 'YScale', 'log')
    
end

    subplot(2,2,[1 2])
    imagesc(cam)
    drawnow
    subplot(2,2,[3 4])
    imagesc([phaseTurb phaseRes])
    axis equal tight
    colorbar()

psf = cam.frame;
% sr = cam.strehl;     