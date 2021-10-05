
%% Experimental script to run on exAO0
function pwfsexperiment(path, mag, c, nFaces, altSlopes, gain, mod, alternative)
%path= '/Users/laurenschatz/Documents/MATLAB/Compass Pull/Data/';
%function experiment(path);
%% Optics Setup
% Wavefront sensor sampling
nSamp=20;

%Telescope Diameter
D=2;
%Full Resolution 
nPx=240;

% Sampling Frequency (exposure time)
freq=300;

%% Source Setup
%On-axis natural guide star in R band 
ngs = source('wavelength',photometry.R);

% Science object in R band
science = source('wavelength',photometry.R);

%% Telescope 
tel= telescope(D,'resolution',nPx,'samplingTime',1/freq);

%% Science Camera
cam=imager('diameter', tel.D);
science.magnitude=0;
science.*tel*cam;
cam.referenceFrame=cam.frame; %reference taken for Strehl
camPerfPsf=cam.frame;

cam.photonNoise=true;
cam.readOutNoise=0.5;
%% Atmosphere

r0=16e-2;
L0=6.25;
v=7.4;

atm=atmosphere(photometry.V0,r0,L0,'windSpeed',v,'windDirection',0);

%% Deformable Mirror

nAct=9;
d=D/(nAct-1);
fmax=1/(2*d); %Max spatial frequency
nMax=ceil(D*fmax/0.37-1); % Max mode corrected by the DM

nModes=sum(1:nMax+1)-1;
%zern=zernike(2:nModes+1, D, 'resolution', nPx); %zernike basis set
bif = influenceFunction('monotonic',30/100);
nActuator = 9;
v=utilities.circle(nActuator,nActuator);
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',logical(v));
%deformable mirror with zernike influence functions
%dm=deformableMirror(nModes,'modes', zern.modes, 'resolution',nPx, 'validActuator', logical(1:nModes));

%% Pyramid Wavefront sensor
wfs = pyramid(nSamp,nPx,'modulation',mod, 'alpha',pi/3, 'nFaces',nFaces, 'src',ngs,'tel',tel, 'c', c, 'rotation',3*pi/2, 'altSlopes',altSlopes,'alternative', alternative);
% AltSlopes==2 for full frame
%propogate
ngs=ngs.*tel*wfs;
wfs.INIT;
%% Save Flux in Pupils

tel=tel-atm;
ngs=ngs.*tel*wfs;
im=wfs.camera.frame;
ave=mean(mean(im));
Mask=im>ave;
ngs.magnitude=mag;
wfs.camera.photonNoise=true;
tel=tel-atm;
ngs=ngs.*tel*wfs;
image=wfs.camera.frame;
image=image.*Mask;
figure; imagesc(image)
flux=sum(sum(image))
save(strcat(path,num2str(mag),'mag',num2str(c), 'c','flux'), 'flux');

wfs.camera.photonNoise=false;
%% Interaction Matrix
ngs=ngs.*tel;
+ngs;




calib=calibration(dm,wfs,ngs, ngs.wavelength/100,10);

% Save Eigen Values
eigenVal=calib.eigenValues;
save(strcat(path,'dmEigenValues'), 'eigenVal');
commandMatrix = calib.M;


%% Noise

wfs.camera.photonNoise=true;
wfs.camera.readOutNoise=0.5;
ngs.magnitude=mag;
ngs=ngs.*tel*wfs;

%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
figure(11)
h = imagesc([turbPhase,ngs.meanRmPhase]);
axis equal tight
colorbar
snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
  set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
  drawnow
      
end

%saving


save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes', num2str(gain),'gain', num2str(mod),'modulation','turbSTD1'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation', 'residSTD1'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase1'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase1'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl1'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl1'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam1'),'saveCam');
%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
%figure(11)
%h = imagesc([turbPhase,ngs.meanRmPhase]);
%axis equal tight
%colorbar
%snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
%   set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
%   drawnow
      
end

%saving

save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','turbSTD2'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','residSTD2'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase2'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase2'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl2'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl2'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam2'),'saveCam');

%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
%figure(11)
%h = imagesc([turbPhase,ngs.meanRmPhase]);
%axis equal tight
%colorbar
%snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
%   set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
%   drawnow
      
end

%saving


save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','turbSTD3'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','residSTD3'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase3'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase3'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl3'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl3'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam3'),'saveCam');
%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
%figure(11)
%h = imagesc([turbPhase,ngs.meanRmPhase]);
%axis equal tight
%colorbar
%snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
%   set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
%   drawnow
      
end

%saving

save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','turbSTD4'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','residSTD4'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase4'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase4'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl4'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl4'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam4'),'saveCam');
%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
%figure(11)
%h = imagesc([turbPhase,ngs.meanRmPhase]);
%axis equal tight
%colorbar
%snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
 %  set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
%   drawnow
      
end

%saving

save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','turbSTD5'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','residSTD5'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase5'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase5'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl5'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl5'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam5'),'saveCam');
%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
%figure(11)
%h = imagesc([turbPhase,ngs.meanRmPhase]);
%axis equal tight
%colorbar
%snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
  % set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
  % drawnow
      
end

%saving

save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','turbSTD6'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','residSTD6'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase6'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase6'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl6'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl6'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam6'),'saveCam');
%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
%figure(11)
%h = imagesc([turbPhase,ngs.meanRmPhase]);
%axis equal tight
%colorbar
%snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
%   set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
%   drawnow
      
end

%saving


save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','turbSTD7'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','residSTD7'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase7'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase7'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl7'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl7'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam7'),'saveCam');
%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
%figure(11)
%h = imagesc([turbPhase,ngs.meanRmPhase]);
%axis equal tight
%colorbar
%snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
%   set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
%   drawnow
      
end

%saving


save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','turbSTD8'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','residSTD8'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase8'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase8'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl8'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl8'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam8'),'saveCam');

%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
%figure(11)
%h = imagesc([turbPhase,ngs.meanRmPhase]);
%axis equal tight
%colorbar
%snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
  % set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
  % drawnow
      
end

%saving
save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','turbSTD9'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','residSTD9'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase9'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase9'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl9'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl9'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam9'),'saveCam');
%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel
draw(tel);
%figure; imagesc(tel)

loopGain=gain; %Loop gain
nIter=100; %number of loop iterations
dm.coefs=0; % 0 The DM\
calib.nThresholded=0; %Cut off low eigen values

%Variables to save
total=zeros(1,nIter);
residual=zeros(1,nIter);

zbasis=zernike(2:2*nModes+1,D, 'resolution', nPx); %zernike modes to fit residual wavefront


SumX2_turb=zeros(2*nModes,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

turbSTD=[];
residSTD=[];
turbPhase = ngs.meanRmPhase;
%figure(11)
%h = imagesc([turbPhase,ngs.meanRmPhase]);
%axis equal tight
%colorbar
%snapnow

% %% Perfect PSF for Strehl
% N=78; %diameter of perfect circle
% perfPupil=utilities.circle(nPx,N);
% perfPsf=fftshift(fft2(fftshift(perfPupil)))/nPx;
% perfPsf=abs(perfPsf).^2;



%%
for n=1:nIter
    
ngs=ngs.*+tel; % propogate through telescope. +tel advances the phase screen one step.    
    
 %statistics of the turbulence
 zbasis.\ngs;
    turbZDcomp(:,n)=zbasis.c;
    SumX2_turb = SumX2_turb + turbZDcomp(:,n).^2;
    SumX_turb = SumX_turb + turbZDcomp(:,n);
    STD_turb= sqrt(SumX2_turb/n-(SumX_turb/n).^2); %standard deviation
    turbSTD(:,n)=STD_turb;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        science=science.*tel*dm*cam;
        psf=cam.frame;
        
%       [strehl,radotf,perfradOtf, r]=utilities.simpsf2strehl(psf,perfPsf);
       [strehl,otf,perfOtf]=utilities.simpsf2strehl(psf,camPerfPsf);
        saveSimStrehl(n-10)=strehl;
        saveStrehl(n-10)=cam.strehl;
        saveCam(:,:,n-10)=cam.frame;
    end
    
    saveResidPhase(:,:,n)=ngs.meanRmOpd; %units in meters
    
    % Mean DM residual coefficients
    
    residualDmCoefs=commandMatrix*wfs.slopes;
   
    %Integrating the DM coefficients
    dm.coefs=dm.coefs-loopGain*residualDmCoefs;

    
    % Display of turbulence and residual phase
%   set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
%   drawnow
      
end

%saving

save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','turbSTD10'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','residSTD10'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase10'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase10'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl10'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl10'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam10'),'saveCam');
end