%% Experimental script to run on exAO0
%function pwfsexperiment(path, mag, c, nFaces, altSlopes, gain)
%path= '/Users/laurenschatz/Documents/MATLAB/Compass Pull/Data/';
%function experiment(path);


%%
mag=0;
c=2;
nFaces=3;
altSlopes=0;
gain=0.8;
% if altSlopes==0
%     alpha=pi/2;
% end
% if altSlopes==2
%     alpha=pi/3;
% end
alpha=pi/3;
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

r0=6e-2;
L0=6.25;
v=7.4;

atm=atmosphere(photometry.V0,r0,L0,'windSpeed',v,'windDirection',0);
 %define the amplitude
lambda=atm.wavelength;
Amp=(2*pi/lambda)*0.5*lambda;

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

%% Zernike Modes
zern=zernike(2:nModes+1, D, 'resolution', nPx); %zernike basis set
sz=size(zern.modes);
for i=1:sz(2)
    zernMode(:,:,i)=reshape(zern.modes(:,i),nPx,nPx)*Amp;
end
commonPathError=(.9*zernMode(:,:,7).*.9+zernMode(:,:,4)*-.4+zernMode(:,:,8)*0.15-1.2*zernMode(:,:,5)+0*zernMode(:,:,3))*.1;
% 
% 
theCommonPathError={tel.pupil,commonPathError};
ngs=ngs.*tel*theCommonPathError*cam;
figure; imagesc(cam.frame); axis equal; %colormap gray
+ngs;

%% Pyramid Wavefront sensor


wfs = pyramid(nSamp,nPx,'modulation',3, 'alpha',pi/3, 'nFaces',4, 'src',ngs,'tel',tel, 'c', 2, 'rotation',3*pi/2, 'altSlopes',0,'alternative', 0);
% AltSlopes==2 for Raw Intensity signal handeling
%propogate
ngs=ngs.*tel*theCommonPathError*wfs; %UNCOMMENT TO PUT IN COMMON PATH ERROR IN YOU IDIOT
%ngs=ngs.*tel*wfs;
wfs.INIT;
%%
ngs=ngs.*tel*theCommonPathError*wfs;
figure; imagesc(wfs.camera.frame); colormap gray
%x=zeros(40,80);% Slopes Map
x=zeros(80,80); %Raw Intensity
wfs.referenceSlopesMap=x;
%% Interaction Matrix
% ngs=ngs.*tel*theCommonPathError*wfs;
% +ngs;




calib=calibration(dm,wfs,ngs, ngs.wavelength/100,10);

% Save Eigen Values
eigenVal=calib.eigenValues;
%save(strcat(path,'dmEigenValues'), 'eigenVal');
commandMatrix = calib.M;

%% Save Flux in Pupils
tel=tel-atm;
ngs=ngs.*tel*wfs;
im=wfs.camera.frame;
ave=mean(mean(im));
Mask=im>ave;
image=im.*Mask;
flux=sum(sum(image));
%save(strcat(path,num2str(mag),'mag',num2str(c), 'c','flux'), 'flux');

%% Noise

wfs.camera.photonNoise=true;
wfs.camera.readOutNoise=0.5;
ngs.magnitude=mag;
ngs=ngs.*tel*wfs;

%% Closed Loop Parameters
tel=tel+atm; %make random Kol phase screen
ngs=ngs.*tel;
draw(tel);


%phase=load('staticError.mat');
%staticError=phase.phase;
%theStaticError={tel.pupil,staticError};
%ngs=ngs.*tel*theStaticError;
%+ngs;
%science=science.*tel*theStaticError;
%+science;
%science*cam;
%draw(tel);
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
%ngs=ngs.*tel;   
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
    science=science.*tel*dm*cam;
    figure(100); imagesc(cam.frame); axis equal; axis tight
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
        %science=science.*tel*dm*cam;
        psf=cam.frame;
        %figure(100); imagesc(cam.frame); axis equal; axis tight
        
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
strehlSim10=mean(saveSimStrehl);
strehl10=mean(saveStrehl);
m=2;
%%
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD1'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD1'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase1'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase1'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl1'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl1'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam1'),'saveCam');
%% Closed Loop Parameters
%tel=tel+atm; %make random Kol phase screen
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
strehlSim9=mean(saveSimStrehl);
strehl9=mean(saveStrehl);
%saving
% 
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD2'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD2'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase2'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase2'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl2'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl2'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam2'),'saveCam');

%% Closed Loop Parameters
%tel=tel+atm; %make random Kol phase screen
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

strehlSim8=mean(saveSimStrehl);
strehl8=mean(saveStrehl);
%saving
% 
% 
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD3'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD3'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase3'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase3'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl3'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl3'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam3'),'saveCam');
%% Closed Loop Parameters
%tel=tel+atm; %make random Kol phase screen
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
strehlSim7=mean(saveSimStrehl);
strehl7=mean(saveStrehl);
%saving
% 
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD4'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD4'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase4'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase4'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl4'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl4'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam4'),'saveCam');
%% Closed Loop Parameters
%tel=tel+atm; %make random Kol phase screen
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
strehlSim1=mean(saveSimStrehl);
strehl1=mean(saveStrehl);
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD5'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD5'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase5'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase5'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl5'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl5'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam5'),'saveCam');
%% Closed Loop Parameters
%tel=tel+atm; %make random Kol phase screen
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

strehlSim2=mean(saveSimStrehl);
strehl2=mean(saveStrehl);
%saving

% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD6'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD6'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase6'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase6'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl6'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl6'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam6'),'saveCam');
%% Closed Loop Parameters
%tel=tel+atm; %make random Kol phase screen
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
strehlSim3=mean(saveSimStrehl);
strehl3=mean(saveStrehl);

% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD7'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD7'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase7'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase7'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl7'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl7'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam7'),'saveCam');
%% Closed Loop Parameters
%tel=tel+atm; %make random Kol phase screen
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
strehlSim4=mean(saveSimStrehl);
strehl4=mean(saveStrehl);

% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD8'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD8'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase8'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase8'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl8'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl8'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam8'),'saveCam');

%% Closed Loop Parameters
%tel=tel+atm; %make random Kol phase screen
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
strehlSim5=mean(saveSimStrehl);
strehl5=mean(saveStrehl);
%saving
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD9'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD9'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase9'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase9'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl9'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl9'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam9'),'saveCam');
%% Closed Loop Parameters
%tel=tel+atm; %make random Kol phase screen
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
strehlSim6=mean(saveSimStrehl);
strehl6=mean(saveStrehl);
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','turbSTD10'), 'turbSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','residSTD10'), 'residSTD');
% save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','savePhase10'),'savePhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveResidPhase10'),'saveResidPhase');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveStrehl10'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveSimStrehl10'),'saveStrehl');
% save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes','saveCam10'),'saveCam');
