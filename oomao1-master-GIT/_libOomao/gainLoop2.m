%% loop script
function [finalStrehl,gain,M]=gainLoop2(path, mag,c,nFaces,altSlopes, mod, alternative)

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
%% Interaction Matrix
ngs=ngs.*tel;
+ngs;




calib=calibration(dm,wfs,ngs, ngs.wavelength/100,10);

% Save Eigen Values
eigenVal=calib.eigenValues;
commandMatrix = calib.M;

%% Save Flux in Pupils
tel=tel-atm;
ngs=ngs.*tel*wfs;
im=wfs.camera.frame;
ave=mean(mean(im));
Mask=im>ave;
image=im.*Mask;
flux=sum(sum(image));


%% Noise

wfs.camera.photonNoise=true;
wfs.camera.readOutNoise=0.5;
ngs.magnitude=mag;
ngs=ngs.*tel*wfs;

counter=0;
counter2=0;
for i=0.5:0.1:2.5
    
%% Closed Loop Parameters
tel=tel-atm;
    atm=atmosphere(photometry.V0,r0,L0,'windSpeed',v,'windDirection',0);
tel=tel+atm; %make random Kol phase screen
%atm2=atm;
%tel2=tel;
ngs=ngs.*tel;
draw(tel);
%figure; imagesc(tel)

loopGain=i; %Loop gain
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
    turbSTD(:,n)=STD_turb.*10^9;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c; % units of meters
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD.*10^9; % units of nanometers
        
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

% nM=2:131;
% figure; plot(nM, STD.^2);
% title('WF variance vs Mode')
counter=counter+1;
finalStrehl(:,counter)=mean(saveStrehl);
index(counter)=i;


%% Closed Loop Parameters

%atm=atmosphere(photometry.V0,r0,L0,'windSpeed',v,'windDirection',0);
tel=tel+atm; %make random Kol phase screen
%atm2=atm;
%tel2=tel;
ngs=ngs.*tel;
draw(tel);
%figure; imagesc(tel)

loopGain=i; %Loop gain
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
    turbSTD(:,n)=STD_turb.*10^9;

    savePhase(:,:,n)=ngs.meanRmOpd; %units in meters
    turbPhase = ngs.meanRmPhase;
    
    
    %Propogate to DM/WFS
    ngs=ngs.*tel*dm*wfs;
    
  %statistics of the residual wavefront
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c; % units of meters
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD.*10^9; % units of nanometers
        
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

% nM=2:131;
% figure; plot(nM, STD.^2);
% title('WF variance vs Mode')
counter2=counter2+1;
finalStrehl2(:,counter2)=mean(saveStrehl);
index(counter2)=i;
end

[M1,I1]=max(finalStrehl);
M1;
gain1=index(I1);

[M2,I2]=max(finalStrehl2);
M2;
gain2=index(I2);
g=[gain1,gain2];
gain=max(g);

%%

%%

pwfsexperiment(path, mag, c, nFaces, altSlopes, gain, mod, alternative);

end