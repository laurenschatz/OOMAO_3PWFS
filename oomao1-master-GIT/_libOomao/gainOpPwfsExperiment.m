%% loop script
function gainOpPwfsExperiment(path, mag,c,nFaces,altSlopes, mod, alternative)
%'/home/lschatz/MATLAB/lschatz/oomao-remote/Data/gainOptimizationPwfsExperiment/3PWFSSM/'
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

cam.photonNoise=false;
cam.readOutNoise=0.0;
%% Atmosphere

r0=16e-2;
L0=6.25;
v=7.4;

atm=atmosphere(photometry.V0,r0,L0,'windSpeed',v,'windDirection',0);
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
ngs=ngs.*tel*wfs;
im=wfs.camera.frame;
image=im.*Mask;
flux=sum(sum(image));


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


%% Interaction Matrix
ngs=ngs.*tel;
+ngs;




calib=calibration(dm,wfs,ngs, ngs.wavelength/100,10);

% Save Eigen Values
eigenVal=calib.eigenValues;
commandMatrix = calib.M;




%% Noise

wfs.camera.photonNoise=true;
wfs.camera.readOutNoise=0.5;
ngs.magnitude=mag;
ngs=ngs.*tel*wfs;

counter=0;
counter2=0;
%for gain=0.1:0.1:1.5
for gain =0.5:0.1:0.7
    
%% Closed Loop Parameters
for realization=1:15;
%% Closed Loop Parameters
tel=tel-atm;
ngs=ngs.*tel;
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
        science=science.*tel*dm*cam;
        psf=cam.frame;
        figure(12); imagesc(psf);
  
    if n>10
        zbasis.\ngs;
        residZDcomp(:,n-10)=zbasis.c;
        SumX2 = SumX2 + residZDcomp(:,n-10).^2;
        SumX = SumX + residZDcomp(:,n-10);
        STD= sqrt(SumX2/(n-10)-(SumX/(n-10)).^2); %standard deviation
        residSTD(:,n-10)=STD;
        %Strehl
%         science=science.*tel*dm*cam;
%         psf=cam.frame;
%         figure(12); imagesc(psf);
        
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


save(strcat(path, num2str(mag),'mag',num2str(c), 'c', num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes', num2str(gain),'gain', num2str(mod),'modulation','turbSTD', num2str(realization), '.mat'), 'turbSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation', 'residSTD',num2str(realization), '.mat'), 'residSTD');
save(strcat(path, num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','savePhase', num2str(realization), '.mat'),'savePhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveResidPhase',num2str(realization), '.mat'),'saveResidPhase');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveStrehl',num2str(realization), '.mat'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveSimStrehl', num2str(realization), '.mat'),'saveStrehl');
save(strcat(path,num2str(mag),'mag',num2str(c), 'c',num2str(nFaces), 'nFaces', num2str(altSlopes),'altSlopes',num2str(gain),'gain', num2str(mod),'modulation','saveCam',num2str(realization), '.mat'),'saveCam');

end
end
end