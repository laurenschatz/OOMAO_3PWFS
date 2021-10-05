%OOMAO Fourier Modes Test Script

%ngs=source
%ngs=ngs*telescope

%F=FourierModes(n)
%ngs.phase=F(i)
nSamp=20;
mod=5;
c=2;
altSlopes=0;
alternative=0;
nFaces=4;

%Telescope Diameter
D=2;
%Full Resolution 
nPx=240;

% Sampling Frequency (exposure time)
freq=300;
%nPix=252; %should be at least twice the number of modes for Nyquist
nModes=6;
nLenslet=16;
%resolution=6*nLenslet;
lambda=700*10^-9;
error=10000*10^-9;

%% Define the Source 
ngs=source('wavelength', photometry.V, 'magnitude',0); %Plane wave, photometric band in V, 0 mag

%Create the phase
f=fourierModes(nModes,nPx); %generates Nmodes^2 number of fourier modes in a Npix roster
%pull and reshape into a data cube
sz=size(f.modes);

for n=1:sz(2)
    modes(:,:,n)=reshape(f.modes(:,n),nPx,nPx);
    
end

k=f.k; %y-spatial frequency index
l=f.l; %x-spatial frequency index
CS=f.CS; % Positive or Negative for sin/cosine
    
%% Define the Telescope

%Creates a telescope object. Diameter 6.5 meters, field of view 0.35 arcmin
%(21 arcseconds), resolution of 128 pixels across the pupil, and a sampling
%time of 1/10th a second. 
tel=telescope(6.5, 'fieldOfViewInArcMin', 0.35, 'resolution', nPx, 'samplingTime', 1/10);


%% Define Shack Hartmann

%wfs=shackHartmann(nLenslet,npix, 0.75);
wfs = pyramid(nSamp,nPx,'modulation',mod, 'alpha',pi/3, 'nFaces',nFaces, 'src',ngs,'tel',tel, 'c', c, 'rotation',3*pi/2, 'altSlopes',altSlopes,'alternative', alternative);


%% Propogate
% you have to propogate the wave through the telescope before applying a
% phase
ngs=ngs*tel;

%% Apply phase
ngs.phase=(modes(:,:,5)*2*pi*error/lambda);

%% propogate through wfs
ngs=ngs.*tel*wfs;

%% 
wfs.INIT
%%
+wfs
%%
figure; imagesc(wfs.camera)

%%
figure; slopesDisplay(wfs)





%Start Loop: Send in Fourier mode)
% 
% for i=1:size(k);
%         ngs.phase=modes(:,:,i)
% %assign a fourier mode
% ngs.phase=modes(:,:,5);


