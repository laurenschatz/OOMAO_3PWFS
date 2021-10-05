
%% ADAPTIVE OPTICS MODELING WITH OOMAO
% Demonstrate how to build a simple closed-loop single conjugated adaptive
% optics system with both a geometric diffractive Shack-Hartmann WFS or a
% Pyramid WFS

%% Choose which WFS mode to use
wfsModel = 'diff'; % Options: 'diff', 'geom'


atm=atmosphere(photometry.I1, 0.20);

nPx = 60;
tel = telescope(3.6,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/250);

% In the following, an on-axis natural guide star in V band is defined.
ngs = source('wavelength',photometry.V);

% and for later use a science object in H band is instantiated
science = source('wavelength',photometry.H);


%% Definition of the wavefront sensor
% Up to now, only the Shack--Hartmann WFS has been implemented in OOMAO.
% The shackHartmann class constructor has 2 required inputs:
%
% * the number of lenslet on one side of the square lenslet array
% * the resolution of the camera
%
% and  1 optionnal input:
%
% * the minimum light ratio that is the ratio between a partially
% illuminated subaperture and a fully illuminated aperture
nLenslet = 10;
d = tel.D/nLenslet;
%wfs = shackHartmann(nLenslet,nPx,0.5);
%wfs = pyramid(nLenslet,nPx,'modulation',5, 'alpha',pi*.3314725, 'rotation',0, 'nFaces',3, 'src',ngs,'tel',tel, 'c', 12, 'rotation',pi/2);
wfs = pyramid(nLenslet,nPx,'modulation',5, 'alpha',pi/2, 'rotation',0, 'nFaces',4, 'src',ngs,'tel',tel, 'c', 4, 'rotation',0);% increase modulation to avoid loss of performance due to small linear range
%wfs = pyramid(nLenslet,nPx,'modulation',5, 'src', ngs,'tel',tel,'minLightRatio', 0.05);
%%
% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*wfs;
%%
% Selecting the subapertures illuminated at 75% or more with respect to a
% fully illuminated subaperture
% setValidLenslet(wfs)
% %%
% % A new frame read-out and slopes computing:
% +wfs;
% %%
% % Setting the reference slopes to the current slopes that corresponds to a
% % flat wavefront
% wfs.referenceSlopes = wfs.slopes;
wfs.INIT
%%
% A new frame read-out and slopes computing:
+wfs;

%%
%% Step Four Propogate
ngs=ngs.*tel*wfs;


% input a Fourier mode

fmode=fourierModes(3, tel.resolution);

fmPhase=reshape(fmode.modes(:,5),nPx,nPx);
%%
theFourierMode={tel.pupil, fmPhase};

% insrt with multiple levels of aberration

ngs=ngs.*tel*theFourierMode*wfs;
%%
amp=-500:50:500;
sz=size(amp);
for i=1:sz(2)
   FourierMode(:,:,i)=fmPhase.*amp(i); 
        ngs=ngs.*tel*FourierMode(:,:,i)*wfs;
            pupils(:,:,i)=wfs.camera.frame;
   
end
%%




