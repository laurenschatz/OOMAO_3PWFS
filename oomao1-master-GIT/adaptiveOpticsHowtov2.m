
%% ADAPTIVE OPTICS MODELING WITH OOMAO
% Demonstrate how to build a simple closed-loop single conjugated adaptive
% optics system with both a geometric diffractive Shack-Hartmann WFS or a
% Pyramid WFS

%% Choose which WFS mode to use
wfsModel = 'diff'; % Options: 'diff', 'geom'

%% Definition of the atmosphere 


% The atmosphere class constructor has 2 required input:22
% 
% * the wavelength [m]
% * the Fried parameter for the previously given wavelength [m]
%
% 1 optionnal input: [m]
%
% * the outer scale
%
% and  parameter/value pairs of optional inputs:
%
% * the altitudes of the turbulence layers  [m]
% * the fractionnal r0 which is the ratio of the r0 at altitude h on the
% integrated r0: $(r_0(h)/r_0)^{5/3}$
% * the wind speeds  [m/s]
% * the wind directions [rd]
%
% In the following the atmosphere is given for a r0=15cm in V band and an
% outer scale of 30m with 3 turbulence layers.

r0 = 6e-2; %7 .8e-2
L0 = 6.25;
v = 7.4;

atm = atmosphere(photometry.V0,r0,L0,'windSpeed',v,'windDirection',0);

% atm = atmosphere(photometry.V,0.15,30,...
%     'altitude',[0,4,10]*1e3,...
%     'fractionnalR0',[0.7,0.25,0.05],...
%     'windSpeed',[5,10,20],...
%     'windDirection',[0,pi/4,pi]);

%% Define the camera

%% Definition of the telescope
% The telescope class constructor has 1 required input:
%
% * the telescope diameter [m]
%
% 1 optionnal input:
%
% * the central obscuration ratio
%
% and  parameter/value pairs of optional inputs:
%
% * the field of view either in arcminute or arcsecond
% * the pupil sampling or resolution in pixels
% * the atmopheric layer motion sampling time [s]
nPx = 112;

D=2;
freq=300;
tel = telescope(D,'resolution',nPx,'samplingTime',1/freq);

% tel = telescope(2,...
%     'fieldOfViewInArcMin',2.5,...
%     'resolution',nPx,...
%     'samplingTime',1/250);




cam = imager('diameter',tel.D);

%flush(cam)
%cam.frame = cam.frame*0;
%cam.clockRate    = 1;
%exposureTime     = 200;
%cam.exposureTime = exposureTime;
%startDelay       = 30;
%cam.startDelay   = startDelay;



%% Definition of a calibration source
% The source class constructor has parameters/value pairs of optional inputs:
%
% * the zenith angle [rd] (default: 0rd)
% * the azimuth angle [rd] (default: 0rd)
% * the wavelength [m] (default: V band)
% * the magnitude
%
% In the following, an on-axis natural guide star in V band is defined.
ngs = source('wavelength',photometry.R);

% and for later use a science object in H band is instantiated
science = source('wavelength',photometry.R);
science.*tel*cam;
cam.referenceFrame=cam.frame;

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
nLenslet =56;
d = tel.D/nLenslet;
%wfs = shackHartmann(nLenslet,nPx,0.5);

%pi*.3314725
wfs = pyramid(nLenslet,nPx,'modulation',5, 'alpha',0.825, 'nFaces',4, 'src',ngs,'tel',tel, 'c',2, 'rotation',0,'alternative',0, 'altSlopes',0);
%wfs = pyramid(nLenslet,nPx,'modulation',5, 'alpha',pi/2, 'rotation',0, 'nFaces',4, 'src',ngs,'tel',tel, 'c', 2, 'rotation',0);% increase modulation to avoid loss of performance due to small linear range
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
close all
figure; imagesc(wfs.camera.frame)
% The WFS camera display:
[cent]=utilities.pyrcen(wfs.camera.frame);
figure
subplot(1,2,1)
imagesc(wfs.camera)
%%
% The WFS slopes display:
subplot(1,2,2)
slopesDisplay(wfs)

% GAIN CALIBRATION
if ~strcmp(wfs.tag,'PYRAMID')
    wfs.pointingDirection = zeros(2,1);
    pixelScale = ngs.wavelength/...
        (2*d*wfs.lenslets.nyquistSampling);
    tipStep = pixelScale/2;
    nStep   = floor(nPx/nLenslet/3)*2;
    sx      = zeros(1,nStep+1);
    u       = 0:nStep;
    wfs.camera.frameListener.Enabled = false;
    wfs.slopesListener.Enabled = false;
    warning('off','oomao:shackHartmann:relay')
    for kStep=u
        ngs.zenith = -tipStep*kStep;
        +ngs;
        drawnow
       % sx(kStep+1) = median(wfs.slopes(1:end/2));
    end
    warning('on','oomao:shackHartmann:relay')
    Ox_in  = u*tipStep*constants.radian2arcsec;
    Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
    %figure
    plot(Ox_in,Ox_out)
    grid
    slopesLinCoef = polyfit(Ox_in,Ox_out,1);
    wfs.slopesUnits = 1/slopesLinCoef(1);
    
    ngs.zenith = 0;
    wfs.pointingDirection = [];
end

%% Definition of the deformable mirror
% The influence functions are made of two cubic Bezier curves. This
% parametric influence function allows modeling a large range of influence
% function types. As examples two influence functions are pre--defined, the
% "monotonic" and "overshoot" models. The second parameter is the
% mechanical coupling between two adjacent actuators
bif = influenceFunction('monotonic',30/100);
%%
% Cut of the influence function
%figure
%show(bif)
%%
% The deformableMirror constructor class has 1 required input:
%
% * the number of actuator on one side of the square array of actuators
%
% and  parameter/value pairs of optional inputs:
% 
% * the influence function object (modes)
% * the influence function resolution in pixels
% * the map of valid actuator
%
% Here, the actuators to lenslets registration obeys Fried geometry so the
% map of valid actuator in the square can be retrieved from the WFS object
nActuator = 9;
v=utilities.circle(nActuator,nActuator);

zern=zernike(2:sum(sum(v))+10, tel.D, 'resolution',nPx);





%needs to be modes
% dm = deformableMirror(nModes,...
%     'modes',zern.modes,...
%     'resolution',nPx,...
%     'validActuator',logical(v));
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',logical(v));




%% Interaction matrix: DM/WFS calibration
ngs=ngs.*tel;
dmCalib = calibration(dm,wfs,ngs,ngs.wavelength/100,10);

if strcmp(wfsModel,'geom')
    wfsG = geomGrad(wfs.validLenslet, tel.pupilLogical);
    dmCalibG = calibration(dm,wfsG,ngs,ngs.wavelength/40);
end
% The influence functions are normalized to 1, the actuator are then
% controlled in stroke in meter, here we choose a half a wavelength stroke.
% stroke = ngs.wavelength/2;
% %%
% % The DM actuator commands or coefficients is set to an identity matrix
% % scaled to the required stroke; each column of the matrix is one set of
% % actuator coefficients (commands)
% dm.coefs = eye(dm.nValidActuator)*stroke;
% %%
% % Propagation of the source through the telescope and the DM to the WFS
% ngs=ngs.*tel*dm*wfs;
% %%
% % The source has been propagated through the DM as many times as the number
% % of column in the DM coefficients matrix. As a result, the slopes in the
% % WFs object is also a matrix, each column correspond to one actuactor. The
% % interaction matrix is then easily derived from the slopes:
% interactionMatrix = wfs.slopes./stroke;
% figure(10)
% subplot(1,2,1)
% imagesc(interactionMatrix)
% xlabel('DM actuators')
% ylabel('WFS slopes')
% ylabel(colorbar,'slopes/actuator stroke')
% 
% 
% %% Command matrix derivation
% % The command matrix is obtained by computing first the singular value
% % decomposition of the interaction matrix,
% [U,S,V] = svd(interactionMatrix);
% eigenValues = diag(S);
% subplot(1,2,2)
% semilogy(eigenValues,'.')
% xlabel('Eigen modes')
% ylabel('Eigen values')
% %%
% % the 4 last eigen values are filtered out
% nThresholded = 4;
% iS = diag(1./eigenValues(1:end-nThresholded));
% [nS,nC] = size(interactionMatrix);
% iS(nC,nS) = 0;
% %%
% % and then the command matrix is derived.
% commandMatrix = V*iS*U';

dmCalib.nThresholded = 5;
if strcmp(wfsModel,'geom')
    commandMatrix = dmCalibG.M;
else
    commandMatrix = dmCalib.M;
end

%%
%
%
%
%
%
%
%
%
%% Noise

wfs.camera.photonNoise=true;
wfs.camera.readOutNoise=0.5;
ngs.magnitude=8;
ngs=ngs.*tel*wfs;


%% The closed loop
% Combining the atmosphere and the telescope
tel = tel+atm;
figure(10)
imagesc(tel)
%%
% Resetting the DM command
dm.coefs = 0;
%%
% Propagation throught the atmosphere to the telescope
ngs=ngs.*tel;
%%
% Saving the turbulence aberrated phase
turbPhase = ngs.meanRmPhase;
%%
% Propagation to the WFS
ngs=ngs*dm*wfs;
%%
% Display of turbulence and residual phase
figure(11)
h = imagesc([turbPhase,ngs.meanRmPhase]);
axis equal tight
colorbar
snapnow
%%
% Closed loop integrator gain:
loopGain = 0.5;
gain_pol = 0.9;
%%
% closing the loop with an integrator controller 
nIteration = 300;
total  = zeros(1,nIteration);
residue = zeros(1,nIteration);
tic


zbasis=zernike(2:90, tel.D, 'resolution',nPx);

SumX2_turb=zeros(89,1);
SumX_turb=SumX2_turb;
SumX2=SumX2_turb;
SumX=SumX2_turb;

for kIteration=1:nIteration
    % Propagation throught the atmosphere to the telescope, +tel means that
    % all the layers move of one step based on the sampling time and the
    % wind vectors of the layers
    ngs=ngs.*+tel;
    % Saving the turbulence aberrated phase
    zbasis.\ngs;
    fc(:,kIteration)=zbasis.c;
    SumX2_turb = SumX2_turb + fc(:,kIteration).^2;
    SumX_turb = SumX_turb + fc(:,kIteration);
    STD_turb= sqrt(SumX2_turb/kIteration-(SumX_turb/kIteration).^2);
 
    figure(22); imagesc(cam.frame);
    turbPhase = ngs.meanRmPhase; %CHECK THIS
    % Variance of the atmospheric wavefront
    total(kIteration) = var(ngs);

    savePhase(:,:,kIteration)=turbPhase;
    
    % Propagation to the dm/WFS
   ngs=ngs.*tel*dm*wfs;
  
    
    % Variance of the residual wavefront
    residue(kIteration) = var(ngs);
    saveResidPhase(:,:,kIteration)=ngs.meanRmPhase;
    
    if kIteration>10
    zbasis.\ngs;
    rc(:,kIteration)=zbasis.c;
    SumX2 = SumX2 + rc(:,kIteration).^2;
    SumX = SumX + rc(:,kIteration);
    STD_resid = sqrt(SumX2/kIteration-(SumX/kIteration).^2);
    end
    
    
    % Computing the DM residual coefficients

    residualDmCoefs = commandMatrix*wfs.slopes;
   
    % Integrating the DM coefficients
    dm.coefs = dm.coefs - loopGain*residualDmCoefs;
    
    % Strehl
    science.*tel*dm*cam;
    saveStrehl(:,kIteration)=cam.strehl;
    saveCam(:,:,kIteration)=cam.frame;
   
    
    
    % POL 
    %dm.coefs = (1-gain_pol)*dm.coefs - gain_pol*commandMatrix*(wfs.slopes - dmCalib.D*dm.coefs);

    % Display of turbulence and residual phase
    set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
    drawnow
    
    cam.flush();
    stdPhase(kIteration) = std2(savePhase(:,:,kIteration));
    stdResidPhase(kIteration) = std2(saveResidPhase(:,:,kIteration));
        figure(2)
    subplot(1,2,1)
    plot(stdPhase*1e9*science.wavelength/(2*pi))
    hold on
    plot(stdResidPhase*1e9*ngs.wavelength/(2*pi));
    set(gca, 'YScale', 'log')
    hold off
    


    
end
%% Standard Deviation
    figure;
    plot(STD_turb*1e6,'b-')
    hold on
    plot(STD_resid*1e6,'r-')
    
%% Varience
figure; 
plot(log((STD_turb*1e6).^2),'b-'); hold on
plot(log((STD_resid*1e6).^2),'r-')
xlabel('Zernike Mode Number')
ylabel('Variance (log micro^2)')

%% Strehl
figure; plot(saveStrehl)
mean(saveStrehl)


    
    
%%    
toc
snapnow
u = (0:nIteration-1).*tel.samplingTime;
atm.wavelength = ngs.wavelength;

% %% Power Spectrum Analysis%%
% 
% 
% for i=1:nIteration
%     zbasis.\saveResidPhase(:,:,i)
%     c(:,i)=zbasis.c;
% end
% 
% figure; plot(log(var(c,0,2))); hold on
% 
% for i=1:nIteration
%     zbasis.\savePhase(:,:,i)
%     t(:,i)=zbasis.c;
% end
% 
% plot(log((var(t,0,2))));





% % Piston removed phase variance
% totalTheory = phaseStats.zernikeResidualVariance(1,atm,tel);
% atm.wavelength = photometry.V;
% %%
% % Phase variance to micron rms converter 
% rmsMicron = @(x) 1e6*sqrt(x).*ngs.wavelength/2/pi;
% figure(12)
% plot(u,rmsMicron(total),u([1,end]),rmsMicron(totalTheory)*ones(1,2),u,rmsMicron(residue))
% grid
% legend('Full','Full (theory)','Residue','Location','Best')
% xlabel('Time [s]')
% ylabel('Wavefront rms [\mum]')
% 
% %% THEORETICAL PERFORMANCE ANALYSIS - CASE OF THE SINGLE INTEGRATOR LOOP WITH GAIN
% varFit      = 0.23*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
% varAlias    = 0.07*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
% varTempo    = phaseStats.closedLoopVariance(atm, tel.samplingTime,0*tel.samplingTime,loopGain)*(atm.wavelength/science(1).wavelength)^2;
% marechalStrehl_lsq_theoretical = 100*exp(-varFit-varAlias-varTempo)
% marechalStrehl_lsq_emperical=100*exp(-(mean(rmsMicron(residue(50:end)))*1e-6/ngs.wavelength*2*pi)^2*(atm.wavelength/science(1).wavelength)^2)
% 
% figure(12)
% text(0.5,1,['SR (Marechal approx):' num2str(marechalStrehl_lsq_theoretical) '%'])
% text(0.5,0.8,['Empirical (Marechal approx):' num2str(...
%     100*exp(-(mean(rmsMicron(residue(50:end)))*1e-6/ngs.wavelength*2*pi)^2*(atm.wavelength/science(1).wavelength)^2)) '%'])
% % %ylim([0 1])
% 
% %% WFS noise
% % Noise can be added to the wavefront sensor but first we need to set the
% % star magnitude.
% ngs.magnitude = 10;
% %%
% % It can be useful to know the number of photon per subaperture. To do so,
% % let separate the atmosphere from the telescope
% tel = tel - atm;
% %%
% % re-propagate the source,
% ngs = ngs.*tel*wfs;
% %%
% % and display the subaperture intensity
% figure
% intensityDisplay(wfs)
% %%
% % Now the readout noise in photo-electron per pixel per frame rms is set
% wfs.camera.readOutNoise = 5;
% %% 
% % Photon-noise is enabled.
% wfs.camera.photonNoise = true;
% %%
% % A pixel threshold is defined
% wfs.framePixelThreshold = 5;
% %%
% % The loop is closed again
% nIteration = 200;
% total  = zeros(1,nIteration);
% residue = zeros(1,nIteration);
% dm.coefs = 0;
% tel = tel + atm;
% for kIteration=1:nIteration
%     % Propagation throught the atmosphere to the telescope, +tel means that
%     % all the layers move of one step based on the sampling time and the
%     % wind vectors of the layers
%     ngs=ngs.*+tel; 
%     % Saving the turbulence aberrated phase
%     turbPhase = ngs.meanRmPhase;
%     % Variance of the atmospheric wavefront
%     total(kIteration) = var(ngs);
%     % Propagation to the WFS
%     ngs=ngs*dm*wfs; 
%     % Variance of the residual wavefront
%     residue(kIteration) = var(ngs);
%     % Computing the DM residual coefficients
%     residualDmCoefs = commandMatrix*wfs.slopes;
%     % Integrating the DM coefficients
%     dm.coefs = dm.coefs - loopGain*residualDmCoefs;
%     % Display of turbulence and residual phase
%     set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
%     drawnow
% end
% %%
% % Updating the display
% set(0,'CurrentFigure',12)
% hold on
% plot(u,rmsMicron(total),'b--',u,rmsMicron(residue),'r--')
% legend('Full','Full (theory)','Residue','Full (noise)','Residue (noise)','Location','Best')
