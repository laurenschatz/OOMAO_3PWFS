%% Zernike Modes

nModes=10;
nPx=240;
D=2;
npix=240*4;
zern=zernike(2:nModes+1, D, 'resolution', nPx); %zernike basis set
sz=size(zern.modes);
for i=1:sz(2)
    zernMode(:,:,i)=reshape(zern.modes(:,i),nPx,nPx);
end
%%
fm=fourierModes(8, nPx);

%%

for i=1:64
    fmMode(:,:,i)=reshape(fm.modes(:,i),nPx,nPx);
end
%%
% figure;
% for i=1:64
%     imagesc(fmMode(:,:,i))
%     pause(0.5)
% end
%%
circ=utilities.circle(240,240);
%%
x=1:240;
[X,Y]=meshgrid(x,x);

% %%
mask1=sin((X+Y)/4);

ave1=mean(mean(mask1));
mask1=mask1>ave1;
% 
% % mask1=fmMode(:,:,8);
% % ave1=mean(mean(mask1));
% % mask1=mask1>ave1;
x=circ;
x(nPx/2:end,:)=0;
x(:,1:nPx/2-1)=0;
mask1=mask1.*x;
figure; imagesc(mask1)
% 
% %% FIX THIS
mask2=sin((X-Y)/4);
ave2=mean(mean(mask2));
mask2=mask2<ave2;
% 
% mask2=fmMode(:,:,8);
% ave2=mean(mean(mask2));
% mask2=mask2>ave2;
x=circ;
x(1:nPx/2-1,:)=0;
x(:,1:nPx/2-1)=0;
mask2=mask2.*x;
figure; imagesc(mask2)
% 
% 
% %%
mask3=sin((Y/4));
ave3=mean(mean(mask3));
mask3=mask3>ave3;
% 
% % mask3=fmMode(:,:,5);
% % mask3=rot90(mask3);
% % ave3=mean(mean(mask3));
% % mask3=mask3>ave3;
x=circ;
x(:,120:end)=0;
mask3=mask3.*x;

% %% Mask 1: on axis circle
% 
% mask1=sin((Y/2));
% m=utilities.circle(240,60);
% ave1=mean(mean(mask1));
% mask1=mask1>ave1;
% mask1=mask1.*m;

% %% Mask 2: ring:
% % 
% % mask2=sin((X+Y)/4);
% % m2=utilities.circle(240,240);
% % m2=m2-m;
% % mask2=mask2>ave1;
% % mask2=mask2.*m2;
% % figure; imagesc(mask2)
% 
% %% Mask 2: 
% mask2=sin((X)/8);
% 
% ave2=mean(mean(mask2));
% mask2=mask2>ave2;
% mask2=mask2.*circ;
% 
% mTL=zeros(240,240);
% mTL(1:119,:)=1;
% x=m==0;
% 
% mask2=mask2.*mTL.*x;
% figure; imagesc(mask2)
% 
% %% Mask 3
% 
% mask3=sin((X+Y)/8);
% 
% ave3=mean(mean(mask3));
% mask3=mask3<ave3;
% mask3=mask3.*circ;
% 
% mBR=zeros(240,240);
% mBR(120:end,120:end)=1;
% x=m==0;
% 
% mask3=mask3.*mBR.*x;
% figure; imagesc(mask3)
% %%
% mask4=sin((X-Y)/8);
% 
% ave4=mean(mean(mask4));
% mask4=mask4>ave4;
% mask4=mask4.*circ;
% 
% mBR=zeros(240,240);
% mBR(1:119,120:end)=1;
% x=m==0;
% 
% mask4=mask4.*mBR.*x;
% figure; imagesc(mask4)
% %%
% mask5=sin((X-Y)/8);
% 
% ave5=mean(mean(mask5));
% mask5=mask5<ave5;
% mask5=mask5.*circ;
% 
% mBR=zeros(240,240);
% mBR(120:end,1:119)=1;
% x=m==0;
% 
% mask5=mask5.*mBR.*x;
% figure; imagesc(mask5)
%%

freq=100;
ngs = source('wavelength',photometry.R);
tel= telescope(D,'resolution',nPx,'samplingTime',1/freq);
cam=imager('diameter', tel.D);


phase=zernMode(:,:,5).*.5+zernMode(:,:,3).*.5;
amp= mask1+mask2+mask3;%+mask3+mask4+mask5;%mask2+mask3+mask4+mask5;%+mask1;
figure; imagesc(amp)
tel.pupil=amp;
%%

%ngs.magnitude=0;
%ngs.*tel*cam;
%figure; imagesc(cam.frame)
%%

phase=padarray(phase,[nPx*4, nPx*4],0,'both');
amp=padarray(amp,[nPx*4, nPx*4],0,'both');

pupil=amp.*exp(i.*phase.*0.01);
FTpupil=fftshift(fft2(fftshift(pupil)))/nPx;
psf=abs(FTpupil).^2;
figure; imagesc(sqrt(psf)); axis equal; axis  tight;
