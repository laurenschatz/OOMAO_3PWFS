clear all; close all

%% Step  0 parameters
% npix=60; %Number of pixels across a pupil
nLenslet=30;
resolution=360;
lambda=700*10^-9;%wavelength
%stroke=lambda/100; %Deformable mirror stroke. (OPD error 1/100 wave)
%amp=(2*pi/lambda)*stroke; %Amplitude of error of the fourier modes.
Amp=(2*pi/lambda)*1*10^-9; %wavefront error in radians.
%% Step One define the telecsope
nPx = 360;
tel = telescope(3.6,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/250);
cam = imager('diameter',tel.D);

atm=atmosphere(photometry.I1, 0.20);



%% Step Two define calibration source
ngs = source('wavelength',photometry.V);

%% Step Three Define the wavefront sensor

%pi*.3315
wfs = pyramid(30,360,'modulation',1, 'alpha',pi/2, 'nFaces',4, 'src',ngs,'tel',tel, 'c', 2, 'rotation',3*pi/2, 'altSlopes',0,'alternative', 0);
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
%Display the WFS
figure; imagesc(wfs.camera)
% % WFS slopes display
% subplot(1,2,2)
% slopesDisplay(wfs)


%% Step Four Propogate
ngs=ngs.*tel*wfs;

Frec=fourierReconstructor(wfs,tel,atm);
Frec.measureFilter(ngs,1,1);
%%
Ez=reshape(Frec.E,360,360,[]);
Gxz=reshape(Frec.GX, [],1);
Gyz=reshape(Frec.GY,[],1);
Lz=reshape(Frec.ll, [],1);
Kz=reshape(Frec.kk, [],1);


%%
figure; imagesc(abs(Frec.GX)); axis equal; title('X-direction')

figure; imagesc(abs(Frec.GY)); axis equal; title('Y-direction')

%%
fmode=fourierModes(11, tel.resolution); %Fourier Modes
%%
l=fmode.l;
l=reshape(l,11,11);
k=fmode.k;
k=reshape(k,11,11);
%% Get rid of piston
sz=size(fmode.modes);
% a=find(fmode.l==0);
% b=find(fmode.k==0);
% c=find(a==b);
% d=a(c);

for i=1:sz(2)
    fmPhase(:,:,i)=reshape(fmode.modes(:,i),nPx,nPx)*Amp;
end
% dummy=fmPhase(:,:,1:d-1);
% dummy2=fmPhase(:,:,d+1:sz(2));
% 
% dl=fmode.l(1:d-1);
% dl2=fmode.l(d+1:sz(2));
% 
% dk=fmode.k(1:d-1);
% dk2=fmode.k(d+1:sz(2));
% 
% fmPhase=cat(3,dummy,dummy2);
% l=cat(2,dl,dl2);
% k=cat(2,dk,dk2);


%%

for i=1:sz(2)
theFourierMode={tel.pupil,fmPhase(:,:,i)};
ngs=ngs.*tel*theFourierMode*wfs;

%% Step Seven Pull Out the Slopes Maps PWFS
szsm=size(wfs.slopesMap);
x=wfs.slopesMap(1:szsm(1),1:szsm(2)/2);
y=wfs.slopesMap(1:szsm(1),szsm(2)/2+1:szsm(2));

Sx(:,:,i)=x(szsm(1)/4+1:szsm(1)/2+szsm(1)/4,szsm(1)/4+1:szsm(1)/2+szsm(1)/4);
%szx=size(Sxx);
Sy(:,:,i)=y(szsm(1)/4+1:szsm(1)/2+szsm(1)/4,szsm(1)/4+1:szsm(1)/2+szsm(1)/4);
%szy=size(Syy);
%Sy(:,:,i)=Syy(szx(1)/2-40:szx(1)/2+40,szx(1)/2-40:szx(1)/2+40);

%% Step Seven Pull Out the Slopes SH
% Sx(:,:,i)=wfs.xSlopesMap;
% Sy(:,:,i)=wfs.ySlopesMap;

%% Fourier Transform
sz=size(Sx);
Fx=fftshift(fft2(fftshift(Sx(:,:,i))))/(sz(1));
Fy=fftshift(fft2(fftshift(Sy(:,:,i))))/(sz(1));

%% Step Eight: Take the magnitude.
Mx(:,:,i)=sqrt(Fx.*conj(Fx));
My(:,:,i)=sqrt(Fy.*conj(Fy));

%figure; imagesc(Mx(:,:,i))
%% calc spatial frequency
%K(i)= sqrt(((k(i).^2)/nPx)+((l(i).^2)/nPx)/2);
end
%%
sz=size(Mx);
N=11;
counter=1;
for i=1:N
    for j=1:N 
Gx(i,j)= Mx(sz(1)/2+1-l(i,j),sz(1)/2+1+k(i,j),counter);        
Gy(i,j)=My(sz(1)/2+1-l(i,j),sz(1)/2+1+k(i,j),counter);
counter=counter+1;
    end
end

figure; imagesc(Gx)
figure; imagesc(Gy)


% Ak=sqrt((maxx.^2+maxy.^2)/2);
% %combine like spatial frequences
% uniq=unique(K);
% ssz=size(uniq);
% for i=1:ssz(2)
%     I=find(K==uniq(i));
%     sz=size(I);
%     for j=1:sz(2)
%         values(j)=Ak(I(j));
%     end
%     sv=size(values);
%     AK(i)=sqrt(sum(values.^2)/sv(2));
% end
% %%
% figure; scatter(K,Ak/Amp)
% 
% %%
% c=polyfit(uniq,AK,1);
% %disp(['Equation is y= ' num2str(c(1)) '*x^2+' num2str(c(2)) '*x+' num2str(c(3))]);
% disp(['Equation is y= ' num2str(c(1)) '*x+' num2str(c(2))])
% Ak_est=polyval(c,uniq);
% 
% figure; scatter(uniq,AK); hold on
% %plot(uniq,Ak_est)
% %hold off
% 
% % %
% % 
% % sz=size(My);
% % v=VideoWriter('MySlopesP.avi');
% % v.FrameRate=10;
% % open(v)
% % figure;
% % for i=1:sz(3)
% %     imagesc(My(:,:,i))
% %     set(gca,'nextplot','replacechildren')
% %     frame=getframe(gcf);
% %     writeVideo(v,frame)
% % end
% % close(v)
% % 
% % %%
% % %
% % 
% % sz=size(Mx);
% % v=VideoWriter('MxSlopesP.avi');
% % v.FrameRate=10;
% % open(v)
% % figure;
% % for i=1:sz(3)
% %     imagesc(Mx(:,:,i))
% %     set(gca,'nextplot','replacechildren')
% %     frame=getframe(gcf);
% %     writeVideo(v,frame)
% % end
% % close(v)
% % %
% % sz=size(Sy);
% % v=VideoWriter('SySlopesP.avi');
% % v.FrameRate=10;
% % open(v)
% % figure;
% % for i=1:sz(3)
% %     imagesc(Sy(:,:,i)); caxis([-1,1]) 
% %     set(gca,'nextplot','replacechildren')
% %     frame=getframe(gcf);
% %     writeVideo(v,frame)
% % end
% % close(v)
% % %
% % sz=size(Sx);
% % v=VideoWriter('SxSlopesP.avi');
% % v.FrameRate=10;
% % open(v)
% % figure;
% % for i=1:sz(3)
% %     imagesc(Sx(:,:,i)); caxis([-1,1]) 
% %     set(gca,'nextplot','replacechildren')
% %     frame=getframe(gcf);
% %     writeVideo(v,frame)
% % end
% % close(v)
% % %