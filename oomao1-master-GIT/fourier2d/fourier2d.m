
%% Optics Setup
% Wavefront sensor sampling
nSamp=50;

%Telescope Diameter
D=2;
%Full Resolution 
nPx=240;

% Sampling Frequency (exposure time)
freq=300;
mod=0;
c=2;
altSlopes=0;
alternative=0;
nFaces=4;
Amp=0.5;


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

nActuator=50;

fmode=fourierModes(9, tel.resolution); %Fourier Modes
%%
l=fmode.l;
l=reshape(l,9,9);
k=fmode.k;
k=reshape(k,9,9);
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

figure; imagesc(Mx(:,:,i))
%% calc spatial frequency
%K(i)= sqrt(((k(i).^2)/nPx)+((l(i).^2)/nPx)/2);
end
%%
sz=size(Mx);
N=9;
counter=1;
for i=1:N
    for j=1:N 
Gx(i,j)= Mx(sz(1)/2+1-l(i,j),sz(1)/2+1+k(i,j),counter);        
Gy(i,j)=My(sz(1)/2+1-l(i,j),sz(1)/2+1+k(i,j),counter);
counter=counter+1;
    end
end




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