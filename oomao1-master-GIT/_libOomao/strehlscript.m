
%[cubepsf,nFrame]=utilities.readMPTiff('CL_3PWFS_FF_ND-flat-15-29-03-19(0).tif');
%[cubepsf,nFrame]=utilities.readMPTiff('CL_3PWFS_FF_ND00-29-03-19(0).tif');
%[cubepsf,nFrame]=utilities.readMPTiff('CL_3PWFS_FF_ND-flat-15-29-03-19(0).tif');
%[cubepsf,nFrame]=utilities.readMPTiff('CL_4PWFS_FF_ND00-flat(0).tif');
%[cubepsf,nFrame]=utilities.readMPTiff('CL_3PWFS_FF_ND10=5-29-03-19(9).tif');
%[cubepsf,nFrame]=utilities.readMPTiff('CL_4PWFS_FF_ND00-flat-29-03-19(3).tif');
%[cubepsf,nFrame]=utilities.readMPTiff('CL_4PWFS_FF_ND00-flatbest-29-03-19(0).tif');

%%
%[cubepsf,nFrame]=utilities.readMPTiff('CL_4PWFS_FF_ND00-flat-r2-best-29-03-19(1).tif');
%[cubepsf,nFrame]=utilities.readMPTiff('CL_3PWFS_SM_ND00-flat-r2(0).tif');


[cubepsf,nFrame]=utilities.readMPTiff('CL_4PWFS_SM_ND20(9).tif');
%[cubepsf,nFrame]=utilities.readMPTiff('CL_4PWFS_best_flat_openloop[(1).tif');

%%
%cubepsf2=cubepsf(:,:,1:2000);
sz=size(cubepsf);
N=78; %diameter of perfect circle
perfPupil=utilities.circle(sz(1),N);
perfPsf=fftshift(fft2(fftshift(perfPupil)))/sz(1);
perfPsf=abs(perfPsf).^2;
%%
%experimental
for i=1:sz(3)
[s,o,p,r]=utilities.experimentalpsf2strehl(cubepsf(:,:,i),expbest);
%[s,o,p,r]=utilities.psf2strehl(cubepsf(:,:,i),perfPsf);
strehl(:,i)=s;
otf(:,i)=o;
perfotf(:,i)=p;
x(:,i)=r;
end

%%
m=mean(strehl)
dev=std(strehl)/sqrt(sz(3))
maxs=max(strehl)
mins=min(strehl)

