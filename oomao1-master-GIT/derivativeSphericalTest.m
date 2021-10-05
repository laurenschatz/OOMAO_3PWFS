
x= 0:0.1:2*pi;
y=x;

zmode = zernike(1:15, 2,'resolution', 500);
%fmode=fourierModes(9, tel.resolution);
mode=zmode.modes(:,11);

sz=size(mode);
fmPhase=reshape(mode,sqrt(sz(1)),sqrt(sz(1)));
figure; imagesc(fmPhase); axis equal

%%
[dx,dy]=gradient(fmPhase.*1000);
figure; imagesc(dy)

%%
c=utilities.circle(500, 496);
figure; imagesc(c)

%%
figure; imagesc(dx.*c); axis equal
title('X Gradient of Spherical Aberration', 'FontSize',24)
%%
figure; imagesc(dy.*c); axis equal
title('Y Gradient of Spherical Aberration', 'FontSize',24)
figure; imagesc(fmPhase.*c)
title('Spherical Aberration', 'FontSize',24)




