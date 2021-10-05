%%% Psf maker

fmode=fourierModes(15, 400);

sz=size(fmode.modes);
for i=1:225
    fmPhase(:,:,i)=reshape(fmode.modes(:,i),sqrt(sz(1)),sqrt(sz(1)));
    %figure; imagesc(fmPhase(:,:,i)); 
   
end
%%
Wi=fmPhase(:,:,4);
B=padarray(Wi, [200,200],0,'both');
E = exp(1i*2*pi*B*0.1);
mask=utilities.circle(800,200);
E=E.*mask;
psf = abs(fftshift(fft2(ifftshift(E)))).^2;
figure; imagesc(psf)