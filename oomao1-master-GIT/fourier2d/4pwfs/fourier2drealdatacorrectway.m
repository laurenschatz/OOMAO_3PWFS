master=load('FrameAndSlopes3WithFluxFourierModes1to81.mat');
S=master.slopes;
Flux=master.Flux;
l=reshape(fmode.l, 9,9);
k=reshape(fmode.k,9,9);

sz=size(S);
for i=1:sz(3)
Fx=S(:,1:20,i);
Fy=S(:,21:40,i);

Fx=fftshift(fft2(fftshift(Fx)))/80;
Fy=fftshift(fft2(fftshift(Fy)))/80;

mx(:,:,i)=sqrt(Fx.*conj(Fx));
my(:,:,i)=sqrt(Fy.*conj(Fy));
end
%%
N=9;
counter=1;
for i=1:N
    for j=1:N 
Gx(i,j)= mx(11-l(i,j),11+k(i,j),counter)/Flux(1,counter);
Gy(i,j)=my(11-l(i,j),11+k(i,j),counter)/Flux(1,counter);
counter=counter+1;
    end
end
%% Plots
figure; imagesc(Gx)
title({'Measured X response- Fourier Domain'},'FontSize',18);
xlabel({'K-mode number'},'FontSize',18);
yticklabels({-4,-3,-2,-1,0,1,2,3,4})
xticklabels({-4,-3,-2,-1,0,1,2,3,4})
ylabel({'L-mode numbler'},'FontSize',18);
colorbar
caxis([10e-10 10e-7])

figure; imagesc(Gy)
title({'Measured Y response- Fourier Domain'},'FontSize',18);
xlabel({'K-mode number'},'FontSize',18);
yticklabels({-4,-3,-2,-1,0,1,2,3,4})
xticklabels({-4,-3,-2,-1,0,1,2,3,4})
ylabel({'L-mode numbler'},'FontSize',18);
colorbar
caxis([10e-10 10e-7])



