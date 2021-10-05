master=load('FrameAndSlopesFourierModes16to81.mat');
l=reshape(fmode.l, 9,9);
k=reshape(fmode.k,9,9);
S=master.slopes;
sz=size(S);
for i=1:sz(3)
Fx=S(:,1:80,i);
Fy=S(:,81:160,i);

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
Gx(i,j)= mx(41-l(i,j),41+k(i,j),counter);        
Gy(i,j)=my(41-l(i,j),41+k(i,j),counter);
counter=counter+1;
    end
end
