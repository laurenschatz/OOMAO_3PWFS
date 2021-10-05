master=load('FrameAndSlopesFourierModes16to81.mat');
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
%% x direction
for i=1:sz(3)
maxx=max(max(mx));
maxx=squeeze(maxx);

[x,y]=find(mx==maxx(i));
kx(i)=(x(2)-x(1))/2;
lx(i)=(y(2)-y(1))/2;   
end

%%
lx(1:36)=lx(1:36)*-1;
x=lx(:,:)+abs(min(lx))+1;
y=kx(:,:)+abs(min(kx))+1;

%%
szx=size(x);
modesx=zeros(9,9);
for i=1:szx(2)
        modesx(x(i),y(i))=maxx(i);
end

%% repeat for y direction
for i=1:sz(3)
maxy=max(max(my));
maxy=squeeze(maxy);

[x,y]=find(my==maxy(i));
sx=size(x)
sy=size(y)
if sx(1) ==2 && sy(1)==2;
ky(i)=(x(2)-x(1))/2;
ly(i)=(y(2)-y(1))/2;  

else
    ky(i)=0;
    ly(i)=0;
end
end

%%
ly(1:36)=ly(1:36)*-1;
x=ly(:,:)+abs(min(ly))+1;
y=ky(:,:)+abs(min(ky))+1;

%%
szy=size(x);
modesy=zeros(9,9);
for i=1:szy(2)
        modesy(x(i),y(i))=maxy(i);
end
