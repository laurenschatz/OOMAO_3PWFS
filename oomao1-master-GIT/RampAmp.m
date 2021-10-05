function [S,plus,minus,multAmp]=RampAmp(ngs,wfs,tel,Amp,fmPhase,rmatrix)
%%
mult=1:100:1000;
multAmp=mult.*Amp;
sz=size(multAmp);
sz2=size(fmPhase);

for j=1:sz2(3)
for i=1:sz(2)
plusfmPhase=fmPhase(:,:,j)*multAmp(i);
    plusFourierMode={tel.pupil, plusfmPhase};
        ngs=ngs.*tel*plusFourierMode*wfs;
            plusslopes=wfs.slopes;

minusfmPhase=fmPhase(:,:,j)*-multAmp(i);
       minusFourierMode={tel.pupil, minusfmPhase};
         ngs=ngs.*tel*minusFourierMode*wfs;
            minusslopes=wfs.slopes;
            

plusAmp=rmatrix*plusslopes;
    plus(j,i)=plusAmp(j);
minusAmp=rmatrix*minusslopes;
    minus(j,i)=minusAmp(j);
    
S(j,i)=(plus(j,i)-minus(j,i))/(multAmp(i)-(-multAmp(i)));    
end
end
%%

figure; plot(multAmp,S(2,:)); xlabel('Radians of error'); ylabel('Slope')
%%
figure; plot(multAmp,plus(2,:)); xlabel('Radians of error'); ylabel('measured error'); hold on
plot(-multAmp,minus(2,:));

end