    
mag='8';
fileID=fopen(strcat('/home/lschatz/MATLAB/lschatz/oomao-remote/Data/resultsHighReadNoiseExp/4PWFSRIresults/Mag',num2str(mag),'.txt'), 'w');

for g=0.1:0.1:0.8
    
    gain=g;
    nFaces='4';

    z=strcat(mag,'mag2c', nFaces,'nFaces2altSlopes',num2str(gain),'gain5modulation40readNoisesaveStrehl');
    for i=1:15
        name=strcat(z,num2str(i), '.mat');
        x=load(name);
        s(i)=mean(x.saveStrehl);
    end
    x=mean(s)
%fname1=strcat(mag,'mag',gain,'gainMeanStrehl.mat');
%fname2=strcat(mag,'mag',gain,'gainStrehl.mat');
%save(fname1, 'x')
%save(fname2,'s')
s=s*100;

    fprintf(fileID, '%6s\r\n',strcat('Strehl gain ', num2str(g)));
    fprintf(fileID, '%6.2f\r\n', s)
    %for i=1:size(s)
        %fprintf(fileID ,'%6.2s',strcat(num2str(s(i)),','));
    %end



end

fclose(fileID)
