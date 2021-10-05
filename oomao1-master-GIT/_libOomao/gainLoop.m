%% loop script
function [finalStrehl,gain,M]=gainLoop(path, mag,c,nFaces,altSlopes, mod, alternative)
counter=0;
for i=0.5:0.1:2.5
    counter=counter+1;
    Strehl=pwfsGains(mag, c, nFaces, altSlopes, i, mod, alternative);
    finalStrehl(:,counter)=Strehl;
    index(counter)=i;
end

[M,I]=max(finalStrehl);
M;
gain=index(I);

pwfsexperiment(path, mag, c, nFaces, altSlopes, gain, mod, alternative);

end