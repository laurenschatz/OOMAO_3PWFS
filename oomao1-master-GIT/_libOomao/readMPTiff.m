function [imageCube,nImage] = readMPTiff(filename)

% reads multiple tiff file, written by Hamamatsu camera

% reads number of frames

structure = imfinfo(filename);

sizeX = structure(1).Width;

sizeY = structure(1).Height;

nImage = length(structure);

%imageCube = zeros(sizeX, sizeY, nImage);
imageCube = zeros(sizeY, sizeX, nImage);

w=waitbar(0/nImage,'Reading ...');
for ind=1:nImage

imageCube(:,:,ind) = imread(filename, ind);
waitbar(ind/nImage,w)
end

end
