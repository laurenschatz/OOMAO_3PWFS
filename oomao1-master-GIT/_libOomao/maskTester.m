%Pupil tester

% set up
nLenslet=20;
nPx=6*nLenslet;
c=10;


pupil=utilities.circle(c*nPx,nPx);


%% Make pyramid Mask

            nFaces_ =4;
            nx      = c*nPx;
            ny      = nx;
            alpha=pi/2;
            rotation=3*pi/2;
            %pwfs.rotation=0;
            
            cx = nx/2+0.5;
            cy = cx;
            
            %x = linspace(-cx+1, nx-cx, nx);%*pi/sqrt(2);
            %y = linspace(-cy+1, ny-cy, ny);%*pi/sqrt(2);
            %[xGrid, yGrid] = meshgrid(x, y);
            [xGrid,yGrid]=freqspace(nx,'meshgrid');
            xGrid = xGrid.*floor(nx/2);
            xGrid=xGrid.*alpha.*2/sqrt(2);
            yGrid = yGrid.*floor(ny/2);
            yGrid=yGrid.*alpha.*2/sqrt(2);

            % DEFINE THE SLOPE GRID for the mask
            angleGrid = atan2(yGrid*sin(rotation) + ...
                xGrid*cos(rotation), ...
                yGrid*cos(rotation) - ...
                xGrid*sin(rotation));

            
            % INITIALIZE PYRAMID MASK
            pyp = zeros(nx, ny);
            for kFaces=0:nFaces_-1
                theta = (pi*(1/nFaces_ - 1) + kFaces*2*pi/nFaces_ + rotation);
                slope = (sin(theta)*xGrid + cos(theta)*yGrid);
                %Take into account the last tile of the pyramid mask
                if kFaces == nFaces_-1
                 
                    slope((-pi+kFaces*2*pi/nFaces_ <= angleGrid) &...
                        (angleGrid <= (-pi+(kFaces + 1)*2*pi/nFaces_))) = 1;                    
                else
                   
                    slope((-pi+kFaces*2*pi/nFaces_ <= angleGrid) &...
                        (angleGrid < (-pi+(kFaces + 1)*2*pi/nFaces_))) = 1;
                end  
                  
                
                pyp = pyp+slope;
            end
           
            if nFaces_==3
                pyp=-pyp;
            end
            pyramidmask = exp(1i*pyp);

     





%%
FTpupil=fftshift(fft2(fftshift(pupil)))/nx;
pyramid=FTpupil.*pyramidmask;
Pupilpyramids=abs((ifftshift(fft2(ifftshift(pyramid))))/nx).^2;

