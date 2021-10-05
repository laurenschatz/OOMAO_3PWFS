%function obj = pyramidSLM(obj, nFaces, angle, centrePos, rotation, rooftop)
centrePos= [512.5 640.5];
nFace=3;
angle=-1.25*pi/8;
rotation=0;
            
            nx = 1024;
            ny = 1280;
            
            cx = centrePos(1);
            cy = centrePos(2);
            
            x = linspace(-cx+1, nx-cx, nx);
            y = linspace(-cy+1, ny-cy, ny);
            [xGrid, yGrid] = meshgrid(x, y);
            xGrid = xGrid';
            yGrid = yGrid';
            
            % DEFINE THE SLOPE GRID
            angleGrid = atan2(yGrid*sin(-rotation) + ...
                xGrid*cos(-rotation), ...
                yGrid*cos(-rotation) - ...
                xGrid*sin(-rotation));
            
 %%  Make Binary 3 mask         

s1=zeros(nx,ny);
s2=zeros(nx,ny);
s3=zeros(nx,ny);
            
s1((pi/6 <= angleGrid) & (angleGrid < pi/2)) = 1;
s2((pi/2 <= angleGrid) & (angleGrid < 5*pi/6)) = 1; 
s=s1+s2;
s3=s<1;

%% Apply 4PWFS tilts

                    
            nFaces_ = 4;            
            % INITIALIZE PYRAMID MASK
            pyp = zeros(nx, ny);

            for kFaces=0:nFaces_-1
                theta = (pi*(1/nFaces_ - 1) + kFaces*2*pi/nFaces_ + rotation);
                slope = (sin(theta)*xGrid + cos(theta)*yGrid);
                %Take into account the last tile of the pyramid mask
                if kFaces == 2
                    slope=-xGrid.*s3; 
                    pyp = pyp+angle*slope;
                else
                    if kFaces==0   
                    slope=slope.*s1;
                    pyp = pyp+angle*slope;
                    end
                    if kFaces==1
                    slope=slope.*s2;
                    pyp = pyp+angle*slope;
                    pyp=-pyp;
                    end
                end


            end
            
            
            
            
               
   
