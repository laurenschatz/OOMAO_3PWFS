        function pyp = pyramidSLM(resolution,nFaces, angle, centrePos, rotation)

            nx = resolution(1);
            ny = resolution(2);
            
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
            
            % INITIALIZE PYRAMID MASK
            pyp = zeros(nx, ny);
            s2=zeros(nx,ny);
            for kFaces=0:nFaces-1
                theta = pi*(1/nFaces - 1) + kFaces*2*pi/nFaces + rotation;
                slope = sin(theta)*xGrid + cos(theta)*yGrid;

                % Take into account the last tile of the pyramid mask
                if kFaces == nFaces-1
                    slope((-pi+kFaces*2*pi/nFaces <= angleGrid) &...
                        (angleGrid <= (-pi+(kFaces + 1)*2*pi/nFaces))) = 1; 
                else
                    slope((-pi+kFaces*2*pi/nFaces <= angleGrid) &...
                        (angleGrid <= (-pi+(kFaces + 1)*2*pi/nFaces))) = 1;

                end

                pyp = pyp + angle*slope;
            end
            
%             obj.phaseMask = pyp;
%             pym = exp(1i*pyp);
%             obj.theMask = fftshift(pym)./sum(abs(pym(:)));
        end
