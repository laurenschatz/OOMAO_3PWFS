classdef lamTools
    % lamTools with shared functionality :: static class
    
    methods (Static)
        %% Estimate compressed profile using mean-weighted compression
        function [Cn2eq altEq] = eqLayers(Cn2, altitudes, nEqLayers, power)
            %{
            Cn2         ::  The input Cn2 profile (vector)
            altitudes   ::  The input altitudes (vector)
            nEqLayers   ::  The number of output equivalent layers (scalar)
            power       ::  the exponent of the turbulence (default 5/3)
            
            See: Saxenhuber17: Comparison of methods for the reduction of
            reconstructed layers in atmospheric tomography, App Op, Vol. 56, No. 10 / April 1 2017 
            %}
            nCn2 = numel(Cn2);
            nAltitudes = numel(altitudes);
            if nargin ~= 4
                power = 5/3;
            end
            
            % if nargin ~= 5
            nSlab = floor(round(nCn2)/fix(nEqLayers));
            ppp = 1;
            posSlab =  round((linspace(0, nEqLayers-1, nEqLayers))*nSlab)+1;
            for iii = 1:nEqLayers-1
                if posSlab(iii) >= posSlab(iii+1)
                    posSlab(iii+1) = posSlab(iii)+1;
                end
            end
            posSlab = [posSlab, nAltitudes+1];
            Cn2eq = zeros(1,nEqLayers);
            altEq = zeros(1,nEqLayers);
            for ii = 1:nEqLayers
                Cn2eq(ii) = sum(Cn2(posSlab(ii):posSlab(ii+1)-1));
                altEq(ii) =  (sum(altitudes(posSlab(ii):posSlab(ii+1)-1) .^ (power) .* Cn2(posSlab(ii):posSlab(ii+1)-1))/Cn2eq(ii)) .^ (1./power);
            end
        end
        
        
        %%
        % PSD: Power Spectral Density
        % Calculates PSD for circular aperture (i.e. not spiders, no central obscuration).
        % Uses a window for measuring PSD of variable [input] defined over 2D grid
        % (i.e. to smooth edges of the data) and data goes to the edges of
        % the 2D grid.
        function [psd, psdCirc, kx] = powerSpecrum(input, m, n)
            % Check that phase is square
            nRes = size(input, 1);
            w = fourierReconstructor.window('hann', nRes);
            psd = (abs(fftshift(fft2(input.*w))).^2).*sum(w(:).^2);
            
            % Calculer circular PSD
            if size(psd, 1) == size(psd, 2); % Input matrix must be square.
                if nargin < 2
                    [m,~] = size(psd);
                    centerm = ceil(m/2+1); %matrix is square , so use m or n
                    centern = centerm;
                else
                    centerm = n;
                    centern = m;
                    %[m,n] = size(psd);
                end
                [psdCirc, kx] = radial(psd, centerm, centern);
            else
                psdCirc = 0;
                kx = 0;
                disp('Not a square array, cannot calculate circular PSD');
            end
        end 
        
        %% distances
        function rho = dist(l,u, nPts)
            [x,y] = meshgrid(linspace(l,u,nPts));
            z = complex(x,y);   
            rho  = abs(bsxfun(@minus,z(:),z(:).'));
        end
        %% crops the central portion of a square matrix, zero pads a square matrix to extend it
        function out = crop(input, ncrop)
            nFrames = size(input,3);
            dim = size(input, 1);
            out = zeros(ncrop, ncrop, nFrames);
            for iFrame = 1:nFrames
                if ncrop < dim
                    deb = round((dim - ncrop) / 2 + 1);
                    fin = round((dim + ncrop) / 2);
                    out(:,:,iFrame) = input(deb:fin, deb:fin, iFrame);
                else
                    deb = round((ncrop-dim) / 2 + 1);
                    fin = round((ncrop+dim) / 2);
                    out(deb:fin, deb:fin, iFrame) = input(:,:,iFrame);
                end
            end
        end
        %% subpixel shift by known amount
        function out = shift(im,x,y)
            %             subx = x-round(x);
            %             suby = y-round(y);
            %             im = circshift(im, round([x, y]));
            %             im = conv2(im, [subx, 1-subx], 'same');
            %             out = conv2(im, [suby, 1-suby], 'same');
            [nPx, nPy] = size(im);
            [X, Y] = meshgrid(1:nPx, 1:nPy);
            out = interp2(X, Y, im, X-x, Y-y, 'cubic', 0);
            
            %             [nPx, nPy] = size(im);
            %             amp = fftshift(abs(fft2(im,nPx*2,nPy*2)));
            %             [X, Y] = meshgrid(linspace(-1, 1, nPx*2), linspace(-1, 1, nPy*2));
            %             eField = amp .* exp(1i * (-X*pi*x-Y*pi*y));
            %             out = fftshift(abs(fft2(eField)));
            %             out = out(nPx/2+1:nPx/2+nPx, nPy/2+1:nPy/2+nPy);
            %             out = out*sum(im(:))/sum(out(:));
        end
        %%
        function iCn = invertMeasurementNoiseCovarianceMatrix(input)
            if iscell(input)
                nGs = size(input,1);
                out = zeros([size(input{1}), nGs]);
                for iGs = 1:nGs
                    out(:,:,iGs) = input{iGs};
                end
            else
                nGs = size(input,3);
                out = input;
            end
            iCn = cell(nGs,1);
            nSlopes = size(out,1);
            for iGs = 1:nGs
                Cn = out(:,:,iGs);
                
                % Extract diagonal and cross-terms -> create sparse matrix
                noiseCovarDiag = diag(Cn);
                noiseCovarDiagP1 = diag(Cn,size(Cn,1)/2);
                B = zeros(size(Cn,1),3);
                B(:,1) = noiseCovarDiag;
                B(1:1:end/2,2) = noiseCovarDiagP1;
                B(end/2+1:1:end,3) = noiseCovarDiagP1;
                CnE = spdiags(B,[0,-nSlopes/2,nSlopes/2],nSlopes,nSlopes);
                iCn{iGs} = pinv(full(CnE));
                %iCn{iGs}(abs(iCn{iGs})< 1e-10) = 0;
                
                % extract only meaningful values on the main and two
                % cross-covar diagonals
                B(:,1) = diag(iCn{iGs});
                B(1:1:end/2,2) = diag(iCn{iGs},size(Cn,1)/2);
                B(end/2+1:1:end,3) = diag(iCn{iGs},size(Cn,1)/2);
                iCn{iGs} = spdiags(B,[0,-nSlopes/2,nSlopes/2],nSlopes,nSlopes);
            end
            iCn = blkdiag(iCn{:});
        end
        function [out, out2,out3] = applySpiders(D, pupil,width,angle)
            x = linspace(-D/2,D/2,size(pupil,1));
            [X Y] = meshgrid(x,x);
            bandIdx = abs(Y)< width/2;
            band = zeros(size(pupil));
            band(bandIdx) = 1;
            band = logical(band + band');
            out = (pupil - double(band)).*pupil;
            if nargout >1 %outputs the petals
                out2 = zeros([size(pupil),4]);
                out3 = zeros([size(pupil),4]);
                petal = zeros(size(pupil));
                petalIdx = (X >0 & Y > 0);
                petal(petalIdx) = 1;
                out2(:,:,1) = petal.*pupil;
                out2(:,:,2) = fliplr(petal).*pupil;
                out2(:,:,3) = flipud(fliplr(petal)).*pupil;
                out2(:,:,4) = flipud(petal).*pupil;
            end
            if exist('angle','var')
                out = imrotate(out,angle,'bicubic','crop');
                out(find(out)) = 1;
                out = out.*pupil;
                out(~pupil) = 0; %force same pixels as in input
                for i = 1:4
                    tmp = imrotate(out2(:,:,i),angle,'bicubic','crop').*out;
                    tmp(find(tmp)) = 1;
                    tmp(~pupil) = 0;
                    piston = mean(tmp(find(out)));
                    tmp = tmp - piston*out;
                    out2(:,:,i) = tmp;
                    out3(:,:,i) = tmp + piston*out;
                end
                pause(0)
                %out = sum(out2,3);
                %out(~pupil) = 0; %force same pixels as in input
            end
        end
        
        
        %% Rotate the DM actuator position by rotAngle in radian
        function [pxx, pyy] = rotateDM(px, py, rotAngle)            
            % function [pxx, pyy] = rotateDM(px,py, rotAngle) 
            % This function rotate the DM actuator positions.
            %
            % px (pxx)  :: The old (new) X actuator positions.
            % py (pyy)  :: The odl (new) Y actuator positions.
            % rotAngle  :: The rotation angle in radian.
            %
            % Created   :: N. Schwartz, Dec 2016
            
            % EXAMPLE:
            % px = real(bifM4.actuatorCoord);
            % py = imag(bifM4.actuatorCoord);
            % [pxx, pyy] = lamTools.rotateDM(px, py, rotAngle);
            % bifM4.actuatorCoord = pxx + 1j*pyy;
            pxx = px*cos(rotAngle) - py*sin(rotAngle);
            pyy = py*cos(rotAngle) + px*sin(rotAngle);
        end
%         %% rotate DM coordinates 
%         function [px, py] = rotateDM(px,py, theta)
%                 res = [cos(theta), sin(theta); -sin(theta), cos(theta)] * [px;py]';
%                 px = res(1,:)';
%                 py = res(2,:)';
%         end
        
        
        %% saveIFCube
        function saveIFCube(dmModes, pupil, directory, fileSuffix)
            % function saveIFCube(dmModes, pupil, directory, fileSuffix)
            % This function saves the influence functions (IFs) in a data
            % cube, and saves the pupil. This function is useful to then
            % create the Karhunen-Loeve (see IDL function).
            %
            % dmModes     :: Data cube of NxMxk where k is the number of images/IFs/Modes.
            % pupil       :: The telescope pupil.
            % directory   :: The directory of where the data is saved.
            % fileSuffix  :: The 2 saved files will have this suffix.
            %
            %
            % Created     :: N. Schwartz, Dec 2016
            
            % EXAMPLE:
            %    dmModes = dm.modes.modes;
            %    pupil = tel.pupil;
            %    directory = '/result/SCAO-H/INPUT/KarhunenLoeve';
            %    fileSuffix = strcat(num2str(tel.resolution), 'pix_', num2str(size(dmModes, 2)),'modes')
            %    lamTools.saveIFCube(dmModes, pupil, directory, fileSuffix);
            
            % Reshape data into 3D data cube
            m = size(dmModes, 2);
            nRes = size(pupil,1);
            tmp = dmModes;
            if size(dmModes,3) ==1; tmp = (reshape(full(dmModes),nRes, nRes, m)); end;
            
            % Save IF cube
            filename = strcat(directory, '/IF_', fileSuffix,'.fits');
            fitswrite(tmp, filename);
            
            % Save pupil
            filename = strcat(directory, '/pupilMask_', fileSuffix,'.fits');
            fitswrite(double(pupil), filename);
        end %End of saveIFCube
        
        
        %% groupIF
        function [newDmModes, newActuatorCoord] = groupIF(dmModes, actuatorCoord, rotAngle)
            % groupIF(dmModes, actuatorCoord, rotAngle)
            % This function groups Edge IFs 2-by-2 (i.e. slaved DM)
            %
            % dmModes       :: Data cube of NMxk where k is the number of images/IFs/modes.
            % actuatorCoord :: DM actuators coordinates
            % rotAngle      :: DM rotation angle in radians.
            %
            % Created       :: N. Schwartz, May 2017
            
            % EXAMPLE:
            % [newDmModes, newActuatorCoord] = lamTools.groupIF(dm.modes.modes, bifM4.actuatorCoord, 0);
            % m4.nb_act = size(newDmModes, 2);
            % bifM4.actuatorCoord = newActuatorCoord;
            % dm = deformableMirror(m4.nb_act,'modes',bifM4,'resolution',tel.resolution,...
            %                       'validActuator',true(1,m4.nb_act),'diameter',tel.D);
            % dm.modes.modes = newDmModes; clear newDmModes
            
            % Default parameters
            if nargin < 3; rotAngle = 0; end;
            
            % Reshape data into 3D data cube
            tmp = dmModes;
            
            % Select and group Edge IFs
            idx = lamTools.whichEdge(actuatorCoord, 0.5, rotAngle);
            newActuatorCoord = actuatorCoord;
            idxEdge = zeros(12,27);
            for ii=1:12; idxEdge(ii,:) = find(idx(:,ii)); end
            idxEdge = [idxEdge(2:12,:); idxEdge(1,:)];
            
            for ii=1:27
                tmp(:,idxEdge(1,ii)) = tmp(:,idxEdge(1,ii)) +  tmp(:,idxEdge(2,ii));
                tmp(:,idxEdge(3,ii)) = tmp(:,idxEdge(3,ii)) +  tmp(:,idxEdge(4,ii));
                tmp(:,idxEdge(5,ii)) = tmp(:,idxEdge(5,ii)) +  tmp(:,idxEdge(6,ii));
                tmp(:,idxEdge(7,ii)) = tmp(:,idxEdge(7,ii)) +  tmp(:,idxEdge(8,ii));
                tmp(:,idxEdge(9,ii)) = tmp(:,idxEdge(9,ii)) +  tmp(:,idxEdge(10,ii));
                tmp(:,idxEdge(11,ii)) = tmp(:,idxEdge(11,ii)) +  tmp(:,idxEdge(12,ii));
            end
            tmp(:,idxEdge(2:2:12,:)) = [];
            newActuatorCoord(idxEdge(2:2:12,:)) = [];
            newDmModes = tmp;
            %             oldDmModes = reshape(full(dmModes),nRes, nRes, m);
            %             figure; imagesc(sum(oldDmModes,3)); colorbar
            %             figure; imagesc(sum(newDmModes, 3)); colorbar
            %             figure; imagesc(sum(oldDmModes,3) - sum(newDmModes, 3)); colorbar
            %             figure; imagesc(newDmModes(:,:,16))
        end %End of groupIF
        
        
        %% groupIF
        function [newDmModes, newActuatorCoord] = groupIFScaledown(dmModes, actuatorCoord, rotAngle)
            % groupIF(dmModes, actuatorCoord, rotAngle)
            % This function groups Edge IFs 2-by-2 (i.e. slaved DM)
            %
            % dmModes       :: Data cube of NMxk where k is the number of images/IFs/modes.
            % actuatorCoord :: DM actuators coordinates
            % rotAngle      :: DM rotation angle in radians.
            %
            % Created       :: N. Schwartz, May 2017
            
            % EXAMPLE:
            % [newDmModes, newActuatorCoord] = lamTools.groupIF(dm.modes.modes, bifM4.actuatorCoord, 0);
            % m4.nb_act = size(newDmModes, 2);
            % bifM4.actuatorCoord = newActuatorCoord;
            % dm = deformableMirror(m4.nb_act,'modes',bifM4,'resolution',tel.resolution,...
            %                       'validActuator',true(1,m4.nb_act),'diameter',tel.D);
            % dm.modes.modes = newDmModes; clear newDmModes
            
            % Default parameters
            if nargin < 3; rotAngle = 0; end;
            
            % Reshape data into 3D data cube
            tmp = dmModes;
            
            % Select and group Edge IFs
            idx = lamTools.whichEdge(actuatorCoord, 0.4, rotAngle);
            n=sum(idx);
            newActuatorCoord = actuatorCoord;
            for ii=1:12; idxEdge(ii,:) = find(idx(:,ii)); end
            idxEdge = [idxEdge(2:12,:); idxEdge(1,:)];
  
            for ii=1:n(1)
                tmp(:,idxEdge(1,ii)) = tmp(:,idxEdge(1,ii)) +  tmp(:,idxEdge(2,ii));
                tmp(:,idxEdge(3,ii)) = tmp(:,idxEdge(3,ii)) +  tmp(:,idxEdge(4,ii));
                tmp(:,idxEdge(5,ii)) = tmp(:,idxEdge(5,ii)) +  tmp(:,idxEdge(6,ii));
                tmp(:,idxEdge(7,ii)) = tmp(:,idxEdge(7,ii)) +  tmp(:,idxEdge(8,ii));
                tmp(:,idxEdge(9,ii)) = tmp(:,idxEdge(9,ii)) +  tmp(:,idxEdge(10,ii));
                tmp(:,idxEdge(11,ii)) = tmp(:,idxEdge(11,ii)) +  tmp(:,idxEdge(12,ii));
            end
            tmp(:,idxEdge(2:2:12,:)) = [];
            newActuatorCoord(idxEdge(2:2:12,:)) = [];
            newDmModes = tmp;
            %             oldDmModes = reshape(full(dmModes),nRes, nRes, m);
            %             figure; imagesc(sum(oldDmModes,3)); colorbar
            %             figure; imagesc(sum(newDmModes, 3)); colorbar
            %             figure; imagesc(sum(oldDmModes,3) - sum(newDmModes, 3)); colorbar
            %             figure; imagesc(newDmModes(:,:,16))
        end %End of groupIF
        
        
        %% groupIF2
        function [newDmModes, newActuatorCoord] = groupIF2(dmModes, actuatorCoord, rotAngle)
            % function groupIF2(dmModes, actuatorCoord, rotAngle)
            % This function groups Edge IFs 4-by-4 (i.e. slaved DM using 2 actuators on either side of the spider)
            % 
            % dmModes       :: Data cube of NMxk where k is the number of images/frames.
            % actuatorCoord :: DM actuators coordinates
            % rotAngle      :: DM rotation angle in radians.
            %
            % Created       :: N. Schwartz, June 2017
            
            % EXAMPLE:
            %    dmModes = dm.modes.modes;
            %    actuatorCoord = bifM4.actuatorCoord;
            %    rotAngle = pi/12;
            %    lamTools.groupIF(dmModes, actuatorCoord, rotAngle);
            
            % Default parameters
            if nargin < 3; rotAngle = 0; end;
            
            % Reshape data into 3D data cube
            tmp = dmModes;
            
            % Select and group Edge IFs
            [row1, row2, ~, ~] = lamTools.whichEdgeExtended(actuatorCoord, rotAngle);
            %             row1 = lamTools.whichEdge(actuatorCoord, 0.5, rotAngle);
            idxEdge = zeros(12,27);
            idxEdge2 = zeros(12,28);
            for ii=1:12;
                idxEdge(ii,:) = find(row1(:,ii));
                idxEdge2(ii,:) = find(row2(:,ii));
            end
            idxEdge = [idxEdge(2:12,:); idxEdge(1,:)];
            idxEdge2 = [idxEdge2(2:12,1:27); idxEdge2(1,1:27)];
            
            for ii=1:27
                tmp(:,idxEdge(1,ii)) = tmp(:,idxEdge(1,ii)) +  tmp(:,idxEdge(2,ii)) + tmp(:,idxEdge2(1,ii)) +  tmp(:,idxEdge2(2,ii));
                tmp(:,idxEdge(3,ii)) = tmp(:,idxEdge(3,ii)) +  tmp(:,idxEdge(4,ii)) + tmp(:,idxEdge2(3,ii)) +  tmp(:,idxEdge2(4,ii));
                tmp(:,idxEdge(5,ii)) = tmp(:,idxEdge(5,ii)) +  tmp(:,idxEdge(6,ii)) + tmp(:,idxEdge2(5,ii)) +  tmp(:,idxEdge2(6,ii));
                tmp(:,idxEdge(7,ii)) = tmp(:,idxEdge(7,ii)) +  tmp(:,idxEdge(8,ii)) + tmp(:,idxEdge2(7,ii)) +  tmp(:,idxEdge2(8,ii));
                tmp(:,idxEdge(9,ii)) = tmp(:,idxEdge(9,ii)) +  tmp(:,idxEdge(10,ii))+ tmp(:,idxEdge2(9,ii)) +  tmp(:,idxEdge2(10,ii));
                tmp(:,idxEdge(11,ii))= tmp(:,idxEdge(11,ii))+  tmp(:,idxEdge(12,ii))+ tmp(:,idxEdge2(11,ii))+  tmp(:,idxEdge2(12,ii));
            end
            
            p1 = idxEdge(2:2:12,:);
            p2 = idxEdge2(1:12,:);
            tmp(:,[p1(:); p2(:)]) = [];
            newActuatorCoord = actuatorCoord;
            newActuatorCoord([p1(:); p2(:)]) = [];
            newDmModes = tmp;
            
            %A = reshape(idxEdge(2:2:12,:),6*27,1);
            %B = reshape(idxEdge2(1:12,:), 12*27,1);
            %tmp(:,:, reshape(cat(1,A,B),12*27+6*27,1)) = [];
            %newDmModes = tmp;
            
            % oldDmModes = reshape(full(dmModes),nRes, nRes, m);
            % figure; imagesc(sum(oldDmModes,3)); colorbar
            % figure; imagesc(sum(newDmModes, 3)); colorbar
            % figure; imagesc(sum(oldDmModes,3) - sum(newDmModes, 3)); colorbar
            % figure; imagesc(newDmModes(:,:,16))
        end %End of groupIF2
        
        
        function [PX,PY]=M4downScaled(D) % coded by Y.O
            %% M4SCALEDOWN computes M4-like DM with given diameter
            
            % First Petal -----------------------------------
            % First line
            n=ceil(D/2/0.5067)+3;
            p1x=ones(n,1)*0.161;
            p1y=(0:n-1).'*0.5067-0.08;
            px = p1x;
            py = p1y;
            % Second Line
            p2x = p1x + 0.51;
            p2y = p1y - 0.25;
            px = [px; p2x];
            py = [py; p2y];
            % Other lines
            t=1;
            pnx = p2x;
            pny = p2y;
            for k=2:n-1
                if k==2
                    pnx = pnx + 0.45;
                else
                    pnx = pnx + 0.434;
                end
                pny = pny + t*0.25;
                px = [px; pnx];
                py = [py; pny];
                t=t*-1;
            end
            % Last line
            Yedge = @(x) tan(30*pi/180)*x+0.161/cos(30*pi/180);
            newpx = px(Yedge(px)+0.3<py);
            newpy = py(Yedge(px)+0.3<py);
            px = newpx;
            py = newpy;
            p1x=ones(n,1)*0.161;
            p1y=(0:n-1).'*0.5035-0.08;
            rp1x = (p1x-0.34)*cos(59.9*pi/180) + p1y*sin(59.9*pi/180);
            rp1y = -(p1x-0.34)*sin(59.9*pi/180) + p1y*cos(59.9*pi/180);
            PX1 = [px; rp1x];
            PY1 = [py; rp1y];
            
            % Second Petals -----------------------------------
            % First line
            n=ceil(D/2/0.5031)+3;
            p1x=ones(n,1)*0.163;
            p1y=(0:n-1).'*0.5031-0.08;
            px = p1x;
            py = p1y;
            % Second Line
            p2x = p1x + 0.51;
            p2y = p1y - 0.25;
            px = [px; p2x];
            py = [py; p2y];
            % Other lines
            t=1;
            pnx = p2x;
            pny = p2y;
            for k=2:n-1
                if k==2
                    pnx = pnx + 0.452;
                else
                    pnx = pnx + 0.438;
                end
                
                pny = pny + t*0.25 - 0.0025;
                
                px = [px; pnx];
                py = [py; pny];
                t=t*-1;
            end
            % Last line
            Yedge = @(x) tan(30*pi/180)*x+0.161/cos(30*pi/180);
            newpx = px(Yedge(px)+0.3<py);
            newpy = py(Yedge(px)+0.3<py);
            px = newpx;
            py = newpy;
            p1x=ones(n,1)*0.163;
            p1y=(0:n-1).'*0.5035-0.08;
            rp1x = (p1x-0.32)*cos(60.4*pi/180) + p1y*sin(60.4*pi/180);
            rp1y = -(p1x-0.32)*sin(60.4*pi/180) + p1y*cos(60.4*pi/180);
            px = [px; rp1x];
            py = [py; rp1y];
            PX2=px*cos(59.8*pi/180)+py*sin(59.8*pi/180);
            PY2=-px*sin(59.8*pi/180)+py*cos(59.8*pi/180);
            
            
            % Third Petals -----------------------------------
            % First line
            n=ceil(D/2/0.5031)+3;
            p1x=ones(n,1)*0.1667;
            p1y=(0:n-1).'*0.5031-0.08;
            px = p1x;
            py = p1y;
            % Second Line
            p2x = p1x + 0.51;
            p2y = p1y - 0.25;
            px = [px; p2x];
            py = [py; p2y];
            % Other lines
            t=1;
            pnx = p2x;
            pny = p2y;
            for k=2:n-1
                if k==2
                    pnx = pnx + 0.452;
                else
                    pnx = pnx + 0.438;
                end
                
                pny = pny + t*0.25+0.0035;
                
                px = [px; pnx];
                py = [py; pny];
                t=t*-1;
            end
            % Last line
            Yedge = @(x) tan(30*pi/180)*x+0.161/cos(30*pi/180);
            newpx = px(Yedge(px)+0.3<py);
            newpy = py(Yedge(px)+0.3<py);
            px = newpx;
            py = newpy;
            p1x=ones(n,1)*0.1667;
            p1y=(0:n-1).'*0.5066-0.08;
            rp1x = (p1x-0.32)*cos(59.8*pi/180) + p1y*sin(59.8*pi/180);
            rp1y = -(p1x-0.32)*sin(59.8*pi/180) + p1y*cos(59.8*pi/180);
            px = [px; rp1x];
            py = [py; rp1y];
            PX3=px*cos(120.22*pi/180)+py*sin(120.22*pi/180);
            PY3=-px*sin(120.22*pi/180)+py*cos(120.22*pi/180);
            
            
            %--------------------------------------------
            PX = [PX1; PX1*cos(pi)+PY1*sin(pi);
                PX2; PX2*cos(pi)+PY2*sin(pi);
                PX3; PX3*cos(pi)+PY3*sin(pi);];
            PY = [PY1; -PX1*sin(pi)+PY1*cos(pi);
                PY2; -PX2*sin(pi)+PY2*cos(pi);
                PY3; -PX3*sin(pi)+PY3*cos(pi);];
            
            r = sqrt(PX.^2+PY.^2);
            PX(r>D/2+1) = [];
            PY(r>D/2+1) = [];
            
            % Check
            %plot(PX,PY,'o');
            %axis equal tight;
            
        end
        
        %% create a scaled-down version of the E-ELT M4 DM (Old codes)
        %{
        function [px,py] = M4downScaled(tel,option,pitch,inPupil)
            %% M4DOWNSCALED
            % [px,py] = M4downScaled(tel,option,pitch,inPupil) computes the
            % actuator postion for different DM  models with 'option' as
            %    LBT
            %    square or Fried
            %    random
            %    hexagon
            %    triangle (scale downed M4)
            %  
            
            D = tel.D + 2*pitch;
            obscurationRatio = tel.obstructionRatio;
            if ~exist('inPupil','var') || isempty(inPupil)
                inPupil = 1;
            end
            [px,py] = lamTools.actuatorCoordinates(tel,option,pitch,inPupil);
            % add double rim of actuators at 0, pi/3, 2pi/3 angles
            %cx = -D/2:pitch:D/2;
            cx = linspace(-1,1,ceil(D/pitch)+4)*((ceil(D/pitch)+4)*pitch)/2;
            %cx = pitch/2*ones(size(cy));
            cy = 0.163*ones(size(cx));
            
            out1 = [];
            out2 = [];
            for iTheta = 1:3
                theta = (iTheta-1)*pi/3;
                out1 = [out1 [cos(theta), sin(theta); -sin(theta), cos(theta)] * [cx;cy]];
                out2 = [out2 [cos(theta), sin(theta); -sin(theta), cos(theta)] * [cx;-cy]];
            end
            critDistance = pitch*0.5;
            for kAct = 1:length(out1)
                dist = sqrt((out1(1,kAct) - px).^2 + (out1(2,kAct) - py).^2);
                removeIdx = (dist < critDistance);
                px(removeIdx) = [];
                py(removeIdx) = [];
                dist = sqrt((out2(1,kAct) - px).^2 + (out2(2,kAct) - py).^2);
                removeIdx = (dist < critDistance);
                px(removeIdx) = [];
                py(removeIdx) = [];
            end
            %hold on
            px = [px,out1(1,:),out2(1,:)];
            py = [py,out1(2,:),out2(2,:)];
            %scatter(px ,py, 'r.')
            % remove those outside the pupil
            dist = sqrt(px.^2 + py.^2);
            px(dist>D) = [];py(dist>D) = [];
            %px(dist<D/2*obscurationRatio-pitch/2) = [];
            %py(dist<D/2*obscurationRatio-pitch/2) = [];
            % scatter plot
            %hold off;
            %figure
            %scatter(px,py,'k.')
            %hold off;
            %box on
            %title('DM actuator locations')
            %xlabel('meters')
            %ylabel('meters')
            %axis equal tight;
        end
        %}
        %% create a DM with pre-defined actuator layout
        function [px,py] = actuatorCoordinates(tel,option,pitch,inPupil)
            %% actuatorCoordinates for different DM models
            % [px py] = actuatorCoordinates(tel,option,pitch,inPupil)
            % Computes the actuator coordinates for different DM models for a
            % given telescope as
            %    LBT
            %    square or Fried
            %    random
            %    hexagon
            %    triangle (scale downed M4)
            % a pitch in meters and either across the whole squared
            % computational domain set by tel.D or within the pupil
            % tel.pupil
            
            
            D = tel.D + 2*pitch;
            obscurationRatio = tel.obstructionRatio;
            if ~exist('inPupil','var') || isempty(inPupil)
                inPupil = 1;
            end
            switch option
                % LBT
                case 'LBT'
                    N = 12;     % Initial number of actuators in first ring
                    n = 3;      % Additional actuators per rings
                    nRings = 9; % Number of rings
                    xcentre = 0;
                    ycentre = 0;
                    px = 0;
                    py = 0;
                    
                    for ind=1:nRings
                        theta = 0:pi/(N/2 + ((ind-1)*n)):2*pi;
                        px = [px, pitch*ind*cos(theta) + xcentre];
                        py = [py, pitch*ind*sin(theta) + ycentre];
                    end
                    
                    if inPupil==1
                        % Select in Pupil
                        dist = sqrt(px.^2 + py.^2);
                        valid = dist <= D/2+pitch/2 & dist >= obscurationRatio * D/2  -pitch/2;
                        px = px(valid);
                        py = py(valid);
                    else
                        vx = find(abs(px) <= D/2);
                        vy = find(abs(py(vx)) <= D/2);
                        px = px(vx(vy));
                        py = py(vx(vy));
                    end
                    % Square
                case {'square', 'fried'}
                    [px, py] = meshgrid(-D/2:pitch:D/2);
                    px = reshape(px, 1, numel(px));
                    py= reshape(py, 1, numel(py));
                    
                    if inPupil==1
                        % Select in Pupil
                        dist = sqrt(px.^2 + py.^2);
                        valid = dist <= D/2+pitch/2 & dist >= obscurationRatio * D/2  -pitch/2;
                        px = px(valid);
                        py = py(valid);
                        %scatter(px, py)
                    else
                        vx = find(abs(px) <= D/2);
                        vy = find(abs(py(vx)) <= D/2);
                        px = px(vx(vy));
                        py = py(vx(vy));
                    end
                    
                    % Random
                case 'random'
                    oc = tel.obstructionRatio * tel.D/2;
                    pr = oc + (tel.D/2-oc).*rand(250,1);
                    ptheta = 2*pi*rand(250,1);
                    px = pr.*cos(ptheta);
                    py = pr.*sin(ptheta);
                    
                    % % Hexagon
                case {'hexagon', 'triangle'}
                    % Hexagon
                    if strcmp(option, 'hexagon')
                        x_hexagon=[-1 -0.5 0.5 1];
                        y_hexagon=[0 -sqrt(3)/2 -sqrt(3)/2 0];
                    end
                    
                    % Hexagon centre
                    if strcmp(option, 'triangle')
                        x_hexagon=[0, 1.5, -1, -0.5,      0.5,        1];
                        %y_hexagon=[0, -0.5,   0, -0.5, -0.5, 0];
                        y_hexagon=[0, -.9,   0, -.9, -.9, 0]; % for M4
                    end
                    
                    N = tel.D/pitch/2+2;
                    M = tel.D/(pitch*0.9)/2+2;
                    px = [];
                    py = [];
                    
                    for nn=0:N-1
                        for mm=0:M-1
                            px = [px, x_hexagon+1.5+3*nn;];
                            py = [py, y_hexagon+.9+.9*2*mm];
                        end
                    end
                    
                    % If number of actuators fitting into square is odd make sure y
                    % coordinate has centre at 0.  I.e. if pitch is exactly and EVEN
                    % division of D.
                    px = px - max(px(:))/2;
                    if (tel.D/pitch)/2==round((tel.D/pitch)/2)
                        py = py - max(py(:))/2 - 1/2;
                    else
                        py = py - max(py(:))/2;
                    end
                    px = px * pitch;
                    py = py * pitch;
                    
                    if inPupil==1
                        % Select in Pupil
                        dist = sqrt(px.^2 + py.^2);
                        valid = dist <= D/2+pitch/2 & dist >= obscurationRatio * D/2  -pitch/2;
                        px = px(valid);
                        py = py(valid);
                    else
                        vx = find(abs(px) <= D/2);
                        vy = find(abs(py(vx)) <= D/2);
                        px = px(vx(vy));
                        py = py(vx(vy));
                    end
            end
            % scatter plot
            %scatter(px,py)
            %box on
            %title('DM actuator locations')
            %xlabel('meters')
            %ylabel('meters')
            %axis equal tight;
            
        end
        
        %% m1 pupil
        function [m1Pupil, m1PupilAllGlass,m1PupilAllGlassNoSpider] = m1Pupil(nOutPoints,m1DiamOut,view,rotAngle)
            %{
            function [m1Pupil, m1PupilAllGlass] = m1Pupil(nOutPoints,m1DiamOut,view)
            INPUT:
                nOutPoints      ::  number of pixels for the output pupil
                m1DiamOut       ::  [mm] diameter of the output grid (not mirror) over the nOutPoints
                view            ::  flag to visualise the pupil on screen
                rotAngle        ::  Pupil rotation in radian
            OUTPUT
                m1Pupil         ::  the resampled M1 pupil at the nOutPoints
                                    resolution
                m1PupilAllGlass ::  the M1 pupil all-glass with spiders
                m1PupilAllGlassNoSpider ::  the M1 pupil all-glass no
                                            spiders
            %}
            if ~exist('rotAngle','var')
                rotAngle = 0; % rotation angle in radian (YO:Added 22/09/2017)
            end
            pwdir = pwd;
            cd /home/matuser/mfiles/_FASTF/TELESCOPE_DATA/
            fileName = 'M1Pupil_10ArcminInnerDiam_10ArcsecOuterDiam.fits';
            M1 = fitsread(fileName);
            nPoints = size(M1,1);
            m1Diam = 38542; % the actual M1 diameter for 10'' unvignetted pupil
            if ~exist('m1DiamOut','var') || isempty(m1DiamOut)
                m1DiamOut = m1Diam; % m4 output diameter is 41609mm projected onto M1 = 2386*9417/540
            end
            %pixelSize = m1Diam/nPoints; % in mm
            xind = linspace(-m1Diam/2, m1Diam/2, nPoints);
            [XXin, YYin] = meshgrid(xind, xind);
            xindOut = linspace(-m1DiamOut/2, m1DiamOut/2, nOutPoints);
            [XXout, YYout] = meshgrid(xindOut, xindOut);
            if abs(rotAngle) > 0
                rXXout = XXout*cos(rotAngle)-YYout*sin(rotAngle);
                rYYout = XXout*sin(rotAngle)+YYout*cos(rotAngle);
                XXout = rXXout;
                YYout = rYYout;
            end
            m1Pupil = interp2(XXin, YYin, M1, XXout, YYout,'cubic');
            m1Pupil(isnan(m1Pupil)) = 0;
            m1Pupil = double(logical(m1Pupil > 0.1)); % NS: Added 19/12/2016 to have pixels 0/1
            if ~exist('view','var') || isempty(m1DiamOut)
                view = 0; % m4 output diameter is 41609mm projected onto M1 = 2386*9417/540
            end
            if view
                imagesc(xind/1000, xind/1000,m1Pupil);
                xlabel('meters'), ylabel('meters')
                title(['M1 pupil'],'fontweight','bold')
                text(-10,-10,fileName)
            end
            if nargout > 1
                M1 = fitsread('M1Pupil_10ArcminInnerDiam_AllGlassOuterDiam.fits');
                nPoints = size(M1,1);
                m1Diam = 36903; %mm
                %m1DiamOut = m1Diam;
                xind = linspace(-m1Diam/2, m1Diam/2, nPoints);
                [XXin, YYin] = meshgrid(xind, xind);
                xindOut = linspace(-m1DiamOut/2, m1DiamOut/2, nOutPoints);
                [XXout, YYout] = meshgrid(xindOut, xindOut);
                m1PupilAllGlass = interp2(XXin, YYin, M1, YYout, XXout,'cubic');
                m1PupilAllGlass(isnan(m1PupilAllGlass)) = 0;
                m1PupilAllGlass = double(logical(m1PupilAllGlass > 0.1)); % NS: Added 19/12/2016 to have pixels 0/1
            end
            if nargout > 2
                [~, RHO] = cart2pol(XXout,YYout);
                m1PupilAllGlassNoSpider = (RHO <= m1Diam/2+50 & RHO > m1Diam/2*0.3);
            end
            cd(pwdir) % move back to where we were
        end
        %%
        function out = m1PupilScaledown(D, pupil, width, rotation)
            
            if nargin < 4; rotation = 0; end; % Rotation in degrees
            if rotation ~=0;
                disp('!WARNING: lamTools.m1PupilScaledown')
                disp('Rotation parameter may introduce artifacts!')
            end
            
            % Create 1 horizontal spider
            x = linspace(-D/2,D/2,size(pupil,1));
            X = meshgrid(x,x);
            bandIdx = abs(X)<= width/2;
            band = zeros(size(pupil));
            band(bandIdx) = 1;
            band = logical(band);
            
%             spider1 = (pupil - double(band)).*pupil;
%             spider1 = imrotate(spider1,60,'bicubic','crop');
%             spider1(find(spider1)) = 1;
%             spider1 = spider1.*pupil;
%             spider1(~pupil) = 0;
%             spider2 = (pupil - double(band)).*pupil;
%             spider2 = imrotate(spider2,-60,'bicubic','crop');
%             spider2(find(spider2)) = 1;
%             spider2 = spider2.*pupil;
%             spider2(~pupil) = 0;
%             spider3 = (pupil - double(band)).*pupil;     
%             out = spider1.*spider2.*spider3;
            
            % Change (16 Oct 2017 NSc, to make sure spiders' width is the
            % same in all 3 directions
            pupilOut = band; 
            tmp  = imrotate(pupilOut,rotation, 'bicubic', 'crop');
            tmp1 = imrotate(tmp, 60, 'bicubic', 'crop');
            tmp2 = imrotate(tmp,-60, 'bicubic', 'crop');
            pupilOut = tmp + tmp1 + tmp2;
            
            pupilOut = (pupil - pupilOut).*pupil;
            threshold = 0.5;
            pupilOut(pupilOut<threshold)=0; 
            pupilOut(pupilOut>=threshold)=1;
            out = pupilOut;
%             % Verification: Check surface is the same for all spiders
%             figure;
%             subplot(2,2,1); imagesc(pupilOut); axis square; title(strcat('Full: ', num2str(sum(pupilOut(:)))));
%             temp1 = (pupil - tmp).*pupil; temp1(temp1<threshold)=0; temp1(temp1>=threshold)=1;
%             subplot(2,2,2); imagesc(temp1); axis square; title(strcat('Tot pix: ', num2str(sum(temp1(:)))));
%             temp2 = (pupil - tmp1).*pupil; temp2(temp2<threshold)=0; temp2(temp2>=threshold)=1;
%             subplot(2,2,3); imagesc(temp2); axis square; title(strcat('Tot pix: ', num2str(sum(temp2(:)))));
%             temp3 = (pupil - tmp2).*pupil; temp3(temp3<threshold)=0; temp3(temp3>=threshold)=1;
%             subplot(2,2,4); imagesc(temp3); axis square; title(strcat('Tot pix: ', num2str(sum(temp3(:)))));
        end
        %% m4Petals
        function petals = m4Petals(pupil,rotAngle,isVLT)
            % function petals = m4Petals(tel)
            %
            % Create the 6 petals/segments of M4
            % pupil   :: The pupil of the telescope (e.g. tel.pupil).
            % rotAngle:: Define the rotation angle [radians] of the DM;
            % petals  :: A data cube [size: tel.resolution, tel.resolution, 6] containing the 6 petals/segments.
            % Petals are arranged in trigonometrical order.
            if ~exist('rotAngle','var') || isempty(rotAngle)
                rotAngle = 0;
            end
            if ~exist('isVLT','var') || isempty(isVLT)
                isVLT = 0;
            end
            sz = size(pupil);
            % Calculate the polar coordinates
            [tmpX, tmpY] = meshgrid(-sz(1)/2:sz(1)/2-1,-sz(2)/2:sz(2)/2-1);
            [tmpTheta, ~] = cart2pol(tmpX,tmpY);
            
            % Create 6 petals
            if ~isVLT
                petals = zeros(sz(1), sz(2), 6);
                
                % Select the 6 petals matching the pupil shape
                tmpPetal = zeros(sz(1), sz(2));
                tmpPetal(tmpTheta>=-pi/6-rotAngle & tmpTheta <pi/6-rotAngle) = 1; petals(:,:,1) = tmpPetal.*pupil;
                tmpPetal = zeros(sz(1), sz(2));
                tmpPetal(tmpTheta>=pi/6-rotAngle & tmpTheta <3*pi/6-rotAngle) = 1; petals(:,:,6) = tmpPetal.*pupil;
                tmpPetal = zeros(sz(1), sz(2));
                tmpPetal(tmpTheta>=3*pi/6-rotAngle & tmpTheta <5*pi/6-rotAngle) = 1; petals(:,:,5) = tmpPetal.*pupil;
                tmpPetal = zeros(sz(1), sz(2));
                tmpPetal(tmpTheta>=5*pi/6-rotAngle | tmpTheta<-5*pi/6-rotAngle) = 1; petals(:,:,4) = tmpPetal.*pupil;
                tmpPetal = zeros(sz(1), sz(2));
                tmpPetal(tmpTheta>=-5*pi/6-rotAngle & tmpTheta <-3*pi/6-rotAngle) = 1; petals(:,:,3) = tmpPetal.*pupil;
                tmpPetal = zeros(sz(1), sz(2));
                tmpPetal(tmpTheta>=-3*pi/6-rotAngle & tmpTheta <-pi/6-rotAngle) = 1; petals(:,:,2) = tmpPetal.*pupil;
            % Create 4 petals
            else         
                petals = zeros(sz(1), sz(2), 4);
                
                % Select the 6 petals matching the pupil shape
                tmpPetal = zeros(sz(1), sz(2));
                tmpPetal(tmpTheta>pi/2-rotAngle & tmpTheta <=pi-rotAngle) = 1; petals(:,:,2) = tmpPetal.*pupil;
                petals(:,:,3) = rot90(petals(:,:,2),-1).*pupil;
                petals(:,:,1) = rot90(petals(:,:,2)).*pupil;
                petals(:,:,4) = rot90(petals(:,:,1)).*pupil;%flipud(petals(:,:,1)).*pupil;
            end
        
        end
        
        %% whichPetal
        function idx = whichPetal(actuatorCoord, rotAngle, nPetals)
            % function idx = whichPetal(actuatorCoord)
            %
            % actuatorCoord :: The actuators coordinates (complex). e.g.: bifM4.actuatorCoord
            % rotAngle      :: Define the rotation angle [radians] of the DM;
            % idx           :: The petal number of where the actuator is located.
            % Petals are arranged in trigonometrical order.
            if ~exist('rotAngle','var') || isempty(rotAngle)
                rotAngle = pi/4;
            end
            if ~exist('nPetals','var') || isempty(nPetals)
                nPetals = 6;
            end
            % petals = reshape(petals, tel.resolution*tel.resolution, 6);
            nAct = size(actuatorCoord);
            idx = zeros(nAct);
           
            if nPetals == 6
                angleDM = angle(actuatorCoord);
                for kActuator=1:nAct
                    %     whichPetal = [sum(full(dm.modes.modes(:,kActuator)).*petals(:,1)),...
                    %         sum(full(dm.modes.modes(:,kActuator)).*petals(:,2)),...
                    %         sum(full(dm.modes.modes(:,kActuator)).*petals(:,3)),...
                    %         sum(full(dm.modes.modes(:,kActuator)).*petals(:,4)),...
                    %         sum(full(dm.modes.modes(:,kActuator)).*petals(:,5)),...
                    %         sum(full(dm.modes.modes(:,kActuator)).*petals(:,6))];
                    %     idx = max(find(whichPetal == (max(whichPetal))));
                    %     idxTab(kActuator,1) = idx;
                    if angleDM(kActuator,1)>-pi/6+rotAngle && angleDM(kActuator,1) <pi/6+rotAngle; idxTmp =1; end;
                    if angleDM(kActuator,1)>pi/6+rotAngle && angleDM(kActuator,1) <3*pi/6+rotAngle; idxTmp =2; end;
                    if angleDM(kActuator,1)>3*pi/6+rotAngle && angleDM(kActuator,1) <5*pi/6+rotAngle; idxTmp =3; end;
                    if angleDM(kActuator,1)>5*pi/6+rotAngle || angleDM(kActuator,1) <-5*pi/6+rotAngle; idxTmp =4; end;
                    if angleDM(kActuator,1)>-5*pi/6+rotAngle && angleDM(kActuator,1) <-3*pi/6+rotAngle; idxTmp =5; end;
                    if angleDM(kActuator,1)>-3*pi/6+rotAngle && angleDM(kActuator,1) <-pi/6+rotAngle; idxTmp =6; end;
                    idx(kActuator,1) = idxTmp;
                end
            elseif nPetals == 4
%                  angleDM = unwrap(angle(actuatorCoord)+2*pi);
                 angleDM = unwrap(angle(actuatorCoord) +pi*2);
                 angleDM = angleDM - rotAngle;
                 angleDM(angleDM <0) = angleDM(angleDM <0)+2*pi;
                 rotAngle = 0;
                for kActuator=1:nAct
                    if angleDM(kActuator,1)>=0+rotAngle && angleDM(kActuator,1) <=pi/2+rotAngle; idxTmp =4; end;
                    if angleDM(kActuator,1)>=pi/2+rotAngle && angleDM(kActuator,1) <=pi+rotAngle; idxTmp =3; end;
                    if angleDM(kActuator,1)>=pi+rotAngle && angleDM(kActuator,1) <=3/2*pi+rotAngle; idxTmp =2; end;
                    if angleDM(kActuator,1)>=3/2*pi+rotAngle && angleDM(kActuator,1) <=2*pi+rotAngle; idxTmp =1; end;
                    idx(kActuator,1) = idxTmp;
                end
            end
%             plot(1:sum(idx == 1), angleDM(idx == 1),'.k'); hold on;
%             plot((1:sum(idx == 2))+sum(idx == 1), angleDM(idx == 2),'.r');
%             plot((1:sum(idx == 3))+sum(idx <= 2), angleDM(idx == 3),'.y'); 
%             plot((1:sum(idx == 4))+sum(idx <=3), angleDM(idx == 4),'.b'); 
%             plot((1:sum(idx == 5))+sum(idx <= 4), angleDM(idx == 5),'.m'); 
%             plot((1:sum(idx == 6))+sum(idx <= 5), angleDM(idx == 6),'.c'); hold off;
%             legend(strcat('1- ', num2str(sum(idx == 1))),...
%                 strcat('2- ', num2str(sum(idx == 2))),...
%                 strcat('3- ', num2str(sum(idx == 3))),...
%                 strcat('4- ', num2str(sum(idx == 4))),...
%                 strcat('5- ', num2str(sum(idx == 5))),...
%                 strcat('6- ', num2str(sum(idx == 6)))); 
        end
        
        %% whichEdge
        function idx = whichEdge(actuatorCoord, pitch, rotAngle, ~)
            % idx = whichEdge(actuatorCoord, view)
            %
            % actuatorCoord :: The actuators coordinates (complex). e.g.: bifM4.actuatorCoord
            % pitch         :: The actuator pitch in meters (optional, default: 0.5).
            % rotAngle      :: Define the rotation angle [radians] of the DM;
            % idx           :: A matrix [nAct, 12] containing the location of the edges.
            % Edges are arranged in trigonometrical order, starting at 4 o'clock.
            
            % Examples:
            %    idx = lamTools.whichEdge(bifM4.actuatorCoord)
            %    idx = lamTools.whichEdge(bifM4.actuatorCoord, 0.5)                
            %    idx = lamTools.whichEdge(bifM4.actuatorCoord, 0.5, pi/12)        
            %    idx = lamTools.whichEdge(bifM4.actuatorCoord, 0.5, pi/12, 'view') 
% % %             bifM4.actuatorCoord ????
            % Parse input
            if nargin == 1
                pitch = 0.5;
                rotAngle = 0;
                view = 0;
            elseif nargin == 2
                rotAngle = 0;
                view = 0;
            elseif nargin == 3
                rotAngle = -rotAngle; % To match the other convention
                view = 0;
            else 
                rotAngle = -rotAngle; % To match the other convention
                view = 1;
            end
            p = pitch;
            
            % Algorithm:
            %   1. Select actuator in petal N (N from 1 to 6),
            %   2. Shift actuators by one pitch (+/-X or +/-Y depending on the location of the petal),
            %   3. Select the edge actuators on the shifted coordinates that are over the Xpi/6 line.
            idx = [...
                angle(actuatorCoord) > -pi/6-rotAngle & angle(actuatorCoord) < pi/6-rotAngle & ...
                angle(real(actuatorCoord)+1i*imag(actuatorCoord)-1i*p) < -pi/6-rotAngle,...% 3 o'clock petal bottom row
                angle(actuatorCoord) > -pi/6-rotAngle & angle(actuatorCoord) < pi/6-rotAngle & ...
                angle(real(actuatorCoord)+1i*imag(actuatorCoord)+1i*p) > pi/6-rotAngle,...% 3 o'clock petal bottom row
                angle(actuatorCoord) > pi/6-rotAngle & angle(actuatorCoord) < 3*pi/6-rotAngle & ...
                angle(real(actuatorCoord)+1i*imag(actuatorCoord)-1i*p) < pi/6-rotAngle,...% 3 o'clock petal bottom row
                angle(actuatorCoord) > pi/6-rotAngle & angle(actuatorCoord) < 3*pi/6-rotAngle & ...
                angle(real(actuatorCoord)-p+1i*imag(actuatorCoord)) > 3*pi/6-rotAngle,...% 3 o'clock petal bottom row
                angle(actuatorCoord) > 3*pi/6-rotAngle & angle(actuatorCoord) < 5*pi/6-rotAngle & ...
                angle(real(actuatorCoord)+p+1i*imag(actuatorCoord)) < 3*pi/6-rotAngle,...% 3 o'clock petal bottom row
                angle(actuatorCoord) > 3*pi/6-rotAngle & angle(actuatorCoord) < 5*pi/6-rotAngle & ...
                angle(real(actuatorCoord)+1i*imag(actuatorCoord)-1i*p) > 5*pi/6-rotAngle,...% 3 o'clock petal bottom row
                (angle(actuatorCoord) > 5*pi/6-rotAngle | angle(actuatorCoord) < -5*pi/6-rotAngle) & ...
                (angle(real(actuatorCoord)+1i*imag(actuatorCoord)+1i*p) < 5*pi/6-rotAngle & angle(actuatorCoord) > 0) ,...% 3 o'clock petal bottom row
                (angle(actuatorCoord) > 5*pi/6-rotAngle | angle(actuatorCoord) < -5*pi/6-rotAngle) & ...
                (angle(real(actuatorCoord)+1i*imag(actuatorCoord)-1i*p) > -5*pi/6-rotAngle & angle(actuatorCoord) < 0),...% 3 o'clock petal bottom row
                angle(actuatorCoord) > -5*pi/6-rotAngle & angle(actuatorCoord) < -3*pi/6-rotAngle & ...
                angle(real(actuatorCoord)+1i*imag(actuatorCoord)+1i*p) < -5*pi/6-rotAngle,...% 3 o'clock petal bottom row
                angle(actuatorCoord) > -5*pi/6-rotAngle & angle(actuatorCoord) < -3*pi/6-rotAngle & ...
                angle(real(actuatorCoord)+p+1i*imag(actuatorCoord)) > -3*pi/6-rotAngle,...% 3 o'clock petal bottom row
                angle(actuatorCoord) > -3*pi/6-rotAngle & angle(actuatorCoord) < -pi/6-rotAngle & ...% here
                angle(real(actuatorCoord)-p+1i*imag(actuatorCoord)) < -3*pi/6-rotAngle,...% 3 o'clock petal bottom row
                angle(actuatorCoord) > -3*pi/6-rotAngle & angle(actuatorCoord) < -pi/6-rotAngle & ...
                angle(real(actuatorCoord)+1i*imag(actuatorCoord)+1i*p) > -pi/6-rotAngle];% 5 o'clock petal top-kright row
            
            
            % Plot DM petal edges if verbose
            if view
                nbEdgeAct = {num2str(sum(idx(:,1))), num2str(sum(idx(:,2))),...
                    num2str(sum(idx(:,3))), num2str(sum(idx(:,4))),...
                    num2str(sum(idx(:,5))), num2str(sum(idx(:,6))),...
                    num2str(sum(idx(:,7))), num2str(sum(idx(:,8))),...
                    num2str(sum(idx(:,9))), num2str(sum(idx(:,10))),...
                    num2str(sum(idx(:,11))), num2str(sum(idx(:,12)))};
                figure
                scatter(real(actuatorCoord), imag(actuatorCoord)); hold on;
                scatter(real(actuatorCoord(idx(:,1))), imag(actuatorCoord(idx(:,1))), 'm*');
                scatter(real(actuatorCoord(idx(:,2))), imag(actuatorCoord(idx(:,2))), 'c*');
                scatter(real(actuatorCoord(idx(:,3))), imag(actuatorCoord(idx(:,3))), 'r*');
                scatter(real(actuatorCoord(idx(:,4))), imag(actuatorCoord(idx(:,4))), 'g*');
                scatter(real(actuatorCoord(idx(:,5))), imag(actuatorCoord(idx(:,5))), 'b*');
                scatter(real(actuatorCoord(idx(:,6))), imag(actuatorCoord(idx(:,6))), 'k*');
                scatter(real(actuatorCoord(idx(:,7))), imag(actuatorCoord(idx(:,7))), 'm+');
                scatter(real(actuatorCoord(idx(:,8))), imag(actuatorCoord(idx(:,8))), 'c+');
                scatter(real(actuatorCoord(idx(:,9))), imag(actuatorCoord(idx(:,9))), 'r+');
                scatter(real(actuatorCoord(idx(:,10))), imag(actuatorCoord(idx(:,10))), 'g+');
                scatter(real(actuatorCoord(idx(:,11))), imag(actuatorCoord(idx(:,11))), 'b+');
                scatter(real(actuatorCoord(idx(:,12))), imag(actuatorCoord(idx(:,12))), 'k+');
                hold off; axis square; grid on;
                title('Location of the 12 DM Petal Edges')
                legend('Actuator Location', strcat('Edge #1 (',nbEdgeAct{1}, ' act)'),...
                    strcat('Edge #2 (',nbEdgeAct{2}, ' act)'), strcat('Edge #3 (',nbEdgeAct{3}, ' act)'),...
                    strcat('Edge #4 (',nbEdgeAct{4}, ' act)'), strcat('Edge #5 (',nbEdgeAct{5}, ' act)'),...
                    strcat('Edge #6 (',nbEdgeAct{6}, ' act)'), strcat('Edge #7 (',nbEdgeAct{7}, ' act)'),...
                    strcat('Edge #8 (',nbEdgeAct{8}, ' act)'), strcat('Edge #9 (',nbEdgeAct{9}, ' act)'),...
                    strcat('Edge #10 (',nbEdgeAct{10}, ' act)'), strcat('Edge #11 (',nbEdgeAct{11}, ' act)'),...
                    strcat('Edge #12 (',nbEdgeAct{12}, ' act)'), 'location', 'eastoutside');
            end
        end
        
        %% Find edge actuators up to 4 rows appart from the edge
        function [row1, row2, row3, row4] = whichEdgeExtended(actuatorCoord, rotAngle)
            % [row1, row2, row3, row4] = whichEdgeExtended(actuatorCoord, rotAngle)
            % This function returns the indexes of the rows parallel to the
            % spiders (4 in total). Also see lamTools.whichEdge.
            %
            % actuatorCoord  :: The actuators coordinates (complex). e.g.: bifM4.actuatorCoord
            % rotAngle       :: Define the rotation angle [radians] of the DM;
            % row1 (row2...) :: A matrix [nAct, 12] containing the location of the edges actuators.
            % Edges are arranged in trigonometrical order.
            
            % Examples:
            %    [row1, row2, row3, row4] = lamTools.edgeActuators(actuatorCoord, rotAngle);
            % Noah Schwartz, March 2017
            
            % Calculate the 1 row of actuators near the DM edge
            row1 = lamTools.whichEdge(actuatorCoord, 0.5, rotAngle);
            
            % Second row index
            row2 = zeros(size(actuatorCoord, 1), 12);
            row2(:,1) = angle(actuatorCoord) > -pi/6 & angle(actuatorCoord-0.5-0.5j) < -pi/6& angle(actuatorCoord)<0;row2(row1(:,1),1) = 0;
            row2(:,2) = angle(actuatorCoord) < pi/6 & angle(actuatorCoord+1j) > pi/6& angle(actuatorCoord)>0;        row2(row1(:,2),2) = 0;
            row2(:,3) = angle(actuatorCoord) > pi/6 & angle(actuatorCoord-1j) < pi/6& angle(actuatorCoord-1j)>0;     row2(row1(:,3),3) = 0;
            row2(:,4) = angle(actuatorCoord) < pi/2 & angle(actuatorCoord-1) > pi/2;                                 row2(row1(:,4),4) = 0;
            row2(:,5) = angle(actuatorCoord) > pi/2 & angle(actuatorCoord+1) < pi/2;                                 row2(row1(:,5),5) = 0;
            row2(:,6) = angle(actuatorCoord) < 5*pi/6 & angle(actuatorCoord-1j) > 5*pi/6;                            row2(row1(:,6),6) = 0;
            row2(:,7) = angle(actuatorCoord) > 5*pi/6 & angle(actuatorCoord+1j) < 5*pi/6;                            row2(row1(:,7),7) = 0;
            row2(:,8) = angle(actuatorCoord) < -5*pi/6 & angle(actuatorCoord-1j) > -5*pi/6;                          row2(row1(:,8),8) = 0;
            row2(:,9) = angle(actuatorCoord) > -5*pi/6 & angle(actuatorCoord+1j) < -5*pi/6;                          row2(row1(:,9),9) = 0;
            row2(:,10) = angle(actuatorCoord) < -pi/2 & angle(actuatorCoord+1) > -pi/2;                              row2(row1(:,10),10) = 0;
            row2(:,11) = angle(actuatorCoord) > -pi/2 & angle(actuatorCoord-1) < -pi/2;                              row2(row1(:,11),11) = 0;
            row2(:,12) = angle(actuatorCoord) < -pi/6 & angle(actuatorCoord+1j) > -pi/6 & angle(actuatorCoord+1j)<0; row2(row1(:,12),12) = 0;
            row2 = logical(row2);
            
            % Third row index
            row3 = zeros(size(actuatorCoord, 1), 12);
            row3(:,1) = angle(actuatorCoord) > -pi/6 & angle(actuatorCoord-1-1j) < -pi/6& angle(actuatorCoord)<0;    row3(row1(:,1),1) = 0;  row3(row2(:,1),1) = 0;
            row3(:,2) = angle(actuatorCoord) < pi/6 & angle(actuatorCoord+1.5j) > pi/6& angle(actuatorCoord)>0;        row3(row1(:,2),2) = 0;  row3(row2(:,2),2) = 0;
            row3(:,3) = angle(actuatorCoord) > pi/6 & angle(actuatorCoord-1.5j) < pi/6& angle(actuatorCoord-1.5j)>0;     row3(row1(:,3),3) = 0;  row3(row2(:,3),3) = 0;
            row3(:,4) = angle(actuatorCoord) < pi/2 & angle(actuatorCoord-1.5) > pi/2;                                 row3(row1(:,4),4) = 0;  row3(row2(:,4),4) = 0;
            row3(:,5) = angle(actuatorCoord) > pi/2 & angle(actuatorCoord+1.5) < pi/2;                                 row3(row1(:,5),5) = 0;  row3(row2(:,5),5) = 0;
            row3(:,6) = angle(actuatorCoord) < 5*pi/6 & angle(actuatorCoord-1.5j) > 5*pi/6;                            row3(row1(:,6),6) = 0;  row3(row2(:,6),6) = 0;
            row3(:,7) = angle(actuatorCoord) > 5*pi/6 & angle(actuatorCoord+1.5j) < 5*pi/6;                            row3(row1(:,7),7) = 0;  row3(row2(:,7),7) = 0;
            row3(:,8) = angle(actuatorCoord) < -5*pi/6 & angle(actuatorCoord-1.5j) > -5*pi/6;                          row3(row1(:,8),8) = 0;  row3(row2(:,8),8) = 0;
            row3(:,9) = angle(actuatorCoord) > -5*pi/6 & angle(actuatorCoord+1.5j) < -5*pi/6;                          row3(row1(:,9),9) = 0;  row3(row2(:,9),9) = 0;
            row3(:,10) = angle(actuatorCoord) < -pi/2 & angle(actuatorCoord+1.5) > -pi/2;                              row3(row1(:,10),10) = 0;row3(row2(:,10),10) = 0;
            row3(:,11) = angle(actuatorCoord) > -pi/2 & angle(actuatorCoord-1.5) < -pi/2;                              row3(row1(:,11),11) = 0;row3(row2(:,11),11) = 0;
            row3(:,12) = angle(actuatorCoord) < -pi/6 & angle(actuatorCoord+1.5j) > -pi/6 & angle(actuatorCoord+1.5j)<0; row3(row1(:,12),12) = 0;row3(row2(:,12),12) = 0;
            row3 = logical(row3);
            
            % Fourth row index
            row4 = zeros(size(actuatorCoord, 1), 12);
            row4(:,1) = angle(actuatorCoord) > -pi/6 & angle(actuatorCoord-1.5-1.0j) < -pi/6& angle(actuatorCoord)<0;row4(row1(:,1),1) = 0;  row4(row2(:,1),1) = 0;  row4(row3(:,1),1) = 0;
            row4(:,2) = angle(actuatorCoord) < pi/6 & angle(actuatorCoord+2j) > pi/6& angle(actuatorCoord)>0;        row4(row1(:,2),2) = 0;  row4(row2(:,2),2) = 0;  row4(row3(:,2),2) = 0;
            row4(:,3) = angle(actuatorCoord) > pi/6 & angle(actuatorCoord-2j) < pi/6& angle(actuatorCoord-2j)>0;     row4(row1(:,3),3) = 0;  row4(row2(:,3),3) = 0;  row4(row3(:,3),3) = 0;
            row4(:,4) = angle(actuatorCoord) < pi/2 & angle(actuatorCoord-1.7) > pi/2;                                 row4(row1(:,4),4) = 0;  row4(row2(:,4),4) = 0;  row4(row3(:,4),4) = 0;
            row4(:,5) = angle(actuatorCoord) > pi/2 & angle(actuatorCoord+1.7) < pi/2;                                 row4(row1(:,5),5) = 0;  row4(row2(:,5),5) = 0;  row4(row3(:,5),5) = 0;
            row4(:,6) = angle(actuatorCoord) < 5*pi/6 & angle(actuatorCoord-2j) > 5*pi/6;                            row4(row1(:,6),6) = 0;  row4(row2(:,6),6) = 0;  row4(row3(:,6),6) = 0;
            row4(:,7) = angle(actuatorCoord) > 5*pi/6 & angle(actuatorCoord+2j) < 5*pi/6;                            row4(row1(:,7),7) = 0;  row4(row2(:,7),7) = 0;  row4(row3(:,7),7) = 0;
            row4(:,8) = angle(actuatorCoord) < -5*pi/6 & angle(actuatorCoord-2j) > -5*pi/6;                          row4(row1(:,8),8) = 0;  row4(row2(:,8),8) = 0;  row4(row3(:,8),8) = 0;
            row4(:,9) = angle(actuatorCoord) > -5*pi/6 & angle(actuatorCoord+2j) < -5*pi/6;                          row4(row1(:,9),9) = 0;  row4(row2(:,9),9) = 0;  row4(row3(:,9),9) = 0;
            row4(:,10) = angle(actuatorCoord) < -pi/2 & angle(actuatorCoord+1.7) > -pi/2;                              row4(row1(:,10),10) = 0;row4(row2(:,10),10) = 0;row4(row3(:,10),10) = 0;
            row4(:,11) = angle(actuatorCoord) > -pi/2 & angle(actuatorCoord-1.7) < -pi/2;                              row4(row1(:,11),11) = 0;row4(row2(:,11),11) = 0;row4(row3(:,11),11) = 0;
            row4(:,12) = angle(actuatorCoord) < -pi/6 & angle(actuatorCoord+2j) > -pi/6 & angle(actuatorCoord+2j)<0; row4(row1(:,12),12) = 0;row4(row2(:,12),12) = 0;row4(row3(:,12),12) = 0;
            row4 = logical(row4);
            % % Plot
            % figure
            % scatter(real(actuatorCoord), imag(actuatorCoord));grid on
            % hold on
            % for ii=1:12; scatter(real(actuatorCoord(row1(:,ii))), imag(actuatorCoord(row1(:,ii))), 'c'); pause(0.0);end;
            % for ii=1:12; scatter(real(actuatorCoord(row2(:,ii))), imag(actuatorCoord(row2(:,ii))), 'r'); pause(0.0);end;
            % for ii=1:12; scatter(real(actuatorCoord(row3(:,ii))), imag(actuatorCoord(row3(:,ii))), 'k'); pause(0.0);end; pause(0);
            % for ii=1:12; scatter(real(actuatorCoord(row4(:,ii))), imag(actuatorCoord(row4(:,ii))), 'm'); pause(0.0);end;
            % axis square; hold off
        end
        
        %% applyPetals
        function applyPetals(dm,actuatorCoord,petals,rotAngle)
            idx = lamTools.whichPetal(actuatorCoord, rotAngle,4);    % Find where the actuators belong (i.e. which petal/segment)
            for kActuator=1:dm.nActuator        % Apply petals to actuators influence functions
                dm.modes.modes(:,kActuator) = dm.modes.modes(:,kActuator).*petals(:,idx(kActuator));
            end
        end
        %% applyM4petals
        function [out,petals] = applyM4petals(dm,pupil,actuatorCoord,flagPupilLayout,rotAngle)
            % function out = applyM4petals(dm,pupil,actuatorCoord)
            %
            % dm            :: The DM object.
            % actuatorCoord :: The actuators coordinates (complex).
            % rotAngle      :: Define the rotation angle [radians] of the DM;
            % pupil         :: The pupil of the telescope (e.g. tel.pupil).
            % out--> TBD.
            if ~exist('flagPupilLayout','var') || isempty(flagPupilLayout)
                flagPupilLayout = 'allGlass';
            end
            if ~exist('rotAngle','var') || isempty(rotAngle)
                rotAngle = 0;
            end
            if strcmp(flagPupilLayout,'allGlassNoSpiders') % Apply the M1 pupil to the DM influence functions. ATTENTION: M4 pupil taken to be M1's. To be updated later.
                    dm.modes.modes = bsxfun(@times,dm.modes.modes,pupil(:));
                    idx = 1;
            else
                petals = lamTools.m4Petals(pupil, rotAngle);           % Define the 6 petals
                %imagesc(petals(:,:,1)+petals(:,:,2)*2+petals(:,:,3)*3+petals(:,:,4)*4+petals(:,:,5)*5+petals(:,:,6)*6); axis equal tight xy;
                sz = size(pupil);
                petals = reshape(petals, sz(1)*sz(2), 6);
                idx = lamTools.whichPetal(actuatorCoord, rotAngle,6);    % Find where the actuators belong (i.e. which petal/segment)
                
                %for kActuator=1:dm.nActuator        % Apply petals to actuators influence functions
                %    dm.modes.modes(:,kActuator) = dm.modes.modes(:,kActuator).*petals(:,idx(kActuator));
                %end
                
%                 % Temp
%                 scatter(real(actuatorCoord(idx ==1)), imag(actuatorCoord(idx == 1))); hold on
%                 scatter(real(actuatorCoord(idx ==2)), imag(actuatorCoord(idx == 2)), 'r');
%                 scatter(real(actuatorCoord(idx ==3)), imag(actuatorCoord(idx == 3)), 'k');
%                 scatter(real(actuatorCoord(idx ==4)), imag(actuatorCoord(idx == 4)), 'y');
%                 scatter(real(actuatorCoord(idx ==5)), imag(actuatorCoord(idx == 5)), 'm');
%                 scatter(real(actuatorCoord(idx ==6)), imag(actuatorCoord(idx == 6))', 'g'); hold off
%                 
%                 figure; imagesc((sum(reshape(full(dm.modes.modes),108,108,358),3))); colorbar;
%              for ii=1:size(dm.modes.modes,2)
% %                  imagesc(reshape(dm.modes.modes(:,ii), 108,108)); drawnow
% %                  scatter(real(actuatorCoord(ii)), imag(actuatorCoord(ii)));
% %                  drawnow; title(idx(ii))
%              end
%              % temp
             
                petals=sparse(petals);
                dm.modes.modes = dm.modes.modes .* petals(:,idx);
            end
            out = idx;
        end
        %         function out = applyM4petals(dm,pupil,range)
        %             if ~exist('range','var')
        %                 range = 1.5; % range in units of meters beyond which an actuator must be not to be affected by cropping at the edges
        %             end
        %             nRes = sqrt(size(dm.modes.modes,1));
        %             petals = zeros(nRes,nRes,6);
        %             x = linspace(-1,1,nRes);
        %             u = x ;
        %             v = x ;
        %
        %             [x,y,r,o] = utilities.cartAndPol(u,v);
        %             for kPetal = 1:6
        %                 angleEnd(kPetal) = lamTools.wraptopi(pi/6 - (kPetal-1)*2*pi/6);
        %                 angleInit(kPetal)   = lamTools.wraptopi(pi/6 - (kPetal)*2*pi/6);
        %                 petals(:,:,kPetal) = double(r <= 1+1e-3 & r >= 0.2 & o > angleInit(kPetal) & o < angleEnd(kPetal));%.*pupil;
        %                 if kPetal == 4
        %                     petals(:,:,kPetal) = double(r <= 1+1e-3 & r >= 0.2 & (o > angleInit(kPetal) | o < angleEnd(kPetal)));%.*pupil;
        %                 end
        %             end
        % %             petals(:,:,4) = pupil - sum(petals,3);
        %
        %             petals = reshape(petals,nRes^2,6);
        %             angleActCoord = unwrap(angle(dm.modes.actuatorCoord));
        %             angleInit = -pi/6:pi/3:2*pi;
        %             angleEnd = angleInit + pi/3;
        %             isCloseToBoundary = [];
        %             rMin = 5;
        %             rMax = 19;
        %             for kActuator = 1:dm.nActuator
        %                 whichSector(kActuator) = find((angleActCoord(kActuator) > angleInit & angleActCoord(kActuator) < angleEnd));
        %                 radius = abs(dm.modes.actuatorCoord(kActuator));
        %                 theta = angle(dm.modes.actuatorCoord(kActuator));
        %                 rth = [radius*exp(1i*angleInit), rMin*exp(1i*theta), rMax*exp(1i*theta)];
        %                 deltaX = real(rth) - real(dm.modes.actuatorCoord(kActuator));
        %                 deltaY = imag(rth) - imag(dm.modes.actuatorCoord(kActuator));
        %                 dist = sqrt(deltaX.^2 + deltaY.^2);
        %                 if any(dist < range)
        %                     isCloseToBoundary = [isCloseToBoundary kActuator];
        %                 end
        %             end
        %
        %             for kActuator = isCloseToBoundary
        %                 dm.modes.modes(:,kActuator) = dm.modes.modes(:,kActuator).*petals(:,whichSector(kActuator));
        %             end
        %             figure
        %             scatter(real(dm.modes.actuatorCoord), imag(dm.modes.actuatorCoord)), hold
        %             scatter(real(dm.modes.actuatorCoord(isCloseToBoundary)), imag(dm.modes.actuatorCoord(isCloseToBoundary)))
        %             box on, axis square
        %             title('M4 actuator locations; red: actuators affected by edge cropping')
        %             xlabel('position, [m]')
        %             ylabel('position, [m]')
        %             out = isCloseToBoundary;
        %         end
        %% wraptopi
        function a = wraptopi(a, a_center )
            % function a = wrap(a,a_center)
            %
            % Wraps angles to a range of 2*pi.
            % Inverse of Matlab's "unwrap", and better than wrapToPi ( which has
            % redundant [-pi,pi])
            % Optional input "a_center" defines the center angle.  Default is 0, giving
            % angles from (-pi,pi], chosen to match angle(complex(-1,0)).  Maximum
            % possible value is pi.
            
            % T.Hilmer, UH
            % 2010.10.18 version 2
            %   removed code from version 1. Have not bug-checked second input
            %   "a_center"
            
            if nargin < 2, a_center = 0; end
            
            % new way
            a = mod(a,2*pi); % [0 2pi)
            
            % shift
            j = a > pi - a_center;
            a(j) = a(j) - 2*pi;
            j = a < a_center - pi;
            a(j) = a(j) + 2*pi;
        end
        
        %% multipleSpectra
        function out = multiZernikeTemporalSpectra(nu,atm,zern,tel)
            nMode = zern.nMode;
            out = zeros(nMode, length(nu),atm.nLayer);
            parfor kmode=1:nMode
                out(kmode,:,:)  = lamTools.temporalSpectrum(nu,atm,tel,kmode);
            end
            
        end
        %% type I OMGI optimiser
        function g = omgi(nu, PSD, samplingFreq, pureDelay, noiseVar,verbose)
            if ~exist('verbose','var')
                verbose = 0;
            end
            T = 1/samplingFreq;
            s = @(x) 2*1i*pi*x;
            
            if size(PSD,1) > size(PSD,2)
                PSD = PSD';
            end
            
            % TFs
            hWfs = @(x) (1-exp(-s(x)*T))./(s(x)*T);
            hDac = hWfs; % Different meaning, same value
            hLag = @(x) exp(-pureDelay*s(x));
            hInt = @(x) 1./(1-exp(-s(x)*T));
            hDm  = @(x) 1;
            G = @(x,g) hWfs(x) .* hDac(x) .* hLag(x) .* g .* hInt(x) .* hDm(x);
            Gn = @(x,g) hDac(x) .* hLag(x) .* g .* hInt(x) .* hDm(x);
            % rejection transfer function
            E = @(x,g) abs(1./(1+G(x,g)));
            %Noise transfer
            hN = @(x,g) Gn(x,g) .* E(x,g);%./hWfs(x);
            %myfun = @(g) sum(hN(nu)*noiseVar + RTF*PSD);
            g = fminsearch(@(g)lamTools.integrate(g,nu,E,hN,PSD,noiseVar,T), 0.1);
            
            if verbose
                gains = linspace(0.001,0.5,100);
                for kGain = 1:length(gains)
                    outN(kGain) = trapz(nu, abs(hN(nu,gains(kGain))).^2*noiseVar*2*T);
                    outS(kGain) = trapz(nu, abs(E(nu,gains(kGain))).^2.*PSD);
                end
                figure
                semilogy(gains, outS,'r')
                hold on
                semilogy(gains, outN,'b')
                semilogy(gains, outS + outN,'k')
                xlabel('gain')
                ylabel('residual')
                title('single integrator optimal modal gian')
                legend('signal residual','noise residual','total residual')
                outNopt = trapz(nu, abs(hN(nu,g)).^2*noiseVar*2*T);
                outSopt = trapz(nu, abs(E(nu,g)).^2.*PSD);
                plot(g,outNopt + outSopt,'ro')
            end
        end
        
        function out = integrate(g,nu,E,hN, PSD, noiseVar,T)
            outN = trapz(nu, abs(hN(nu,g)).^2*noiseVar*2*T);
            outS = trapz(nu, abs(E(nu,g)).^2.*PSD);
            out = outS + outN;
        end
        
        %%
        function out = doubleIntParamOptim(nu, PSD, samplingFreq, pureDelay, varNoise)
            %% DOUBLE INTEGRATOR PARAMETER OPTIMISATION
            
            % out = doubleIntParamOptim(nu, PSD, samplingFreq, pureDelay,
            % varNoise) computes the double integrator gain and lead-filter
            % parameters for a 45 phase margin stability.
            % nu            :: the temporal frequency vector
            % PSD           :: the temporal PSD vector in units^2
            % samplingFreq  :: the loop sampling frequency in Hz
            % pureDelay     :: the loop pure delay
            % varNoise      :: noise variance in units^2
            % Created       :: C. Correia, Dec'15
            % Comment: based on do_optim_typeII from HIA and TMT,
            % developed originally by JPVeran
            
            T = 1/samplingFreq;
            s = @(x) 2*1i*pi*x;
            
            g0 = 1e-3;
            varTotalOld = inf;
            
            % OPEN-LOOP TRANSFER FUNCTION
            
            %             G = @(x) ((1-exp(-s(x)*T))./(s(x)*T)).^2.*...   % hWfs = @(x) (1-exp(-s(x)*T))./(s(x)*T); hDac = hWfs; % Different meaning, same value
            %                 exp(-tau*s(x)).*...                         % hLag = @(x) exp(-tau*s(x));
            %                 (1./(1-exp(-s(x)*T)).^2).*...               % hInt = @(x) 1./(1-exp(-s(x)*T)); squared for the double integrator
            %                 1;                                          % hDm  = @(x) 1;
            
            hWfs = @(x) (1-exp(-s(x)*T))./(s(x)*T);
            hDac = hWfs; % Different meaning, same value
            hLag = @(x) exp(-pureDelay*s(x));
            hInt = @(x) 1./(1-exp(-s(x)*T));
            hDm  = @(x) 1;
            G = @(x) hWfs(x) .* hDac(x) .* hLag(x) .* hInt(x) .* hInt(x) .* hDm(x);
            gcc = [];
            
            while 1
                % Open loop gain without phase lead
                %hOl = @(x) g0*G(x);
                
                [phaseMargin, fcross] = calcPhiMargin(g0*G(nu), nu);
                
                % Phase lead needed
                additionalPhaseNeeded=pi/4-phaseMargin;
                
                % calculating the lead filter parameters to achieve the required phase lead
                a=(1-sin(additionalPhaseNeeded))/(1+sin(additionalPhaseNeeded));
                f0=fcross*sqrt(a); %fprintf('fs=%g, fcross=%f\n',fs,fcross);
                Tlead=1/(2*pi*f0);
                
                % gain 1/g created by phase lead filter
                g=sqrt(a);
                
                % complete Hol. g is here to adjust g0 to remove the scaling caused by lead filter.
                %hOl1 = @(x) g * hOl(x) .* (1+Tlead*s(x))./(1+a*Tlead*s(x));
                hOl1 = @(x) g * g0*G(x) .* (1+Tlead*s(x))./(1+a*Tlead*s(x));
                % Rejection transfer function
                E = @(x) abs(1./(1+hOl1(x)));
                
                % Closed-loop transfer function
                %Ecl = @(x) hOl1(x) .* E(x);
                %Noise transfer
                hN = @(x) hOl1(x) .* E(x);%./hWfs(x);
                
                % RESIDUAL SIGNAL VARIANCE
                %rms=trapz(nus, PSD.*abs(Hrej.*Hrej));
                %varSignal = 2*quadgk( @(nu) PSD.*abs(E(nu)).^2 , 0 , Inf);
                varSignal = trapz(nu, PSD.*abs(E(nu)).^2);
                % NOISE GAIN
                %noiseGain=(trapz(nus, abs(Hn.*Hn))./trapz(nus, ones(length(nus),1)));
                %noiseGain = 2*quadgk( @(nu) PSD.*abs(hN(nu)).^2 , 0 , Inf)
                noiseGain = trapz(nu, abs(hN(nu)).^2)/(1/2/T);
                %Tot Error.
                varPropNoise=noiseGain*varNoise;
                varTotal=varSignal+varPropNoise;
                
                if varTotal < varTotalOld
                    %increase the gain.
                    g0=g0*1.01;
                    gcc = [gcc g0];
                    varTotalOld = varTotal;
                    out{1} = g0*g;
                    out{2} = a;
                    out{3} = Tlead;
                    out{4} = varSignal;
                    out{5} = varPropNoise;
                    out{6} = noiseGain;
                    out{7} = G;
                    out{8} = hN;
                    
                else
                    break
                end
            end
            function [margin, fcross]=calcPhiMargin(Hol, nus)
                %Computes Phase margin of the open loop transfer function.
                %Defined as the phase of Hol when abs(Hol) crosses 1.
                ind=abs(Hol)>1;
                indl=find(ind);
                indl=indl(end);
                if indl==length(Hol)
                    ph2=angle(Hol(end));
                    fcross=nus(end);
                else
                    abs1=log10(abs(Hol(indl:indl+1)));%interpolate the logrithm which is close to linear.
                    frac1=abs1(1)/(abs1(1)-abs1(2));%linear interpolation.
                    ph1=(angle(Hol(indl:indl+1)));
                    %Reduce distance
                    diff=(ph1(2)-ph1(1))/(2*pi);
                    diff=diff-fix(diff);
                    ph1(2)=ph1(1)+2*pi*diff;
                    
                    ph2=ph1(1)-(ph1(1)-ph1(2))*frac1;
                    
                    nu1=nus(indl:indl+1);
                    fcross=nu1(1)-(nu1(1)-nu1(2))*frac1;
                end
                %bring [-pi;pi] to [-2pi;0];
                normph=ph2/(2*pi);
                normph=normph-floor(normph)-1;
                ph2=normph*2*pi;
                if ph2>0
                    ph2=ph2-2*pi;
                end
                margin=pi+ph2;
            end
            
        end
        
        function out = ZernikeTemporalSpectrum(nu,atm,tel,kmode)
            %% SPECTRUM Phase power spectrum density
            %
            % out = phaseStats.spectrum(f,atm) computes the phase power
            % spectrum density from the spatial frequency f and an
            % atmosphere object
            %
            % See also atmosphere
            
            out = zeros(size(nu,2),atm.nLayer);
            for kLayer = 1:atm.nLayer
                nPx = size(atm.layer(kLayer).phase,1);
                pupil = logical(utilities.piston(nPx));
                D = tel.diameterAt(atm.layer(kLayer).altitude);
                
                zern = zernike(kmode+1,'resolution',nPx,'pupil',pupil,'D',D);
                
                atmSlab = slab(atm,kLayer);
                [vx,vy] = pol2cart(atmSlab.layer.windDirection,atmSlab.layer.windSpeed);
                for k=1:numel(nu)
                    if abs(vx)>eps(atmSlab.layer.windSpeed)
                        out(k,kLayer) = out(k,kLayer) + quadgk( @integrandFy , -Inf, Inf);
                    else
                        out(k,kLayer) = out(k,kLayer) + quadgk( @integrandFx , -Inf, Inf);
                    end
                end
                if vx == 0;
                    signvx = 1;
                else
                    signvx = sign(vx);
                end
                if vy == 0;
                    signvy = 1;
                else
                    signvy = sign(vy);
                end
                out(:,kLayer) = out(:,kLayer)*signvx*signvy;
                %%%%
                % --- apply wind direction using Conan95, Eq. 30
            end
            
            function int = integrandFy(fy)
                fx = (nu(k) -fy*vy)/vx;
                %fx = (nu(k))/vx;
                int = zernikeStats.spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab , zern)/vx;
            end
            
            function int = integrandFx(fx)
                fy = (nu(k) -fx*vx)/vy;
                int = zernikeStats.spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab , zern)/vy;
            end
        end
        
        function plotDMWFSgrid(tel,dm,wfs)
           
            figure;
            
            % plot dm actuator
            plot(real(dm.modes.actuatorCoord),imag(dm.modes.actuatorCoord),'.k');

            %
            n = size(wfs.validLenslet,1);
            [sx,sy] = meshgrid(linspace(-1,1,n)*(tel.D/2-tel.D/n/2));
            hold on;
            plot(sx(wfs.validLenslet),sy(wfs.validLenslet),'or','MarkerSize',1);
            plot(sx(~wfs.validLenslet),sy(~wfs.validLenslet),'ob','MarkerSize',1);
            hold off;

        end
    end
end