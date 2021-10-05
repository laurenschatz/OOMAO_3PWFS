classdef m4dmScaledown < deformableMirror
    %M4DMSCALEDOWN Create a deformable mirror object with the scaled-downed M4 configuration
    %
    % Examples:
    %
    % % Continuous scale-downed M4
    % dm = elt.m4dmScaledown(dmCrossCouplingCoeff,tel,...
    %    'resolution',tel.resolution,...
    %    'pupil',tel.pupil,...
    %    'isSegmented',false,...
    %    'isSlaving',false);
    %
    % % Segmented scale-downed M4
    % dm = elt.m4dmScaledown(dmCrossCouplingCoeff,tel,...
    %    'resolution',tel.resolution,...
    %    'pupil',tel.pupil,...
    %    'isSegmented',true,...
    %    'isSlaving',false);
    %
    % % Segmented scale-downed M4 with actuator slaving
    % dm = elt.m4dmScaledown(dmCrossCouplingCoeff,tel,...
    %    'resolution',tel.resolution,...
    %    'pupil',tel.pupil,...
    %    'isSegmented',true,...
    %    'isSlaving',true);
    %
    % See also deformableMirror
    % See also influenceFunction
    properties
        % influence function model: Gaussian, Bezier, ...
        influenceFunctionType
    end
    
    properties (SetAccess=protected)
        % flag for M4 segmentation
        isSegmented;
        % flag for slaving actuators
        isSlaving;
        %
        log;
    end
    
    methods
        
        %% Constructor
        function obj = m4dmScaledown(dmCrossCouplingCoeff,tel,varargin)
            
            p = inputParser;
            p.addRequired('dmCrossCouplingCoeff', @isnumeric);
            p.addRequired('telescope', @(x) isa(x,'telescope'));
            p.addParameter('resolution', [], @isnumeric);
            p.addParameter('zLocation', 0, @isnumeric);
            p.addParameter('offset', [0,0], @isnumeric);
            p.addParameter('influenceFunctionType', 'gaussian', @ischar)
            p.addParameter('isSegmented',false,@islogical);
            p.addParameter('pupil', @isnumeric);
            p.addParameter('isSlaving',false,@islogical);
            
            p.parse(dmCrossCouplingCoeff, tel, varargin{:});
            
            % compute actuator position of the scale downed M4
            dmPitch                 = 0.5; % normalisation : directly given on M1 coordinates
            [px,py]=lamTools.M4downScaled(tel.D);
            
            % influenceFunction
            if strcmp(p.Results.influenceFunctionType,'gaussian') || strcmp(p.Results.influenceFunctionType,'Gaussian')
                bifM4               = gaussianInfluenceFunction(dmCrossCouplingCoeff, dmPitch);
            elseif strcmp(p.Results.influenceFunctionType,'bezier') || strcmp(p.Results.influenceFunctionType,'Bezier')
                bifM4               = influenceFunction(dmCrossCouplingCoeff, dmPitch);
            end
            % Select "Valid" actuators
            bifM4.actuatorCoord = (px(:) + 1j*py(:));
            idx                 = (abs(bifM4.actuatorCoord) <= (tel.D+dmPitch*2)/2) &...
                (abs(bifM4.actuatorCoord) >= (tel.D*tel.obstructionRatio-dmPitch*2)/2); % keep actuators within a distance below 1 dm pitch off of the border
            bifM4.actuatorCoord = (px(idx) + 1j*py(idx));
            nAct                = nnz(idx);
            
            % test with 37m telescope
            %bifM4.actuatorCoord = (px(:) + 1j*py(:));
            %idx                 = (abs(bifM4.actuatorCoord) <= (tel.D+2*dmPitch)/2) &...
            %    (abs(bifM4.actuatorCoord) >= tel.D*tel.obstructionRatio/2-dmPitch); % keep actuators within a distance below 1 dm pitch off of the border
            %bifM4.actuatorCoord = (px(idx) + 1j*py(idx));
            %nAct                = nnz(idx);
            
            % deformable mirror object instantiation
            obj = obj@deformableMirror(nAct,...
                'modes',bifM4,...
                'resolution',p.Results.resolution,...
                'validActuator',true(1,nAct),...
                'diameter', tel.D);
            
            % Apply segmentation
            if p.Results.isSegmented
                lamTools.applyM4petals(obj, tel.pupil, bifM4.actuatorCoord, 'allGlass', 0);
            else
                lamTools.applyM4petals(obj, tel.pupil,bifM4.actuatorCoord, 'allGlassNoSpiders', 0);
            end
            
            % Apply slaving
            
            if p.Results.isSlaving
                %apply slaving
                [newDmModes,newActuatorCoord] = lamTools.groupIFScaledown(obj.modes.modes,bifM4.actuatorCoord,0);
                nAct = size(newDmModes,2);
                bifM4.actuatorCoord = newActuatorCoord;
                vA = obj.validActuator;
                vA(nAct+1:end)=false;
                obj.validActuator = vA;
                obj.modes = bifM4;
                obj.modes.modes = newDmModes;
                obj.nActuator = nAct;
            end
            
            obj.log = logBook.checkIn(obj);
            display(obj)
            
            obj.surfaceListener     = addlistener(obj,'surface','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.surfaceListener.Enabled = false;
        end
    end

end