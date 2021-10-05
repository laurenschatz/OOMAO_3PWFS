classdef m4dm < deformableMirror
    %M4DM Create a deformable mirror object with the real ELT M4 configuration
    %
    % Examples:
    %
    % % Continuous M4
    % dm = elt.m4dm(dmCrossCouplingCoeff,...
    %    'resolution',tel.resolution,...
    %    'pupil',tel.pupil,...
    %    'isSegmented',false,...
    %    'isSlaving',false);
    %
    % % Segmented M4
    % dm = elt.m4dm(dmCrossCouplingCoeff,...
    %    'resolution',tel.resolution,...
    %    'pupil',tel.pupil,...
    %    'isSegmented',true,...
    %    'isSlaving',false);
    %
    % % Segmented M4 with actuator slaving
    % dm = elt.m4dm(dmCrossCouplingCoeff,...
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
        function obj = m4dm(dmCrossCouplingCoeff,varargin)
            
            p = inputParser;
            p.addRequired('dmCrossCouplingCoeff', @isnumeric);
            p.addParameter('resolution', [], @isnumeric);
            p.addParameter('zLocation', 0, @isnumeric);
            p.addParameter('offset', [0,0], @isnumeric);
            p.addParameter('diameter', 37, @isnumeric)
            p.addParameter('influenceFunctionType', 'gaussian', @ischar)
            p.addParameter('isSegmented',false,@islogical);
            p.addParameter('pupil', @isnumeric);
            p.addParameter('isSlaving',false,@islogical);
            
            p.parse(dmCrossCouplingCoeff, varargin{:});
            
            % load actuator model
            actCentrePos            = fitsread('/result/SCAO-H/INPUT/m4true_pos_m.fits'); % m4 actuator coord : estimated from m4 IF COG
            dmPitch                 = 0.5;                            % normalisation : directly given on M1 coordinates

            % influenceFunction
            if strcmp(p.Results.influenceFunctionType,'gaussian') || strcmp(p.Results.influenceFunctionType,'Gaussian')
                bifM4               = gaussianInfluenceFunction(dmCrossCouplingCoeff, dmPitch);
            elseif strcmp(p.Results.influenceFunctionType,'bezier') || strcmp(p.Results.influenceFunctionType,'Bezier')
                bifM4               = nfluenceFunction(dmCrossCouplingCoeff, dmPitch);
            end
            bifM4.actuatorCoord     =  (actCentrePos(:,1) + 1j*actCentrePos(:,2));
            idx                     = (abs(bifM4.actuatorCoord) <= 38/2) & (abs(bifM4.actuatorCoord) >= 10.1/2); % keep actuators within a distance below 1 dm pitch off of the border
            bifM4.actuatorCoord     =  (actCentrePos(idx,1) + 1j*actCentrePos(idx,2));
            nAct                    = nnz(idx);
            
            
            % deformable mirror object instantiation
            obj = obj@deformableMirror(nAct,...
                'modes',bifM4,...
                'resolution',p.Results.resolution,...
                'validActuator',true(1,nAct),...
                'diameter', p.Results.diameter);
            
            % apply the segmentation by 6 petals
            if p.Results.isSegmented
                lamTools.applyM4petals(obj, p.Results.pupil, bifM4.actuatorCoord, 'allGlass', 0);
            end
            
            if p.Results.isSlaving
                %apply slaving
                [newDmModes,newActuatorCoord] = lamTools.groupIF(obj.modes.modes,bifM4.actuatorCoord,0);
                nAct = size(newDmModes,2);
                bifM4.actuatorCoord = newActuatorCoord;
                %tmpDM = deformableMirror(nAct,...
                %    'modes',bifM4,...
                %    'resolution',p.Results.resolution,...
                %    'validActuator',true(1,nAct),...
                %    'diameter', p.Results.diameter);
                %obj = tmpDM;
                vA = obj.validActuator;
                vA(nAct+1:end)=false;
                obj.validActuator = vA;
                obj.modes = bifM4;
                obj.modes.modes = newDmModes;
            end

            
            obj.log = logBook.checkIn(obj);
            display(obj)
            
            obj.surfaceListener     = addlistener(obj,'surface','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.surfaceListener.Enabled = false;
        end
    end
    
end