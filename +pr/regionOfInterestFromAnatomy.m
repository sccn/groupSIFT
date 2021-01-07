classdef regionOfInterestFromAnatomy < pr.regionOfInterest
    % holds information and methods (functions) to work with a (optionally) probabilistic region of interest (ROI).
    % Projected measure, Dipole desnity or a combination of these.
    
    properties (SetAccess = 'protected')
        keyword
        anatomicalLabel % cell array containing anatomical labels
        atlas
    end; % properties that cannot be changed outside the condtructor
    
    properties (Dependent = true)
        label % a string label for the ROI, obtained by concatenating all anatomical labels
    end
    
    methods(Static)
        function list = getAllAnatomicalLabels(varargin)
            
            inputOptions = finputcheck(varargin, ...
                {'atlas'       'string'   {'loni', 'aal'}     'aal';...
                'onlyEEGSources'    'boolean'   [] true;...                
                'includeHeschl'  'boolean'  [] false;... % whether to inclide the very time (few 2-3 discontineous voxels) Heschl area.
                });
            
            
            if strcmpi(inputOptions.atlas, 'loni') % from LONI
                anatomicalDataForHeadGrid = load([pr.getToolboxFolder '+pr/precompute/' 'anatomical_information_for_headGrid_8_mm.mat']);
                if inputOptions.onlyEEGSources
                    anatomicalDataForHeadGrid.brainArealabel = anatomicalDataForHeadGrid.brainArealabel(1:end-2); % exclude brain stem
                    anatomicalDataForHeadGrid.brainArealabel([7 34]) = []; % exclude Hippocampus (?)
                end;
                list = anatomicalDataForHeadGrid.brainArealabel;
            else  % from AAL
                headGrid = pr.headGrid;
                [labelIdCube areaLabel] = headGrid.getAALData;
                if inputOptions.onlyEEGSources
                    idToRemove  = 91:length(areaLabel); % these are Cerebelum and Vermis (around cerebellum), so they are not likely source of EEG.
                    
                    if ~inputOptions.includeHeschl
                        idToRemove = [idToRemove 79 80]; % Heschl_L and Heschl_R
                    end;
                    
                    areaLabel(idToRemove) = []; 
                end;
                list = areaLabel;
            end;
        end;
    end;
    methods
        function obj = regionOfInterestFromAnatomy(headGrid, keyword, varargin)
            superClassArgs = {};
            
            if nargin>1
                superClassArgs{1} = headGrid;
            end;
            
            if nargin>2
                for i=1:length(varargin) % pass additional parameters in varargin to the parent class
                    superClassArgs{1+i} = varargin{i};
                end;
            end;
            
            obj = obj@pr.regionOfInterest(superClassArgs{:}); % call the super class constructor
            
            obj.keyword = keyword;

              inputOptions = finputcheck(varargin, ...
                {'atlas'       'string'   {'loni', 'aal'}     'aal';...
                'partialMatch' 'boolean'  [] false;...  % find all anatomical labels that their names partially match the keyboard, for ex
                });
            
            obj.atlas = inputOptions.atlas;
            
            if strcmpi(obj.atlas, 'loni') % from LONI
                
                anatomicalDataForHeadGrid = headGrid.getAnatomicalData;
                
                if inputOptions.partialMatch
                    regionIds = [];
                    for i=1:length(anatomicalDataForHeadGrid.brainArealabel)
                        if ~isempty(strfind(lower(anatomicalDataForHeadGrid.brainArealabel{i}), lower(keyword)))
                            regionIds = [regionIds i];
                        end;
                    end;
                else
                     regionIds = find(strcmpi(keyword, anatomicalDataForHeadGrid.brainArealabel), 1);
                end;
                
                obj.anatomicalLabel = anatomicalDataForHeadGrid.brainArealabel(regionIds);
                
                obj.membershipProbabilityCube = double(headGrid.insideBrainCube * 0);
                for iRegion = 1:length(regionIds)
                    obj.membershipProbabilityCube(headGrid.insideBrainCube(:)) = obj.membershipProbabilityCube(headGrid.insideBrainCube(:)) + anatomicalDataForHeadGrid.probabilityOfEachLocationAndBrainArea(:,regionIds(iRegion));
                end;
            else % from AAL
                [labelIdCube areaLabel] = headGrid.getAALData;
                
                if inputOptions.partialMatch
                regionIds = [];
                for i=1:length(areaLabel)
                    if ~isempty(strfind(lower(areaLabel{i}), lower(keyword)))
                        regionIds = [regionIds i];
                    end;
                end;
                else
                    regionIds = find(strcmpi(keyword, areaLabel), 1);
                end;
                
                obj.anatomicalLabel = areaLabel(regionIds);
                
                obj.membershipProbabilityCube = double(headGrid.insideBrainCube * 0);
                for iRegion = 1:length(regionIds)
                    obj.membershipProbabilityCube = obj.membershipProbabilityCube +  double(labelIdCube == regionIds(iRegion));
                end;                
                
                obj.membershipProbabilityCube(:) = min(1, obj.membershipProbabilityCube(:));
            end;
        end;
        
         function varargout = plotVolume(obj, varargin)
             plotVolume@pr.regionOfInterest(obj, varargin{:});
             set(gcf, 'name', cell2mat(obj.anatomicalLabel));
         end;
         
         function roiLabel = get.label(obj)
             roiLabel = strjoin(',', obj.anatomicalLabel);
         end;
    end
end