classdef headGrid % a grid used in Measure Projection
    properties
        spacing  = 8; % space between grid points (in mm).
        minGridLocation = [-88 -121  -77];
        maxGridLocation = [89   90   99];
        xCube, yCube, zCube % 3d array containing x, y or z values of grid points, meshgrid() output style.
        insideBrainCube % for each point in the 3D grid, it specifies whetehr it is inside brain volume or not.
        cubeSize
    end % properties
    methods
        function obj = headGrid(spacing, minGridLocation, maxGridLocation) % constructor
            
            if nargin > 0
                obj.spacing = spacing;
            else
                obj.spacing = 8; % mm
            end;                        
            
            if nargin > 1
                obj.minGridLocation = minGridLocation;
            end;
            
            if nargin > 2
                obj.maxGridLocation = maxGridLocation;
            end;
            
            obj.spacing = round(obj.spacing * 10) / 10; % round spaicng to one point after zero (up to 0.1). For exmaple 8.22 becomes 8.2. This is to be able toi use precomputed anatomical and Brodmann values.

            [obj.xCube obj.yCube obj.zCube] = meshgrid(obj.minGridLocation(1):obj.spacing:obj.maxGridLocation(1), obj.minGridLocation(2):obj.spacing:obj.maxGridLocation(2), obj.minGridLocation(3):obj.spacing:obj.maxGridLocation(3));
            
            obj.cubeSize = size(obj.xCube);
            
            % first try loading insideBrainCube from a computed file (to increase robustness by not
            % using some mex files in fieldtrip), if the file did not exist, compute and save into a
            % file for future use.
            
            fullPath = which('pr.dipole');
            prFolderPath = fullPath(1:end - length('dipole.m'));
            precomputePath = [prFolderPath  'precompute/'];
            
            filename = [precomputePath 'inside_brain_information_for_headGrid_' num2str(obj.spacing) '_mm.mat'];
            if exist(filename, 'file') % see if the file exist to load it from
                loadedData = load(filename);
                obj.insideBrainCube = loadedData.insideBrainCube;
            else % make it, since it takes some time or it may not work on some Macs and such, we will save the results in /+pr/precompute
                
                % find grid locations that are inside brain
                gridPosition = [obj.xCube(:) obj.yCube(:) obj.zCube(:)];
                
                dipfitdefs; % get the location of standard BEM volume file
                tmp = load('-mat',DIPOLEDENSITY_STDBEM); % load MNI mesh
                
                insideBrainCube = reshape(inside_vol(gridPosition, tmp.vol), size(obj.xCube));
                obj.insideBrainCube = insideBrainCube;
                
                % save for future use               
                try % try saving, but prevent errors if save privilages were not available
                    save(filename, 'insideBrainCube');
                catch err
                    fprintf([err.message '\n']);
                    fprintf('Measure Projection: saving precomputed inside-brain information was not successful.\n');
                end;
            end;
        end;
        function plotScatter(obj) % plot a general scatter overview of the grid and brain
            figure;
            scatter3(obj.xCube(:), obj.yCube(:), obj.zCube(:), 1+obj.insideBrainCube(:) * 50, 'filled');
            axis equal;
        end;
        function plot(obj) % plot a general volumetric overview of the grid and brain
            figure;
            plot_dipplot_with_cortex;
            pr.plot_head_region(obj, obj.insideBrainCube, 'regionColor', 'g');
        end;
        function position = getPosition(obj, id)
            % position = getPosition(obj, id)
            % return an N x 3 array containing grid positions.
            
            if nargin < 2
                id = obj.insideBrainCube;
            end;
            
            if ischar(id) && strcmpi(id, 'all')
                id = (1:numel(obj.xCube))';
            end;
            
            position = [obj.xCube(id) obj.yCube(id) obj.zCube(id)];
        end;
        function anatomicalDataForHeadGrid = getAnatomicalData(obj)
            % anatomicalDataForHeadGrid = getAnatomicalData;
            %
            % obtain anatomical label probabilities based on LPBA40 probabilistic atlas of an
            % cortical structures provided by LONI project [http://www.loni.ucla.edu/Atlases/Atlas_Detail.jsp?atlas_id=12]
            % returns a structure with two field:
            %
            %   probabilityOfEachLocationAndBrainArea: a Nx56 matrix that contains probabilities for
            %       each of the ~56 anatomical label at each location inside brain volume.
            %
            %   brainArealabel: a cell array which contains text labels for 56 anatomical labels.
                                    
            fullPath = which('pr.dipole');
            prFolderPath = fullPath(1:end - length('dipole.m'));
            precomputePath = [prFolderPath  'precompute/'];
            
            filename = [precomputePath 'anatomical_information_for_headGrid_' num2str(obj.spacing) '_mm.mat'];
            if exist(filename, 'file')
                anatomicalDataForHeadGrid = load(filename);
            else % make it, since it takes some time we will save the results in /+pr/precompute
                try
                    fprintf('Calculating anatomical label probabilities...\n');
                    [probabilityOfEachLocationAndBrainArea, brainArealabel] = label_dipoles(obj.getPosition(obj.insideBrainCube));
                    
                    anatomicalDataForHeadGrid.probabilityOfEachLocationAndBrainArea = probabilityOfEachLocationAndBrainArea;
                    anatomicalDataForHeadGrid.brainArealabel = brainArealabel;
                    
                    % save for future use                    
                    try
                        save(filename, 'probabilityOfEachLocationAndBrainArea', 'brainArealabel');
                    catch err
                        fprintf([err.message '\n']);
                        fprintf('Measure Projection: saving precomputed anatomical information was not successful.\n');
                    end;
                    
                catch err
                    fprintf([err.message '\n']);
                    error('Measure Projection: You need to install (and add to path) LONI probabilitisc atlas LPBA40 from http://www.loni.ucla.edu/Atlases/Atlas_Detail.jsp?atlas_id=12 in order to get anatomical regions.');
                end;
            end;
        end;
        function insideBrainGridLocationBrodmannAreaCount = getBrodmannData(obj)
           
            fullPath = which('pr.dipole');
            prFolderPath = fullPath(1:end - length('dipole.m'));
            precomputePath = [prFolderPath  'precompute/'];
            
            filename = [precomputePath 'brodmann_information_for_headGrid_' num2str(obj.spacing) '_mm.mat'];
            if exist(filename, 'file')
                brodmannDataForHeadGrid = load(filename);
                insideBrainGridLocationBrodmannAreaCount = brodmannDataForHeadGrid.insideBrainGridLocationBrodmannAreaCount;
            else % make it, since it takes some time we will save the results in /+pr/precompute
                try
                    fprintf('Calculating Brodmann area probabilities...\n');
                    
                    javaaddpath(which('talairach.jar'));
                    % load Talairach java object and data
                    db = org.talairach.Database;
                    db.load(which('talairach.nii'));
                    
                    maxNumberOfBrodMannAreas  = 52;
                    locationInMNI = obj.getPosition(obj.insideBrainCube);
                    locationInTalarrach = icbm_spm2tal(locationInMNI);
                    
                    % produce a pool of labels by searching spacing / 2 (e..g 4 mm for 8 m mspacing) around each location. Each location produces multiple
                    % labels.
                    insideBrainGridLocationBrodmannAreaCount = zeros(size(locationInTalarrach,1), maxNumberOfBrodMannAreas);
                    
                    for locationId =  1: size(locationInTalarrach,1)
                        labelsForLocation = db.search_range(locationInTalarrach(locationId,1), locationInTalarrach(locationId,2), locationInTalarrach(locationId,3), obj.spacing / 2);
                        
                        for i=1: length(labelsForLocation)
                            splitText = hlp_split(char(labelsForLocation(i)), ',');
                            potentialBrodmannAreaName = splitText{end};
                            
                            
                            if strfind(potentialBrodmannAreaName, 'Brodmann area')
                                areaNumber =str2num(potentialBrodmannAreaName(length('Brodmann area '):end));
                                insideBrainGridLocationBrodmannAreaCount(locationId, areaNumber) = insideBrainGridLocationBrodmannAreaCount(locationId, areaNumber) + 1;
                            end;
                        end;
                    end;
                    
                   % there should not be a normalization since the probabilities of belonging to bordmann areas do not necessarily sum tup to one. 

                    % save for future use
                 
                    
                    try
                        save(filename, 'insideBrainGridLocationBrodmannAreaCount');
                    catch err
                        fprintf([err.message '\n']);
                        fprintf('Measure Projection: saving precomputed Brodman area information was not successful.\n');
                    end;
                    
                catch err
                    
                    % to do: add to static java path here.
                    fprintf([err.message '\n']);
                    error('Measure Projection: Talairach java interface could not be initiated. It is now added to your statics java path and you need to restart Matlab to be able to use it.');
                end;
            end;
            
            
        end;
        function [labelIdCube label] = getAALData(obj);
         
            fullPath = which('pr.dipole');
            prFolderPath = fullPath(1:end - length('dipole.m'));
            precomputePath = [prFolderPath  'precompute/'];
            
            filename = [precomputePath 'AALData.mat'];
            if exist(filename, 'file')
                AALInputData = load(filename);
                
                
                brainLocation = obj.getPosition('all');
                
                xDifference = abs(repmat(vec(AALInputData.xCube(:,1,1))', [size(brainLocation,1) 1]) - repmat(brainLocation(:,1), [1 size(AALInputData.xCube,1)]));
                [dummy brainLocationXIndex] = min(xDifference, [], 2);
                
                
                yDifference = abs(repmat(vec(AALInputData.yCube(1,:,1))', [size(brainLocation,1) 1]) - repmat(brainLocation(:,2), [1 size(AALInputData.yCube,2)]));
                [dummy brainLocationYIndex] = min(yDifference, [], 2);
                
                zDifference = abs(repmat(vec(AALInputData.zCube(1,1,:))', [size(brainLocation,1) 1]) - repmat(brainLocation(:,3), [1 size(AALInputData.zCube,3)]));
                [dummy brainLocationZIndex] = min(zDifference, [], 2);
                
                brainLocationLabelId = zeros(size(brainLocation,1),1, 'uint8');
                
                for i=1:size(brainLocation, 1)
                    brainLocationLabelId(i) = round(AALInputData.atlasLabelIdCube(brainLocationXIndex(i), brainLocationYIndex(i), brainLocationZIndex(i)));
                end;
               
                labelIdCube = obj.xCube * 0;
                labelIdCube(:) = brainLocationLabelId(:);
                label = AALInputData.labelString;
            else
                fprintf('Measure Projection Toolbox: AAL data file cannot be found in %s.\n.', filename);
            end;
        end;
    end;
end