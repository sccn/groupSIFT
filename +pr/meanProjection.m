% 12/09/2019 Makoto. Line 66. If gaussianWeightMatrix is all zero, include the closest voxel. 

classdef meanProjection
    % holds information about an individual mean (across epochs) projection from (multiple) dipole space to head (brain) grid
    % locations.
    properties %(SetAccess = 'immutable')
        headGrid;
        projectionParameter % encaplustae all projection parameters, like number of Std., Gaussian Width and whether in-brain density should be normalized
    end; % properties that cannot be changed outside the condtructor
    
    properties
        similarity;
        convergence = [];
        convergenceSignificance = [];
        domain = [];% array containing domains (areas with similar projections).
        numberOfPermutations = 2500;
        % linearizedProjectedMeasure; % removed since it can consume a lot of memoery when many
        % locations are significant and measure vector is long (e.g. ERSP or ITC).
        dipoleDensity;
        domainNumberCube = [];
        significanceLevelForDomain
        similarityThresholdForDomain
    end % properties
    
    methods(Static)
        function [projectionMatrix totalDipoleDenisty gaussianWeightMatrix]= getProjectionMatrix(dipole, headGrid, projectionParameter, regionOfInterestCube)
            % gaussianWeightMatrix is dipoles x (requested) grid points. It contains dipole density
            % at each grid point for each dipole.
            
            if nargin < 4 || isempty(regionOfInterestCube)
                regionOfInterestCube = headGrid.insideBrainCube;
            end;
            
            if ischar(regionOfInterestCube) && strcmpi(regionOfInterestCube, 'all')
                regionOfInterestCube = true(headGrid.cubeSize);
            end;
            
            standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo = projectionParameter.standardDeviationOfEstimatedDipoleLocation ^ 2;
            
            % projection matrix is number of dipoles x number of grid point inside brain volume
            numberOfPointsInTheResgionOfInterest = sum(regionOfInterestCube(:));
            
            projectionMatrix = zeros(size(dipole.location , 1), numberOfPointsInTheResgionOfInterest);
            totalDipoleDenisty = zeros(1, numberOfPointsInTheResgionOfInterest);
            gaussianWeightMatrix = zeros(size(dipole.location , 1), numberOfPointsInTheResgionOfInterest);
            distanceFromDipoleToGridLocationMatrix = zeros(size(dipole.location , 1), numberOfPointsInTheResgionOfInterest);
            
            
            % a N x 3 matrix (N is the number of grid points inside brain volume
            gridPosition = [headGrid.xCube(regionOfInterestCube) headGrid.yCube(regionOfInterestCube) headGrid.zCube(regionOfInterestCube);];
            
            if projectionParameter.normalizeInBrainDipoleDenisty
                dipoleInBrainDensityNormalizationFactor = pr.meanProjection.calculateDipoleInBrainDensityNormalizationFactor(dipole, headGrid, projectionParameter);
            end;
            
            for dipoleNumber = 1:size(dipole.location , 1)
                % first place distance in the array
                distanceFromDipoleToGridLocationMatrix(dipoleNumber,:) = sum( (gridPosition - repmat(dipole.location(dipoleNumber,:), size(gridPosition,1), 1)) .^2,2 ) .^ 0.5;
                
                normalizationFactor = 1 / (projectionParameter.standardDeviationOfEstimatedDipoleLocation ^3 * sqrt(8 * (pi^3))) ; 
                gaussianWeightMatrix(dipoleNumber,:) = normalizationFactor * exp(-distanceFromDipoleToGridLocationMatrix(dipoleNumber,:).^2 / (2 * standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo));
                
                % truncate the dipole denisty Gaussian at ~3 standard deviation
                gaussianWeightMatrix(dipoleNumber, distanceFromDipoleToGridLocationMatrix(dipoleNumber,:) > (projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian * projectionParameter.standardDeviationOfEstimatedDipoleLocation) ) = 0;
                
                % If gaussianWeightMatrix is all zero, include the closest voxel. 12/09/2019 Makoto.
                if sum(gaussianWeightMatrix(dipoleNumber,:)) == 0
                    [~, minIdx] = min(distanceFromDipoleToGridLocationMatrix(dipoleNumber,:));
                    currentGWM = normalizationFactor * exp(-distanceFromDipoleToGridLocationMatrix(dipoleNumber,:).^2 / (2 * standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo));
                    gaussianWeightMatrix(dipoleNumber, minIdx) = currentGWM(minIdx);
                end
                
                % normalize the dipole in-brain denisty (make it sum up to one)
                if projectionParameter.normalizeInBrainDipoleDenisty
                    gaussianWeightMatrix(dipoleNumber,:) = gaussianWeightMatrix(dipoleNumber,:) * dipoleInBrainDensityNormalizationFactor(dipoleNumber);
                end;
            end;
            
            % normalize gaussian weights to have the sum of 1 at each grid location
            for gridId = 1:size(gaussianWeightMatrix, 2)
                totalDipoleDenisty(gridId) = sum(gaussianWeightMatrix(:, gridId));
                if totalDipoleDenisty(gridId) > 0
                    projectionMatrix(:, gridId) = gaussianWeightMatrix(:, gridId) / totalDipoleDenisty(gridId);
                end;
            end;
        end
        
        function [projectionMatrix totalDipoleDenisty gaussianWeightMatrix]= getProjectionMatrixForArbitraryLocation(dipole, projectionParameter, location, headGrid, varargin)
            % [projectionMatrix totalDipoleDenisty gaussianWeightMatrix]= getProjectionMatrixForArbitraryLocation(dipole, projectionParameter, location, headGrid)
            %
            % location is an N x 3 matrix (N is the number of points to be projected into.
            % gaussianWeightMatrix is dipoles x (requested) grid points. It contains dipole density
            % at each grid point for each dipole.
            
            if nargin < 4
                headGrid = pr.headGrid(8); % by default use 8 mm headgrid to calculate dipole normalization factors.
            end;
            
            standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo = projectionParameter.standardDeviationOfEstimatedDipoleLocation ^ 2;
            
            % projection matrix is number of dipoles x number of grid point inside brain volume
            numberOfPoints = size(location, 1);
            
            projectionMatrix = zeros(size(dipole.location , 1), numberOfPoints);
            totalDipoleDenisty = zeros(1, numberOfPoints);
            gaussianWeightMatrix = zeros(size(dipole.location , 1), numberOfPoints);
            distanceFromDipoleToGridLocationMatrix = zeros(size(dipole.location , 1), numberOfPoints);
            
            if projectionParameter.normalizeInBrainDipoleDenisty
                dipoleInBrainDensityNormalizationFactor = pr.meanProjection.calculateDipoleInBrainDensityNormalizationFactor(dipole, headGrid, projectionParameter);
            end;
            
            for dipoleNumber = 1:size(dipole.location , 1)
                % first place distance in the array
                distanceFromDipoleToGridLocationMatrix(dipoleNumber,:) = sum( (location - repmat(dipole.location(dipoleNumber,:), size(location,1), 1)) .^2,2 ) .^ 0.5;
                                
                normalizationFactor = 1 / (projectionParameter.standardDeviationOfEstimatedDipoleLocation ^3 * sqrt(8 * (pi^3))) ; 
                gaussianWeightMatrix(dipoleNumber,:) = normalizationFactor * exp(-distanceFromDipoleToGridLocationMatrix(dipoleNumber,:).^2 / (2 * standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo));
                
                % truncate the dipole denisty Gaussian at ~3 standard deviation
                gaussianWeightMatrix(dipoleNumber, distanceFromDipoleToGridLocationMatrix(dipoleNumber,:) > (projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian * projectionParameter.standardDeviationOfEstimatedDipoleLocation) ) = 0;
                
                % normalize the dipole in-brain density (make it sum up to one)
                if projectionParameter.normalizeInBrainDipoleDenisty
                    gaussianWeightMatrix(dipoleNumber,:) = gaussianWeightMatrix(dipoleNumber,:) * dipoleInBrainDensityNormalizationFactor(dipoleNumber);
                end;
            end;
            
            % normalize gaussian weights to have the sum of 1 at each grid location
            for gridId = 1:size(gaussianWeightMatrix, 2)
                totalDipoleDenisty(gridId) = sum(gaussianWeightMatrix(:, gridId));
                if totalDipoleDenisty(gridId) > 0
                    projectionMatrix(:, gridId) = gaussianWeightMatrix(:, gridId) / totalDipoleDenisty(gridId);
                end;
            end;
        end
        
        function dipoleInBrainDensityNormalizationFactor = calculateDipoleInBrainDensityNormalizationFactor(dipole, headGrid, projectionParameter)
            % calculate the factors (scalar values, each for a dipole) that normalize the dipole
            % projected density inside brain volume (makes its sum to be equal to one).
            
            % create a parameter set but withough normalization
            newProjectionParemeter = pr.projectionParameter(projectionParameter.standardDeviationOfEstimatedDipoleLocation, projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian, false);
            
            
            % get each dipole density from gaussianWeightMatrix
            [projectionMatrix totalDipoleDenisty gaussianWeightMatrix]= pr.meanProjection.getProjectionMatrix(dipole, headGrid, newProjectionParemeter, headGrid.insideBrainCube);
            
            % calculate its sum inside brain volume, 1/this gives the factor which should be
            % multiplied by dipole denisty in projection (because we assume dipole should be
            % somewhere in the brain volume).
            dipoleInBrainDensityNormalizationFactor = repmat(1, size(gaussianWeightMatrix,1),1) ./ sum(gaussianWeightMatrix, 2);
        end
        
        function [sumWeightedSimilarityPvalue sumWeightedSimilarity] = calculateMeasureConvergenceSignificance(similarity, dipole, headGrid, projectionParameter, numberOfPermutations)
            
            if nargin < 5
                numberOfPermutations = 2000;
            end;            
            
            % if it is not float type, like uint32, then 1/numberOfPermutations becomes zero, so we
            % make sure here that it is float (double)
            numberOfPermutations = double(numberOfPermutations);
            
            sumWeightedSimilarity = headGrid.xCube * 0;
            
            standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo = projectionParameter.standardDeviationOfEstimatedDipoleLocation ^ 2;
            
            similarity = similarity - diag(diag(similarity)); % remove ones on the diagonal
            dipoleLocation = dipole.location;
            
            surrogateSumWeightedSimilarity = zeros(1, numberOfPermutations);
            sumWeightedSimilarityPvalue = ones(size(headGrid.xCube));
            
            randomPermutation = zeros(numberOfPermutations, size(similarity,1));
            
            % To make p values deterministic (not change based on random values generated)
            % we always set the seed to zero.
            currentRandomStream = RandStream('mt19937ar','Seed',0);
            
            for permutationNumber = 1:numberOfPermutations
                % randomPermutation(permutationNumber,:) = randperm(size(similarity,1)); % without repetition
                randomPermutation(permutationNumber,:) = currentRandomStream.randi(size(similarity,1), 1, size(similarity,1)); % with repetition
            end;
            
            
            pr.progress('init'); % start the text based progress bar
            
            if projectionParameter.normalizeInBrainDipoleDenisty
                dipoleInBrainDensityNormalizationFactor = pr.meanProjection.calculateDipoleInBrainDensityNormalizationFactor(dipole, headGrid, projectionParameter);
            end;
            
            for i=1:numel(headGrid.xCube)
                if headGrid.insideBrainCube(i)
                    
                    if mod(i,10) ==0
                        %fprintf('Percent done = %d\n', round(100 * i / numel(headGrid.xCube)));
                        pr.progress(i / numel(headGrid.xCube), sprintf('\npercent done %d/100',round(100*i / numel(headGrid.xCube))));
                    end;
                    
                    pos = [headGrid.xCube(i) headGrid.yCube(i) headGrid.zCube(i)];
                    distanceToDipoles = sum((dipoleLocation - repmat(pos, size(dipoleLocation,1), 1))' .^2) .^ 0.5;
                    
                    % truncate the dipole denisty Gaussian at ~3 standard deviation
                    dipoleWithNonzeroWeightIds = distanceToDipoles < (projectionParameter.standardDeviationOfEstimatedDipoleLocation * projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian);
                    
                    if ~isempty(dipoleWithNonzeroWeightIds)
                        
                        % pass distance to dipoles through a gaussian kernel with specified standard deviation.
                        gaussianPassedDistanceToDipoles = sqrt(1/(2 * pi * standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo)) * exp(-distanceToDipoles(dipoleWithNonzeroWeightIds).^2 / (2 * standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo));
                        
                        if projectionParameter.normalizeInBrainDipoleDenisty
                            % make dipole denisty sum inside brain volume equal to 1
                            gaussianPassedDistanceToDipoles = gaussianPassedDistanceToDipoles .* dipoleInBrainDensityNormalizationFactor(dipoleWithNonzeroWeightIds)';
                        end;
                        
                        gaussianWeightMatrix = repmat(gaussianPassedDistanceToDipoles,length(gaussianPassedDistanceToDipoles),1);
                        
                        normalizationMatrix = gaussianWeightMatrix .* gaussianWeightMatrix';
                        sumPairwiseWeights = sum(normalizationMatrix(:));
                        
                        normalizationMatrix = normalizationMatrix / sumPairwiseWeights;
                        
                        similarityWeightedByGauissian = similarity(dipoleWithNonzeroWeightIds, dipoleWithNonzeroWeightIds) .* normalizationMatrix;
                        sumWeightedSimilarity(i) = sum(similarityWeightedByGauissian(:));
                        
                        % bootstapping with permutation
                        surrogateSumWeightedSimilarity = zeros(1, numberOfPermutations);
                        for permutationNumber = 1:numberOfPermutations
                            surrogateSimilarityWeightedByGauissian = similarity(randomPermutation(permutationNumber, dipoleWithNonzeroWeightIds), randomPermutation(permutationNumber, dipoleWithNonzeroWeightIds)) .* normalizationMatrix;
                            surrogateSumWeightedSimilarity(permutationNumber) = sum(surrogateSimilarityWeightedByGauissian(:));                                                        
                        end;
                        
                        sumWeightedSimilarityPvalue(i) = sum(surrogateSumWeightedSimilarity >= sumWeightedSimilarity(i)) / numberOfPermutations;
                    end;
                end;
            end;
            
            pause(.1);
            pr.progress('close'); % duo to some bug need a pause() before
            fprintf('\n');            
                        
            % we can only be sure that the significance is less than 1/numberOfPermutations
            sumWeightedSimilarityPvalue = max(sumWeightedSimilarityPvalue, (1 / double(numberOfPermutations)) - eps);
        end
    end;
    methods (Access = 'protected') % may cause version compatibility issue with older Matlab versions..?
        function plotColoredByMeasure(obj, plotType, significanceLevel, dipoleAndMeasure, varargin)
            % plotColoredByMeasure(significanceLevel, dipoleAndMeasure, plotType, 'key', 'value'..)
            % plotType has to be either 'voxel' or 'volume'
            %
            %    'colormap': ['hsv' or 'ycrcb'] default: 'hsv'
            %    'whiteness': 0<real<1,  default: 0.45  (only used for 'ycrcb' color space)
            %    'projection' ['on' or 'off'] deafult: 'off'. show projction of voxels on side MRIs
            %
            
            
            inputOptions = finputcheck(varargin, ...
                {'colormap'         'string'  {'hsv' 'ycrcb'} 'hsv'; ...
                'whiteness'         'real'    [0 1] 0.45;...
                'projection'        'boolean' [] true;...
                'colorbar'          'boolean' []  true;...
                'newFigure'        'boolean'  []   true;...
                });            
            
            if nargin < 3
                error('Measure Projection: please provide a variable of type dipoleAndMeasure in the second argument.');
            end;
            
			if nargin<3 % if significanceLevel has not been provided as an input.
				if isempty(obj.significanceLevelForDomain)
					significanceLevel = 0.002;
				else
					significanceLevel = obj.significanceLevelForDomain;
				end;
			end;
            
            significantId = obj.convergenceSignificance <= significanceLevel;
            
            if ~any(significantId(:))
                fprintf('Measure Projection: there are no significant portions to be plotted.\nConsider increasing the significance value (or disabling FDR correction).\n');
                return;
            end;
            
            [projectionMatrix totalDipoleDenisty]= obj.getProjectionMatrix(dipoleAndMeasure, obj.headGrid, obj.projectionParameter, significantId);
            linearizedProjectedMeasure = dipoleAndMeasure.linearizedMeasure * projectionMatrix;
            
            % create a 1-D multi-dimensional scaling (MDS) from projected similarities.
            projectedMeasureSimilarity = squareform(pdist(linearizedProjectedMeasure','correlation'));
            
            % colorbased on a multi-dimensional scaling of the projected measure
            
            
            if strcmpi(inputOptions.colormap, 'ycrcb')
                % use YCrCb color mapping
                [pos stress]= robust_mdscale(projectedMeasureSimilarity, 2);
                
                pos = pos /  max(range(pos));
                pos(:,1) = pos(:,1) - min(pos(:,1));
                pos(:,2) = pos(:,2) - min(pos(:,2));
                Y = ones(size(pos,1), 1) * inputOptions.whiteness; % Y value in ycbcr, related to luminescence (between 0 for dark and 1 for all white)
                
                voxelColor = ycbcr2rgb([Y pos]);
            else % hsv color mapping
                [pos stress]= robust_mdscale(projectedMeasureSimilarity, 1);
                
                % map the 1-D representation of projected measures into color
                hsvFromRedToBlue = hsv(64);
                hsvFromRedToBlue = hsvFromRedToBlue(1:44,:);
                
                voxelColor = value2color(pos, hsvFromRedToBlue);
            end;
            
            fprintf('Multi-Dimensional Scaling stress (deformation) for representing projected measures by color = %%%d\n', round(100 * stress));
            
            if inputOptions.newFigure
                figure;
                set(gcf, 'name', ['P <=' num2str(significanceLevel)]);
            end;
                        
            plot_dipplot_with_cortex;
            
            if strcmp(plotType , 'voxel') % colored voxel plot
                idNumbers = find(significantId);
                zeroCube  = false(obj.headGrid.cubeSize);
                for i=1:length(pos)
                    membershipCube = zeroCube;
                    membershipCube(idNumbers(i)) = 1;
                    pr.plot_head_region(obj.headGrid, membershipCube, 'regionColor', voxelColor(i,:), 'showProjectedOnMrs', inputOptions.projection, 'projectionAlpha', 0.1);
                end;
            else   % colored volume (surafce) plot
                pr.plot_head_surface(obj.headGrid, significantId, 'surfaceColor', voxelColor, 'showProjectedOnMrs', inputOptions.projection, 'projectionAlpha', 0.1, 'reductionFactor', 1/2);
            end;
            
            % add a colobar
            if inputOptions.colorbar
                colormap(hsv);
                handle = cbar('vert', 1:44);
                
                % remove on y axis text and make the colorbar shorter.
                set(handle, 'ytick', []);
                set(handle, 'position', [ 0.9250    0.1100    0.0310    0.4])
            end;
        end;
        
        function plotColoredBySignificance(obj, plotType, varargin)
            % plotColoredBySignificance(plotType, significanceLevel, varargin)
            % plotType has to be either 'voxel' or 'volume'
            
            inputOptions = finputcheck(varargin, ...
                {'colormap'         'string'  {'hot' 'jet'} 'hot'; ...                
                'projection'        'boolean' [] true;...
                'colorbar'          'boolean' []  true;...
                'newFigure'        'boolean'  []   true;...
                });
            
            
            switch inputOptions.colormap
                case 'hot'
                    inputColormap = hot();
                case  'jet'
                    inputColormap = jet();
            end;
            
            if inputOptions.newFigure
                figure;
                plot_dipplot_with_cortex;
            end;
            
            significantId = obj.convergenceSignificance < 0.05;
            statisticalPower = -log(obj.convergenceSignificance(significantId));
            voxelColor = value2color(statisticalPower, inputColormap);                
            
            inputOptions.projection = true;
            inputOptions.colorbar = true;
            
            if strcmp(plotType , 'voxel') % colored voxel plot
                idNumbers = find(significantId);
                zeroCube  = false(obj.headGrid.cubeSize);
                for i=1:length(idNumbers)
                    membershipCube = zeroCube;
                    membershipCube(idNumbers(i)) = 1;
                    pr.plot_head_region(obj.headGrid, membershipCube, 'regionColor', voxelColor(i,:), 'showProjectedOnMrs', inputOptions.projection, 'projectionAlpha', 0.1);
                end;
            else   % colored volume (surafce) plot
                pr.plot_head_surface(obj.headGrid, significantId, 'surfaceColor', voxelColor, 'showProjectedOnMrs', inputOptions.projection, 'projectionAlpha', 0.1, 'reductionFactor', 1/2);
            end;
            
            % add a colobar
            if inputOptions.colorbar
                colormap(inputColormap);
                handle = cbar('vert', 1:size(inputColormap,1), [0.05 min(obj.convergenceSignificance(significantId(:)))]);
                
                cbarLabelCharArray = get(handle, 'yticklabel');
                cbarLabel = {};
                for i=1:size(cbarLabelCharArray,1)
                    cbarLabel{i} = strtrim(cbarLabelCharArray(i,:));
                end;
                
                % we need fliplr since the order of values in cbar label is actually in decreasing order
                actualPvalues = fliplr(logspace(log10(min(obj.convergenceSignificance(significantId(:)))), log10(0.05), length(cbarLabel)));
                for i=2:length(cbarLabel)
                    cbarLabel{i} = num2str(actualPvalues(i),  '%2.1e');
                end;
                
                set(handle, 'yticklabel',cbarLabel )
                set(handle, 'ycolor', [1 1 1])
                % remove on y axis text and make the colorbar shorter.
                % set(handle, 'ytick', []);
                set(handle, 'position', [ 0.8800    0.1100    0.0310    0.4])
            end;            
        end;
        
        function domainColor = makeDomainColor(obj, inputOptions)
            % assign colors to Domains based on their similarities or other parameters.
            if isempty(obj.domain)
                error('No domains are present. Please first create domain(s) with createDomain() method.');
            end;
            
            % if only one domain exist, just plot it (no need for different colors with mds()...)
            if numel(obj.domain) == 1
                domainColor = [0.1 0.7 0.1];
                return;
            end;
            
            for i=1:length(obj.domain)
                linearizedProjectedMeasure(i,:) =  obj.domain(i).exemplarLinearMeasure;
            end;
            % color based on maximally distinct colors.
            if inputOptions.maximallyDistinctColors
                domainColor = pr.get_distinct_colors(numel(obj.domain));
            else
                
                % create a 1-D multi-dimensional scaling (MDS) from projected similarities.
                projectedMeasureDissimilarity = squareform(pdist(linearizedProjectedMeasure,'correlation'));
                
                if strcmpi(inputOptions.colormap, 'ycrcb')
                    % use YCrCb color mapping
                    [pos stress] = robust_mdscale(projectedMeasureDissimilarity, 2);
                    
                    pos = pos /  max(range(pos));
                    pos(:,1) = pos(:,1) - min(pos(:,1));
                    pos(:,2) = pos(:,2) - min(pos(:,2));
                    Y = ones(size(pos,1), 1) * inputOptions.whiteness; % Y value in ycbcr, related to luminescence (between 0 for dark and 1 for all white)
                    
                    domainColor = ycbcr2rgb([Y pos]);
                elseif strcmpi(inputOptions.colormap, 'hsv')% hsv color mapping
                    
                    [pos stress] = robust_mdscale(projectedMeasureDissimilarity, 1);
                    
                    % make sure that domain colors are not exactly the same or very similar to each
                    % other.
                    for i=1:100
                        d = squareform(pdist(pos));
                        for j=1:size(d,1)
                            for k=1:(j-1)
                                if d(j,k) < 0.15
                                    middle = (pos(j) + pos(k)) / 2;
                                    pos(j) = middle - 0.05;
                                    pos(k) = middle + 0.05;
                                end;
                            end;
                        end;
                    end;
                    
                    % try to normalize the colors (1-D MDS) so they are more or less the same across
                    % different runs of the function. This is done by keeping the mean after
                    % normalization less than 0.5 (it is not perfect though).
                    normalizedPos = pos - min(pos);
                    normalizedPos = normalizedPos / max(normalizedPos);
                    if mean(normalizedPos) > 0.5
                        pos = -pos;
                    end;
                    
                    % map the 1-D representation of projected measures into color
                    hsvFromRedToBlue = hsv(64);
                    hsvFromRedToBlue = hsvFromRedToBlue(1:44,:);
                    
                    domainColor = value2color(pos, hsvFromRedToBlue);
                end;
                
                fprintf('Multi-Dimensional Scaling stress (deformation) for representing domain projected measure exemplars by color = %%%d\n', round(100 * stress));
            end;
        end;
    end
    
    methods
        function obj = meanProjection(dipoleAndMeasure, similarity, headGrid, varargin)
            % meanProjection(dipoleAndMeasure, similarity, headGrid, [optional key,values])
            % options:
            %
            % stdOfDipoleGaussian {default = 12 mm}
            % numberOfPermutations
            % numberOfStdsToTruncateGaussian {default = 3}
            % normalizeInBrainDipoleDenisty {string = 'on' or 'off'}
            
            inputOptions = finputcheck(varargin, ...
                {'stdOfDipoleGaussian'         'real'  [] 12; ...
                'numberOfPermutations'        'real'    [] obj.numberOfPermutations;...
                'numberOfStdsToTruncateGaussian'   'real'    []  3;...
                'normalizeInBrainDipoleDenisty'   'string'    {'on', 'off'}  'on';...
                });
            
            if nargin < 3
                obj.headGrid = pr.headGrid;
            end;
            
            obj.numberOfPermutations = double(inputOptions.numberOfPermutations);
            obj.projectionParameter = pr.projectionParameter(inputOptions.stdOfDipoleGaussian, inputOptions.numberOfStdsToTruncateGaussian, strcmpi(inputOptions.normalizeInBrainDipoleDenisty, 'on'));
            
            obj.headGrid = headGrid;
            obj.similarity = similarity;
            
            [obj.convergenceSignificance obj.convergence] = pr.meanProjection.calculateMeasureConvergenceSignificance(obj.similarity, dipoleAndMeasure, obj.headGrid, obj.projectionParameter , obj.numberOfPermutations);
        end;
        
        function plotScatter(obj, significanceLevel)
            
            if nargin<2
                significanceLevel = 0.002;
            end;
            
            figure;
            
            t = obj.convergence - min(obj.convergence(:)) + 0.001;
            t = t / max(t(:));
            t(~obj.headGrid.insideBrainCube) = eps;
            
            t(obj.convergenceSignificance > significanceLevel) = eps;
            
            color = value2color(t(:), jet);
            
            scatter3(obj.headGrid.xCube(obj.headGrid.insideBrainCube), obj.headGrid.yCube(obj.headGrid.insideBrainCube), obj.headGrid.zCube(obj.headGrid.insideBrainCube), (t(obj.headGrid.insideBrainCube) * 40)+5, color(obj.headGrid.insideBrainCube), 'filled');
            axis equal;
            cbar;
        end;
        
        function [valueOnFineGridOutput mri] = plotMri(obj, significanceLevel, mri3dplotOptions)
            
            if nargin<2
                significanceLevel = 0.002;
            end;
            
            if nargin<3
                mri3dplotOptions = {'mriview' , 'top','mrislices', [-50 -30 -20 -15 -10 -5 0 5 10 15 20 25 30 40 50]};
            end;
            
            valueOnCoarseGrid = obj.convergence;
            valueOnCoarseGrid(obj.convergenceSignificance > significanceLevel) = 0;
            
            valueOnCoarseGrid(valueOnCoarseGrid<0) = 0;
            
            [valueOnFineGrid mri] = convert_coarse_grid_to_mri3d(valueOnCoarseGrid, obj.headGrid.xCube, obj.headGrid.yCube, obj.headGrid.zCube);
            mri3dplot(valueOnFineGrid, mri, mri3dplotOptions{:}); % for some reason, this function hicjacks a currently open figure. even if a new figur is just created.
            
            if nargout > 0
                valueOnFineGridOutput = valueOnFineGrid;
            end;
        end;
        
        function [obj position]= createDomain(obj, dipoleAndMeasure, similarityThreshold, significanceLevel, minNumberOfClusters, outlierThreshold, varargin)
            % [obj position]= createDomain(obj, dipoleAndMeasure, similarityThreshold, significanceLevel)
            
            if nargin<3
                similarityThreshold = 0.8;
            end;
            
            if nargin<4
                significanceLevel = 0.002;
            end;
            
            if nargin<5
                minNumberOfClusters = 2;
            end;
            if nargin<6
                outlierThreshold = 0.6; % or 0.65
            end;
            
            inputOptions = finputcheck(varargin, ...
                {'plotMds'    'boolean'  [] false;... % whether to show Multi-dimensional scaling of projected domain values.
                'topologicalSegmentation' 'boolean' [] true;...
                'sortDomainBy'  'string' {'dipolemass' 'volume'  'sessionentropy' 'subjectentropy'} 'dipolemass';...
                'minDomainSize'  'integer' [1 Inf] 4;... % smalled domain size permitted, smaller domains will be removed.
                'clusteringMethod' 'string' {'affinity' 'connectedComponent'} 'affinity';...
                });
            
            obj.significanceLevelForDomain = significanceLevel;
            
            significantId = obj.convergenceSignificance <= significanceLevel;
            
            if any(significantId(:))
                
                %position = [obj.headGrid.xCube(significantId) obj.headGrid.yCube(significantId) obj.headGrid.zCube(significantId)];
                %[obj.linearizedProjectedMeasure obj.dipoleDensity] = dipoleAndMeasure.projectTo(position, obj.standardDeviationOfEstimatedDipoleLocation);
                
                [projectionMatrix obj.dipoleDensity]= obj.getProjectionMatrix(dipoleAndMeasure, obj.headGrid, obj.projectionParameter, significantId);
                linearizedProjectedMeasure = dipoleAndMeasure.linearizedMeasure * projectionMatrix;
                
                
				if size(linearizedProjectedMeasure,1) == 1 % if the measure is just a sungle scalar value
					projectedMeasureSimilarity = -squareform(pdist(linearizedProjectedMeasure,'euclidean'));
				else % the measure is a vector
					projectedMeasureSimilarity = 1-squareform(pdist(linearizedProjectedMeasure','correlation'));
				end;
				
                % use the signed mutual infromation instead of correlation as a similarity measure.
                % projectedMeasureSimilarity = pr.estimate_mutual_information_from_correlation(projectedMeasureSimilarity);
                % similarityThreshold = pr.estimate_mutual_information_from_correlation(similarityThreshold);
                
                fprintf('Mastering your domains...\n');
                if strcmpi(inputOptions.clusteringMethod, 'connectedComponent')
                    clusterInfo.IDX = pr.connected_components(projectedMeasureSimilarity > similarityThreshold);
                    clusterInfo.numberOfOutliers = 0;
                    clusterInfo.numberOfClusters = length(unique(clusterInfo.IDX));
                    
                    for domainNumber = 1:clusterInfo.numberOfClusters% for test only , we should fix this
                        
                        % find the point that has the maximum sum of similaritie to other points.
                        similarity = projectedMeasureSimilarity;
                        similarity = similarity - diag(diag(similarity));
                        similarity(clusterInfo.IDX ~= domainNumber, clusterInfo.IDX ~= domainNumber) = 0;
                        [dummy maxId] = max(sum(similarity));
                        
                        clusterInfo.examplarId(domainNumber) = maxId;
                    end;
                else
                    clusterInfo = cluster_based_on_similarities(projectedMeasureSimilarity, outlierThreshold, similarityThreshold, minNumberOfClusters, inputOptions.plotMds);
                end;
                
                obj.domainNumberCube = obj.headGrid.xCube * 0;
                obj.domainNumberCube(significantId) = clusterInfo.IDX;
                
                obj.domain = pr.domain;
                domainLabel = {};
                
                if clusterInfo.numberOfOutliers > 0
                    numberOfNonOutlierClusters = clusterInfo.numberOfClusters - 1;
                else
                    numberOfNonOutlierClusters = clusterInfo.numberOfClusters;
                end;
                
                for domainNumber = 1:numberOfNonOutlierClusters % one cluster is for the outliers
                    significantIndiceNumbers = find(significantId);
                    exemplarIdInHeadGrid = significantIndiceNumbers(clusterInfo.examplarId(domainNumber));
                    meanLinearizedProjectedMeasure = linearizedProjectedMeasure(:, clusterInfo.IDX == domainNumber) * obj.dipoleDensity(:, clusterInfo.IDX == domainNumber)';
                    obj.domain(domainNumber) = pr.domain(obj.domainNumberCube == domainNumber, obj.headGrid, meanLinearizedProjectedMeasure, exemplarIdInHeadGrid, linearizedProjectedMeasure(:,clusterInfo.examplarId(domainNumber)), ['Domain ' num2str(domainNumber)], dipoleAndMeasure, obj.projectionParameter);
                    domainLabel{domainNumber} = ['Domain ' num2str(domainNumber)];
                end;
                
                % place legend labels on the plot created by cluster_based_on_similarities()
                % add an Outlier label if needed.
%                 if clusterInfo.numberOfOutliers > 0
%                     legend(['Outlier' domainLabel]);
%                 else
%                     legend(domainLabel);
%                 end;
            else
                fprintf('There is no location with specified significance to create domain from. \nYou may choose a less strict significance criteria (higher p value) and try again.\n');
            end;
            
            
            % separate non-contigious segments of domain into diffreent domains (Topological Segmentation)
            if inputOptions.topologicalSegmentation
                numberOfDomainBeforeSegmentation = length(obj.domain);
                obj = obj.performTopographicSegmentationOnAllDomains(dipoleAndMeasure);
                fprintf('Topological segmentation created %d domains from the original %d.\n',  length(obj.domain), numberOfDomainBeforeSegmentation);
            end;
            
            % remove very small domains.
            obj = obj.removeVerySmallDomains(inputOptions.minDomainSize);
            
            % sort domains by dipole mass or other quantities.
            obj = obj.sortDomain(inputOptions.sortDomainBy);
            
            % assign different colors to domains.
            obj = obj.assignDomainColor;
            
            % label the domains.
            for domainNumber = 1:numel(obj.domain)
                obj.domain(domainNumber).label = ['Domain ' num2str(domainNumber)];
            end;
            
        end;
        
        function plotDomainScatter(obj, domainNumber)
            
            if nargin<2
                domainNumber = 1:length(obj.domain);
            end;
            
            figure;
            domainId = ismember(obj.domainNumberCube, domainNumber);
            pointSize = ones(size(obj.headGrid.xCube));
            pointSize(domainId) = 50;
            
            valueforColor = pointSize * 0;
            valueforColor(domainId) = obj.domainNumberCube(domainId);
            
            color = value2color(valueforColor(:), hsv(length(domainNumber)+1));
            
            % set the color of points that belong to no domains as black
            color(valueforColor(:) == 0,:) = 0;
            
            scatter3(obj.headGrid.xCube(obj.headGrid.insideBrainCube), obj.headGrid.yCube(obj.headGrid.insideBrainCube), obj.headGrid.zCube(obj.headGrid.insideBrainCube), pointSize(obj.headGrid.insideBrainCube) +5, color(obj.headGrid.insideBrainCube,:), 'filled');
            axis equal;
        end;
        
        function difference = findGroupDifference(obj, groupNameOrId, dipoleAndMeasure, similarity, headGrid, standardDeviationOfEstimatedDipoleLocation, numberOfPermutations)
            % groupId can be only either one(1) or two(2)
            
            groupId1 = ismemer(dipoleAndMeasure.groupName, groupNameOrId{1});
            groupId2 = ismemer(dipoleAndMeasure.groupName, groupNameOrId{2});
            
            groupId = groupId1 + 2 * groupId2;
            
            if nargin<7
                numberOfPermutations = 300;
            end;
            
            dipoleLocation = dipoleAndMeasure.location;
            
            standardDeviationOfEstimatedDipoleLocationPowerTwo = standardDeviationOfEstimatedDipoleLocation ^ 2;
            
            group1SignificanceInNeighborhood = ones(headGrid.cubeSize);
            group2SignificanceInNeighborhood = ones(headGrid.cubeSize);
            minGroupSignificanceInNeighborhood = ones(headGrid.cubeSize);
            
            groupDensityDifferenceSignificance =  ones(headGrid.cubeSize);
            normalizedGroup1Density = ones(headGrid.cubeSize);
            normalizedGroup2Density = ones(headGrid.cubeSize);
            
            surrogateGroupIds = create_surrogate_group_ids(groupId, numberOfPermutations, subjectId);
            
            for i=1:numel(gridX)
                if insideBrain(i) %% & sumWeightedSimilarityPvalue(i) < 0.02
                    if mod(i,10) ==0
                        fprintf('Percent done = %d\n', round(100 * i / numel(gridX)));
                    end;
                    
                    pos = [gridX(i) gridY(i) gridZ(i)];
                    distanceToDipoles = sum((dipoleLocation - repmat(pos, size(dipoleLocation,1), 1))' .^2) .^ 0.5;
                    
                    % pass distance to dipoles through a gaussian kernel with specified standard deviation.
                    gaussianPassedDistanceToDipoles = sqrt(1/(2 * pi * standardDeviationOfEstimatedDipoleLocationPowerTwo)) * exp(-distanceToDipoles.^2 / (2 * standardDeviationOfEstimatedDipoleLocationPowerTwo));
                    
                    gaussianWeightMatrix = repmat(gaussianPassedDistanceToDipoles,length(gaussianPassedDistanceToDipoles),1);
                    
                    % a matrix that weights each pair of dipoles according to their gaussians
                    neighborhoodWeightMatrix = gaussianWeightMatrix .* gaussianWeightMatrix';
                    
                    groupSignificanceInNeighborhood = calculate_group_significance_based_on_pairwise_similarity(groupId, similarity, surrogateGroupIds, neighborhoodWeightMatrix);
                    group1SignificanceInNeighborhood(i) = groupSignificanceInNeighborhood(1);
                    group2SignificanceInNeighborhood(i) = groupSignificanceInNeighborhood(2);
                    minGroupSignificanceInNeighborhood(i) = min(groupSignificanceInNeighborhood);
                    
                    [groupDensityDifferenceSignificance(i) normalizedGroup1Densities]= calculate_group_density_difference_significance(groupId, gaussianPassedDistanceToDipoles, surrogateGroupIds);
                    normalizedGroup1Density(i) = normalizedGroup1Densities(1);
                    normalizedGroup2Density(i) = normalizedGroup1Densities(2);
                end;
            end;
            
            
            difference.numberOfPermutations = numberOfPermutations;
            difference.group1SignificanceInNeighborhood = group1SignificanceInNeighborhood;
            difference.group2SignificanceInNeighborhood = group2SignificanceInNeighborhood;
            difference.minGroupSignificanceInNeighborhood = minGroupSignificanceInNeighborhood;
            difference.groupDensityDifferenceSignificance = groupDensityDifferenceSignificance;
            difference.normalizedGroup1Density = normalizedGroup1Density;
            difference.normalizedGroup2Density = normalizedGroup2Density;
            
        end
        
        function plotVoxel(obj, significanceLevel, varargin)
            % plots grid locations, as cubic voxels, with significant convergence.
            
            if nargin<2
                if isempty(obj.significanceLevelForDomain)
                    significanceLevel = 0.002;
                else
                    significanceLevel = obj.significanceLevelForDomain;
                end;
            end;
            
            inputOptions = finputcheck(varargin, ...
                {'projection'    'boolean'  [] true;...
                'newFigure'        'boolean'  []   true;...
                });
            
            if inputOptions.newFigure
                figure;
                set(gcf, 'name', ['P <=' num2str(significanceLevel)]);
            end;
            
            plot_dipplot_with_cortex;
            
            membershipCube = obj.convergenceSignificance <= significanceLevel;
            
            if ~any(membershipCube(:))
                fprintf('Measure Projection: there are no significant portions to be plotted.\nYou may consider increasing the significance value (or disabling FDR correction).\n');
                return;
            end;
            
            pr.plot_head_region(obj.headGrid, membershipCube);
        end;
        
        function plotVolume(obj, significanceLevel, varargin)
            % plots the volume surface with significant convergence.
            
            if nargin<2
                if isempty(obj.significanceLevelForDomain)
                    significanceLevel = 0.002;
                else
                    significanceLevel = obj.significanceLevelForDomain;
                end;
            end;
            
            inputOptions = finputcheck(varargin, ...
                {'projection'    'boolean'  [] true;...
                'surfaceOptions' 'cell' {} {} ;... % options to be aplied to crerated patch objects.
                'plotOptions' 'cell' {} {} ;... % options directly passed to plot_head_surface()
                'newFigure'        'boolean'  []   true;...
                });
            
            if inputOptions.newFigure
                figure;
                set(gcf, 'name', ['P <=' num2str(significanceLevel)]);
            end;
            
            plot_dipplot_with_cortex;
            
            membershipCube = obj.convergenceSignificance <= significanceLevel;
            
            if ~any(membershipCube(:))
                fprintf('Measure Projection: there are no significant portions to be plotted.\nYou may consider increasing the significance value (or disabling FDR correction).\n');
                return;
            end;
            
            pr.plot_head_surface(obj.headGrid, membershipCube,  'showProjectedOnMrs', inputOptions.projection, 'surfaceOptions', inputOptions.surfaceOptions, inputOptions.plotOptions{:});
        end;
        
        function plotCortex(obj, significanceLevel, varargin)
            % plotCortex(significanceLevel, {key, value pairs})
            % plot significant locations on MNI example cortical surface.
            
            if nargin<2
                if isempty(obj.significanceLevelForDomain)
                    significanceLevel = 0.002;
                else
                    significanceLevel = obj.significanceLevelForDomain;
                end;
            end;
            
            inputOptions = finputcheck(varargin, ...
                {'color'    'real'  [] [0 1 1];...
                'newFigure'          'boolean' []  true;...
                });
            
            membershipCube = obj.convergenceSignificance <= significanceLevel;
            significantLocation = obj.headGrid.getPosition(membershipCube);
            
            fsf = load('MNImesh_dipfit.mat');
            cortexVertices = fsf.vertices;
            
            headGridSpacing =  obj.headGrid.spacing;
            
            cortexPointDomainDenisty = pr.project_domain_to_cortex(significantLocation, cortexVertices, headGridSpacing);
            
            pr.plot_cortex(cortexPointDomainDenisty, inputOptions.color, 'newFigure', inputOptions.newFigure,  'minValueNormalizedTo' , 0.3);
            
            set(gcf, 'name', ['P <=' num2str(significanceLevel)]);
        end;
        
        function plotVoxelColoredByMeasure(obj,  varargin)
            plotColoredByMeasure(obj, 'voxel', varargin{:})
        end;
        
        function plotVolumeColoredByMeasure(obj,  varargin)
            plotColoredByMeasure(obj, 'volume', varargin{:})
        end;
        
        function plotVoxelColoredByDomain(obj, varargin)
            % plotVoxelColoredByDomain(dipoleAndMeasure, significanceLevel, 'key', 'value'..)
            
            inputOptions = finputcheck(varargin, ...
                {'colormap'         'string'  {'hsv' 'ycrcb' } 'hsv'; ...
                'whiteness'          'real'    [0 1] 0.45;...
                'newFigure'          'boolean' []  true;...
                'maximallyDistinctColors'    'boolean' [] true;...
                'showProjectedOnMrs'   'boolean' [] true;
                });
            
            allDomainsHaveColor = true;
            for i=1:numel(obj.domain)
                if isempty(obj.domain(i).color)
                    allDomainsHaveColor = false;
                    break;
                end;
            end;
            
            if allDomainsHaveColor % if domains already have color, use these colors
                domainColor = zeros(numel(obj.domain), 3);
                
                for i=1:numel(obj.domain)
                    domainColor(i,:) = obj.domain(i).color;
                end;
            else
                domainColor = makeDomainColor(obj, inputOptions);
            end;
            
            if inputOptions.newFigure
                figureHandle = figure;
                set(figureHandle, 'name', ['P <=' num2str(obj.significanceLevelForDomain)]);
            end;
            
            plot_dipplot_with_cortex;
            
            % remove all legends (set them to 'off')
            pr.remove_all_legends_from_figure;
            
            for i=1:length(obj.domain)
                
                % artifact from MPA paper
                %  if i==2
                %      voxelColor(i,:) = [0.0 255 231] / 255;
                %  end;
                
                domainHgGroup(i) = hggroup('DisplayName', ['Domain ' num2str(i)]);
                set(get(get(domainHgGroup(i), 'annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                pr.plot_head_region(obj.headGrid, obj.domain(i).membershipCube, 'regionColor', domainColor(i,:), 'regionOptions', {'parent', domainHgGroup(i)}, 'showProjectedOnMrs', inputOptions.showProjectedOnMrs);
            end;
            
            
            legendHandle = legend('show');
            set(legendHandle,  'textcolor', [1 1 1]);
        end;
        
        function plotVolumeColoredByDomain(obj, varargin)
            % plotVolumeColoredByMeasure(dipoleAndMeasure, significanceLevel, 'key', 'value'..)
            
            inputOptions = finputcheck(varargin, ...
                {'colormap'         'string'  {'hsv' 'ycrcb' } 'hsv'; ...
                'whiteness'          'real'    [0 1] 0.45;...
                'newFigure'          'boolean' []  true;...
                'showProjectedOnMrs' 'boolean' [] true;
                'surfaceOptions'     'cell'    {} {};... % options to be aplied to created patch objects.
                'plotOptions'        'cell'    {} {};...  % options directly passed to plot_head_surface()
                'plotSeparate'       'boolean' [] true;...
                'maximallyDistinctColors'    'boolean' [] true;...
                });
            
            allDomainsHaveColor = true;
            for i=1:numel(obj.domain)
                if isempty(obj.domain(i).color)
                    allDomainsHaveColor = false;
                    break;
                end;
            end;
            
            if allDomainsHaveColor % if domains already have color, use these colors
                domainColor = zeros(numel(obj.domain), 3);
                
                for i=1:numel(obj.domain)
                    domainColor(i,:) = obj.domain(i).color;
                end;
            else
                domainColor = makeDomainColor(obj, inputOptions);
            end;
            
            if inputOptions.newFigure
                figure;
            end;
            
            set(gcf, 'name', ['P <=' num2str(obj.significanceLevelForDomain)]);
            plot_dipplot_with_cortex;
            
            % remove all legends (set them to 'off')
            pr.remove_all_legends_from_figure;
            
            if inputOptions.plotSeparate     % plot each domain as a separate surface.
                for i=1:length(obj.domain)
                    domainHgGroup(i) = hggroup('DisplayName', ['Domain ' num2str(i)]);
                    set(get(get(domainHgGroup(i), 'annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                    pr.plot_head_surface(obj.headGrid, obj.domain(i).membershipCube, 'surfaceColor', domainColor(i,:), 'surfaceOptions', {'parent', domainHgGroup(i), inputOptions.surfaceOptions{:}}, 'showProjectedOnMrs', inputOptions.showProjectedOnMrs, inputOptions.plotOptions{:});
                end;
            else % plot all domains in a single surface with interpolated domain colors
                voxelColor = domainColor(obj.domainNumberCube(find(obj.domainNumberCube)),:);
                
                domainHgGroupWithNoDisplay = hggroup('DisplayName', ['']);
                pr.plot_head_surface(obj.headGrid, logical(obj.domainNumberCube), 'surfaceColor', voxelColor, 'surfaceOptions', {'parent', domainHgGroupWithNoDisplay, inputOptions.surfaceOptions{:}}, 'showProjectedOnMrs', inputOptions.showProjectedOnMrs, inputOptions.plotOptions{:});
                
                % adding legends by adding dummy patches and assigning them legen colors and text.
                for i=1:length(obj.domain)
                    domainHgGroup(i) = hggroup('DisplayName', ['Domain ' num2str(i)]);
                    set(get(get(domainHgGroup(i), 'annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                    t = zeros(3,3);
                    patch(t,t,t,  domainColor(i,:), 'parent', domainHgGroup(i))
                end;
            end;
            
            legendHandle = legend('show');
            set(legendHandle,  'textcolor', [1 1 1]);
        end;
        
        function plotCortexColoredByDomain(obj, varargin)
            % plotCortexColoredByDomain({key, value pairs})
            
            inputOptions = finputcheck(varargin, ...
                {'colormap'         'string'  {'hsv' 'ycrcb' } 'hsv'; ...
                'newFigure'          'boolean' []  true;...
                'plotOptions'        'cell'    {} {};...  % options directly passed to plot_cortex()
                'maximallyDistinctColors'    'boolean' [] true;...
                });
            
            allDomainsHaveColor = true;
            for i=1:numel(obj.domain)
                if isempty(obj.domain(i).color)
                    allDomainsHaveColor = false;
                    break;
                end;
            end;
            
            if allDomainsHaveColor % if domains already have color, use these colors
                domainColor = zeros(numel(obj.domain), 3);
                
                for i=1:numel(obj.domain)
                    domainColor(i,:) = obj.domain(i).color;
                end;
            else
                domainColor = makeDomainColor(obj, inputOptions);
            end;
            
            fsf = load('MNImesh_dipfit.mat');
            cortexVertices = fsf.vertices;
            
            for domainNumber = 1:length(obj.domain)
                [dummy, dipoleDensity]= pr.meanProjection.getProjectionMatrix(obj.domain(domainNumber).dipoleAndMeasure, obj.headGrid, obj.projectionParameter, obj.domain(domainNumber).membershipCube);
                dipoleDensity = dipoleDensity / sum(dipoleDensity);
               
                
                domainLocation = obj.headGrid.getPosition(obj.domain(domainNumber).membershipCube);
                headGridSpacing =  obj.headGrid.spacing;
                
                cortexPointDomainDenisty(:,domainNumber) = pr.project_domain_to_cortex(domainLocation, cortexVertices, headGridSpacing, dipoleDensity);
            end;
                        
            pr.plot_cortex(cortexPointDomainDenisty, domainColor, 'newFigure', inputOptions.newFigure, inputOptions.plotOptions{:});            
            
            set(gcf, 'name', ['P <=' num2str(obj.significanceLevelForDomain)]);
            
            % remove all legends (set them to 'off')
            pr.remove_all_legends_from_figure;
            
            % adding legends by adding dummy patches and assigning them legen colors and text.
            for i=1:length(obj.domain)
                domainHgGroup(i) = hggroup('DisplayName', ['Domain ' num2str(i)]);
                set(get(get(domainHgGroup(i), 'annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                t = zeros(3,3);
                patch(t,t,t,  domainColor(i,:), 'parent', domainHgGroup(i))
            end;
            
            legendHandle = legend('show');
            set(legendHandle, 'textcolor', [1 1 1]);
            set(legendHandle, 'color', [0 0 0]);
            set(legendHandle, 'box', 'off');
        end;
        
        function plotVoxelColoredBySignificance(obj,  varargin)
            plotColoredBySignificance(obj, 'voxel', varargin{:})
        end;
        
        function plotVolumeColoredBySignificance(obj,  varargin)
            plotColoredBySignificance(obj, 'volume', varargin{:})
        end;
        
        function obj = assignDomainColor(obj, varargin)
            % colors Domains with either distinct colors or an MDS of their exemplar projections.
            
            inputOptions = finputcheck(varargin, ...
                {'colormap'         'string'  {'hsv' 'ycrcb' } 'hsv'; ...
                'whiteness'          'real'    [0 1] 0.45;...
                'maximallyDistinctColors'  'boolean' [] true;...
                });
            
            domainColor = makeDomainColor(obj, inputOptions);
            
            numberOfDomains = length(obj.domain);
            for i=1:numberOfDomains
                obj.domain(i).color = domainColor(i,:);
            end;
        end;
        
        function [obj newOrder valueToSortBy] = sortDomain(obj, measureToSortby, varargin)
            % [obj newOrder]= sortDomain(measureToSortby)
            % Sort Domain based on the method provided in 'measureToSortby' input parameter and
            % returns a new pr.meanProjection object with sorted Domains (descending order)
            %
            % measureToSortby    string, any of 'dipoleMass', 'volume', 'sessionEntropy' or 'subjectEntropy'. The default is 'dipoleMass'.
            
            if nargin < 2
                measureToSortby = 'dipoleMass';
            end;
            
            numberOfDomains = numel(obj.domain);
            valueToSortBy = zeros(numberOfDomains, 1);
            
            switch lower(measureToSortby)
                case 'dipolemass'
                    for i=1:numberOfDomains
                        valueToSortBy(i) = obj.domain(i).getTotalDipoleMass;
                    end;
                case 'volume'
                    for i=1:numberOfDomains
                        valueToSortBy(i) = obj.domain(i).getTotalVolume;
                    end;
                case 'sessionentropy'
                    for i=1:numberOfDomains
                        valueToSortBy(i) = obj.domain(i).getSessionEntropy;
                    end;
                case 'subjectentropy'
                    for i=1:numberOfDomains
                        valueToSortBy(i) = obj.domain(i).getSubjectEntropy;
                    end;
                otherwise
                    error('Measure Projection: provided sorting method is not supported');
            end;
            
            [dummy, newOrder] = sort(valueToSortBy, 'descend');
            
            obj.domain = obj.domain(newOrder);
        end;
        
        function obj = performTopographicSegmentationOnAllDomains(obj, dipoleAndMeasure)
            newDomain = [];
            
            for i=1:length(obj.domain)
                newDomain = cat(2, newDomain, obj.domain(i).segmentByDomainTopography(dipoleAndMeasure));
            end;
            
            obj.domain = newDomain;
        end;
        
        function obj = removeVerySmallDomains(obj, sizeThreshold)
            if nargin < 2
                sizeThreshold = 4;
            end;
            
            for i=1:length(obj.domain)
                domainSize(i) = sum(obj.domain(i).membershipCube(:));
            end;
            
            obj.domain(domainSize < sizeThreshold) = [];
			fprintf('%d vrey small domains removed.\n', sum(domainSize < sizeThreshold));
        end;
    end
end