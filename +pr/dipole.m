classdef dipole % contains positions of dipoles for Measure Projection
    properties
        location = []; % an N x 3 array containing N dipole positions.
        direction = []; % an N x 3 array containing N dipole directions.
        residualVariance = []; % an N x 1 array containing N dipole residual variances.
        coordinateFormat = 'mni';
    end % properties
    methods
        function plotDipoleDensity(obj, varargin)
            
            if isempty(varargin)
                varargin = {'plotargs',{'mriview' , 'top','mrislices', [-50 -30 -20 -15 -10 -5 0 5 10 15 20 25 30 40 50 60 70]}};
            end;
            
            inputOptions = finputcheck(varargin, ...
                { 'plotargs'         'cell'   {}       {}; ...
                  'normalizeInBrainDensity'   'boolean'   []       true; ...
                  'gaussianWidth'    'real'     [0 Inf]    12}); % use 12 mm Gaussians by default
            
              if inputOptions.normalizeInBrainDensity
                  dipoleWeight = pr.meanProjection.calculateDipoleInBrainDensityNormalizationFactor(obj, pr.headGrid, pr.projectionParameter(inputOptions.gaussianWidth));
              else % no normalization
                  dipoleWeight = ones(1, obj.numberOfDipoles);
              end;
              
            pr.dipoledensity(obj.location, 'coordformat', obj.coordinateFormat, 'plot','on',  'methodparam', inputOptions.gaussianWidth,  'plotargs' ,inputOptions.plotargs, 'weight', dipoleWeight); 
            
            set(gcf, 'InvertHardcopy', 'off');
        end;
        
        function plotDipole(obj, withCortex, createNewFigure, varargin)
            % plotDipole(obj, withCortex, createNewFigure, varargin)
            
            if nargin < 2
                withCortex = true;
            end;
            
            if nargin < 3
                createNewFigure = true;
            end;
            
            if createNewFigure
                figure;
            end;
                  
            % only shows 3D spheres for dipoles when their number if less
            % than a certain limit. This is to prevent verylong plot time.
            
            if obj.numberOfDipoles < 600
                plot3dSpheres = 'on';
            else
                plot3dSpheres = 'off';
            end;
            
            plot_dipplot_with_cortex(obj.location, withCortex, 'coordformat', obj.coordinateFormat, 'spheres',plot3dSpheres, 'gui', 'off', varargin{:});
            %dipplot(source, 'coordformat', obj.coordinateFormat, 'spheres',plot3dSpheres 'gui', 'off', varargin{:});
        end;
        
        function [dipoleId sortedDipoleDensity orderOfDipoles dipoleDenisty dipoleDenistyInRegion] = getDipoleDensityContributionToRegionOfInterest(obj, regionOfInterest, projection, cutoffRatio)
            % [dipoleId sortedDipoleDensity orderOfDipoles dipoleDenistyInRegion dipoleDenistyInRegion] = getDipoleDensityContributionToRegionOfInterest(obj, projection, regionOfInterest, cutoffRatio)
            % projection is a type of meanProjection (or compatible).
            % return dipoles ids (in dipoleId)  which cumulatively contribute at least 'cutoffRatio' to the dipole
            % denisty over the region.
            %
            % cutoffRatio is either a scalar or a 2-vector. The first number contains the percent of
            % region (domain) dipole mass explained by selected dipoles after which there will be a no 
            % other dipole selected. The second is the miminum dipole mass ratio contribution to the 
            % region (dipoles with a  contribution less than this value will not be selected).
            % For example cutoffRatio = 0.98 requires selected dipoles to at least explain 98% of
            % dipoles mass in the region.
            % cutoffRatio = [1 0.05] means that all dipoles that at least contribute %5 of their
            % mass to the region will be selected.
            % default cutoffRatio is [1 0.05].
                                     
            if nargin < 4
                cutoffRatio = [1 0.05];
            end;
            
            if length(cutoffRatio) == 1
                cutoffRatio = [cutoffRatio 0.05];
            end;
            
            [projectionMatrix dipoleDenisty gaussianWeightMatrix]= pr.meanProjection.getProjectionMatrix(obj, projection.headGrid, projection.projectionParameter, regionOfInterest);
                       
            dipoleDenistyInRegion = sum(gaussianWeightMatrix,2);
            [sortedDipoleDensity orderOfDipoles] = sort(dipoleDenistyInRegion,'descend');
            
            cumulativePercentProjected = cumsum(sortedDipoleDensity) / sum(sortedDipoleDensity);
            numberOfSortedDipolesToReturnBasedOnCumulativeDensity = find(cumulativePercentProjected >= cutoffRatio(1),1,  'first');
            numberOfSortedDipolesToReturnBasedOnDipoleDensityInRegion = find(sortedDipoleDensity < cutoffRatio(2),1,  'first');
            
            numberOfSortedDipolesToReturn = min([numberOfSortedDipolesToReturnBasedOnCumulativeDensity, numberOfSortedDipolesToReturnBasedOnDipoleDensityInRegion]);
            
            dipoleId = orderOfDipoles(1:numberOfSortedDipolesToReturn);
        end;
        
        function obj = convertDipoleCoordinatesFromSphericalToMni(obj)                      
            fprintf('converting dipoles from spherical to MNI coordinates...\n');
            for i=1:size(obj.location,1)
                locationInMniCoordinate = sph2spm * [obj.location(i,:) 1]';
                obj.location(i,:) = locationInMniCoordinate(1:3);
                
                directionInMniCoordinate = sph2spm * [obj.direction(i,:) 1]';               
                obj.direction(i,:) = directionInMniCoordinate(1:3);
            end;
        end;
        
        function theNumber = numberOfDipoles(obj)
            theNumber = size(obj.location, 1);
        end;
        
        function plotDipoleDensityOnCortex(obj, minValueNormalizedTo, projectionParameter, varargin)
            
            if nargin < 2
                minValueNormalizedTo = 0.06;
            end;
            
            if nargin < 3
                projectionParameter = pr.projectionParameter;
            end;
            
            headGrid = pr.headGrid(8);
            [dummy, dipoleDensity]= pr.meanProjection.getProjectionMatrix(obj, headGrid, projectionParameter, headGrid.insideBrainCube);
            
            fsf = load('MNImesh_dipfit.mat');
            cortexVertices = fsf.vertices;
            
            cortexPointDomainDenisty = pr.project_domain_to_cortex(headGrid.getPosition, cortexVertices, headGrid.spacing, dipoleDensity);
            pr.plot_cortex(cortexPointDomainDenisty, [0 1 0], 'minValueNormalizedTo', minValueNormalizedTo);
        end;
        
        function [obj dipoleModelNumber]= addFromDipfitStructure(obj, dipfit, residualvarianceThreshold)
            % [obj dipoleModelNumber]= addFromDipfitStructure(obj, dipfit, residualvarianceThreshold)
            
            if nargin < 3
                residualvarianceThreshold = 0.15;
            end;
            
            dipoleModelNumber = [];           
            for dipoleNumber = 1:length(dipfit.model)
                  if dipfit.model(dipoleNumber).rv <= residualvarianceThreshold                         
                        dipoleModelNumber = [dipoleModelNumber dipoleNumber];
                        
                        obj.location(end+1, :) = dipfit.model(dipoleNumber).posxyz(1,:);
                        obj.direction(end+1, :) = dipfit.model(dipoleNumber).momxyz(1,:);
                        obj.residualVariance(end+1) = dipfit.model(dipoleNumber).rv;
                        
                        if size(dipfit.model(dipoleNumber).posxyz, 1) == 2 % bilateral dipoles                           
                            dipoleModelNumber = [dipoleModelNumber dipoleNumber];
                            
                            obj.location(end+1, :) = dipfit.model(dipoleNumber).posxyz(2,:);
                            obj.direction(end+1, :) = dipfit.model(dipoleNumber).momxyz(2,:);
                            obj.residualVariance(end+1) = dipfit.model(dipoleNumber).rv;
                        end;
                  end;
             end;
        end;
    end;
end