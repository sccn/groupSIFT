classdef regionOfInterest
    % holds information and methods (functions) to work with a (optionally) probabilistic region of interest (ROI).
    % Projected measure, Dipole desnity or a combination of these.
    
    properties %(SetAccess = 'immutable')
        headGrid;
        membershipCube
        probabilityThreshold
        membershipProbabilityCube
    end;
    
    methods
        function obj = regionOfInterest(headGrid, varargin)
            obj.probabilityThreshold = 0.001;
            
            % by default select the whole brain as the regionOfInterest
            obj.membershipCube = headGrid.insideBrainCube;
            
            if nargin > 0
                obj.headGrid = headGrid;
            end;
        end;
    
        
        function value = get.membershipCube(obj)
            value = obj.membershipProbabilityCube >= obj.probabilityThreshold;
        end;
        
        % set it also for the robabilityCube
        function obj = set.membershipCube(obj, value)
            obj.membershipProbabilityCube = value;
            obj.membershipCube = value;
        end;
        
        function varargout = plotVolume(obj, varargin)          
            
             inputOptions = finputcheck(varargin, ...
                {'surfaceColor'      'real'  []  [0.4 1 0.4];...         
                 'surfaceOptions'    'cell' [] {};...
                 'newFigure'        'boolean'  []   true;...                
                });
            
            if inputOptions.newFigure
                figure;
            end;
            
            set(gcf, 'renderer', 'opengl');
            plot_dipplot_with_cortex;
            regionColor = [0.4 1 0.4];
            pr.plot_head_surface(obj.headGrid, obj.membershipCube, 'surfaceColor', inputOptions.surfaceColor, 'surfaceOptions', inputOptions.surfaceOptions);
        end;
        
        function varargout = plotVoxel(obj, varargin)
            figure;
            set(gcf, 'renderer', 'opengl');
            plot_dipplot_with_cortex;
            regionColor = [0.4 1 0.4];
            pr.plot_head_region(obj.headGrid, obj.membershipCube, 'regionColor', regionColor);
        end;
        
        function varargout = plotCortex(obj, varargin)
            figure;
            set(gcf, 'renderer', 'opengl');            
            regionColor = [0.4 1 0.4];
            
            fsf = load('MNImesh_dipfit.mat');
            cortexVertices = fsf.vertices;
            
            regionLocation = obj.headGrid.getPosition(obj.membershipCube);
            headGridSpacing =  obj.headGrid.spacing;
            
            cortexPointDomainDenisty = pr.project_domain_to_cortex(regionLocation, cortexVertices, headGridSpacing, obj.membershipProbabilityCube(obj.membershipCube(:)));
            pr.plot_cortex(cortexPointDomainDenisty, regionColor, 'newFigure', false);
            
            
            %pr.plot_head_region(obj.headGrid, obj.membershipCube, 'regionColor', regionColor);
        end;
    end
end