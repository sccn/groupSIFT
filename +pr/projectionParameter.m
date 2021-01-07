classdef projectionParameter
    % contains all parameters related to a projection
    % encaplustae all projection parameters, like number of Std., Gaussian Width and whether in-brain density should be normalized
    properties % (SetAccess = 'immutable'); % immutable properties are not supported on at least 2009a
        standardDeviationOfEstimatedDipoleLocation = 12; % in mm, default std. for gaussian cloud representing each dipole.
        numberOfStandardDeviationsToTruncatedGaussaian = 3;
        normalizeInBrainDipoleDenisty = true;
    end;
    
    methods
        function obj = projectionParameter(standardDeviationOfEstimatedDipoleLocation, numberOfStandardDeviationsToTruncatedGaussaian, normalizeInBrainDipoleDenisty)
            % obj = projectionParameter(standardDeviationOfEstimatedDipoleLocation, numberOfStandardDeviationsToTruncatedGaussaian, normalizeInBrainDipoleDenisty)
            if nargin >= 1
                obj.standardDeviationOfEstimatedDipoleLocation = standardDeviationOfEstimatedDipoleLocation;
            end;
            if nargin >= 2
                obj.numberOfStandardDeviationsToTruncatedGaussaian = numberOfStandardDeviationsToTruncatedGaussaian;
            end;
            if nargin >= 3
                obj.normalizeInBrainDipoleDenisty = normalizeInBrainDipoleDenisty;
            end;
        end;
    end;
end