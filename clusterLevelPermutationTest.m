% clusterLevelPermutationTest() - Performs permutation test between two samples. For multiple comparison
%                                 correction, cluster-level correction is applied [1].
%                              
% Reference: [1] Groppe, Urbach, Kutas, 2011. Mass univariate analysis of
%                event-related brain potentials/fields I: A critical tutorial
%                review. Psychophysiology, xx, 1-15. (See also Korn et al.,
%                2004)
%
% Usage
%   >> [mask, tScore, pValue] = clusterLevelPermutationTest(input1, input2, repeatedMeasuresFlag, pValForPreselection, numIterations)
% 
% Input
%              input1, input2: Data matrices. For example, the first dimension is ERP or
%                              vectorized time-frequency measure, and the second dimension is
%                              ICs or subjects.
%        repeatedMeasuresFlag: 1, paired t-test; 0, two-sample t-test.
%         pValForPreselection: This determines the cluster size.
%                numIteration: Number of bootstrap iteration (recommended: 10000) 
% 
% Output
%           mask  : logical mask for significant
%
%           tScore: The tScore is for input1-input2. Positive result means input1 > input2
%                   (same as ttest2). This is computed by Zhou-Gao-Hui
%                   bootstrap method to compute tScore.
%           pValue: This is computed by standard bootstrap test.

% History
% 06/28/2020 Makoto. 
% 05/23/2020 Makoto. Cleaned code.
% 03/02/2017 Makoto. Created.
%
% Copyright (C) 2017 Makoto Miyakoshi, SCCN,INC,UCSD;
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [mask, tScore, pValue, surroMassOfCluster] = clusterLevelPermutationTest(input1, input2, repeatedMeasuresFlag, pValForPreselection, numIterations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform t-test for true difference %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input1_2D = reshape(input1, [size(input1,1)*size(input1,2) size(input1,3)]);
input2_2D = reshape(input2, [size(input2,1)*size(input2,2) size(input2,3)]);
if     repeatedMeasuresFlag == 1
    [~, pValues, ~, stats] = ttest( input1_2D', input2_2D');
else
    [~, pValues, ~, stats] = ttest2(input1_2D', input2_2D');
end

% Compute cluster statistics.
pVal_2D   = reshape(pValues,     [size(input1,1) size(input1,2)]);
tScore_2D = reshape(stats.tstat, [size(input2,1) size(input2,2)]);
pvalMask = pVal_2D< pValForPreselection;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% If no significant result in the uncorrected result, exit (to save time). %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(pvalMask(:))==0
    disp('No significant result.')
    mask               = zeros(size(input1,1), size(input1,2));
    tScore             = zeros(size(input1,1), size(input1,2));
    pValue             = ones( size(input1,1), size(input1,2));
    surroMassOfCluster = zeros(numIterations, 2);
    return
end



% Extract clusters of significant pixels.
connectedComponentLabels = bwlabeln(pvalMask); % This requires image processing toolbox
[entryCount, blobId]  = hist(connectedComponentLabels(:), unique(connectedComponentLabels(:)));
massOfCluster = zeros(length(blobId),1);
for n = 1:length(blobId)
    currentMask = connectedComponentLabels==blobId(n);
    massOfCluster(n) = sum(sum(currentMask.*tScore_2D));
end
massOfCluster = massOfCluster(2:end);

% Prepare outputs.
tScore = tScore_2D;
pValue = pVal_2D;
mask   = connectedComponentLabels;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform surrogate test %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surroMassOfCluster = zeros(numIterations,2);
combinedData = [input1_2D input2_2D];
permSize     = size(combinedData,2);
for iterationIdx = 1:numIterations
    
%     % Report progress.
%     if mod(iterationIdx, 500) == 0
%         disp(sprintf('%.0f/%.0f...', iterationIdx, numIteration))
%     end

    % Perform t-test
    permIdx = randperm(permSize);
    surro1 = combinedData(:,permIdx(1:size(input1_2D,2)));
    surro2 = combinedData(:,permIdx(size(input1_2D,2)+1:end));
    
    if     repeatedMeasuresFlag == 1
        [~, pValuesSurro, ~, statsSurro] = ttest(surro1', surro2');
    else
        [~, pValuesSurro, ~, statsSurro] = ttest2(surro1', surro2');
    end
    
    % Compute cluster statistics
    pValSurro_2D   = reshape(pValuesSurro,     [size(input1,1) size(input1,2)]);
    tScoreSurro_2D = reshape(statsSurro.tstat, [size(input2,1) size(input2,2)]);
    pvalSurroMask = pValSurro_2D< pValForPreselection;
    
    % If no significant result in the uncorrected result, exit (to save time).
    if any(pvalSurroMask(:))==0
        % disp('No significant result.')
        continue
    end
    
    % Extract clusters of significant pixels
    connectedComponentLabels = bwlabeln(pvalSurroMask); % This requires image processing toolbox
    [entryCount, blobId]  = hist(connectedComponentLabels(:), unique(connectedComponentLabels(:)));
    massOfClusterSurro = zeros(length(blobId),1);
    for n = 1:length(blobId)
        currentMask = connectedComponentLabels==blobId(n);
        massOfClusterSurro(n) = sum(sum(currentMask.*tScore_2D));
    end
    
    if sum(pvalSurroMask(:)) ~= size(pvalSurroMask,1)*size(pvalSurroMask,2) % if equal, it means all the points/pixels are significant so no 'background'.
        massOfClusterSurro = massOfClusterSurro(2:end);
    end
    
    % Store the minimum and maximum of massOfClusterSurro
    surroMassOfCluster(iterationIdx,:) = [min(massOfClusterSurro) max(massOfClusterSurro)];
end