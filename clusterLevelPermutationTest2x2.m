% clusterLevelPermutationTest2x2() - Performs permutation test between 2x2 samples. For multiple comparison
%                                    correction, cluster-level correction is applied [1].
%                              
% Reference: [1] Groppe, Urbach, Kutas, 2011. Mass univariate analysis of
%                event-related brain potentials/fields I: A critical tutorial
%                review. Psychophysiology, xx, 1-15. (See also Korn et al.,
%                2004)
%
% Usage
%   >> [mask, tScore, pValue, surroMassOfCluster] = clusterLevelPermutationTest2x2(input1, input2, input3, input4, repeatedMeasuresFlag, pValForPreselection, numIterations)
% 
% Input
%  input1, input2, input3, input4: Data matrices. For example, the first dimension is ERP or
%                                  vectorized time-frequency measure, and the second dimension is
%                                  ICs or subjects.
%            repeatedMeasuresFlag: 1-paired (aka repeated measures) test; 2-mixed-design test, 1-2 and 3-4 must be paired; 3-two two-sample tests.
%             pValForPreselection: This determines the cluster size.
%                    numIteration: Number of bootstrap iteration (recommended: 10000) 
% 
% Output   
%                mask  : logical mask for significant
%                tScore: The tScore is for input1-input2. Positive result means input1 > input2
%                       (same as ttest2). This is computed by Zhou-Gao-Hui
%                       bootstrap method to compute tScore.
%               pValue: This is computed by standard bootstrap test.
%   surroMassOfCluster: [minSurroStats maxSurroStats] The surrogate maximum
%                       statistics. The data length is 2 x numSorro, in
%                       which the first half is min-stats and the latter
%                       half is the max-stats.

% History
% 03/29/2023 Makoto. Bug fixed. The mixed-design test was using the paird t-test inputs. Annotations corrected.
% 05/29/2020 Makoto. Bug fixed. tcdf() only left tails was tested, but now both tails.
% 05/27/2020 Makoto. Sample size added.
% 05/25/2020 Makoto. Subtraction of subtraction supported.
% 03/02/2017 Makoto. Created.
%
% Copyright (C) 2020, Makoto Miyakoshi (mmiyakoshi@ucsd.edu) , SCCN,INC,UCSD
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [mask, tScore, pValue, surroMassOfCluster] = clusterLevelPermutationTest2x2(input1, input2, input3, input4, repeatedMeasuresFlag, pValForPreselection, numIterations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform t-test for true difference %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input1_2D = reshape(input1, [size(input1,1)*size(input1,2) size(input1,3)]);
input2_2D = reshape(input2, [size(input2,1)*size(input2,2) size(input2,3)]);
input3_2D = reshape(input3, [size(input3,1)*size(input3,2) size(input3,3)]);
input4_2D = reshape(input4, [size(input4,1)*size(input4,2) size(input4,3)]);

% For a paired (aka repeated-measure) test.
if     repeatedMeasuresFlag == 1
    [~, pValues, ~, stats] = ttest((input1_2D-input2_2D)-(input3_2D-input4_2D));

% For a mixed-design test.
elseif repeatedMeasuresFlag == 2
    [~, pValues, ~, stats] = ttest2(input1_2D-input2_2D, input3_2D-input4_2D);

% For two two-sample tests. 05/25/2020 Makoto.
elseif repeatedMeasuresFlag == 3
    
    % Compute mean.
    X1_bar = mean(input1_2D,2);
    X2_bar = mean(input2_2D,2);
    X3_bar = mean(input3_2D,2);
    X4_bar = mean(input4_2D,2);
    
    % Compute variance.
    sSq1 = var(input1_2D,1,2);
    sSq2 = var(input2_2D,1,2);
    sSq3 = var(input3_2D,1,2);
    sSq4 = var(input4_2D,1,2);
    
    % Compute mean, variance, and sample size for (X1-X2) and (X3-X4)
    X12_bar = X1_bar - X2_bar;
    sSq12   = sSq1   + sSq2;
    X34_bar = X3_bar - X4_bar;
    sSq34   = sSq3   + sSq4;    
    % N12     = (size(input1_2D,2) + size(input2_2D,2))/2; % This was not found in a textbook, but this should be fine.
    % N34     = (size(input3_2D,2) + size(input4_2D,2))/2;
    N12     = size(input1_2D,2) + size(input2_2D,2); % On the second thought, I thought sample sizes needs to be added.
    N34     = size(input3_2D,2) + size(input4_2D,2);
    
    % Perform Welch's t-test. See Wikipedia Welch's t-test.
    stats.tstat = (X12_bar - X34_bar) ./ sqrt(sSq12/N12 + sSq34/N34);
    
    % Compute degrees of freedom using Welch-Satterthwaite equation. This
    % equasion is used to obtain effective degrees of freedom in the
    % calse of linear combination of independent samples.
    nu               = (sSq12/N12 + sSq34/N34).^2 ./ (sSq12.^2/(N12.^2*(N12-1)) + sSq34.^2/(N34.^2*(N34-1)));
    percentileValues = tcdf(stats.tstat, nu);
    combinedPercentileValues = [percentileValues 1-percentileValues];
    pValues = min(combinedPercentileValues,[],2);

        % % Visualization for debugging.
        % plotData = reshape(stats.tstat, [40 103]);
        % figure
        % imagesc(plotData); axis xy
        % hold on
        % significantMask = bwlabeln(reshape(pValues<0.05, [40 103]));
        % contour(logical(significantMask), 'color', [0 0 0])
end

% Compute cluster statistics.
pVal_2D   = reshape(pValues,     [size(input1,1) size(input1,2)]);
tScore_2D = reshape(stats.tstat, [size(input2,1) size(input2,2)]);
pvalMask  = pVal_2D < pValForPreselection;



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
combinedData = [input1_2D input2_2D input3_2D input4_2D];
input1Length = size(input1_2D,2);
input2Length = size(input2_2D,2);
input3Length = size(input3_2D,2);
input4Length = size(input4_2D,2);
permSize     = size(combinedData,2);
for iterationIdx = 1:numIterations
    
    % Generate surrogate data using permutation.
    permIdx = randperm(permSize);
    surro1_2D = combinedData(:,permIdx(1:input1Length));
    surro2_2D = combinedData(:,permIdx(input1Length+1:(input1Length+input2Length)));
    surro3_2D = combinedData(:,permIdx((input1Length+input2Length)+1:(input1Length+input2Length+input3Length)));
    surro4_2D = combinedData(:,permIdx((input1Length+input2Length+input3Length)+1:end));
    
    % For fixed-effect design test.
    if     repeatedMeasuresFlag == 1
        [~, pValues, ~, stats] = ttest((surro1_2D-surro2_2D)-(surro3_2D-surro4_2D));
        
        % For mixed-effects design test.
    elseif repeatedMeasuresFlag == 2
        [~, pValues, ~, stats] = ttest2((surro1_2D-surro2_2D)-(surro3_2D-surro4_2D));
        
        % For random-effect design test. 05/25/2020 Makoto.
    elseif repeatedMeasuresFlag == 3
        
        % Compute mean.
        X1_bar = mean(surro1_2D,2);
        X2_bar = mean(surro2_2D,2);
        X3_bar = mean(surro3_2D,2);
        X4_bar = mean(surro4_2D,2);
        
        % Compute variance.
        sSq1 = var(surro1_2D,1,2);
        sSq2 = var(surro2_2D,1,2);
        sSq3 = var(surro3_2D,1,2);
        sSq4 = var(surro4_2D,1,2);
        
        % Compute mean, variance, and sample size for (X1-X2) and (X3-X4)
        X12_bar = X1_bar - X2_bar;
        sSq12   = sSq1   + sSq2;
        X34_bar = X3_bar - X4_bar;
        sSq34   = sSq3   + sSq4;
        % N12     = (size(surro1_2D,2) + size(surro2_2D,2))/2; % This was not found in a textbook, but this should be fine.
        % N34     = (size(surro3_2D,2) + size(surro4_2D,2))/2;
        N12     = size(surro1_2D,2) + size(surro2_2D,2); % On the second thought, I thought sample sizes needs to be added.
        N34     = size(surro3_2D,2) + size(surro4_2D,2);
        
        % Perform Welch's t-test. See Wikipedia Welch's t-test.
        statsSurro.tstat = (X12_bar - X34_bar) ./ sqrt(sSq12/N12 + sSq34/N34);
        
        % Compute degrees of freedom using Welch-Satterthwaite equation. This
        % equasion is used to obtain effective degrees of freedom in the
        % calse of linear combination of independent samples.
        nuSurro      = (sSq12/N12 + sSq34/N34).^2 ./ (sSq12.^2/(N12.^2*(N12-1)) + sSq34.^2/(N34.^2*(N34-1)));
        percentileValuesSurro = tcdf(statsSurro.tstat, nuSurro);
        combinedPercentileValuesSurro = [percentileValuesSurro 1-percentileValuesSurro];
        pValuesSurro = min(combinedPercentileValuesSurro,[],2);
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
    massOfClusterSurro = massOfClusterSurro(2:end);
    
    % Store the minimum and maximum of massOfClusterSurro
    surroMassOfCluster(iterationIdx,:) = [min(massOfClusterSurro) max(massOfClusterSurro)];
end