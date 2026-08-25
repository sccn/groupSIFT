function [mask, tScore, pValue, surroMassOfCluster] = clusterLevelPermutationTest(input1, input2, repeatedMeasuresFlag, pValForPreselection, numIterations, randomSeed)
% clusterLevelPermutationTest() - Two-sided cluster-mass permutation test.
%
% Usage:
%   [mask, tScore, pValue, surroMassOfCluster] = ...
%       clusterLevelPermutationTest(input1, input2, repeatedMeasuresFlag, ...
%       pValForPreselection, numIterations [, randomSeed])
%
% Inputs:
%   input1, input2       Frequency x time x subject data arrays.
%   repeatedMeasuresFlag 1 for paired samples; 0 for independent samples.
%   pValForPreselection Two-sided pointwise cluster-forming threshold.
%   numIterations       Number of random permutations.
%   randomSeed          Optional scalar seed or randomization-plan
%                       structure. Paired plans use signs; independent
%                       plans use uniformScores. The caller's RNG state is
%                       restored after scalar-seed use.
%
% Outputs:
%   mask                Integer cluster labels. Positive and negative
%                       clusters are formed separately.
%   tScore              Observed input1-minus-input2 t statistics.
%   pValue              Observed two-sided pointwise p values.
%   surroMassOfCluster  numIterations x 2 matrix containing the most
%                       negative and most positive surrogate cluster mass.
%
% Within-subject permutations independently swap the two condition labels
% within every subject (equivalently, sign-flip paired differences).

if nargin < 6
    randomSeed = [];
end
validateInputs(input1, input2, repeatedMeasuresFlag, pValForPreselection, numIterations, randomSeed);

randomizationPlan = [];
if isstruct(randomSeed)
    randomizationPlan = randomSeed;
elseif ~isempty(randomSeed)
    previousRngState = rng;
    restoreRng = onCleanup(@() rng(previousRngState)); %#ok<NASGU>
    rng(randomSeed, 'twister');
end

spatialSize = [size(input1, 1), size(input1, 2)];
input1_2D = reshape(input1, prod(spatialSize), size(input1, 3));
input2_2D = reshape(input2, prod(spatialSize), size(input2, 3));

if repeatedMeasuresFlag == 1
    [~, pValues, ~, stats] = ttest(input1_2D', input2_2D');
else
    [~, pValues, ~, stats] = ttest2(input1_2D', input2_2D', 'Vartype', 'unequal');
end
observedEffect = mean(input1_2D, 2)' - mean(input2_2D, 2)';
[pValues, stats.tstat] = normalizeDegenerateStatistics(pValues, stats.tstat, observedEffect);

pValue = reshape(pValues, spatialSize);
tScore = reshape(stats.tstat, spatialSize);
[mask, ~] = groupSIFT_labelClusters(pValue, tScore, pValForPreselection);

surroMassOfCluster = zeros(numIterations, 2);
if repeatedMeasuresFlag == 1
    pairedDifference = input1_2D - input2_2D;
    if ~isempty(randomizationPlan)
        validateRandomizationPlan(randomizationPlan, 'signs', numIterations, size(pairedDifference, 2));
    end
else
    combinedData = [input1_2D, input2_2D];
    numberInInput1 = size(input1_2D, 2);
    if ~isempty(randomizationPlan)
        validateRandomizationPlan(randomizationPlan, 'uniformScores', numIterations, size(combinedData, 2));
    end
end

for iterationIdx = 1:numIterations
    if repeatedMeasuresFlag == 1
        if isempty(randomizationPlan)
            randomSigns = ones(1, size(pairedDifference, 2));
            randomSigns(rand(1, size(pairedDifference, 2)) < 0.5) = -1;
        else
            randomSigns = randomizationPlan.signs(iterationIdx, :);
        end
        surrogateDifference = bsxfun(@times, pairedDifference, randomSigns);
        [~, pValuesSurro, ~, statsSurro] = ttest(surrogateDifference');
        surrogateEffect = mean(surrogateDifference, 2)';
    else
        if isempty(randomizationPlan)
            permutationIdx = randperm(size(combinedData, 2));
        else
            [~, permutationIdx] = sort(randomizationPlan.uniformScores(iterationIdx, :));
        end
        surrogate1 = combinedData(:, permutationIdx(1:numberInInput1));
        surrogate2 = combinedData(:, permutationIdx(numberInInput1 + 1:end));
        [~, pValuesSurro, ~, statsSurro] = ttest2( ...
            surrogate1', surrogate2', 'Vartype', 'unequal');
        surrogateEffect = mean(surrogate1, 2)' - mean(surrogate2, 2)';
    end
    [pValuesSurro, statsSurro.tstat] = normalizeDegenerateStatistics( ...
        pValuesSurro, statsSurro.tstat, surrogateEffect);

    pValueSurro = reshape(pValuesSurro, spatialSize);
    tScoreSurro = reshape(statsSurro.tstat, spatialSize);
    [~, clusterMasses] = groupSIFT_labelClusters( ...
        pValueSurro, tScoreSurro, pValForPreselection);
    surroMassOfCluster(iterationIdx, :) = extremeMasses(clusterMasses);
end
end

function [pValues, tStatistics] = normalizeDegenerateStatistics(pValues, tStatistics, effect)
invalidMask = isnan(pValues) | isnan(tStatistics);
zeroEffectMask = invalidMask & (effect == 0);
nonzeroEffectMask = invalidMask & (effect ~= 0);
tStatistics(zeroEffectMask) = 0;
pValues(zeroEffectMask) = 1;
tStatistics(nonzeroEffectMask) = sign(effect(nonzeroEffectMask)) .* Inf;
pValues(nonzeroEffectMask) = 0;
end

function extremes = extremeMasses(clusterMasses)
negativeMasses = clusterMasses(clusterMasses < 0);
positiveMasses = clusterMasses(clusterMasses > 0);
extremes = [0, 0];
if ~isempty(negativeMasses)
    extremes(1) = min(negativeMasses);
end
if ~isempty(positiveMasses)
    extremes(2) = max(positiveMasses);
end
end

function validateInputs(input1, input2, repeatedMeasuresFlag, pThreshold, numIterations, randomSeed)
if ~isnumeric(input1) || ~isreal(input1) || isempty(input1) || ...
        ~isnumeric(input2) || ~isreal(input2) || isempty(input2)
    error('groupSIFT:InvalidInput', 'Inputs must be non-empty real numeric arrays.');
end
if any(~isfinite(input1(:))) || any(~isfinite(input2(:)))
    error('groupSIFT:NonfiniteInput', 'Inputs must not contain NaN or Inf values.');
end
if size(input1, 1) ~= size(input2, 1) || size(input1, 2) ~= size(input2, 2)
    error('groupSIFT:SpatialSizeMismatch', ...
        'The first two input dimensions must match.');
end
if ~isscalar(repeatedMeasuresFlag) || ~ismember(repeatedMeasuresFlag, [0, 1])
    error('groupSIFT:InvalidDesignFlag', ...
        'repeatedMeasuresFlag must be 1 (paired) or 0 (independent).');
end
if repeatedMeasuresFlag == 1 && size(input1, 3) ~= size(input2, 3)
    error('groupSIFT:PairedSampleSizeMismatch', ...
        'Paired inputs must contain the same number of subjects.');
end
if size(input1, 3) < 2 || size(input2, 3) < 2
    error('groupSIFT:InsufficientSampleSize', ...
        'Each input must contain at least two subjects.');
end
if ~isscalar(pThreshold) || ~isfinite(pThreshold) || pThreshold <= 0 || pThreshold >= 1
    error('groupSIFT:InvalidClusterThreshold', ...
        'pValForPreselection must be a scalar strictly between 0 and 1.');
end
if ~isscalar(numIterations) || ~isfinite(numIterations) || ...
        numIterations < 1 || numIterations ~= floor(numIterations)
    error('groupSIFT:InvalidIterationCount', ...
        'numIterations must be a positive integer.');
end
if isstruct(randomSeed) && ~isscalar(randomSeed)
    error('groupSIFT:InvalidRandomizationPlan', ...
        'The randomization plan must be a scalar structure.');
end
if ~isempty(randomSeed) && ~isstruct(randomSeed) && ...
        (~isscalar(randomSeed) || ~isfinite(randomSeed) || ...
        randomSeed < 0 || randomSeed ~= floor(randomSeed))
    error('groupSIFT:InvalidRandomSeed', ...
        'randomSeed must be empty, a nonnegative integer scalar, or a randomization-plan structure.');
end
end

function validateRandomizationPlan(plan, fieldName, numberOfIterations, numberOfSubjects)
if ~isfield(plan, fieldName)
    error('groupSIFT:InvalidRandomizationPlan', ...
        'The randomization plan must contain a %s field.', fieldName);
end
values = plan.(fieldName);
if ~isnumeric(values) || ~isequal(size(values), [numberOfIterations, numberOfSubjects]) || ...
        any(~isfinite(values(:)))
    error('groupSIFT:InvalidRandomizationPlan', ...
        '%s must be a finite numIterations-by-subject matrix.', fieldName);
end
if strcmp(fieldName, 'signs') && any(~ismember(values(:), [-1, 1]))
    error('groupSIFT:InvalidRandomizationPlan', ...
        'Every entry in signs must be -1 or 1.');
end
end
