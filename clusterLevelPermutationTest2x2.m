function [mask, tScore, pValue, surroMassOfCluster] = clusterLevelPermutationTest2x2(input1, input2, input3, input4, repeatedMeasuresFlag, pValForPreselection, numIterations, randomSeed)
% clusterLevelPermutationTest2x2() - Cluster test for a 2-by-2 interaction.
%
% repeatedMeasuresFlag:
%   1 - all four cells contain the same subjects.
%   2 - cells 1/2 and 3/4 are paired within two independent groups.
%   3 - all four cells are independent. Permutation tests the strong null
%       that observations are exchangeable across all four cells.
%
% The tested contrast is (input1-input2) - (input3-input4). Positive and
% negative clusters are formed separately. The optional final input may be
% a scalar RNG seed or a randomization-plan structure. A repeated design
% uses plan.signs; mixed and independent designs use plan.uniformScores.
% The caller's RNG state is restored after scalar-seed use.

if nargin < 8
    randomSeed = [];
end
inputs = {input1, input2, input3, input4};
validateInputs(inputs, repeatedMeasuresFlag, pValForPreselection, numIterations, randomSeed);

randomizationPlan = [];
if isstruct(randomSeed)
    randomizationPlan = randomSeed;
elseif ~isempty(randomSeed)
    previousRngState = rng;
    restoreRng = onCleanup(@() rng(previousRngState)); %#ok<NASGU>
    rng(randomSeed, 'twister');
end

spatialSize = [size(input1, 1), size(input1, 2)];
input2D = cell(1, 4);
for inputIdx = 1:4
    input2D{inputIdx} = reshape(inputs{inputIdx}, prod(spatialSize), size(inputs{inputIdx}, 3));
end

[pValues, observedT] = interactionStatistic(input2D, repeatedMeasuresFlag);
pValue = reshape(pValues, spatialSize);
tScore = reshape(observedT, spatialSize);
[mask, ~] = groupSIFT_labelClusters(pValue, tScore, pValForPreselection);

surroMassOfCluster = zeros(numIterations, 2);
if repeatedMeasuresFlag == 1
    interactionDifference = (input2D{1} - input2D{2}) - (input2D{3} - input2D{4});
    if ~isempty(randomizationPlan)
        validateRandomizationPlan(randomizationPlan, 'signs', ...
            numIterations, size(interactionDifference, 2));
    end
elseif repeatedMeasuresFlag == 2
    group1Difference = input2D{1} - input2D{2};
    group2Difference = input2D{3} - input2D{4};
    combinedDifferences = [group1Difference, group2Difference];
    numberInGroup1 = size(group1Difference, 2);
    if ~isempty(randomizationPlan)
        validateRandomizationPlan(randomizationPlan, 'uniformScores', ...
            numIterations, size(combinedDifferences, 2));
    end
else
    combinedData = [input2D{:}];
    cellSizes = cellfun(@(x) size(x, 2), input2D);
    cumulativeCellSizes = cumsum(cellSizes);
    if ~isempty(randomizationPlan)
        validateRandomizationPlan(randomizationPlan, 'uniformScores', ...
            numIterations, size(combinedData, 2));
    end
end

for iterationIdx = 1:numIterations
    if repeatedMeasuresFlag == 1
        if isempty(randomizationPlan)
            randomSigns = ones(1, size(interactionDifference, 2));
            randomSigns(rand(1, size(interactionDifference, 2)) < 0.5) = -1;
        else
            randomSigns = randomizationPlan.signs(iterationIdx, :);
        end
        surrogateDifference = bsxfun(@times, interactionDifference, randomSigns);
        [~, pValuesSurro, ~, statsSurro] = ttest(surrogateDifference');
        surrogateT = statsSurro.tstat;
        surrogateEffect = mean(surrogateDifference, 2)';
        [pValuesSurro, surrogateT] = normalizeDegenerateStatistics( ...
            pValuesSurro, surrogateT, surrogateEffect);
    elseif repeatedMeasuresFlag == 2
        if isempty(randomizationPlan)
            permutationIdx = randperm(size(combinedDifferences, 2));
        else
            [~, permutationIdx] = sort(randomizationPlan.uniformScores(iterationIdx, :));
        end
        surrogateGroup1 = combinedDifferences(:, permutationIdx(1:numberInGroup1));
        surrogateGroup2 = combinedDifferences(:, permutationIdx(numberInGroup1 + 1:end));
        [~, pValuesSurro, ~, statsSurro] = ttest2( ...
            surrogateGroup1', surrogateGroup2', 'Vartype', 'unequal');
        surrogateT = statsSurro.tstat;
        surrogateEffect = mean(surrogateGroup1, 2)' - mean(surrogateGroup2, 2)';
        [pValuesSurro, surrogateT] = normalizeDegenerateStatistics( ...
            pValuesSurro, surrogateT, surrogateEffect);
    else
        if isempty(randomizationPlan)
            permutationIdx = randperm(size(combinedData, 2));
        else
            [~, permutationIdx] = sort(randomizationPlan.uniformScores(iterationIdx, :));
        end
        surrogateInputs = cell(1, 4);
        startIdx = 1;
        for cellIdx = 1:4
            stopIdx = cumulativeCellSizes(cellIdx);
            surrogateInputs{cellIdx} = combinedData(:, permutationIdx(startIdx:stopIdx));
            startIdx = stopIdx + 1;
        end
        [pValuesSurro, surrogateT] = welchInteraction(surrogateInputs);
    end

    pValueSurro = reshape(pValuesSurro, spatialSize);
    tScoreSurro = reshape(surrogateT, spatialSize);
    [~, clusterMasses] = groupSIFT_labelClusters( ...
        pValueSurro, tScoreSurro, pValForPreselection);
    surroMassOfCluster(iterationIdx, :) = extremeMasses(clusterMasses);
end
end

function [pValues, tStatistics] = interactionStatistic(input2D, designFlag)
if designFlag == 1
    interactionDifference = (input2D{1} - input2D{2}) - (input2D{3} - input2D{4});
    [~, pValues, ~, stats] = ttest(interactionDifference');
    tStatistics = stats.tstat;
    effect = mean(interactionDifference, 2)';
    [pValues, tStatistics] = normalizeDegenerateStatistics(pValues, tStatistics, effect);
elseif designFlag == 2
    group1Difference = input2D{1} - input2D{2};
    group2Difference = input2D{3} - input2D{4};
    [~, pValues, ~, stats] = ttest2( ...
        group1Difference', group2Difference', 'Vartype', 'unequal');
    tStatistics = stats.tstat;
    effect = mean(group1Difference, 2)' - mean(group2Difference, 2)';
    [pValues, tStatistics] = normalizeDegenerateStatistics(pValues, tStatistics, effect);
else
    [pValues, tStatistics] = welchInteraction(input2D);
end
end

function [pValues, tStatistics] = welchInteraction(input2D)
means = cellfun(@(x) mean(x, 2), input2D, 'UniformOutput', false);
variances = cellfun(@(x) var(x, 0, 2), input2D, 'UniformOutput', false);
sampleSizes = cellfun(@(x) size(x, 2), input2D);

contrast = means{1} - means{2} - means{3} + means{4};
varianceTerms = cell(1, 4);
for cellIdx = 1:4
    varianceTerms{cellIdx} = variances{cellIdx} ./ sampleSizes(cellIdx);
end
standardErrorSquared = varianceTerms{1} + varianceTerms{2} + ...
    varianceTerms{3} + varianceTerms{4};
tStatistics = contrast ./ sqrt(standardErrorSquared);

degreesDenominator = zeros(size(standardErrorSquared));
for cellIdx = 1:4
    degreesDenominator = degreesDenominator + ...
        varianceTerms{cellIdx}.^2 ./ (sampleSizes(cellIdx) - 1);
end
degreesOfFreedom = standardErrorSquared.^2 ./ degreesDenominator;
pValues = 2 .* tcdf(-abs(tStatistics), degreesOfFreedom);

zeroVarianceMask = standardErrorSquared == 0;
tStatistics(zeroVarianceMask & contrast == 0) = 0;
pValues(zeroVarianceMask & contrast == 0) = 1;
tStatistics(zeroVarianceMask & contrast ~= 0) = sign(contrast(zeroVarianceMask & contrast ~= 0)) .* Inf;
pValues(zeroVarianceMask & contrast ~= 0) = 0;
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

function validateInputs(inputs, designFlag, pThreshold, numIterations, randomSeed)
spatialSize = [size(inputs{1}, 1), size(inputs{1}, 2)];
sampleSizes = zeros(1, 4);
for inputIdx = 1:4
    currentInput = inputs{inputIdx};
    if ~isnumeric(currentInput) || ~isreal(currentInput) || isempty(currentInput)
        error('groupSIFT:InvalidInput', 'All inputs must be non-empty real numeric arrays.');
    end
    if any(~isfinite(currentInput(:)))
        error('groupSIFT:NonfiniteInput', 'Inputs must not contain NaN or Inf values.');
    end
    if size(currentInput, 1) ~= spatialSize(1) || size(currentInput, 2) ~= spatialSize(2)
        error('groupSIFT:SpatialSizeMismatch', ...
            'The first two dimensions of all four inputs must match.');
    end
    sampleSizes(inputIdx) = size(currentInput, 3);
    if sampleSizes(inputIdx) < 2
        error('groupSIFT:InsufficientSampleSize', ...
            'Every cell must contain at least two subjects.');
    end
end
if ~isscalar(designFlag) || ~ismember(designFlag, [1, 2, 3])
    error('groupSIFT:InvalidDesignFlag', ...
        'repeatedMeasuresFlag must be 1, 2, or 3.');
end
if designFlag == 1 && any(sampleSizes ~= sampleSizes(1))
    error('groupSIFT:RepeatedSampleSizeMismatch', ...
        'A fully repeated design requires equal subject counts in all four cells.');
end
if designFlag == 2 && (sampleSizes(1) ~= sampleSizes(2) || sampleSizes(3) ~= sampleSizes(4))
    error('groupSIFT:MixedSampleSizeMismatch', ...
        'A mixed design requires equal subject counts within cells 1/2 and 3/4.');
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
