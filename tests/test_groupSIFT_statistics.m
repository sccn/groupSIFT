function tests = test_groupSIFT_statistics
tests = functiontests(localfunctions);
end


function testSubjectMatchingReordersByFileName(testCase)
list1 = {'/conditionA/S01.set'; '/conditionA/S02.set'; '/conditionA/S03.set'};
list2 = {'/conditionB/s03.set'; '/conditionB/s01.set'; '/conditionB/s02.set'};

[subjectIds, indexByList, unmatchedIds] = groupSIFT_matchSubjects({list1, list2});

verifyEqual(testCase, subjectIds, {'S01'; 'S02'; 'S03'});
verifyEqual(testCase, indexByList{1}, [1; 2; 3]);
verifyEqual(testCase, indexByList{2}, [2; 3; 1]);
verifyEmpty(testCase, unmatchedIds{1});
verifyEmpty(testCase, unmatchedIds{2});
end


function testSubjectMatchingRejectsDuplicates(testCase)
list1 = {'S01.set'; 's01.SET'};
verifyError(testCase, @() groupSIFT_matchSubjects({list1}), ...
    'groupSIFT:DuplicateSubjectId');
end


function testPairedPermutationUsesSurrogateTStatistics(testCase)
input1 = reshape([1.2, 2.1, 2.9, 4.4, 5.2, 6.3], [1, 1, 6]);
input2 = zeros(1, 1, 6);
numberOfIterations = 20;
seed = 2501;

rngBefore = rng;
[~, observedT, observedP, surrogateMass] = clusterLevelPermutationTest( ...
    input1, input2, 1, 0.999, numberOfIterations, seed);
verifyEqual(testCase, rng, rngBefore);

[~, expectedP, ~, expectedStats] = ttest(reshape(input1, [], 1));
verifyEqual(testCase, observedT, expectedStats.tstat, 'AbsTol', 1e-12);
verifyEqual(testCase, observedP, expectedP, 'AbsTol', 1e-12);

pairedDifference = reshape(input1 - input2, 1, []);
expectedMass = zeros(numberOfIterations, 2);
rng(seed, 'twister');
for iterationIdx = 1:numberOfIterations
    randomSigns = ones(1, numel(pairedDifference));
    randomSigns(rand(1, numel(pairedDifference)) < 0.5) = -1;
    [~, surrogateP, ~, surrogateStats] = ttest((pairedDifference .* randomSigns)');
    if surrogateP < 0.999
        if surrogateStats.tstat < 0
            expectedMass(iterationIdx, 1) = surrogateStats.tstat;
        elseif surrogateStats.tstat > 0
            expectedMass(iterationIdx, 2) = surrogateStats.tstat;
        end
    end
end
rng(rngBefore);
verifyEqual(testCase, surrogateMass, expectedMass, 'AbsTol', 1e-12);
end


function testPositiveAndNegativeClustersAreSeparate(testCase)
positiveDifference = [2.1, 2.4, 2.9, 3.2, 3.8, 4.1, 4.7, 5.0];
negativeDifference = -[2.0, 2.5, 2.7, 3.4, 3.6, 4.2, 4.4, 5.1];
input1 = zeros(1, 2, numel(positiveDifference));
input1(1, 1, :) = positiveDifference;
input1(1, 2, :) = negativeDifference;
input2 = zeros(size(input1));

[mask, tScore] = clusterLevelPermutationTest(input1, input2, 1, 0.05, 2, 17);

verifyGreaterThan(testCase, tScore(1), 0);
verifyLessThan(testCase, tScore(2), 0);
verifyNotEqual(testCase, mask(1), 0);
verifyNotEqual(testCase, mask(2), 0);
verifyNotEqual(testCase, mask(1), mask(2));
end


function testIndependentTwoByTwoUsesWelchContrast(testCase)
input1 = reshape([4.0, 4.4, 5.1, 5.3, 6.2], [1, 1, 5]);
input2 = reshape([1.0, 1.4, 2.2, 2.5, 3.0, 3.4], [1, 1, 6]);
input3 = reshape([3.1, 3.5, 3.7, 4.2], [1, 1, 4]);
input4 = reshape([2.0, 2.4, 2.8, 3.2, 3.6, 3.9, 4.1], [1, 1, 7]);

[~, observedT, observedP] = clusterLevelPermutationTest2x2( ...
    input1, input2, input3, input4, 3, 0.05, 2, 31);

groups = {reshape(input1, 1, []), reshape(input2, 1, []), ...
    reshape(input3, 1, []), reshape(input4, 1, [])};
sampleSizes = cellfun(@numel, groups);
means = cellfun(@mean, groups);
variances = cellfun(@(x) var(x, 0), groups);
varianceTerms = variances ./ sampleSizes;
contrast = means(1) - means(2) - means(3) + means(4);
standardErrorSquared = sum(varianceTerms);
expectedT = contrast / sqrt(standardErrorSquared);
degreesOfFreedom = standardErrorSquared^2 / ...
    sum(varianceTerms.^2 ./ (sampleSizes - 1));
expectedP = 2 * tcdf(-abs(expectedT), degreesOfFreedom);

verifyEqual(testCase, observedT, expectedT, 'AbsTol', 1e-12);
verifyEqual(testCase, observedP, expectedP, 'AbsTol', 1e-12);
end


function testRandomSeedIsReproducible(testCase)
input1 = reshape(1:12, [1, 2, 6]);
input2 = zeros(size(input1));
[~, ~, ~, firstMass] = clusterLevelPermutationTest(input1, input2, 1, 0.5, 25, 99);
[~, ~, ~, secondMass] = clusterLevelPermutationTest(input1, input2, 1, 0.5, 25, 99);
verifyEqual(testCase, firstMass, secondMass);
end


function testExplicitPairedRandomizationPlanIsUsedExactly(testCase)
input1 = reshape([1.1, 1.8, 2.7, 4.2], [1, 1, 4]);
input2 = zeros(size(input1));
plan.signs = [1, 1, 1, 1; -1, 1, -1, 1; 1, -1, -1, 1];

rngBefore = rng;
[~, ~, ~, surrogateMass] = clusterLevelPermutationTest( ...
    input1, input2, 1, 0.999, size(plan.signs, 1), plan);
verifyEqual(testCase, rng, rngBefore);

pairedDifference = reshape(input1 - input2, 1, []);
expectedMass = zeros(size(plan.signs, 1), 2);
for iterationIdx = 1:size(plan.signs, 1)
    surrogateDifference = pairedDifference .* plan.signs(iterationIdx, :);
    [~, surrogateP, ~, surrogateStats] = ttest(surrogateDifference');
    if surrogateP < 0.999
        if surrogateStats.tstat < 0
            expectedMass(iterationIdx, 1) = surrogateStats.tstat;
        elseif surrogateStats.tstat > 0
            expectedMass(iterationIdx, 2) = surrogateStats.tstat;
        end
    end
end
verifyEqual(testCase, surrogateMass, expectedMass, 'AbsTol', 1e-12);
end


function testMixedTwoByTwoAcceptsSharedPermutationPlan(testCase)
input1 = reshape([4.2, 5.1, 6.0, 7.3], [1, 1, 4]);
input2 = reshape([1.1, 1.7, 2.3, 2.8], [1, 1, 4]);
input3 = reshape([3.0, 3.8, 4.9], [1, 1, 3]);
input4 = reshape([1.2, 1.5, 2.1], [1, 1, 3]);
plan.uniformScores = [0.7, 0.1, 0.4, 0.2, 0.8, 0.3, 0.6; ...
                      0.2, 0.9, 0.1, 0.6, 0.5, 0.7, 0.3];

rngBefore = rng;
restoreRng = onCleanup(@() rng(rngBefore)); %#ok<NASGU>
rng(11, 'twister');
[~, ~, ~, firstMass] = clusterLevelPermutationTest2x2( ...
    input1, input2, input3, input4, 2, 0.999, 2, plan);
rng(987, 'twister');
[~, ~, ~, secondMass] = clusterLevelPermutationTest2x2( ...
    input1, input2, input3, input4, 2, 0.999, 2, plan);
verifyEqual(testCase, firstMass, secondMass, 'AbsTol', 1e-12);
end


function testRandomizationPlanDimensionsAreValidated(testCase)
input1 = reshape(1:4, [1, 1, 4]);
input2 = zeros(size(input1));
plan.signs = ones(2, 3);
verifyError(testCase, @() clusterLevelPermutationTest( ...
    input1, input2, 1, 0.05, 2, plan), ...
    'groupSIFT:InvalidRandomizationPlan');
end


function testDegenerateNullDataReturnFiniteNeutralStatistics(testCase)
input1 = zeros(2, 3, 5);
input2 = zeros(2, 3, 5);
[mask, tScore, pValue, surrogateMass] = clusterLevelPermutationTest( ...
    input1, input2, 1, 0.05, 4, 7);
verifyEqual(testCase, mask, zeros(2, 3));
verifyEqual(testCase, tScore, zeros(2, 3));
verifyEqual(testCase, pValue, ones(2, 3));
verifyEqual(testCase, surrogateMass, zeros(4, 2));
end


function testPairedSampleSizeMismatchIsRejected(testCase)
input1 = zeros(2, 2, 4);
input2 = zeros(2, 2, 5);
verifyError(testCase, @() clusterLevelPermutationTest(input1, input2, 1, 0.05, 10), ...
    'groupSIFT:PairedSampleSizeMismatch');
end
