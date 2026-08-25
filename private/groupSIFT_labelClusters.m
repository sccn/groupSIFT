function [labelMatrix, clusterMasses] = groupSIFT_labelClusters(pValueMatrix, tScoreMatrix, pThreshold)
% groupSIFT_labelClusters() - Label positive and negative clusters separately.
%
% Positive and negative samples must not be joined into the same cluster in
% a two-sided cluster-mass test. Labels are positive integers; positive
% clusters are numbered first, followed by negative clusters.

positiveLabels = bwlabeln((pValueMatrix < pThreshold) & (tScoreMatrix > 0));
negativeLabels = bwlabeln((pValueMatrix < pThreshold) & (tScoreMatrix < 0));

numberOfPositiveClusters = max(positiveLabels(:));
numberOfNegativeClusters = max(negativeLabels(:));

labelMatrix = positiveLabels;
negativeMask = negativeLabels > 0;
labelMatrix(negativeMask) = negativeLabels(negativeMask) + numberOfPositiveClusters;

clusterMasses = zeros(numberOfPositiveClusters + numberOfNegativeClusters, 1);
for clusterIdx = 1:numberOfPositiveClusters
    currentMask = positiveLabels == clusterIdx;
    clusterMasses(clusterIdx) = sum(tScoreMatrix(currentMask));
end
for clusterIdx = 1:numberOfNegativeClusters
    currentMask = negativeLabels == clusterIdx;
    outputIdx = numberOfPositiveClusters + clusterIdx;
    clusterMasses(outputIdx) = sum(tScoreMatrix(currentMask));
end
end
