% 07/03/2020 Makoto. Created.

addpath('/data/mobi/Daisuke/code')

% Determine p-value for cluster-level correction.
clusterLevelPvalue = 0.0001;
latencyWindowLimit = [];
frequencyBandLimit = [30 50];

saveName = 'Gamma';

%%%%%%%%%%%%%%%%%%%%%%
%%% Load .mat data %%%
%%%%%%%%%%%%%%%%%%%%%%
load('/data/mobi/Daisuke/p5100_groupSeparated/SzCt/Sz_minus_Ct_tStatistics.mat');
load('/data/mobi/Daisuke/p5100_groupSeparated/SzCt/Sz_minus_Ct_dipolePairDensity.mat');

% store the loaded data to axis data
connectivityData = struct();
connectivityData.dimensionLabels        = dimensionLabels;
connectivityData.frequencies            = frequencies;
connectivityData.latencies              = latencies;
connectivityData.connectivityType       = connectivityType;
connectivityData.fileNameList           = fileNameList;
connectivityData.baselineIdx            = baselineIdx;
connectivityData.roiLabels              = roiLabels;
connectivityData.symmetricRoiCentroids  = symmetricRoiCentroids;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform cluster-level correction %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the critical Mass of Cluster from pooled surrogate results.
criticalMassOfCluster = prctile(surroMassOfCluster(:), [clusterLevelPvalue*100/2 100-(clusterLevelPvalue*100/2)]); % Blob sizes are nicely bounded!

% Build the mask.
significantClusterMask = zeros(size(clusterMask));
for edgeIdx = 1:size(clusterMask,1)
    
    % Identify all blobs that survive the cluster-level threthold (pooled across edges)
    tmpBlobMask    = squeeze(clusterMask(edgeIdx,:,:));
    tmpTStatistics = squeeze(tStatistics(edgeIdx,:,:));
    [entryCount, blobId]  = hist(tmpBlobMask(:), unique(tmpBlobMask(:)));
    tmpMassOfCluster = zeros(length(blobId)-1,1);
    for m = 2:length(blobId)
        currentMask = tmpBlobMask==blobId(m);
        tmpMassOfCluster(m-1) = sum(sum(currentMask.*tmpTStatistics));
    end
    survivedClusterMaskIdx = find(tmpMassOfCluster < criticalMassOfCluster(1) | tmpMassOfCluster > criticalMassOfCluster(2));
    
    % Combine the survived blobs.
    if ~isempty(survivedClusterMaskIdx)
        combinedMask = zeros([size(clusterMask,2) size(clusterMask,3)]);
        for clusterMaskIdx = 1:length(survivedClusterMaskIdx)
            currentMask = tmpBlobMask==survivedClusterMaskIdx(clusterMaskIdx);
            
            % If currentMask is a row vector, transpose it. (06/28/2020 Makoto)
            if size(currentMask,1)==1 
                currentMask = currentMask';
            end
            
            combinedMask = logical(combinedMask + currentMask);
        end
    else
        % disp('Nothing survived.')
        continue
    end
    
    % Store the result.
    significantClusterMask(edgeIdx,:,:) = combinedMask;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Overwrite imported data by remove non-significant edges %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
significantEdgeIdx     = find(sum(sum(significantClusterMask,2),3));
finallySelectedEdgeIdx = finallySelectedEdgeIdx(significantEdgeIdx);
significantClusterMask = significantClusterMask(significantEdgeIdx,:,:);
tStatistics            = tStatistics(significantEdgeIdx,:,:);
pValues                = pValues(significantEdgeIdx,:,:);
if exist('tStatistics_beforeSubtraction1', 'var')
    tStatistics_beforeSubtraction1 = tStatistics_beforeSubtraction1(significantEdgeIdx,:,:);
    tStatistics_beforeSubtraction2 = tStatistics_beforeSubtraction2(significantEdgeIdx,:,:);
end
if exist('tStatistics_beforeSubtraction3', 'var')
    tStatistics_beforeSubtraction3 = tStatistics_beforeSubtraction3(significantEdgeIdx,:,:);
    tStatistics_beforeSubtraction4 = tStatistics_beforeSubtraction4(significantEdgeIdx,:,:);
end


%%%%%%%%%%%%%%%%%%%%
%%% Make a mask. %%%
%%%%%%%%%%%%%%%%%%%%
    % % Limit latency window if necessary
    % latencyWindowMask  = logical(zeros(size(significantClusterMask)));
    % if any(latencyWindowLimit)
    %     latencyIdx = find(latencies>=latencyWindowLimit(1) & latencies<=latencyWindowLimit(2));
    %     latencyWindowMask(:,:,latencyIdx) = 1;
    %     significantClusterMask = logical(significantClusterMask.*latencyWindowMask);
    % else
    %     latencyIdx = 1:length(latencies);
    % end

% Limit frequency band if necessary
freqRangeMask = logical(zeros(size(significantClusterMask)));
if any(frequencyBandLimit)
    freqIdx = find(frequencies>=frequencyBandLimit(1) & frequencies<=frequencyBandLimit(2));
    freqRangeMask(:,freqIdx,:) = 1;
    significantClusterMask = logical(significantClusterMask.*freqRangeMask);
else
    freqIdx = 1:length(frequencies);
end

    % % Obtain positive and negative masks.
    % positiveMask = maskedT>0;
    % massOfClusterPosi = sum(sum(maskedT.*positiveMask,3),2);
    % negativeMask = maskedT<0;
    % massOfClusterNega = sum(sum(maskedT.*negativeMask,3),2);
    % 
    % % Combine the masks.
    % massOfClusterCombined = massOfClusterPosi + massOfClusterNega;

% Apply the inclusive mask --> this could be used for movie!
maskedT_diff = tStatistics.*significantClusterMask;

% Obtain the same for Sz and Ct.
maskedT_Sz = tStatistics_beforeSubtraction1.*significantClusterMask;
maskedT_Ct = tStatistics_beforeSubtraction2.*significantClusterMask;

% Compute mass of cluser (MOC)
mocDiff = sum(maskedT_diff,2);
mocSz   = sum(maskedT_Sz,  2);
mocCt   = sum(maskedT_Ct,  2);

% Skip zeros.
nonzeroIdx = find(mocDiff~=0);
mocDiff = mocDiff(nonzeroIdx);
mocSz = mocSz(nonzeroIdx);
mocCt = mocCt(nonzeroIdx);
finallySelectedEdgeIdx = finallySelectedEdgeIdx(nonzeroIdx);

% Obtain the edge label.
[toIdx, fromIdx] = ind2sub(76, finallySelectedEdgeIdx);
edgeList = cell(length(finallySelectedEdgeIdx),1);
for edgeIdx = 1:length(finallySelectedEdgeIdx)
    edgeList{edgeIdx} = sprintf('From %s to %s', roiLabels{fromIdx(edgeIdx)}, roiLabels{toIdx(edgeIdx)});
end
    
% Sort the results.
[mocDiffSorted, sortingIdx] = sort(mocDiff, 'descend');
mocSzSorted = mocSz(sortingIdx);
mocCtSorted = mocCt(sortingIdx);
edgeListSorted = edgeList(sortingIdx);

% Rename ROI labels. 12/13/2019 Makoto.
for labelIdx = 1:length(edgeListSorted)
    currentLabel = edgeListSorted{labelIdx};
    underScoreIdx = strfind(currentLabel, '_');
    currentLabel(underScoreIdx) = deal('-');
    edgeListSorted{labelIdx} = currentLabel;
end

% % Plot figures.
% figure
% bar(mocDiffSorted)
% set(gca, 'xtick', 1:length(mocDiffSorted), 'xticklabel', edgeListSorted)
% rotateXLabels(gca, 90)

% Create a table.

% Prepare the table format.
tableContents = table(mocSzSorted, mocCtSorted, mocDiffSorted, ...
                      'VariableNames', {'Sz' 'Ct' 'SzCt'}, ...
                      'RowNames', edgeListSorted);

% Save as .xlsx.
writetable(tableContents, ['/data/mobi/Daisuke/p5110_extractDataFromConnMatrix/' saveName '.xlsx'], 'FileType','spreadsheet', 'WriteRowNames', true);
