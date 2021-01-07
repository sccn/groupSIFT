% 01/26/2016 Makoto. Created.

% define physiologically valid ROI agreed by Scott and Makoto in Dec 2015
% excludeRoiIdx = [21 22 37:42 71:78]; % 29 30 are insula... better to include?
% includeRoiIdx = setdiff(1:88, excludeRoiIdx);
% 
% includeRoiIdx = [37 87]; % Hippocampus_L + Temporal_Inf_L *** This looks a good combination
% excludeRoiIdx = setdiff(1:88, includeRoiIdx);
% 
% includeRoiIdx = [39 87]; % Prahippocampus_L + Temporal_Inf_L *** This looks a good combination
% excludeRoiIdx = setdiff(1:88, includeRoiIdx);
% 
% includeRoiIdx = [33 35]; % Cingulum_Mid_L + Cingulum_Post_L
% excludeRoiIdx = setdiff(1:88, includeRoiIdx);
% 
% includeRoiIdx = [25 27]; % Frontal_Med_Orbital_L + Rectus_L *** This looks a good combination
% excludeRoiIdx = setdiff(1:88, includeRoiIdx);
% 
% includeRoiIdx = [33 71]; % Cingulum_Mid_L + Caudate_L
% excludeRoiIdx = setdiff(1:88, includeRoiIdx);
% 
% includeRoiIdx = [31 71]; % Cingulum_Ant_L + Caudate_L *** This looks a good combination
% excludeRoiIdx = setdiff(1:88, includeRoiIdx);
% 
% includeRoiIdx = [33 77]; % Cingulum_Mid_L + Thalamus_L
% excludeRoiIdx = setdiff(1:88, includeRoiIdx);
% 
% includeRoiIdx = [21 22]; % Ocfactory
% excludeRoiIdx = setdiff(1:88, includeRoiIdx);
% 
% includeRoiIdx = [35 67]; % Cingulum_Post_L + Precuneus_L
% excludeRoiIdx = setdiff(1:88, includeRoiIdx);

includeRoiIdx = 1:88; % All brain regions
excludeRoiIdx = setdiff(1:88, includeRoiIdx);

% obtain ROI labels
roiLabels = pr.regionOfInterestFromAnatomy.getAllAnatomicalLabels;
includeRoiLabels = roiLabels(includeRoiIdx);
excludeRoiLabels = roiLabels(excludeRoiIdx);

% visualize color-labeled including regions
customColors = rand(length(includeRoiLabels),3);
figure;
set(gcf, 'numberTitle', 'off', 'name', 'Anatomical ROIs to include')
plot_dipplot_with_cortex;
for n = 1:length(includeRoiLabels)
    tmpRoi = pr.regionOfInterestFromAnatomy(pr.headGrid, includeRoiLabels{n});
    %tmp = tmpRoi.membershipCube;
    pr.plot_head_surface(tmpRoi.headGrid,...
                         tmpRoi.membershipCube,...
                         'showProjectedOnMrs', 0,...
                         'surfaceColor', customColors(n,:),...
                         'surfaceOptions', {'facealpha', 1});
end

azimuthValue   = 180;
elevationValue = 0;

delete(findobj(gca, 'type','light'));
camlight('left');
camlight('right');
camlight('headlight');
camlight('headlight', 'infinite');
camlight(azimuthValue, elevationValue)

view([40 34]); % default
view([90 0]);
view([80 10]);
view([70 20]);
view([60 30]);
view([50 40]);

view([32 7]);

delete(findall(gcf, 'Tag', 'img'))
camlight(-90,0)
camlight(-90,0)
camlight(-90,0)



%% visualize color-labeled excluding regions
excludeRoiIdx = [21 22 37:42 71:78];

excludeRoiIdx = [21 22 71:78]; % Upper basal: Olfactory, Caudate, Putamen, Pallidum, and Thalamus

excludeRoiIdx = [37:42]; % Lower basal: Hippocampus, Parahippocampus, Amygdala

includeRoiIdx = setdiff(1:88, excludeRoiIdx);

roiLabels = pr.regionOfInterestFromAnatomy.getAllAnatomicalLabels;
includeRoiLabels = roiLabels(includeRoiIdx);
excludeRoiLabels = roiLabels(excludeRoiIdx);

customColors = rand(length(excludeRoiLabels),3);
figure;
set(gcf, 'numberTitle', 'off', 'name', 'Anatomical ROIs to exclude')
plot_dipplot_with_cortex;
for n = 2:2:length(excludeRoiLabels)
    tmpRoi = pr.regionOfInterestFromAnatomy(pr.headGrid, excludeRoiLabels{n});
    pr.plot_head_surface(tmpRoi.headGrid,...
                         tmpRoi.membershipCube,...
                         'showProjectedOnMrs', 0,...
                         'surfaceColor', customColors(n,:),...
                         'surfaceOptions', {'facealpha', 1});
end
delete(findall(gcf, 'Tag', 'img'))
camlight(-90,0)
camlight(-90,0)
camlight(-90,0)

%% Compute ROI volues
numberOfRegionsOfInterest = length(roiLabels);
roiCentroids = zeros(length(roiLabels),3);
firstROI = pr.regionOfInterestFromAnatomy(pr.headGrid, roiLabels{1});
voxelSizeInCm = firstROI.headGrid.spacing/10; % mm
roiVolumesInCc = zeros(numberOfRegionsOfInterest,1);
numberedRoiLabels = cell(numberOfRegionsOfInterest,1);
for i = 1:numberOfRegionsOfInterest
    disp(sprintf('%.0f/%.0f ROI', i, numberOfRegionsOfInterest));
    regionOfInterest = pr.regionOfInterestFromAnatomy(pr.headGrid, roiLabels{i});
    roiVolumesInCc(i,1) = sum(regionOfInterest.membershipCube(:)).*voxelSizeInCm^3;
    numberedRoiLabels{i} = [sprintf('%2.0f. ',i) roiLabels{i}];
end

% Mean volume: 14.3 cc (SD 9.3)
mean(roiVolumesInCc)
std(roiVolumesInCc)

% Total volume: 1258 cc.
% round(sum(roiVolumesInCc))

% Total excluded: 92 cc (7.31 %)
% round(sum(roiVolumesInCc([21 22 37:42 71:78])))

% Upper basal: 57 cc
% round(sum(roiVolumesInCc([21 22 71:78])))

% Lower basal: 35 cc
% round(sum(roiVolumesInCc([37:42])))


figure; set(gcf, 'color', [0.66 0.76 1])
bar(roiVolumesInCc)
xlim([0.5 length(numberedRoiLabels)+0.5])
set(gca, 'XTick', 1:size(numberedRoiLabels), 'XTickLabel', numberedRoiLabels, 'position', [0.0456 0.2007 0.9451 0.7832])
rotateXLabels(gca, 90)
set(findall(gca, '-property', 'interpreter'), 'interpreter', 'none')
set(findall(gca, '-property', 'fontsize'), 'fontsize', 12)
set(get(gca, 'ylabel'), 'string', 'Volume of ROI (cc)', 'fontsize', 18)
