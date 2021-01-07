% pop_groupSIFT_convertToGroupAnatomicalRois()
%
% History: 12/09/2019 Makoto. selectDipWithLargerMoment() is supported.

% Copyright (C) 2016, Makoto Miyakoshi (mmiyakoshi@ucsd.edu) , SCCN,INC,UCSD
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

function varargout = pop_groupSIFT_convertToGroupAnatomicalRois(varargin)
% POP_GROUPSIFT_CONVERTTOGROUPANATOMICALROIS MATLAB code for pop_groupSIFT_convertToGroupAnatomicalRois.fig
%      POP_GROUPSIFT_CONVERTTOGROUPANATOMICALROIS, by itself, creates a new POP_GROUPSIFT_CONVERTTOGROUPANATOMICALROIS or raises the existing
%      singleton*.
%
%      H = POP_GROUPSIFT_CONVERTTOGROUPANATOMICALROIS returns the handle to a new POP_GROUPSIFT_CONVERTTOGROUPANATOMICALROIS or the handle to
%      the existing singleton*.
%
%      POP_GROUPSIFT_CONVERTTOGROUPANATOMICALROIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_GROUPSIFT_CONVERTTOGROUPANATOMICALROIS.M with the given input arguments.
%
%      POP_GROUPSIFT_CONVERTTOGROUPANATOMICALROIS('Property','Value',...) creates a new POP_GROUPSIFT_CONVERTTOGROUPANATOMICALROIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_groupSIFT_convertToGroupAnatomicalRois_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_groupSIFT_convertToGroupAnatomicalRois_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_groupSIFT_convertToGroupAnatomicalRois

% Last Modified by GUIDE v2.5 26-Jun-2020 17:10:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_groupSIFT_convertToGroupAnatomicalRois_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_groupSIFT_convertToGroupAnatomicalRois_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pop_groupSIFT_convertToGroupAnatomicalRois is made visible.
function pop_groupSIFT_convertToGroupAnatomicalRois_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_groupSIFT_convertToGroupAnatomicalRois (see VARARGIN)

% Choose default command line output for pop_groupSIFT_convertToGroupAnatomicalRois
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_groupSIFT_convertToGroupAnatomicalRois wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pop_groupSIFT_convertToGroupAnatomicalRois_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




function gaussianSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to gaussianSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gaussianSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of gaussianSizeEdit as a double


% --- Executes during object creation, after setting all properties.
function gaussianSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gaussianSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fileNameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to fileNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileNameEdit as text
%        str2double(get(hObject,'String')) returns contents of fileNameEdit as a double


% --- Executes during object creation, after setting all properties.
function fileNameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in selectFilesButton.
function selectFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain user input file name
userInputFilename = get(handles.fileNameEdit, 'String');
if isempty(userInputFilename)
    error('Enter file name.')
end

% Obtain multiple .set files
[allFiles, workingFolder] = uigetfile('*.set', 'MultiSelect', 'on');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end

% Display process start
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''3.Convert to group anatomical ROIs'' started. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))

% Move to the working folder
cd(workingFolder)

% Load empty dipolePairAndMeasure from eeglabroot/plugins/groupSIFT/+pr: Thanks Clement!
dipolePairAndMeasureObj = pr.dipolePairAndMeasure;

% Load all .set data located under the working folder
ALLEEG = [];
for n = 1:length(allFiles)
    loadName = allFiles{n};
    EEG = pop_loadset('filename', loadName,'filepath',workingFolder, 'loadmode', 'info');
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
end

% Compute how many IC combinations in all subjects
numICs = zeros(length(ALLEEG),1);
for n = 1:length(ALLEEG)
    numICs(n) = size(ALLEEG(1,n).icaweights,1);
end
numIcCombination = sum(numICs.^2);

% Assign unique IDs to dipoles
if get(handles.rpdcButton,'value')==1
    timeFreqSize = size(squeeze(ALLEEG(1).CAT.Conn.RPDC(  1,1,:,:)));
else
    timeFreqSize = size(squeeze(ALLEEG(1).CAT.Conn.dDTF08(1,1,:,:)));
end
latencies = ALLEEG(1).CAT.Conn.erWinCenterTimes;
frequencies = ALLEEG(1).CAT.Conn.freqs;
dimensionLabels= ALLEEG(1).CAT.Conn.dims;
dimensionLabels{1,5} = 'subjects';
dipolePairAndMeasureObj.linearizedMeasure = zeros([numIcCombination timeFreqSize(1)*timeFreqSize(2)]);
dipoleCounter = 0;
counter = 1;
for n = 1:length(ALLEEG)
    numberOfSessiondipoles = length(ALLEEG(n).dipfit.model);
    dipoleLocation = [];
    dipoleResidualVariance = [];
    for m = 1:numberOfSessiondipoles
        %dipoleLocation(m,:) = ALLEEG(n).dipfit.model(m).posxyz(1,:); % 12/09/2019 Makoto.
        dipoleLocation(m,:) = selectDipWithLargerMoment(ALLEEG(n).dipfit.model(m).posxyz, ALLEEG(n).dipfit.model(m).momxyz);
        dipoleResidualVariance(m) = ALLEEG(n).dipfit.model(m).rv;
    end;
    dipoleId = (dipoleCounter+1):(dipoleCounter+length(dipoleLocation));
    dipoleCounter = dipoleCounter + length(dipoleLocation);
    
    for m = 1:length(dipoleLocation) % from location: confirm it with ALLEEG(1,1).CAT.Conn: dims: {'var_to'  'var_from'  'freq'  'time'}
        for k = 1:length(dipoleLocation) % to location
           %dipolePairAndMeasureObj.from.location = cat(1, dipolePairAndMeasureObj.from.location, ALLEEG(n).dipfit.model(m).posxyz(1,:));
            dipolePairAndMeasureObj.from.location = cat(1, dipolePairAndMeasureObj.from.location, selectDipWithLargerMoment(ALLEEG(n).dipfit.model(m).posxyz, ALLEEG(n).dipfit.model(m).momxyz));
            dipolePairAndMeasureObj.from.residualVariance(end+1,1) = ALLEEG(n).dipfit.model(m).rv;
           %dipolePairAndMeasureObj.to.location = cat(1, dipolePairAndMeasureObj.to.location, ALLEEG(n).dipfit.model(k).posxyz(1,:));
            dipolePairAndMeasureObj.to.location = cat(1, dipolePairAndMeasureObj.to.location, selectDipWithLargerMoment(ALLEEG(n).dipfit.model(m).posxyz, ALLEEG(n).dipfit.model(m).momxyz));
            if m == k
                dipolePairAndMeasureObj.linearizedMeasure(counter,:) = zeros(timeFreqSize(1)*timeFreqSize(2),1);
            else
                if get(handles.rpdcButton,'value')==1 % Again, ALLEEG(1,1).CAT.Conn: dims: {'var_to'  'var_from'  'freq'  'time'}
                    tmpTimeFreqMatrix = squeeze(ALLEEG(n).CAT.Conn.RPDC(  k,m,:,:));
                else
                    tmpTimeFreqMatrix = squeeze(ALLEEG(n).CAT.Conn.dDTF08(k,m,:,:));
                end
                dipolePairAndMeasureObj.linearizedMeasure(counter,:) = tmpTimeFreqMatrix(:);
            end
            dipolePairAndMeasureObj.sessionId(counter,1) = n;
            dipolePairAndMeasureObj.fromDipoleId(counter,1) = dipoleId(m);
            dipolePairAndMeasureObj.toDipoleId(counter,1)   = dipoleId(k);
            counter = counter + 1;
        end
    end
end

% Find unique dipoles
uniqueDipole = pr.dipole;
[~,fromId,fromIdReverse] = unique(dipolePairAndMeasureObj.fromDipoleId, 'stable');

% Extract coordinates + residual variance of the unique dipoles
uniqueDipole.location = dipolePairAndMeasureObj.from.location(fromId,:);
uniqueDipole.residualVariance = dipolePairAndMeasureObj.from.residualVariance(fromId);

% Load headGrid cubes
headGrid = pr.headGrid;

% Set Gaussian smoothing kernel
    % FWHM = 2.355*sigma See https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    % 4.2  (FWHM==10mm, 12/23/2015)
    % 8.5  (FWHM==20mm, 06/24/2015)
    % 12.8 (FWHM==30mm, 09/05/2014)
    % 17.1 (FWHM==40mm, 01/06/2015)
    % Note that Gaussian is NOT truncated until reaching 150mm apart from the center
userInputKernelSize = str2num(get(handles.gaussianSizeEdit, 'String'));
standardDeviationOfEstimatedDipoleLocation = userInputKernelSize/2.355; % this calculates sigma in Gaussian equation
projectionParameter = pr.projectionParameter(standardDeviationOfEstimatedDipoleLocation);
%projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian = 150/standardDeviationOfEstimatedDipoleLocation;
projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian = 3;
[~,~,gaussianWeightMatrix] = pr.meanProjection.getProjectionMatrix(uniqueDipole, headGrid, projectionParameter);

% Define valid anatomical ROIs as EEG sources agreed by Scott and Makoto in Dec 2015; See below for the list
% excludeRoiIdx = [21 22 37:42 71:78]; % 29 30 are insula... better to include?
excludeRoiIdx = [];
includeRoiIdx = setdiff(1:88, excludeRoiIdx);
roiLabels = pr.regionOfInterestFromAnatomy.getAllAnatomicalLabels;

        % % Compute ROI volues
        % numberOfRegionsOfInterest = length(roiLabels);
        % dipoleProbabilityInAllRegion = zeros(uniqueDipole.numberOfDipoles, numberOfRegionsOfInterest);
        % roiCentroids = zeros(length(roiLabels),3);
        % firstROI = pr.regionOfInterestFromAnatomy(pr.headGrid, roiLabels{1});
        % voxelSizeInCm = firstROI.headGrid.spacing/10; % mm
        % roiVolumesInCc = zeros(numberOfRegionsOfInterest,1);
        % numberedRoiLabels = cell(numberOfRegionsOfInterest,1);
        % for i = 1:numberOfRegionsOfInterest
        %     disp(sprintf('%.0f/%.0f ROI', i, numberOfRegionsOfInterest));
        %     regionOfInterest = pr.regionOfInterestFromAnatomy(pr.headGrid, roiLabels{i});
        %     roiVolumesInCc(i,1) = sum(regionOfInterest.membershipCube(:)).*voxelSizeInCm^3;
        %     numberedRoiLabels{i} = [sprintf('%2.0f. ',i) roiLabels{i}];
        % end
        % figure; set(gcf, 'color', [0.66 0.76 1])
        % bar(roiVolumesInCc)
        % xlim([0.5 length(numberedRoiLabels)+0.5])
        % set(gca, 'XTick', 1:size(numberedRoiLabels), 'XTickLabel', numberedRoiLabels, 'position', [0.0456 0.2007 0.9451 0.7832])
        % rotateXLabels(gca, 90)
        % set(findall(gca, '-property', 'interpreter'), 'interpreter', 'none')
        % set(findall(gca, '-property', 'fontsize'), 'fontsize', 12)
        % set(get(gca, 'ylabel'), 'string', 'Volume of ROI (cc)', 'fontsize', 18

% Include and exclude ROIs
includedRoiLabels = roiLabels(includeRoiIdx);
excludedRoiLabels = roiLabels(excludeRoiIdx);

% These regions are to be included
%     'Precentral_L'
%     'Precentral_R'
%     'Frontal_Sup_L'
%     'Frontal_Sup_R'
%     'Frontal_Sup_Orb_L'
%     'Frontal_Sup_Orb_R'
%     'Frontal_Mid_L'
%     'Frontal_Mid_R'
%     'Frontal_Mid_Orb_L'
%     'Frontal_Mid_Orb_R'
%     'Frontal_Inf_Oper_L'
%     'Frontal_Inf_Oper_R'
%     'Frontal_Inf_Tri_L'
%     'Frontal_Inf_Tri_R'
%     'Frontal_Inf_Orb_L'
%     'Frontal_Inf_Orb_R'
%     'Rolandic_Oper_L'
%     'Rolandic_Oper_R'
%     'Supp_Motor_Area_L'
%     'Supp_Motor_Area_R'
%     'Frontal_Sup_Medial_L'
%     'Frontal_Sup_Medial_R'
%     'Frontal_Med_Orb_L'
%     'Frontal_Med_Orb_R'
%     'Rectus_L'
%     'Rectus_R'
%     'Insula_L'
%     'Insula_R'
%     'Cingulum_Ant_L'
%     'Cingulum_Ant_R'
%     'Cingulum_Mid_L'
%     'Cingulum_Mid_R'
%     'Cingulum_Post_L'
%     'Cingulum_Post_R'
%     'Calcarine_L'
%     'Calcarine_R'
%     'Cuneus_L'
%     'Cuneus_R'
%     'Lingual_L'
%     'Lingual_R'
%     'Occipital_Sup_L'
%     'Occipital_Sup_R'
%     'Occipital_Mid_L'
%     'Occipital_Mid_R'
%     'Occipital_Inf_L'
%     'Occipital_Inf_R'
%     'Fusiform_L'
%     'Fusiform_R'
%     'Postcentral_L'
%     'Postcentral_R'
%     'Parietal_Sup_L'
%     'Parietal_Sup_R'
%     'Parietal_Inf_L'
%     'Parietal_Inf_R'
%     'SupraMarginal_L'
%     'SupraMarginal_R'
%     'Angular_L'
%     'Angular_R'
%     'Precuneus_L'
%     'Precuneus_R'
%     'Paracentral_Lobule_L'
%     'Paracentral_Lobule_R'
%     'Temporal_Sup_L'
%     'Temporal_Sup_R'
%     'Temporal_Pole_Sup_L'
%     'Temporal_Pole_Sup_R'
%     'Temporal_Mid_L'
%     'Temporal_Mid_R'
%     'Temporal_Pole_Mid_L'
%     'Temporal_Pole_Mid_R'
%     'Temporal_Inf_L'
%     'Temporal_Inf_R'
%
% These regions are to be combined
%     'Hippocampus_L'
%     'Hippocampus_R'
%     'ParaHippocampal_L'
%     'ParaHippocampal_R'
%     'Amygdala_L'
%     'Amygdala_R'
%                --> Lower Basal
%
%     'Olfactory_L'
%     'Olfactory_R'
%     'Caudate_L'
%     'Caudate_R'
%     'Putamen_L'
%     'Putamen_R'
%     'Pallidum_L'
%     'Pallidum_R'
%     'Thalamus_L'
%     'Thalamus_R'
%               --> Upper Basal
%
% One can visualize these regions by running visualizeAnatomicalRoiWithNHimasBlobs.m contained by the groupSIFT folder.

% Define Upper and Lower Basal ROIs
upperBasalLIdx = [21 71:2:77]; % Upper basal: Olfactory, Caudate, Putamen, Pallidum, and Thalamus
upperBasalRIdx = [22 72:2:78]; % Upper basal: Olfactory, Caudate, Putamen, Pallidum, and Thalamus
lowerBasalLIdx = [37:2:41];    % Lower basal: Hippocampus, Parahippocampus, Amygdala
lowerBasalRIdx = [38:2:42];    % Lower basal: Hippocampus, Parahippocampus, Amygdala

% Obtain dipole density and centroids in each ROI
numberOfRegionsOfInterest = length(includedRoiLabels);
dipoleProbabilityInRegion = zeros(uniqueDipole.numberOfDipoles, numberOfRegionsOfInterest);
roiCentroids = zeros(length(includedRoiLabels),3);
tmpROI = pr.regionOfInterestFromAnatomy(pr.headGrid, includedRoiLabels{1});
upperBasalLMembershipCube = zeros(size(tmpROI.membershipCube));
upperBasalRMembershipCube = zeros(size(tmpROI.membershipCube));
lowerBasalLMembershipCube = zeros(size(tmpROI.membershipCube));
lowerBasalRMembershipCube = zeros(size(tmpROI.membershipCube));
for i=1:numberOfRegionsOfInterest
    regionOfInterest(i) = pr.regionOfInterestFromAnatomy(pr.headGrid, includedRoiLabels{i});
    dipoleProbabilityInRegion(:,i) = gaussianWeightMatrix * regionOfInterest(i).membershipProbabilityCube(headGrid.insideBrainCube);
    
    % compute centroids of ROIs
    xCube = regionOfInterest(i).headGrid.xCube;
    yCube = regionOfInterest(i).headGrid.yCube;
    zCube = regionOfInterest(i).headGrid.zCube;
    membershipCube = regionOfInterest(i).membershipCube;
    xCentroid = mean(xCube(membershipCube));
    yCentroid = mean(yCube(membershipCube));
    zCentroid = mean(zCube(membershipCube));
    roiCentroids(i,:) = [xCentroid yCentroid zCentroid];
    
    if     ismember(i, upperBasalLIdx)
        upperBasalLMembershipCube = upperBasalLMembershipCube + membershipCube;
    elseif ismember(i, upperBasalRIdx)
        upperBasalRMembershipCube = upperBasalRMembershipCube + membershipCube;
    elseif ismember(i, lowerBasalLIdx)
        lowerBasalLMembershipCube = lowerBasalLMembershipCube + membershipCube;
    elseif ismember(i, lowerBasalRIdx)
        lowerBasalRMembershipCube = lowerBasalRMembershipCube + membershipCube;
    end
end

% Exclude ROIs included in Upper and Lower Basal
roiLabelsReduced = roiLabels(setdiff(1:88, [upperBasalLIdx upperBasalRIdx lowerBasalLIdx lowerBasalRIdx]));
roiLabelsReduced(73:76) = {'UpperBasal_L' 'UpperBasal_R' 'LowerBasal_L' 'LowerBasal_R'};
dipoleProbabilityInRegionReduced = dipoleProbabilityInRegion(:, setdiff(1:size(dipoleProbabilityInRegion,2), [upperBasalLIdx upperBasalRIdx lowerBasalLIdx lowerBasalRIdx]));
roiCentroidsReduced   = roiCentroids(setdiff(1:size(dipoleProbabilityInRegion,2), [upperBasalLIdx upperBasalRIdx lowerBasalLIdx lowerBasalRIdx]),:);

% Add integrated Upper and Lower Basals for density and centroid
upperBasalLDipDensity = sum(dipoleProbabilityInRegion(:,upperBasalLIdx),2);
upperBasalRDipDensity = sum(dipoleProbabilityInRegion(:,upperBasalRIdx),2);
lowerBasalLDipDensity = sum(dipoleProbabilityInRegion(:,lowerBasalLIdx),2);
lowerBasalRDipDensity = sum(dipoleProbabilityInRegion(:,lowerBasalRIdx),2);
dipoleProbabilityInRegionReduced(:, end+1:end+4) = [upperBasalLDipDensity upperBasalRDipDensity lowerBasalLDipDensity lowerBasalRDipDensity];
upperBasalLCentroid = [mean(xCube(logical(upperBasalLMembershipCube))) mean(yCube(logical(upperBasalLMembershipCube))) mean(zCube(logical(upperBasalLMembershipCube)))];
upperBasalRCentroid = [mean(xCube(logical(upperBasalRMembershipCube))) mean(yCube(logical(upperBasalRMembershipCube))) mean(zCube(logical(upperBasalRMembershipCube)))];
lowerBasalLCentroid = [mean(xCube(logical(lowerBasalLMembershipCube))) mean(yCube(logical(lowerBasalLMembershipCube))) mean(zCube(logical(lowerBasalLMembershipCube)))];
lowerBasalRCentroid = [mean(xCube(logical(lowerBasalRMembershipCube))) mean(yCube(logical(lowerBasalRMembershipCube))) mean(zCube(logical(lowerBasalRMembershipCube)))];
roiCentroidsReduced(end+1:end+4, :) = [upperBasalLCentroid; upperBasalRCentroid; lowerBasalLCentroid; lowerBasalRCentroid];
numberOfRegionsOfInterestReduced = length(roiCentroidsReduced);
includeRoiIdxReduced = 1:length(roiCentroidsReduced);

% Force the ROI centroids to be symmetric
leftHemisphereCentroids = roiCentroidsReduced(1:2:end,:);
xPositiveLeftHemisphereControids = bsxfun(@times, leftHemisphereCentroids, [-1 1 1]);
meanAbsRoiCentroids = round((xPositiveLeftHemisphereControids+roiCentroidsReduced(2:2:end,:))/2);
symmetricRoiCentroids = zeros(size(roiCentroidsReduced));
symmetricRoiCentroids(1:2:end,:) = bsxfun(@times, meanAbsRoiCentroids, [-1 1 1]);
symmetricRoiCentroids(2:2:end,:) = meanAbsRoiCentroids;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make a mask to determine ROIs that have more than X % of subjects. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain preselected ROI index.
icSubjIdList = [];
for subjId = 1:length(ALLEEG)
    tmpNumICs = ALLEEG(1,subjId).CAT.nbchan;
    icSubjIdList = [icSubjIdList; repmat(subjId, [tmpNumICs 1])];
end
subjDipoleDensityTable = zeros(length(ALLEEG), size(dipoleProbabilityInRegionReduced,2));
for subjId = 1:length(ALLEEG)
    tmpSubjIdx = find(icSubjIdList == subjId);
    subjDipoleDensityTable(subjId,:) = sum(dipoleProbabilityInRegionReduced(tmpSubjIdx,:),1);
end
roiNonzeroSubjDetectionVector = sum(subjDipoleDensityTable ~= 0);
userInputPercentage = str2num(get(handles.minSubjPercentEdit, 'String'))/100;
preselectedRoiIdx = find(roiNonzeroSubjDetectionVector>= userInputPercentage*size(subjDipoleDensityTable,1)); 

    %     % Sanity check
    %     sum(dipoleProbabilityInRegionReduced)
    %     subjDipoleDensityTable(:,preselectedRoiIdx)
    %     figure; bar(sort(vec(subjDipoleDensityTable(:,preselectedRoiIdx))))
    
% Report results.
roiDipoleDensiyPct = 100*sum(subjDipoleDensityTable(:,preselectedRoiIdx))/sum(dipoleProbabilityInRegionReduced(:));
[sortedRoiValuesPct, sortIdxPct] = sort(roiDipoleDensiyPct, 'descend');
sortedRoiLabels                  = roiLabelsReduced(preselectedRoiIdx(sortIdxPct));
sumDipoleDensity                 = sum(sortedRoiValuesPct);
reportCell = cell(length(sortedRoiLabels)+1,2);
reportCell(1,1) = {'--ROI labels--'};
reportCell(1,2) = {'--Dipole Density (%)--'};
reportCell(2:end,1) = sortedRoiLabels;
reportCell(2:end,2) = num2cell(sortedRoiValuesPct');
disp(reportCell)
preselectionLog = sprintf('%.0f/%.0f anatomical ROIs passed the selection.\nTotal of %.1f%% dipole density is accounted for.',...
                            length(preselectedRoiIdx), size(subjDipoleDensityTable,2), sumDipoleDensity);
set(handles.resultText, 'String', preselectionLog)
drawnow

% Convert the preselected ROI indices to edge (dipolePairdensity) index.
% Max == 76^2 - 76 (these are self referencing edges) == 5700. 
preselectedRoiIdxRepeatedMatrix = repmat((preselectedRoiIdx'-1)*length(roiLabelsReduced), [1 length(preselectedRoiIdx')]);
preselectedEdgeIdx = bsxfun(@plus, preselectedRoiIdxRepeatedMatrix, preselectedRoiIdx)';
preselectedEdgeIdx(logical(eye(size(preselectedEdgeIdx)))) = 0; % Edge's self connection is excluded.
preselectedEdgeIdx = nonzeros(preselectedEdgeIdx(:));



        % % Obtain dipole density that are excluded
        % numberOfRegionsOfNonInterest = length(excludedRoiLabels);
        % dipoleProbabilityInBasal = zeros(uniqueDipole.numberOfDipoles, numberOfRegionsOfNonInterest);
        % roiCentroidsBasal = zeros(length(excludedRoiLabels),3);
        % for i=1:numberOfRegionsOfNonInterest
        %     regionOfNonInterest(i) = pr.regionOfInterestFromAnatomy(pr.headGrid, excludedRoiLabels{i});
        %     dipoleProbabilityInBasal(:,i) = gaussianWeightMatrix * regionOfNonInterest(i).membershipProbabilityCube(headGrid.insideBrainCube);
        %     
        %     % compute centroids of ROIs
        %     xCubeBasal = regionOfNonInterest(i).headGrid.xCube;
        %     yCubeBasal = regionOfNonInterest(i).headGrid.yCube;
        %     zCubeBasal = regionOfNonInterest(i).headGrid.zCube;
        %     membershipCube = regionOfInterest(i).membershipCube;
        %     xCentroidBasal = mean(xCubeBasal(membershipCube));
        %     yCentroidBasal = mean(yCubeBasal(membershipCube));
        %     zCentroidBasal = mean(zCubeBasal(membershipCube));
        %     roiCentroidsBasal(i,:) = [xCentroidBasal yCentroidBasal zCentroidBasal];
        % end
        % dipDensityRemoved = sum(dipoleProbabilityInBasal(:))/(sum(dipoleProbabilityInRegion(:))+sum(dipoleProbabilityInBasal(:)));
        % disp(sprintf('\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
        % disp(sprintf('%.2f%% dipole density removed due to subcortical localization.', dipDensityRemoved*100))
        % disp(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute dipole pair density (non-normalized) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate subject list
fileNameList = {ALLEEG.filename}';

% Calculate dipole densities in source regions and destination regions. Reuse indices from 'from' dipole for 'to' dipole
fromDipoleProbabilityInRegion = zeros(dipolePairAndMeasureObj.from.numberOfDipoles, numberOfRegionsOfInterestReduced);
toDipoleProbabilityInRegion   = zeros(dipolePairAndMeasureObj.to.numberOfDipoles,   numberOfRegionsOfInterestReduced);
[~,~,toIdReverse] = unique(dipolePairAndMeasureObj.toDipoleId, 'stable');
for i = 1:numberOfRegionsOfInterestReduced
    fromDipoleProbabilityInRegion(:,i) = dipoleProbabilityInRegionReduced(fromIdReverse, i);
    toDipoleProbabilityInRegion(:,  i) = dipoleProbabilityInRegionReduced(toIdReverse,   i);
end

% Prepare dipole pair density by multiplying from-density and to-density. The first dimension is the sum of each subject's IC*IC, the second is ROI^2
% Note that this maps the combination of ICs to ROI*ROI.
dipolePairDensityFromIcSquareToRoiSquare = zeros(dipolePairAndMeasureObj.from.numberOfDipoles, length(includeRoiIdxReduced) * length(includeRoiIdxReduced));
for i = 1:size(dipolePairDensityFromIcSquareToRoiSquare,1)
    tmpDipoleProbability = (fromDipoleProbabilityInRegion(i,:)' *  toDipoleProbabilityInRegion(i,:))';
    tmpDipoleProbability(logical(eye(size(tmpDipoleProbability)))) = 0; % Edge's self connection is excluded.
    dipolePairDensityFromIcSquareToRoiSquare(i, :) = vec(tmpDipoleProbability)';
end

    %     % Check if diagonal dipole density is removed.
    %     traceList = zeros(size(dipolePairDensityFromIcSquareToRoiSquare,1),1);
    %     for n = 1:size(dipolePairDensityFromIcSquareToRoiSquare,1)
    %         tmp = reshape(dipolePairDensityFromIcSquareToRoiSquare(n,:), [76 76]);
    %         traceList(n) = trace(tmp);
    %     end
        
    % % Sanity check--passed.
    % sumDipDensity = sum(dipoleProbabilityInRegionReduced(:)); % Sum of dip density of all subjects
    % unselectedIdx = setdiff(1:size(dipolePairProbabilityOnFromRegionToRegionSubstrate,2), preselectedEdgeIdx(:));
    % mustBeZero = sum(sum(dipolePairProbabilityOnFromRegionToRegionSubstrate(:,unselectedIdx))); % 0 Confirmed, when 100% of dipole density is accounted for while not all ROIs passed the selection.

    % Check if all IC*IC has dipolePairDensity.
    % Note that some of them have zero dipole density in IC^2 * ROI^2 data (why?)
    %     figure; bar(sort(sum(dipolePairProbabilityOnFromRegionToRegionSubstrate,2)));
    %     zeroEdgeDensityIdx = find(sum(dipolePairProbabilityOnFromRegionToRegionSubstrate,2)==0);
    %     zeroEdgeSubjects = dipolePairAndMeasureObj.sessionId(zeroEdgeDensityIdx);
    %     tmp  = dipolePairAndMeasureObj.linearizedMeasure(zeroEdgeDensityIdx(7),:);
    %     tmp2 = reshape(tmp, [30 100]);
    %     figure; imagesc(tmp2); axis xy; colorbar

% Compute dipole pair density for only preselected pairs (actively zero-outs the non-selected pairs; also saves time.)
numberOfsessions = max(dipolePairAndMeasureObj.sessionId); % Session means subjects.
dipolePairDensity = zeros(size(dipolePairDensityFromIcSquareToRoiSquare, 2), numberOfsessions);
for i = 1:length(preselectedEdgeIdx)
    currentEdgeIdx = preselectedEdgeIdx(i);
    currentEdgeIcSquare = dipolePairDensityFromIcSquareToRoiSquare(:,currentEdgeIdx);
    sessionDensity = zeros(1,numberOfsessions);
    for j=1:numberOfsessions
        sessionDensity(j) = sum(currentEdgeIcSquare(dipolePairAndMeasureObj.sessionId == j));
    end

    % Store dipole pair density in ROI^2 x subj.
    dipolePairDensity(currentEdgeIdx,:) = sessionDensity;
end

    %     % Test
    %     tmp = reshape(dipolePairDensity, [76 76 10]);
    %     traceList = zeros(size(tmp,3),1);
    %     for n = 1:size(tmp,3)
    %         traceList(n) = trace(tmp(:,:,n));
    %     end

    % Note that sum(sum(dipolePairDensity,2)>0) == numNonzereoEdges confirmed.
    % sum(sum(dipolePairDensity,2)>0)
    
dipolePairDensity = reshape(dipolePairDensity, [length(includeRoiIdxReduced) length(includeRoiIdxReduced) numberOfsessions]);

        % % Calculate dipole pair density for each unique edge
        % numberOfsessions = max(dipolePairAndMeasureObj.sessionId); % session means subjects
        % dipolePairDensity = zeros(size(dipolePairProbabilityOnFromRegionToRegionSubstrate, 2), numberOfsessions);
        % for i = 1:size(dipolePairDensity,1)
        %     dipolePairProbabilityOnFromRegionToRegion = dipolePairProbabilityOnFromRegionToRegionSubstrate(:,i);
        %     sessionDensity = zeros(1,numberOfsessions);
        %     for j=1:numberOfsessions
        %         sessionDensity(j) = sum(dipolePairProbabilityOnFromRegionToRegion(dipolePairAndMeasureObj.sessionId == j));
        %     end
        %     
        %     % store dipole pair density
        %     dipolePairDensity(i,:) = sessionDensity;
        % end
        % dipolePairDensity = reshape(dipolePairDensity, [length(includeRoiIdxReduced) length(includeRoiIdxReduced) numberOfsessions]);

        % Note: The following 'convergence similarity' thing is disabled due to circular statistics (a.k.a. double-dipping) issue. Sorry Nima!
        % % calculate dipole-pair pairwise similarity (could be further optimized by not calculating pair that do not share dipole-pair density in any ROIs)
        % % Note: correlationSimilarity must have a size of (numberOfDipolePairs)*(numberOfDipolePairs), NOT (linearizedTimeFreqs)*(linearizedTimeFreqs) 05/02/2016 Clement and Makoto.
        % correlationSimilarity = 1-squareform(pdist(dipolePairAndMeasureObj.linearizedMeasure, 'correlation'));
        % connectivitySimilarity = estimate_mutual_information_from_correlation(correlationSimilarity);
        %
        % % save dipolePairDensity and others that is needed for similarity test
        % dipoleProbabilityInRegion = dipoleProbabilityInRegionReduced;
        % roiLabels = roiLabelsReduced;
        % save([workingFolder filesep userInputFilename '_dipolePairDensity'], 'dipolePairDensity', 'dipoleProbabilityInRegion', 'roiLabels', 'symmetricRoiCentroids', 'fileNameList')
        % 
        % % save convergenceData
        % linearizedMeasure = dipolePairAndMeasureObj.linearizedMeasure;
        % save([workingFolder filesep userInputFilename '_measureConvergence'], 'dipolePairProbabilityOnFromRegionToRegionSubstrate', 'linearizedMeasure', 'connectivitySimilarity');

% Save dipolePairDensity and others that is needed for similarity test.
dipoleProbabilityInRegion = dipoleProbabilityInRegionReduced;
linearizedMeasure = dipolePairAndMeasureObj.linearizedMeasure;
roiLabels = roiLabelsReduced;
save([workingFolder userInputFilename '_dipolePairDensity'],...
    'dipolePairDensity', 'dipoleProbabilityInRegion', 'roiLabels', 'symmetricRoiCentroids',...
    'fileNameList', 'dipolePairDensityFromIcSquareToRoiSquare', 'linearizedMeasure',...
    'preselectedRoiIdx', '-v7.3');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate individual subject's connectivity map %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the map
numOfRegionsOfInterest = size(dipolePairDensity,1);
toRegionNumber   = repmat((1:numOfRegionsOfInterest)', [1, numOfRegionsOfInterest]); % first deimension (rows) contains toRegion Numbers.
toRegionNumber   = toRegionNumber(:);
fromRegionNumber = repmat( 1:numOfRegionsOfInterest,   [numOfRegionsOfInterest, 1]); % second deimension (columns) contains fromRegion Numbers.
fromRegionNumber = fromRegionNumber(:);
    % Double check the from and to!!
    % ALLEEG(1,1).CAT.Conn: dims: {'var_to'  'var_from'  'freq'  'time'}
    % figure; imagesc(toRegionNumber); title('To Indices')
    % figure; imagesc(fromRegionNumber); title('From Indices')

if get(handles.rpdcButton,'value')==1
    connectivityType = 'rPDC';
else
    connectivityType = 'dDTF08';
end



allConnectivityStack = single(zeros(size(dipolePairDensity,1), size(dipolePairDensity,2), length(frequencies), length(latencies), length(allFiles)));
for sessionIdx = unique(dipolePairAndMeasureObj.sessionId)'
    effectiveConnectivityTimeFreq = single(zeros(timeFreqSize(1)*timeFreqSize(2), size(dipolePairDensity,1)*size(dipolePairDensity,2)));
    tic
    for fromIdx = preselectedRoiIdx
        for toIdx = preselectedRoiIdx
            
            % Exclude the same ROI to ROI connection.
            if fromIdx == toIdx
                continue
            end
            currentEdgeIdx = find((fromRegionNumber==fromIdx) & (toRegionNumber==toIdx)); % Up to ROI^2
            currentSessionIdx = find(dipolePairAndMeasureObj.sessionId==sessionIdx);      % Obtain indices to extract each subject's IC*IC (vectorized)
            dipolePairDensitiesInIcSquare = dipolePairDensityFromIcSquareToRoiSquare(currentSessionIdx, currentEdgeIdx);
            
                    % %%%%%%%%%%%% % SINGLE IC CAN CREATE INTER-ROI CONNECTIVITY!!
                    % %%% Test %%% % Confirmed that multiple combinations of 'fromIdx' and 'toIdx' points to the same combination of ICs diagonal dipole density with zero inforflow.
                    % %%%%%%%%%%%% % This is a problem of modeling dipole pair density.
                    % tmp = reshape(dipolePairDensitiesInIcSquare, [12 12]);
                    % figure; imagesc(tmp)
                    % 
                    % nonZeroIdx = find(dipolePairDensitiesInIcSquare);
                    % tmp2 = dipolePairAndMeasureObj.linearizedMeasure(nonZeroIdx,:);
                    % figure
                    % for n = 1:length(nonZeroIdx)
                    %     subplot(1,length(nonZeroIdx),n)
                    %     imagesc(reshape(squeeze(tmp2(n,:)), [30 100]))
                    % end
                    % disp('test')

            % 03/04/2017 Makoto. Outsize the 3SD is now zero. Subjects CAN have all zero data, since less than 100% of subject overlap is possible.
            % 02/10/2016 Makoto. Put it back after discussing with Nima and Scott (during and after my presentation at lab meeting)
            % 12/23/2015 Makoto. Disabled this normalization since this will inflate near-zero values
            normFactor = sum(dipolePairDensitiesInIcSquare);
            if normFactor > 0
                dipolePairDensitiesInIcSquareNormalized = dipolePairDensitiesInIcSquare / normFactor;
            else
                continue
            end
            
            % This is the projected measure
            projectedMeasure = bsxfun(@times, dipolePairAndMeasureObj.linearizedMeasure(currentSessionIdx,:), dipolePairDensitiesInIcSquareNormalized);
            effectiveConnectivityTimeFreq(:,currentEdgeIdx) = squeeze(sum(projectedMeasure));
        end
    end
    effectiveConnectivityTimeFreq = single(effectiveConnectivityTimeFreq); % use single to save time and disk space
    effectiveConnectivityTimeFreq = permute(reshape(effectiveConnectivityTimeFreq, timeFreqSize(1), timeFreqSize(2), max(fromRegionNumber), max(fromRegionNumber)), [3 4 1 2]); % Use most consistent dimensions for SIFT and most intuitive for users
    %save(sprintf([workingFolder filesep userInputFilename '_%04.f_connectivity'], sessionIdx), 'effectiveConnectivityTimeFreq', 'connectivityType', 'latencies', 'frequencies', 'dimensionLabels', 'fileNameList');
    
    allConnectivityStack(:,:,:,:,sessionIdx) = effectiveConnectivityTimeFreq;
    
    timeLapse = toc;
    disp(sprintf('%2.0d/%2.0d subjects done (%0.1d sec lapsed)', sessionIdx, length(unique(dipolePairAndMeasureObj.sessionId)), round(timeLapse)));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exclude edges that have less than the specified number of subjects %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate who shall die for fake connections by single dipole ranging across multiple ROIs.
sumEdgeMatrix        = sum(logical(squeeze(sum(sum(allConnectivityStack,3),4))),3);
numSubjThreshold     = size(allConnectivityStack,5)*userInputPercentage;
stillGoodMask        = sumEdgeMatrix>=numSubjThreshold;

% Re-define the results and indices by applying the mask. This kills a lot.
allConnectivityStack = bsxfun(@times, allConnectivityStack, stillGoodMask);
finallySelectedEdgeIdx   = find(stillGoodMask);

% Generate final report.
finalSelectionLog = sprintf('%.0f graph edges have cross-IC AND cross-ROI connections\nThese ones are valid and submitted to the final analyses.', length(finallySelectedEdgeIdx));
set(handles.finalSelectionText, 'String', finalSelectionLog)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save 5-D data for statistics as Matlab version 7.3 compatible to support >2GB data. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([workingFolder userInputFilename '_allSubjStack'],...
    'allConnectivityStack', 'dimensionLabels', 'frequencies', 'latencies',...
    'connectivityType', 'fileNameList', 'finalSelectionLog', 'finallySelectedEdgeIdx', '-v7.3')



% Display process end
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''3.Convert to group anatomical ROIs'' finished. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))





function minSubjPercentEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minSubjPercentEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minSubjPercentEdit as text
%        str2double(get(hObject,'String')) returns contents of minSubjPercentEdit as a double



% --- Executes during object creation, after setting all properties.
function minSubjPercentEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSubjPercentEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in minSubjPercentTestPushbutton.
function minSubjPercentTestPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to minSubjPercentTestPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain multiple .set files
[allFiles, workingFolder] = uigetfile('*.set', 'MultiSelect', 'on');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end

% Move to the working folder
cd(workingFolder)

% Load empty dipolePairAndMeasure from eeglabroot/plugins/groupSIFT/+pr: Thanks Clement!
dipolePairAndMeasureObj = pr.dipolePairAndMeasure;

% Load all .set data located under the working folder
ALLEEG = [];
for n = 1:length(allFiles)
    loadName = allFiles{n};
    EEG = pop_loadset('filename', loadName,'filepath',workingFolder, 'loadmode', 'info');
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
end

% Display process start
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% Compute the number of ROIs and estimate the upper bound of total dipole density accounted for %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))

% Compute how many IC combinations in all subjects
numICs = zeros(length(ALLEEG),1);
for n = 1:length(ALLEEG)
    numICs(n) = size(ALLEEG(1,n).icaweights,1);
end
numIcCombination = sum(numICs.^2);

% Assign unique IDs to dipoles
if get(handles.rpdcButton,'value')==1
    timeFreqSize = size(squeeze(ALLEEG(1).CAT.Conn.RPDC(  1,1,:,:)));
else
    timeFreqSize = size(squeeze(ALLEEG(1).CAT.Conn.dDTF08(1,1,:,:)));
end
latencies = ALLEEG(1).CAT.Conn.erWinCenterTimes;
frequencies = ALLEEG(1).CAT.Conn.freqs;
dimensionLabels= ALLEEG(1).CAT.Conn.dims;
dipolePairAndMeasureObj.linearizedMeasure = zeros([numIcCombination timeFreqSize(1)*timeFreqSize(2)]);
dipoleCounter = 0;
counter = 1;
for n = 1:length(ALLEEG)
    numberOfSessiondipoles = length(ALLEEG(n).dipfit.model);
    dipoleLocation = [];
    dipoleResidualVariance = [];
    for m = 1:numberOfSessiondipoles
       %dipoleLocation(m,:)       = ALLEEG(n).dipfit.model(m).posxyz(1,:);
        dipoleLocation(m,:)       = selectDipWithLargerMoment(ALLEEG(n).dipfit.model(m).posxyz, ALLEEG(n).dipfit.model(m).momxyz);
        dipoleResidualVariance(m) = ALLEEG(n).dipfit.model(m).rv;
    end;
    dipoleId = (dipoleCounter+1):(dipoleCounter+length(dipoleLocation));
    dipoleCounter = dipoleCounter + length(dipoleLocation);
    
    for m = 1:size(dipoleLocation,1) % from location: confirm it with ALLEEG(1,1).CAT.Conn: dims: {'var_to'  'var_from'  'freq'  'time'}
        for k = 1:size(dipoleLocation,1) % to location
           %dipolePairAndMeasureObj.from.location = cat(1, dipolePairAndMeasureObj.from.location, ALLEEG(n).dipfit.model(m).posxyz(1,:));
            dipolePairAndMeasureObj.from.location = cat(1, dipolePairAndMeasureObj.from.location, selectDipWithLargerMoment(ALLEEG(n).dipfit.model(m).posxyz, ALLEEG(n).dipfit.model(m).momxyz));
            dipolePairAndMeasureObj.from.residualVariance(end+1,1) = ALLEEG(n).dipfit.model(m).rv;
           %dipolePairAndMeasureObj.to.location = cat(1, dipolePairAndMeasureObj.to.location, ALLEEG(n).dipfit.model(m).posxyz(1,:));
            dipolePairAndMeasureObj.to.location = cat(1, dipolePairAndMeasureObj.to.location, selectDipWithLargerMoment(ALLEEG(n).dipfit.model(m).posxyz, ALLEEG(n).dipfit.model(m).momxyz));
            if m == k
                dipolePairAndMeasureObj.linearizedMeasure(counter,:) = zeros(timeFreqSize(1)*timeFreqSize(2),1);
            else
                if get(handles.rpdcButton,'value')==1 % Again, ALLEEG(1,1).CAT.Conn: dims: {'var_to'  'var_from'  'freq'  'time'}
                    tmpTimeFreqMatrix = squeeze(ALLEEG(n).CAT.Conn.RPDC(  k,m,:,:));
                else
                    tmpTimeFreqMatrix = squeeze(ALLEEG(n).CAT.Conn.dDTF08(k,m,:,:));
                end
                dipolePairAndMeasureObj.linearizedMeasure(counter,:) = tmpTimeFreqMatrix(:);
            end
            dipolePairAndMeasureObj.sessionId(counter,1) = n;
            dipolePairAndMeasureObj.fromDipoleId(counter,1) = dipoleId(m);
            dipolePairAndMeasureObj.toDipoleId(counter,1)   = dipoleId(k);
            counter = counter + 1;
        end
    end
end

% Find unique dipoles
uniqueDipole = pr.dipole;
[~,fromId,fromIdReverse] = unique(dipolePairAndMeasureObj.fromDipoleId, 'stable');

% Extract coordinates + residual variance of the unique dipoles
uniqueDipole.location = dipolePairAndMeasureObj.from.location(fromId,:);
uniqueDipole.residualVariance = dipolePairAndMeasureObj.from.residualVariance(fromId);

% Load headGrid cubes
headGrid = pr.headGrid;

% Set Gaussian smoothing kernel
userInputKernelSize = str2num(get(handles.gaussianSizeEdit, 'String'));
standardDeviationOfEstimatedDipoleLocation = userInputKernelSize/2.355; % this calculates sigma in Gaussian equation
projectionParameter = pr.projectionParameter(standardDeviationOfEstimatedDipoleLocation);
%projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian = 150/standardDeviationOfEstimatedDipoleLocation;
projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian = 3;
[~,~,gaussianWeightMatrix]= pr.meanProjection.getProjectionMatrix(uniqueDipole, headGrid, projectionParameter);

% Define valid anatomical ROIs as EEG sources agreed by Scott and Makoto in Dec 2015; See below for the list
% excludeRoiIdx = [21 22 37:42 71:78]; % 29 30 are insula... better to include?
excludeRoiIdx = [];
includeRoiIdx = setdiff(1:88, excludeRoiIdx);
roiLabels = pr.regionOfInterestFromAnatomy.getAllAnatomicalLabels;

% Include and exclude ROIs
includedRoiLabels = roiLabels(includeRoiIdx);
excludedRoiLabels = roiLabels(excludeRoiIdx);

% Define Upper and Lower Basal ROIs
upperBasalLIdx = [21 71:2:77]; % Upper basal: Olfactory, Caudate, Putamen, Pallidum, and Thalamus
upperBasalRIdx = [22 72:2:78]; % Upper basal: Olfactory, Caudate, Putamen, Pallidum, and Thalamus
lowerBasalLIdx = [37:2:41];    % Lower basal: Hippocampus, Parahippocampus, Amygdala
lowerBasalRIdx = [38:2:42];    % Lower basal: Hippocampus, Parahippocampus, Amygdala

% Obtain dipole density and centroids in each ROI
numberOfRegionsOfInterest = length(includedRoiLabels);
dipoleProbabilityInRegion = zeros(uniqueDipole.numberOfDipoles, numberOfRegionsOfInterest);
roiCentroids = zeros(length(includedRoiLabels),3);
tmpROI = pr.regionOfInterestFromAnatomy(pr.headGrid, includedRoiLabels{1});
upperBasalLMembershipCube = zeros(size(tmpROI.membershipCube));
upperBasalRMembershipCube = zeros(size(tmpROI.membershipCube));
lowerBasalLMembershipCube = zeros(size(tmpROI.membershipCube));
lowerBasalRMembershipCube = zeros(size(tmpROI.membershipCube));
for i=1:numberOfRegionsOfInterest
    regionOfInterest(i) = pr.regionOfInterestFromAnatomy(pr.headGrid, includedRoiLabels{i});
    dipoleProbabilityInRegion(:,i) = gaussianWeightMatrix * regionOfInterest(i).membershipProbabilityCube(headGrid.insideBrainCube);
    
    % compute centroids of ROIs
    xCube = regionOfInterest(i).headGrid.xCube;
    yCube = regionOfInterest(i).headGrid.yCube;
    zCube = regionOfInterest(i).headGrid.zCube;
    membershipCube = regionOfInterest(i).membershipCube;
    xCentroid = mean(xCube(membershipCube));
    yCentroid = mean(yCube(membershipCube));
    zCentroid = mean(zCube(membershipCube));
    roiCentroids(i,:) = [xCentroid yCentroid zCentroid];
    
    if     ismember(i, upperBasalLIdx)
        upperBasalLMembershipCube = upperBasalLMembershipCube + membershipCube;
    elseif ismember(i, upperBasalRIdx)
        upperBasalRMembershipCube = upperBasalRMembershipCube + membershipCube;
    elseif ismember(i, lowerBasalLIdx)
        lowerBasalLMembershipCube = lowerBasalLMembershipCube + membershipCube;
    elseif ismember(i, lowerBasalRIdx)
        lowerBasalRMembershipCube = lowerBasalRMembershipCube + membershipCube;
    end
end

% Exclude ROIs included in Upper and Lower Basal
roiLabelsReduced = roiLabels(setdiff(1:88, [upperBasalLIdx upperBasalRIdx lowerBasalLIdx lowerBasalRIdx]));
roiLabelsReduced(73:76) = {'UpperBasal_L' 'UpperBasal_R' 'LowerBasal_L' 'LowerBasal_R'};
dipoleProbabilityInRegionReduced = dipoleProbabilityInRegion(:, setdiff(1:size(dipoleProbabilityInRegion,2), [upperBasalLIdx upperBasalRIdx lowerBasalLIdx lowerBasalRIdx]));
roiCentroidsReduced   = roiCentroids(setdiff(1:size(dipoleProbabilityInRegion,2), [upperBasalLIdx upperBasalRIdx lowerBasalLIdx lowerBasalRIdx]),:);

% Add integrated Upper and Lower Basals for density and centroid
upperBasalLDipDensity = sum(dipoleProbabilityInRegion(:,upperBasalLIdx),2);
upperBasalRDipDensity = sum(dipoleProbabilityInRegion(:,upperBasalRIdx),2);
lowerBasalLDipDensity = sum(dipoleProbabilityInRegion(:,lowerBasalLIdx),2);
lowerBasalRDipDensity = sum(dipoleProbabilityInRegion(:,lowerBasalRIdx),2);
dipoleProbabilityInRegionReduced(:, end+1:end+4) = [upperBasalLDipDensity upperBasalRDipDensity lowerBasalLDipDensity lowerBasalRDipDensity];
upperBasalLCentroid = [mean(xCube(logical(upperBasalLMembershipCube))) mean(yCube(logical(upperBasalLMembershipCube))) mean(zCube(logical(upperBasalLMembershipCube)))];
upperBasalRCentroid = [mean(xCube(logical(upperBasalRMembershipCube))) mean(yCube(logical(upperBasalRMembershipCube))) mean(zCube(logical(upperBasalRMembershipCube)))];
lowerBasalLCentroid = [mean(xCube(logical(lowerBasalLMembershipCube))) mean(yCube(logical(lowerBasalLMembershipCube))) mean(zCube(logical(lowerBasalLMembershipCube)))];
lowerBasalRCentroid = [mean(xCube(logical(lowerBasalRMembershipCube))) mean(yCube(logical(lowerBasalRMembershipCube))) mean(zCube(logical(lowerBasalRMembershipCube)))];
roiCentroidsReduced(end+1:end+4, :) = [upperBasalLCentroid; upperBasalRCentroid; lowerBasalLCentroid; lowerBasalRCentroid];
numberOfRegionsOfInterestReduced = length(roiCentroidsReduced);
includeRoiIdxReduced = 1:length(roiCentroidsReduced);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make a mask to determine ROIs that receives user-specified percent of unique subjects. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain preselected ROI index.
icSubjIdList = [];
for subjId = 1:length(ALLEEG)
    tmpNumICs = ALLEEG(1,subjId).CAT.nbchan;
    icSubjIdList = [icSubjIdList; repmat(subjId, [tmpNumICs 1])];
end
subjDipDensityList = zeros(length(ALLEEG), size(dipoleProbabilityInRegionReduced,2));
for subjId = 1:length(ALLEEG)
    tmpSubjIdx = find(icSubjIdList == subjId);
    subjDipDensityList(subjId,:) = sum(dipoleProbabilityInRegionReduced(tmpSubjIdx,:),1);
end
roiNonzeroSubjDetectionVector = sum(subjDipDensityList ~= 0);
userInputPercentage = str2num(get(handles.minSubjPercentEdit, 'String'))/100;
preselectedRoiIdx = find(roiNonzeroSubjDetectionVector>= userInputPercentage*size(subjDipDensityList,1)); 
roiDipoleDensiyPct = 100*sum(subjDipDensityList(:,preselectedRoiIdx))/sum(dipoleProbabilityInRegionReduced(:));

%%%%%%%%%%%%%%%%%%%%%%
%%% Report resutls %%%
%%%%%%%%%%%%%%%%%%%%%%
[sortedRoiValuesPct, sortIdxPct] = sort(roiDipoleDensiyPct, 'descend');
sortedRoiLabels                  = roiLabelsReduced(preselectedRoiIdx(sortIdxPct));
sumDipoleDensity                 = sum(sortedRoiValuesPct);
reportCell = cell(length(sortedRoiLabels)+1,2);
reportCell(1,1) = {'-ROI labels-'};
reportCell(1,2) = {'-Dip Density-'};
reportCell(2:end,1) = sortedRoiLabels;
for n = 1:length(sortedRoiValuesPct')
    reportCell{n+1,2} = sprintf('%.2f', sortedRoiValuesPct(n));
end
disp(reportCell)
preselectionLog = sprintf('%.0f/%.0f anatomical ROIs passed the selection.\nTotal of %.1f%% dipole density is accounted for.',...
                            length(preselectedRoiIdx), size(subjDipDensityList,2), sumDipoleDensity);
set(handles.resultText, 'String', preselectionLog)

% Visualize the ROIs.
roiSum = logical(zeros(size(regionOfInterest(1).membershipCube)));
for n = 1:length(sortedRoiLabels)
    tmpRoi = pr.regionOfInterestFromAnatomy(pr.headGrid, sortedRoiLabels{n});
    roiSum = roiSum + double(tmpRoi.membershipCube)*sortedRoiValuesPct(n);
end

% Adjust roiSum position by 1 slide lower.
roiSumAdjusted = cat(3, roiSum(:,:,2:end), zeros(size(roiSum,1), size(roiSum,2)));
[Xq, Yq, Zq]       = meshgrid(23/91:23/91:23, 27/109:27/109:27, 23/91:23/91:23); % Resize roiSum to 109 x 91 x 91, which is the size brainBlobBrowser takes.
interpolatedRoiSum = permute(interp3(roiSumAdjusted, Xq, Yq, Zq), [2 1 3]);
nanMask            = isnan(interpolatedRoiSum);
interpolatedRoiSum(nanMask) = 0;

% Show dipole density using brainBlobBrowser.
brainBlobBrowserCustom('data', interpolatedRoiSum, 'roiDipoleDensityReport', reportCell)

% Display process End
% Display process start
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% "Compute the number of ROIs and estimate the upper bound of total dipole density accounted for" finished. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))



% 12/09/2019 Makoto. When dual dipoles are detected, compare their moments and pick up the larger one.
function output = selectDipWithLargerMoment(posxyz, momxyz)

if norm(posxyz)> 200
    error('Abnormal dipole location detected.')
end

if size(posxyz, 1) == 2
    mom1 = norm(momxyz(1,:));
    mom2 = norm(momxyz(2,:));
    if mom1 > mom2
        output = posxyz(1,:);
    else
        output = posxyz(2,:);
    end
else
    output = posxyz;
end
