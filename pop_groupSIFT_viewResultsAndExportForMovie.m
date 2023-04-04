% pop_groupSIFT_viewResultsAndExportForMovie(varargin)
%
% History
% 04/04/2023 Makoto. Bug fixed. 'normlen', 'on' added for a workaround for the bug in dipplot().
% 08/20/2020 Makoto. Bug fixed. Edge name highlighted in red fixed.
% 08/15/2020 Makoto. Updated. The relation between 'pooling edges' and graph edge-wise multiple comparison correction was clarified. GUI and recommendations changed accordingly. The nature of the surrogate stat distribution of 'mass of cluster' seems to have an interesting property.
% 07/20/2020 Makoto. Updated. Supported Mike X Cohen's 'Matlab for Brain and Cognitive Scientists' p.245 calculation for determining mass of cluster threshold. Edge-pooled extreme value statistics. GFWER control (u=1), for which see Groppe et al. (2011)
% 05/31/2020 Makoto. Updated. Supported 2x2 test.
% 02/25/2020 Makoto. Updated. Contour plot for the case of subtraction.
% 12/20/2019 Makoto. Updated. The minor graphic details fixed.
% 12/13/2019 Makoto. Updated. Compatible with Matlab 2017b now.

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

function varargout = pop_groupSIFT_viewResultsAndExportForMovie(varargin)
% POP_GROUPSIFT_VIEWRESULTSANDEXPORTFORMOVIE MATLAB code for pop_groupSIFT_viewResultsAndExportForMovie.fig
%      POP_GROUPSIFT_VIEWRESULTSANDEXPORTFORMOVIE, by itself, creates a new POP_GROUPSIFT_VIEWRESULTSANDEXPORTFORMOVIE or raises the existing
%      singleton*.
%
%      H = POP_GROUPSIFT_VIEWRESULTSANDEXPORTFORMOVIE returns the handle to a new POP_GROUPSIFT_VIEWRESULTSANDEXPORTFORMOVIE or the handle to
%      the existing singleton*.
%
%      POP_GROUPSIFT_VIEWRESULTSANDEXPORTFORMOVIE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_GROUPSIFT_VIEWRESULTSANDEXPORTFORMOVIE.M with the given input arguments.
%
%      POP_GROUPSIFT_VIEWRESULTSANDEXPORTFORMOVIE('Property','Value',...) creates a new POP_GROUPSIFT_VIEWRESULTSANDEXPORTFORMOVIE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_groupSIFT_viewResultsAndExportForMovie_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_groupSIFT_viewResultsAndExportForMovie_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_groupSIFT_viewResultsAndExportForMovie

% Last Modified by GUIDE v2.5 15-Aug-2020 10:43:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_groupSIFT_viewResultsAndExportForMovie_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_groupSIFT_viewResultsAndExportForMovie_OutputFcn, ...
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



function [processedLogoData, customColorMap] = prepareLogoData(inputImage, targetAxesHandle)

% Obtain the pixel size of the target region.
set(targetAxesHandle, 'Units', 'pixels');
targetPositionInPixel = get(targetAxesHandle, 'position');
targetSize = round(targetPositionInPixel([4 3])); % In the order of vertical, horizontal.

% Process the input image.
logoData = mean(double(inputImage),3);
shorterEdgeLength = min(targetSize);
longerEdgeLength  = max(targetSize);
lengthRatio = longerEdgeLength/shorterEdgeLength;
outputLongerEdgeLength = round(size(logoData,1)*lengthRatio);
if targetSize(1) == targetSize(2)
    processedLogoData = logoData;
elseif targetSize(1) > targetSize(2)
    processedLogoData = repmat(logoData(1,1), outputLongerEdgeLength, size(logoData,2));
    oneSideMergin = floor((outputLongerEdgeLength-size(logoData,2))/2);
    processedLogoData(oneSideMergin+1:oneSideMergin+size(logoData,2),:) = logoData;
else
    processedLogoData = repmat(logoData(1,1), size(logoData,1), outputLongerEdgeLength);
    oneSideMergin = floor((outputLongerEdgeLength-size(logoData,1))/2);
    processedLogoData(:, oneSideMergin+1:oneSideMergin+size(logoData,2)) = logoData;
end

% Obtain the colormap.
clockOutputs = clock;
currentMonth = clockOutputs(2);
colorRatio = 0.12221946;
if currentMonth>=3 & currentMonth<=5
    originalColorMap = colormap(spring(256));
elseif currentMonth>=6 & currentMonth<=9
    originalColorMap = colormap(summer(256));
elseif currentMonth>=10 & currentMonth<=11
    originalColorMap = colormap(autumn(256)); 
else
    originalColorMap = colormap(winter(256));
end
customColorMap = originalColorMap*colorRatio + repmat(1-colorRatio, size(originalColorMap));



% --- Executes just before pop_groupSIFT_viewResultsAndExportForMovie is made visible.
function pop_groupSIFT_viewResultsAndExportForMovie_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_groupSIFT_viewResultsAndExportForMovie (see VARARGIN)

% Auto-fill the Edit Boxes
try
    tStatisticsPath       = dir('*_tStatistics.mat');
    dipolePairDensityPath = dir('*_dipolePairDensity.mat');
    
    set(handles.loadTstatisticsEdit,       'String', [pwd filesep tStatisticsPath.name]);
    set(handles.loadDipolePairDensityEdit, 'String', [pwd filesep dipolePairDensityPath.name]);
end

% % display no data logo on connectivityAxes
% % Updated 12/13/2019 Makoto.
% axes(handles.connectivityAxes)
% logoPath = which('sccnLogo.jpg'); % Scott mentioned the relation between network analysis and mission of SCCN during tea time
% logoData = imread(logoPath);
% logoData = mean(double(logoData),3);
% clockOutputs = clock;
% currentMonth = clockOutputs(2);
% colorRatio = 0.12221946;
% if currentMonth>=3 & currentMonth<=5
%     originalColorMap = colormap(spring(256));
% elseif currentMonth>=6 & currentMonth<=9
%     originalColorMap = colormap(summer(256));
% elseif currentMonth>=10 & currentMonth<=11
%     originalColorMap = colormap(autumn(256)); 
% else
%     originalColorMap = colormap(winter(256));
% end
% customColorMap   = originalColorMap*colorRatio + repmat(1-colorRatio, size(originalColorMap));
% colormap(customColorMap)
% imagesc(logoData);
% set(gca, 'XTick', [], 'YTick', [])
% drawnow
% 
% axes(handles.timeFrequencyAxes);
% rectangleLogo = repmat(logoData(1,1), [length(logoData) round(length(logoData)*4/3)]);
% startPoint = round(size(rectangleLogo,2)/2)-round(size(logoData,2)/2);
% rectangleLogo(:,startPoint:startPoint+length(logoData)-1) = logoData;
% clockOutputs = clock;
% currentMonth = clockOutputs(2);
% colorRatio = 0.12221946;
% if currentMonth>=3 & currentMonth<=5
%     originalColorMap = colormap(spring(256));
% elseif currentMonth>=6 & currentMonth<=9
%     originalColorMap = colormap(summer(256));
% elseif currentMonth>=10 & currentMonth<=11
%     originalColorMap = colormap(autumn(256)); 
% else
%     originalColorMap = colormap(winter(256));
% end
% customColorMap   = originalColorMap*colorRatio + repmat(1-colorRatio, size(originalColorMap));
% colormap(customColorMap)
% imagesc(rectangleLogo);
% set(gca, 'XTick', [], 'YTick', [])
% drawnow
% 
% % display no data logo on dipole plot axes
% rectangleLogo = repmat(logoData(1,1), [round(length(logoData)*6/5) length(logoData)]);
% startPoint = round(size(rectangleLogo,1)/2)-round(size(logoData,1)/2);
% rectangleLogo(startPoint:startPoint+length(logoData)-1,:) = logoData;
% clockOutputs = clock;
% currentMonth = clockOutputs(2);
% colorRatio = 0.12221946;
% if currentMonth>=3 & currentMonth<=5
%     originalColorMap = colormap(spring(256));
% elseif currentMonth>=6 & currentMonth<=9
%     originalColorMap = colormap(summer(256));
% elseif currentMonth>=10 & currentMonth<=11
%     originalColorMap = colormap(autumn(256)); 
% else
%     originalColorMap = colormap(winter(256));
% end
% customColorMap   = originalColorMap*colorRatio + repmat(1-colorRatio, size(originalColorMap));
% colormap(customColorMap)
% 
% axes(handles.coronalAxes);
% imagesc(rectangleLogo);
% set(gca, 'XTick', [], 'YTick', [])
% freezeColors(handles.coronalAxes)
% 
% axes(handles.sagittalAxes);
% imagesc(rectangleLogo);
% set(gca, 'XTick', [], 'YTick', [])
% freezeColors(handles.sagittalAxes)
% 
% axes(handles.axialAxes);
% imagesc(rectangleLogo);
% set(gca, 'XTick', [], 'YTick', [])
% freezeColors(handles.axialAxes)
% drawnow

% Updated 12/13/2019 Makoto.
% Display no data logo on connectivityAxes
logoPath = which('sccnLogo.jpg');
logoData = imread(logoPath);
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.connectivityAxes);
imagesc(handles.connectivityAxes, processedLogoData);
colormap(handles.connectivityAxes, customColormap);
set(handles.connectivityAxes, 'XTick', [], 'YTick', []);
drawnow

% Display no data logo on timeFrequencyAxes
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.timeFrequencyAxes);
imagesc(handles.timeFrequencyAxes, processedLogoData);
colormap(handles.timeFrequencyAxes, customColormap);
set(handles.timeFrequencyAxes, 'XTick', [], 'YTick', []);
drawnow

% Display no data logo on timeFrequencyAxes
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.coronalAxes);
imagesc(handles.coronalAxes, processedLogoData);
colormap(handles.coronalAxes, customColormap);
set(handles.coronalAxes, 'XTick', [], 'YTick', []);
drawnow

% Display no data logo on timeFrequencyAxes
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.sagittalAxes);
imagesc(handles.sagittalAxes, processedLogoData);
colormap(handles.sagittalAxes, customColormap);
set(handles.sagittalAxes, 'XTick', [], 'YTick', []);
drawnow

% Display no data logo on timeFrequencyAxes
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.axialAxes);
imagesc(handles.axialAxes, processedLogoData);
colormap(handles.axialAxes, customColormap);
set(handles.axialAxes, 'XTick', [], 'YTick', []);
drawnow

% Obtain the initial connectivity matrix position and save it to the figure object
connectivityData.axesOriginalPosition = get(handles.connectivityAxes, 'position');
set(handles.calculateButton, 'UserData', connectivityData);

% Choose default command line output for pop_groupSIFT_viewResultsAndExportForMovie
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Display process start
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''6. View Results and Export for Movie'' start. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))


% UIWAIT makes pop_groupSIFT_viewResultsAndExportForMovie wait for user response (see UIRESUME)
% uiwait(handles.viewResultsFigure);


% --- Outputs from this function are returned to the command line.
function varargout = pop_groupSIFT_viewResultsAndExportForMovie_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadTstatisticsButton.
function loadTstatisticsButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadTstatisticsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% obtain *_tStatistics.mat
[loadFile, workingFolder] = uigetfile('*_tStatistics.mat');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end

% set *_tStatistics.mat to Edit Box
set(handles.loadTstatisticsEdit, 'string', [workingFolder loadFile]);



function loadTstatisticsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to loadTstatisticsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadTstatisticsEdit as text
%        str2double(get(hObject,'String')) returns contents of loadTstatisticsEdit as a double


% --- Executes during object creation, after setting all properties.
function loadTstatisticsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadTstatisticsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in loadDipolePairDensityButton.
function loadDipolePairDensityButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadDipolePairDensityButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% obtain *_dipolePairDensity.mat
[loadFile, workingFolder] = uigetfile('*_dipolePairDensity.mat');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end

% set *_dipolePairDensity.mat to Edit Box
set(handles.loadDipolePairDensityEdit, 'string', [workingFolder loadFile]);


function loadDipolePairDensityEdit_Callback(hObject, eventdata, handles)
% hObject    handle to loadDipolePairDensityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadDipolePairDensityEdit as text
%        str2double(get(hObject,'String')) returns contents of loadDipolePairDensityEdit as a double


% --- Executes during object creation, after setting all properties.
function loadDipolePairDensityEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadDipolePairDensityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function clusterLevelPvalueEdit_Callback(hObject, eventdata, handles)
% hObject    handle to clusterLevelPvalueEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clusterLevelPvalueEdit as text
%        str2double(get(hObject,'String')) returns contents of clusterLevelPvalueEdit as a double



% --- Executes during object creation, after setting all properties.
function clusterLevelPvalueEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterLevelPvalueEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calculateButton.
function calculateButton_Callback(hObject, eventdata, handles)
% hObject    handle to calculateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% display start message
disp('Now loading...')

% uncheck checkboxes
set(handles.applyMaskCheckbox, 'Value', 0);

% if connectivityAxesColorbar is present, delete it
if isfield(handles, 'connectivityAxesColorbar')
    try
        delete(handles.connectivityAxesColorbar);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Display no data logo on connectivityAxes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes(handles.connectivityAxes);
% cla(handles.connectivityAxes, 'reset');
% logoPath = which('sccnLogo.jpg'); % Scott mentioned the relation between network analysis and mission of SCCN during tea time
% logoData = imread(logoPath);
% logoData = mean(double(logoData),3);
% clockOutputs = clock;
% currentMonth = clockOutputs(2);
% colorRatio = 0.12221946;
% if currentMonth>=3 & currentMonth<=5
%     originalColorMap = colormap(spring(256));
% elseif currentMonth>=6 & currentMonth<=9
%     originalColorMap = colormap(summer(256));
% elseif currentMonth>=10 & currentMonth<=11
%     originalColorMap = colormap(autumn(256)); 
% else
%     originalColorMap = colormap(winter(256));
% end
% customColorMap   = originalColorMap*colorRatio + repmat(1-colorRatio, size(originalColorMap));
% colormap(customColorMap)
% imagesc(logoData);
% set(gca, 'XTick', [], 'YTick', [])
% 
% % display no data logo on timeFrequencyAxes
% axes(handles.timeFrequencyAxes);
% cla(handles.timeFrequencyAxes, 'reset');
% rectangleLogo = repmat(logoData(1,1), [length(logoData) round(length(logoData)*4/3)]);
% startPoint = round(size(rectangleLogo,2)/2)-round(size(logoData,2)/2);
% rectangleLogo(:,startPoint:startPoint+length(logoData)-1) = logoData;
% clockOutputs = clock;
% currentMonth = clockOutputs(2);
% colorRatio = 0.12221946;
% if currentMonth>=3 & currentMonth<=5
%     originalColorMap = colormap(spring(256));
% elseif currentMonth>=6 & currentMonth<=9
%     originalColorMap = colormap(summer(256));
% elseif currentMonth>=10 & currentMonth<=11
%     originalColorMap = colormap(autumn(256)); 
% else
%     originalColorMap = colormap(winter(256));
% end
% customColorMap = originalColorMap*colorRatio + repmat(1-colorRatio, size(originalColorMap));
% colormap(customColorMap);
% imagesc(rectangleLogo);
% freezeColors(handles.timeFrequencyAxes)
% set(gca, 'XTick', [], 'YTick', [])
% drawnow
% 
% % display no data logo on dipole plot axes
% cla(handles.axialAxes,    'reset')
% cla(handles.sagittalAxes, 'reset')
% cla(handles.coronalAxes,  'reset')
% rectangleLogo = repmat(logoData(1,1), [round(length(logoData)*6/5) length(logoData)]);
% startPoint = round(size(rectangleLogo,1)/2)-round(size(logoData,1)/2);
% rectangleLogo(startPoint:startPoint+length(logoData)-1,:) = logoData;
% clockOutputs = clock;
% currentMonth = clockOutputs(2);
% colorRatio = 0.12221946;
% if currentMonth>=3 & currentMonth<=5
%     originalColorMap = colormap(spring(256));
% elseif currentMonth>=6 & currentMonth<=9
%     originalColorMap = colormap(summer(256));
% elseif currentMonth>=10 & currentMonth<=11
%     originalColorMap = colormap(autumn(256)); 
% else
%     originalColorMap = colormap(winter(256));
% end
% customColorMap   = originalColorMap*colorRatio + repmat(1-colorRatio, size(originalColorMap));
% colormap(customColorMap)
% 
% axes(handles.coronalAxes);
% imagesc(rectangleLogo);
% set(gca, 'XTick', [], 'YTick', [])
% freezeColors(handles.coronalAxes)
% 
% axes(handles.sagittalAxes);
% imagesc(rectangleLogo);
% set(gca, 'XTick', [], 'YTick', [])
% freezeColors(handles.sagittalAxes)
% 
% axes(handles.axialAxes);
% imagesc(rectangleLogo);
% set(gca, 'XTick', [], 'YTick', [])
% freezeColors(handles.axialAxes)
% drawnow

% Updated 12/13/2019 Makoto.
% Display no data logo on connectivityAxes
logoPath = which('sccnLogo.jpg');
logoData = imread(logoPath);
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.connectivityAxes);
imagesc(handles.connectivityAxes, processedLogoData);
colormap(handles.connectivityAxes, customColormap);
set(handles.connectivityAxes, 'XTick', [], 'YTick', []);
drawnow

% Display no data logo on timeFrequencyAxes
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.timeFrequencyAxes);
imagesc(handles.timeFrequencyAxes, processedLogoData);
colormap(handles.timeFrequencyAxes, customColormap);
set(handles.timeFrequencyAxes, 'XTick', [], 'YTick', []);
drawnow

% Display no data logo on timeFrequencyAxes
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.coronalAxes);
imagesc(handles.coronalAxes, processedLogoData);
colormap(handles.coronalAxes, customColormap);
set(handles.coronalAxes, 'XTick', [], 'YTick', []);
drawnow

% Display no data logo on timeFrequencyAxes
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.sagittalAxes);
imagesc(handles.sagittalAxes, processedLogoData);
colormap(handles.sagittalAxes, customColormap);
set(handles.sagittalAxes, 'XTick', [], 'YTick', []);
drawnow

% Display no data logo on timeFrequencyAxes
[processedLogoData, customColormap] = prepareLogoData(logoData, handles.axialAxes);
imagesc(handles.axialAxes, processedLogoData);
colormap(handles.axialAxes, customColormap);
set(handles.axialAxes, 'XTick', [], 'YTick', []);
drawnow



%%%%%%%%%%%%%%%%%%%%%%
%%% Load .mat data %%%
%%%%%%%%%%%%%%%%%%%%%%
load(get(handles.loadTstatisticsEdit,      'String'));
load(get(handles.loadDipolePairDensityEdit,'String'));

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
% Obtain user-selected P-value for cluster-level correction.
clusterLevelPvalue = str2num(get(handles.clusterLevelPvalueEdit, 'String'));

% Obtain the critical Mass of Cluster from pooled surrogate results.
surroMassOfClusterMinMax = reshape(surroMassOfCluster, [size(surroMassOfCluster,1)/2 2 size(surroMassOfCluster,2)]);
surroMassOfClusterMin    = squeeze(surroMassOfClusterMinMax(:,1,:));
surroMassOfClusterMax    = squeeze(surroMassOfClusterMinMax(:,2,:));
    % 07/19/2020 Makoto.
    % From 'Matlab for Brain and Cognitive scientists' p.245, the lower and 
    % upper thresholds are deteremind by max stats and min stats separately.
    % This approach seems slightly more sensitive than pooling the two
    % distributions and apply a two-tailed test.
    %criticalMassOfCluster = prctile(surroMassOfCluster(:), [clusterLevelPvalue*100/2 100-(clusterLevelPvalue*100/2)]); % Blob sizes are nicely bounded!

% Compute critical mass of cluster.
criticalMassOfCluster    = zeros(1,2);
criticalMassOfCluster_G1 = zeros(1,2);

% ROI mode: no multiple comparison correction for graph edges.
if     get(handles.mccForGraphEdgesCheckbox, 'Value')==0
    criticalMassOfCluster(1,1) = prctile(surroMassOfClusterMin(:), clusterLevelPvalue*100); % Use one-tailed test for each (07/20/2020 Makoto)
    criticalMassOfCluster(1,2) = prctile(surroMassOfClusterMax(:), 100-clusterLevelPvalue*100);

% Omnibus correction mode: multiple comparison correction for graph edges on.
elseif get(handles.mccForGraphEdgesCheckbox, 'Value')==1
    
    % This is the rational--For real omnibus correction across all graph edges,
    % we may want to use the strongest even among the graph edges. This
    % should provide protection for all edges. However, including all the
    % graph edge extreme values should also provide a certain type of
    % protection with some limitation (but I cannot describe exactly what
    % limitation it is.) The default value of this option should be 1.
    
    surroMassOfClusterMinSorted = sort(surroMassOfClusterMin, 2, 'ascend');
    surroMassOfClusterMaxSorted = sort(surroMassOfClusterMax, 2, 'descend');
    
    % Support generalized family-wise error rate (GFWER) control with u = 1 (Groppe et al., 2011).
    if     get(handles.useGfwerCheckbox, 'Value')==0
        surroMassOfClusterMinMin = surroMassOfClusterMinSorted(:,1);
        surroMassOfClusterMaxMax = surroMassOfClusterMaxSorted(:,1);
    elseif get(handles.useGfwerCheckbox, 'Value')==1
        disp('Generalized FWER control option is on; Note at least 1 result is false positive.')
        surroMassOfClusterMinMin = surroMassOfClusterMinSorted(:,2);
        surroMassOfClusterMaxMax = surroMassOfClusterMaxSorted(:,2);
    end
    
    criticalMassOfCluster(1,1) = prctile(min(surroMassOfClusterMinMin, [], 2), clusterLevelPvalue*100);
    criticalMassOfCluster(1,2) = prctile(max(surroMassOfClusterMaxMax, [], 2), 100-clusterLevelPvalue*100);
    clear surroMassOfClusterMinSorted surroMassOfClusterMaxSorted
end

% Build the mask.
significantClusterMask = zeros(size(clusterMask));
for n = 1:size(clusterMask,1)
    
    % Identify all blobs that survive the cluster-level threthold (pooled across edges)
    tmpBlobMask    = squeeze(clusterMask(n,:,:));
    tmpTStatistics = squeeze(tStatistics(n,:,:));
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
    significantClusterMask(n,:,:) = combinedMask;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Overwrite imported data by removing non-significant edges %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Limit latency window if necessary
latencyWindowMask  = logical(zeros(size(significantClusterMask)));
latencyWindowLimit = str2num(get(handles.latencyEdit, 'String'));
if any(latencyWindowLimit)
    latencyIdx = find(latencies>=latencyWindowLimit(1) & latencies<=latencyWindowLimit(2));
    latencyWindowMask(:,:,latencyIdx) = 1;
    significantClusterMask = logical(significantClusterMask.*latencyWindowMask);
else
    latencyIdx = 1:length(latencies);
end

% Limit frequency band if necessary
freqRangeMask = logical(zeros(size(significantClusterMask)));
frequencyBandLimit = str2num(get(handles.freqEdit, 'String'));
if any(frequencyBandLimit)
    freqIdx = find(frequencies>=frequencyBandLimit(1) & frequencies<=frequencyBandLimit(2));
    freqRangeMask(:,freqIdx,:) = 1;
    significantClusterMask = logical(significantClusterMask.*freqRangeMask);
else
    freqIdx = 1:length(frequencies);
end

% Apply the inclusive mask --> this could be used for movie!
maskedT = tStatistics.*significantClusterMask;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute mass of cluster after cluster-level correction (i.e. blob rejection) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Take absolute t-scores
if     get(handles.absPosiNegaTscorePopupmenu, 'value') == 1
    massOfCluster = sum(sum(abs(maskedT),3),2);
    colorBarMinMax = [0 max(abs(massOfCluster))];
    
% Take positive t-scores
elseif get(handles.absPosiNegaTscorePopupmenu, 'value') == 2
    positveMask = maskedT>0;
    massOfCluster = sum(sum(maskedT.*positveMask,3),2);
    colorBarMinMax = [0 max(massOfCluster)];
    
% Take negative t-scores
elseif get(handles.absPosiNegaTscorePopupmenu, 'value') == 3
    negativeMask = maskedT<0;
    massOfCluster = sum(sum(maskedT.*negativeMask,3),2);
    colorBarMinMax = [min(massOfCluster) 0];
end

% Return if no significant result
if ~any(massOfCluster)
    disp('Nothing survived.')
    return
end

% Restore 2-D matrix.
connectivityMatrix1D = zeros(length(roiLabels)*length(roiLabels),1);
connectivityMatrix1D(finallySelectedEdgeIdx) = massOfCluster;
connectivityMatrix = reshape(connectivityMatrix1D, [length(roiLabels) length(roiLabels)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sort data and roiLabels so that 1) L and R are separated, 2) ROIs are sorted from frontal to occipital %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symmetricRoiCentroids_LR    = symmetricRoiCentroids([1:2:size(connectivityMatrix,1) 2:2:size(connectivityMatrix,1)],:,:);
symmetricRoiCentroids_L     = symmetricRoiCentroids_LR(1:length(symmetricRoiCentroids_LR)/2,:,:);
yCoordinate                 = symmetricRoiCentroids_L(:,2);
[~,frontOccipitalSortL_Idx] = sort(yCoordinate, 'descend');
frontOccipitalSortIdx       = [frontOccipitalSortL_Idx; frontOccipitalSortL_Idx+length(symmetricRoiCentroids_LR)/2];
originalOrder      = 1:size(connectivityMatrix,1);
LR_separateOrder   = originalOrder([1:2:size(connectivityMatrix,1) 2:2:size(connectivityMatrix,1)]);
roiOrderConvertIdx = LR_separateOrder(frontOccipitalSortIdx);
connectivityMatrixReordered  = connectivityMatrix(roiOrderConvertIdx,roiOrderConvertIdx);
roiLabelsReordered = roiLabels(roiOrderConvertIdx);

% Rename ROI labels. 12/13/2019 Makoto.
for labelIdx = 1:length(roiLabelsReordered)
    currentLabel = roiLabelsReordered{labelIdx};
    underScoreIdx = strfind(currentLabel, '_');
    currentLabel(underScoreIdx) = deal('-');
    roiLabelsReordered{labelIdx} = currentLabel;
end



    % % IMPORTANT CONFIRMATION!
    % Ground Truth: connectivityMatrix(8,1)
    % Wrong  : connectivityMatrixReordered(roiOrderConvertIdx(8), roiOrderConvertIdx(1))
    % Correct: connectivityMatrixReordered(find(roiOrderConvertIdx==8), find(roiOrderConvertIdx==1))

    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Plot connectivity matrix with abs(massOfCluster) %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes(handles.connectivityAxes)
% cla(handles.connectivityAxes, 'reset');
% imagesc(connectivityMatrixReordered, colorBarMinMax);
%     
% % Display a short colorBar without deforming the connectivity matrix plot
% customColorMap = colormap(jet(256));
% 
% if     get(handles.absPosiNegaTscorePopupmenu, 'value') == 1 | get(handles.absPosiNegaTscorePopupmenu, 'value') == 2
%     customColorMap = customColorMap(129:end, :);
%     customColorMap(1,:) = 0.7; % Gray out zero values
% elseif get(handles.absPosiNegaTscorePopupmenu, 'value') == 3
%     customColorMap = customColorMap(1:128, :);
%     customColorMap(end,:) = 0.7; % Gray out zero values
% end
% plotConnectivityButtonUserData = get(handles.calculateButton, 'UserData');
% colormap(customColorMap);
% colorbar;
% colorbarHandle = cbfreeze(handles.connectivityAxes);
% set(handles.connectivityAxes, 'Position', plotConnectivityButtonUserData.axesOriginalPosition);
% set(colorbarHandle, 'Position', [0.5935    0.0195    0.0160    0.1901], 'YLim', colorBarMinMax)
% handles.connectivityAxesColorbar = colorbarHandle;
% set(get(colorbarHandle, 'Title'), 'String', '   Summed t')
% 
% % Display anatomical labels
% set(gca, 'YTick', 1:size(connectivityMatrix,1), 'YTickLabel', roiLabelsReordered)
% XTickRoiLabels = roiLabelsReordered;
% for n = 1:length(XTickRoiLabels)
%     XTickRoiLabels{n} = [' ' XTickRoiLabels{n}];
% end
% set(gca, 'XTick', 1:size(connectivityMatrix,1), 'XTickLabel', XTickRoiLabels, 'XAxisLocation', 'top')
% rotateXLabels(gca, 90)
% set(findall(gca, '-property', 'interpreter'), 'interpreter', 'none')
% 
% % Overlay quadrant lines
% hold on
% line([0 size(connectivityMatrix,1)+0.5], [size(connectivityMatrix,1)/2+0.5 size(connectivityMatrix,1)/2+0.5], 'color', [0 0 0])
% line([size(connectivityMatrix,1)/2+0.5 size(connectivityMatrix,1)/2+0.5], [0 size(connectivityMatrix,1)+0.5], 'color', [0 0 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot connectivity matrix with abs(massOfCluster) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.connectivityAxes)
cla(handles.connectivityAxes, 'reset');
imagesc(connectivityMatrixReordered, colorBarMinMax);
    
% Display a short colorBar without deforming the connectivity matrix plot
customColorMap = colormap(jet(256));

if     get(handles.absPosiNegaTscorePopupmenu, 'value') == 1 | get(handles.absPosiNegaTscorePopupmenu, 'value') == 2
    customColorMap = customColorMap(129:end, :);
    customColorMap(1,:) = 0.7; % Gray out zero values
elseif get(handles.absPosiNegaTscorePopupmenu, 'value') == 3
    customColorMap = customColorMap(1:128, :);
    customColorMap(end,:) = 0.7; % Gray out zero values
end
plotConnectivityButtonUserData = get(handles.calculateButton, 'UserData');
colormap(handles.connectivityAxes, customColorMap);
colorbarHandle = colorbar;
set(handles.connectivityAxes, 'Position', plotConnectivityButtonUserData.axesOriginalPosition);

    % 12/13/2019 Makoto. Updated for Matlab 2017b.
    colorbarPosition = get(colorbarHandle, 'position');
    set(colorbarHandle, 'Position', [colorbarPosition(1) colorbarPosition(2) colorbarPosition(3) colorbarPosition(4)/2], 'YLim', colorBarMinMax)

handles.connectivityAxesColorbar = colorbarHandle;
set(get(colorbarHandle, 'Title'), 'String', '   Summed t')

% Display anatomical labels
set(handles.connectivityAxes, 'YTick', 1:size(connectivityMatrix,1), 'YTickLabel', roiLabelsReordered, 'FontSize', 7)
XTickRoiLabels = roiLabelsReordered;
for n = 1:length(XTickRoiLabels)
    XTickRoiLabels{n} = [' ' XTickRoiLabels{n}];
end
set(handles.connectivityAxes, 'XTick', 1:size(connectivityMatrix,1), 'XTickLabel', XTickRoiLabels, 'XAxisLocation', 'top' , 'FontSize', 7)
rotateXLabels(handles.connectivityAxes, 90)
%set(findall(gca, '-property', 'interpreter'), 'interpreter', 'none')

% Highlight labels of significant graph edges (John Iversen's idea. SMART.) 08/15/2020 Makoto. 08/20/2020 Makoto bug fixed.
try
    significantGraphEdgeAddress = find(connectivityMatrixReordered);
    [toEdges, fromEdges]        = ind2sub(size(connectivityMatrixReordered), significantGraphEdgeAddress);
    uniqueToEdges   = unique(toEdges);
    uniqueFromEdges = unique(fromEdges);
    for significantFromEdgeIdx = 1:length(uniqueFromEdges)
        currentLabel = handles.connectivityAxes.XTickLabel(uniqueFromEdges(significantFromEdgeIdx));
        handles.connectivityAxes.XTickLabel(uniqueFromEdges(significantFromEdgeIdx)) = {['\color{red}\fontsize{7}\bf' currentLabel{1,1}]};
    end
    for significantToEdgeIdx = 1:length(uniqueToEdges)
        currentLabel = handles.connectivityAxes.YTickLabel(uniqueToEdges(significantToEdgeIdx));
        handles.connectivityAxes.YTickLabel(uniqueToEdges(significantToEdgeIdx)) = {['\color{red}\fontsize{7}\bf ' currentLabel{1,1}]};
    end
catch
    warning('Highlighting significant graph edgee labels failed. Sorry John!')
end

% Overlay quadrant lines
hold on
line([0 size(connectivityMatrix,1)+0.5], [size(connectivityMatrix,1)/2+0.5 size(connectivityMatrix,1)/2+0.5], 'color', [0 0 0])
line([size(connectivityMatrix,1)/2+0.5 size(connectivityMatrix,1)/2+0.5], [0 size(connectivityMatrix,1)+0.5], 'color', [0 0 0])

% Set the colorbar font size. 12/13/2019 Makoto.
set(colorbarHandle, 'fontsize', 7)
set(get(colorbarHandle, 'Title'), 'String', '   Summed t')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Display report: Ranking of positive and negative total information inflow and outflow %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increase of information in/outflow (i.e. positive t-scores)
positiveMask = maskedT>0;
positveT = maskedT.*positiveMask;
positiveROIs = sum(sum(positveT,2),3);
posiConnMatrix1D = zeros(length(roiLabels)*length(roiLabels),1);
posiConnMatrix1D(finallySelectedEdgeIdx) = positiveROIs;
posiConnMatrix   = reshape(posiConnMatrix1D, [length(roiLabels) length(roiLabels)]);
totalPosiOutflow = sum(posiConnMatrix,1);
totalPosiInflow  = sum(posiConnMatrix,2);
[posiOutSorted, posiOutSortOrder] = sort(totalPosiOutflow, 'descend');
[posiInSorted,  posiInSortOrder]  = sort(totalPosiInflow,  'descend');
posiOutSortedLabels = roiLabels(posiOutSortOrder);
posiInSortedLabels  = roiLabels(posiInSortOrder);

% Decrease of information in/outflow (i.e. negative t-scores)
negativeMask = maskedT<0;
negativeT = maskedT.*negativeMask;
negativeROIs = sum(sum(negativeT,2),3);
negaConnMatrix1D = zeros(length(roiLabels)*length(roiLabels),1);
negaConnMatrix1D(finallySelectedEdgeIdx) = negativeROIs;
negaConnMatrix   = reshape(negaConnMatrix1D, [length(roiLabels) length(roiLabels)]);
totalPosiOutflow = sum(negaConnMatrix,1);
totalPosiInflow  = sum(negaConnMatrix,2);
[negaOutSorted, negaOutSortOrder] = sort(totalPosiOutflow, 'ascend');
[negaInSorted,  negaInSortOrder]  = sort(totalPosiInflow,  'ascend');
negaOutSortedLabels = roiLabels(negaOutSortOrder);
negaInSortedLabels  = roiLabels(negaInSortOrder);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate summary report of total information in/outflow %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Report total number of edges.
set(handles.resultReportPanel, 'Title', sprintf('Significant edges: %.0f', length(finallySelectedEdgeIdx)));

% Construct report of total information flows.
tmpStr =  'Increase in Info Outflow';
for n = 1:6
    if any(posiOutSorted(n))
        tmpSumTscore = posiOutSorted(n);
        tmpROI       = posiOutSortedLabels{n};
        tmpStr       = [tmpStr sprintf('\n  %4.0f  %s', tmpSumTscore, tmpROI)];
    else
        tmpStr = [tmpStr sprintf('\n  (none)')];
    end
end
tmpStr = [tmpStr sprintf('\n\n')];

tmpStr = [tmpStr 'Increase in Info Inflow'];
for n = 1:6
    if any(posiInSorted(n))
        tmpSumTscore = posiInSorted(n);
        tmpROI       = posiInSortedLabels{n};
        tmpStr       = [tmpStr sprintf('\n  %4.0f  %s', tmpSumTscore, tmpROI)];
    else
        tmpStr = [tmpStr sprintf('\n  (none)')];
    end
end
tmpStr = [tmpStr sprintf('\n\n')];

tmpStr = [tmpStr 'Decrease in Info Outflow'];
for n = 1:6
    if any(negaOutSorted(n))
        tmpSumTscore = negaOutSorted(n);
        tmpROI       = negaOutSortedLabels{n};
        tmpStr       = [tmpStr sprintf('\n  %4.0f  %s', tmpSumTscore, tmpROI)];
    else
        tmpStr = [tmpStr sprintf('\n  (none)')];
    end
end
tmpStr = [tmpStr sprintf('\n\n')];

tmpStr = [tmpStr 'Decrease in Info Inflow'];
for n = 1:6
    if any(negaInSorted(n))
        tmpSumTscore = negaInSorted(n);
        tmpROI       = negaInSortedLabels{n};
        tmpStr       = [tmpStr sprintf('\n  %4.0f  %s', tmpSumTscore, tmpROI)];
    else
        tmpStr = [tmpStr sprintf('\n  (none)')];
    end
end

% Display the result
set(handles.rankingText, 'String', tmpStr)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update handles structure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('tStatistics_beforeSubtraction1', 'var') & exist('tStatistics_beforeSubtraction2', 'var')
    connectivityData.tStatistics_beforeSubtraction1 = tStatistics_beforeSubtraction1;
    connectivityData.tStatistics_beforeSubtraction2 = tStatistics_beforeSubtraction2;
end
if exist('tStatistics_beforeSubtraction3', 'var') & exist('tStatistics_beforeSubtraction4', 'var')
    connectivityData.tStatistics_beforeSubtraction3 = tStatistics_beforeSubtraction3;
    connectivityData.tStatistics_beforeSubtraction4 = tStatistics_beforeSubtraction4;
end
connectivityData.tStatistics            = tStatistics;
connectivityData.finallySelectedEdgeIdx = finallySelectedEdgeIdx;
connectivityData.significantClusterMask = significantClusterMask;
connectivityData.pValues                = pValues;
connectivityData.roiOrderConvertIdx     = roiOrderConvertIdx;
connectivityData.latencyIdx             = latencyIdx;
connectivityData.freqIdx                = freqIdx;
connectivityData.meanSurvivedStatistics = massOfCluster;
connectivityData.significantClusterMask = significantClusterMask;
connectivityData.posiConnMatrix         = posiConnMatrix;
connectivityData.negaConnMatrix         = negaConnMatrix;

set(handles.connectivityAxes, 'UserData', connectivityData);
guidata(hObject, handles);

% % Lock axis color map
% freezeColors(handles.connectivityAxes);



disp('Done.')



% --- Executes on mouse motion over figure - except title and menu.
function viewResultsFigure_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to viewResultsFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% obtain cursor position within the connectivityAxes
connectivityData = get(handles.connectivityAxes, 'UserData');
if isempty(connectivityData)
    return
end
cursorCoordinate = get(handles.connectivityAxes, 'CurrentPoint');

% % display coordinate
% disp(sprintf('row==%d column=%d', round(cursorCoordinate(1,2)), round(cursorCoordinate(1,1))));

% show source and sink ROI
if cursorCoordinate(1,1)>0.5 & cursorCoordinate(1,1)<length(connectivityData.roiLabels)+0.5 & cursorCoordinate(1,2)>0.5 & cursorCoordinate(1,2)<length(connectivityData.roiLabels)+0.5
% % show source and sink ROI
% if cursorCoordinate(1,1)>0.5 & cursorCoordinate(1,1)<72.5 & cursorCoordinate(1,2)>0.5 & cursorCoordinate(1,2)<72.5
   
    % obtain roi label list (if empty, return)
    connectivityData = get(handles.connectivityAxes, 'UserData');
    if isempty(connectivityData)
        return
    end
    
    % Obtain roi index
    fromIdx = round(cursorCoordinate(1,1));
    toIdx   = round(cursorCoordinate(1,2));
    
    % Reorder roiLabels
    sortedRoiLabels = connectivityData.roiLabels(connectivityData.roiOrderConvertIdx);
    
    % Obtain mass of cluster for positive and negative values separately.
    posiMassOfCluster = connectivityData.posiConnMatrix(connectivityData.roiOrderConvertIdx(toIdx), connectivityData.roiOrderConvertIdx(fromIdx));
    negaMassOfCluster = connectivityData.negaConnMatrix(connectivityData.roiOrderConvertIdx(toIdx), connectivityData.roiOrderConvertIdx(fromIdx));

    % Display anatomical ROIs
    set(handles.fromText, 'String', ['From:' sortedRoiLabels{fromIdx}]);
    set(handles.toText,   'String', ['To:' sortedRoiLabels{toIdx}]);
    if any(posiMassOfCluster) | any(negaMassOfCluster)
        set(handles.zText, 'String', sprintf('Summed t (increase): %.0f\n                (decrease): %.f', posiMassOfCluster, negaMassOfCluster));
    else
        set(handles.zText, 'String', sprintf('t-score: n.s.'));
    end
end



% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function viewResultsFigure_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to viewResultsFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% obtain cursor position within the connectivityAxes
cursorCoordinate = get(handles.connectivityAxes, 'CurrentPoint');
connectivityData = get(handles.connectivityAxes, 'UserData');

% show source and sink ROI
if cursorCoordinate(1,1)>0.5 & cursorCoordinate(1,1)<length(connectivityData.roiLabels)+0.5 & cursorCoordinate(1,2)>0.5 & cursorCoordinate(1,2)<length(connectivityData.roiLabels)+0.5
    
    % obtain roi index
    fromIdxIdx = round(cursorCoordinate(1,1));
    toIdxIdx   = round(cursorCoordinate(1,2));
    fromIdx    = connectivityData.roiOrderConvertIdx(fromIdxIdx);
    toIdx      = connectivityData.roiOrderConvertIdx(toIdxIdx);
    
    disp(sprintf('Selected edge: From %s to %s', connectivityData.roiLabels{fromIdx}, connectivityData.roiLabels{toIdx}))
    
    % Obtain finallySelectedEdgeIdx
    finallySelectedEdgeIdx = connectivityData.finallySelectedEdgeIdx;
    
    % Convert from/toIdx to finallySelectedEdgeIdx
    pointerToOriginalMatrixAddress = sub2ind([length(connectivityData.roiLabels) length(connectivityData.roiLabels)], toIdx, fromIdx);
    ponterIdxToTstatistics         = find(connectivityData.finallySelectedEdgeIdx == pointerToOriginalMatrixAddress);

    % Exit if pointing to masked edge
    if isempty(ponterIdxToTstatistics)
        disp('Reference to a masked edge.')
        return
    end
    
    % Extract edge time-freq data
    tStatsTimeFreq = squeeze(connectivityData.tStatistics(ponterIdxToTstatistics, :, :));
    
    % If it is a row vector, transpose it into a column vector. 06/28/2020 Makoto.
    if size(tStatsTimeFreq,1)==1
        tStatsTimeFreq = tStatsTimeFreq';
    end
        
    % Mask it with frequency and time
    userSpecifiedMask = zeros(size(tStatsTimeFreq));
    userSpecifiedMask(connectivityData.freqIdx, connectivityData.latencyIdx) = 1;
    tStatsTimeFreq = tStatsTimeFreq.*userSpecifiedMask;

    % Create masks
    significantClusterMask = squeeze(connectivityData.significantClusterMask(ponterIdxToTstatistics,:,:));
    
    % If it is a row vector, transpose it into a column vector. 06/28/2020 Makoto.
    if size(significantClusterMask,1)==1
        significantClusterMask = significantClusterMask';
    end
    
    % Apply the inclusive mask --> this could be used for movie!
    maskedStats = tStatsTimeFreq.*significantClusterMask;
    
    % plot results
    axes(handles.timeFrequencyAxes);
    cla(handles.timeFrequencyAxes, 'reset');
    
    if size(tStatsTimeFreq,2) == 1 % For the case of continuous data. 06/28/2020 Makoto.
        plot(tStatsTimeFreq, 'k', 'linewidth', 2)
        xlim([1 length(tStatsTimeFreq)])
        plotRange = abs(max(tStatsTimeFreq)-min(tStatsTimeFreq));
        ylim([min(tStatsTimeFreq)-plotRange*0.05 max(tStatsTimeFreq)+plotRange*0.05])
        
        % Obtain edge information.
        connectivityData = get(handles.connectivityAxes, 'UserData');
        timeFreqTitle = ['From ' connectivityData.roiLabels{fromIdx} ' to ' connectivityData.roiLabels{toIdx}];
        freqs = connectivityData.frequencies;

        % Plot axis labels
        set(get(gca, 'XLabel'), 'string', 'Frequency (Hz)')
        [~,idx1] = min(abs(freqs-2));
        [~,idx2] = min(abs(freqs-4));
        [~,idx3] = min(abs(freqs-8));
        [~,idx4] = min(abs(freqs-13));
        [~,idx5] = min(abs(freqs-20));
        [~,idx6] = min(abs(freqs-30));
        [~,idx7] = min(abs(freqs-50));
        set(gca, 'XTick', [idx1 idx2 idx3 idx4 idx5 idx6 idx7], 'XTickLabel', [2 4 8 13 20 30 50])
        set(get(gca, 'YLabel'), 'string', 't-Statistics')
        set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
        
        % Plot title
        title(timeFreqTitle, 'interpreter', 'none', 'fontsize', 10)
        
        
    else
        colormap(handles.timeFrequencyAxes,jet(256))
        colorScale = [-max(abs(tStatsTimeFreq(:))) max(abs(tStatsTimeFreq(:)))];
        switch get(handles.applyMaskCheckbox, 'Value')
            case 0
                imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), tStatsTimeFreq, colorScale); axis xy
            case 1
                imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), maskedStats, colorScale);   axis xy
        end
    end
else
    disp('Selected point is outside the connectivity matrix.')
    return
end

if size(tStatsTimeFreq,2) > 1 % For the case of continuous data. 06/28/2020 Makoto.
    
    % Edit YTick. One decimal display when < 10 Hz. 12/13/2019 Makoto.
    % Worked on this again. 02/25/2020 Makoto.
    YTickInterval = round(length(connectivityData.frequencies)/7);
    YTickIdx      = [1:YTickInterval:length(connectivityData.frequencies)];
    YTickLabel_original            = connectivityData.frequencies(YTickIdx);
    YTickLabel_rounded_zeroDecimal = round(YTickLabel_original);
    YTickLabel_rounded_oneDecimal  = (round(YTickLabel_original*10))/10;
    largerThanTenIdx = find(YTickLabel_original>10);
    YTickFreq = YTickLabel_rounded_oneDecimal;
    YTickFreq(largerThanTenIdx) = YTickLabel_rounded_zeroDecimal(largerThanTenIdx);
    set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickFreq);
    
    % Add vertical bar for latency zero
    if min(connectivityData.latencies)<0 & max(connectivityData.latencies)>0
        zeroLatencyIdx = find(connectivityData.latencies>0,1);
        hold on
        line([0 0], ylim, 'color', [0 0 0], 'linewidth', 2)
    end
    
    % plot title and labels
    connectivityData = get(handles.connectivityAxes, 'UserData');
    timeFreqTitle = ['From ' connectivityData.roiLabels{fromIdx} ' to ' connectivityData.roiLabels{toIdx}];
    
    % Plot axis labels
    set(get(gca, 'XLabel'), 'string', 'Latency (s)')
    set(get(gca, 'YLabel'), 'string', 'Frequency (Hz)')
    set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
    
    % Plot title
    title(timeFreqTitle, 'interpreter', 'none', 'fontsize', 10)
    
    % set(gcf, 'ToolBar', 'figure')
    
    % Display a colorBar without deforming the connectivity matrix plot. 12/13/2019 Updated. Makoto.
    originalPosition = get(handles.timeFrequencyAxes, 'Position');
    colorbarHandle = colorbar;
    set(handles.timeFrequencyAxes, 'Position', originalPosition);
    colorbarPosition = get(colorbarHandle, 'position');
    set(colorbarHandle, 'position', [colorbarPosition(1) colorbarPosition(2) colorbarPosition(3) colorbarPosition(4)/1.1])
    set(get(colorbarHandle, 'Title'), 'String', 't-score')
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show dipole plot %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% obtain user input
symmetricRoiCentroids = connectivityData.symmetricRoiCentroids;

% clear dipole axes
cla(handles.axialAxes,    'reset')
cla(handles.sagittalAxes, 'reset')
cla(handles.coronalAxes,  'reset')

% set 2 dipoles
sources(1,1).posxyz = symmetricRoiCentroids(fromIdx,:);
sources(1,1).momxyz = [0 0 0];
sources(1,1).rv     = 0;
sources(1,2).posxyz = symmetricRoiCentroids(toIdx,:);
sources(1,2).momxyz = [0 0 0];
sources(1,2).rv     = 0;

% 2 colors names officially supported by W3C specification for HTML
colors{1,1}  = [1.00 0.66 0.76]; % source == red
colors{2,1}  = [0.66 0.76 1.00]; % sink == blue

% define colors used for dipoles
plotColor = {colors{1,1} colors{2,1}};

% plot axial
axes(handles.axialAxes)
dipplot(sources, 'projimg', 'on', 'color', plotColor,  'dipolesize', 40,...
                 'projlines', 'off', 'coordformat', 'MNI', 'spheres', 'on', 'gui', 'off', 'normlen', 'on');
set(findall(gca, '-property', 'linewidth'), 'linewidth', 3);       
view([0 90]);
hold on % add lines that connects dipoles
posxyz(1,:) = sources(1,1).posxyz;
posxyz(2,:) = sources(1,2).posxyz;
arrow3d([posxyz(1,1) posxyz(2,1)], [posxyz(1,2) posxyz(2,2)], [posxyz(1,3) posxyz(2,3)], 0.7, 3, 7, [0.66 0.66 0.34]);

% plot sagittal
axes(handles.sagittalAxes)
dipplot(sources, 'projimg', 'on', 'color', plotColor,  'dipolesize', 40,...
                 'projlines', 'off', 'coordformat', 'MNI', 'spheres', 'on', 'gui', 'off', 'normlen', 'on');
set(findall(gca, '-property', 'linewidth'), 'linewidth', 3);          
view([90 0]);
hold on % add lines that connects dipoles
posxyz(1,:) = sources(1,1).posxyz;
posxyz(2,:) = sources(1,2).posxyz;
arrow3d([posxyz(1,1) posxyz(2,1)], [posxyz(1,2) posxyz(2,2)], [posxyz(1,3) posxyz(2,3)], 0.7, 3, 7, [0.66 0.66 0.34]);

% plot coronal
axes(handles.coronalAxes)
dipplot(sources, 'projimg', 'on', 'color', plotColor,  'dipolesize', 40,...
                 'projlines', 'off', 'coordformat', 'MNI', 'spheres', 'on', 'gui', 'off', 'normlen', 'on');
set(findall(gca, '-property', 'linewidth'), 'linewidth', 3);             
set(gcf, 'color', [0.66 0.76 1]);
view([0 0]);
hold on % add lines that connects dipoles
posxyz(1,:) = sources(1,1).posxyz;
posxyz(2,:) = sources(1,2).posxyz;
arrow3d([posxyz(1,1) posxyz(2,1)], [posxyz(1,2) posxyz(2,2)], [posxyz(1,3) posxyz(2,3)], 0.7, 3, 7, [0.66 0.66 0.34]);

% disable rotate3d
rotate3d(gcf)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show extra time-frequency plots if 2x2 conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(connectivityData, 'tStatistics_beforeSubtraction1') &...
   isfield(connectivityData, 'tStatistics_beforeSubtraction2') &...
   isfield(connectivityData, 'tStatistics_beforeSubtraction3') &...
   isfield(connectivityData, 'tStatistics_beforeSubtraction4')
   
    % Check if previous figure is open and is so, delete it
    figureHandleList = findall(0, 'type', 'figure');
    for n = 1:length(figureHandleList)
        tagList{n} = get(figureHandleList(n), 'Tag');
    end
    if any(strcmp(tagList, 'viewResultsTwoConditionSubtractionFigure'))
        subFigureIdx = find(strcmp(tagList, 'viewResultsTwoConditionSubtractionFigure'));
        close(figureHandleList(subFigureIdx));
    end
    
    % Obtain time-frequency t-statistics for Condition 1 and Condition 2.
    tStatsTimeFreq1 = squeeze(connectivityData.tStatistics_beforeSubtraction1(ponterIdxToTstatistics, :, :));
    tStatsTimeFreq2 = squeeze(connectivityData.tStatistics_beforeSubtraction2(ponterIdxToTstatistics, :, :));
    tStatsTimeFreq3 = squeeze(connectivityData.tStatistics_beforeSubtraction3(ponterIdxToTstatistics, :, :));
    tStatsTimeFreq4 = squeeze(connectivityData.tStatistics_beforeSubtraction4(ponterIdxToTstatistics, :, :));
    
    % Mask it with frequency and time.
    tStatsTimeFreq1 = tStatsTimeFreq1.*userSpecifiedMask;
    tStatsTimeFreq2 = tStatsTimeFreq2.*userSpecifiedMask;
    tStatsTimeFreq3 = tStatsTimeFreq3.*userSpecifiedMask;
    tStatsTimeFreq4 = tStatsTimeFreq4.*userSpecifiedMask;
    
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     %% Extracting the masked data mean values. 06/11/2020. %%%
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     barData1 = tStatsTimeFreq1(tStatsTimeFreq1~=0);
        %     barData2 = tStatsTimeFreq2(tStatsTimeFreq2~=0);
        %     barData3 = tStatsTimeFreq3(tStatsTimeFreq3~=0);
        %     barData4 = tStatsTimeFreq4(tStatsTimeFreq4~=0);
        %     
        %     figure; set(gcf, 'color', [0.93 0.96 1])
        %     hold on
        %     plot([1 2], [mean(barData2) mean(barData1)], '.', 'markersize', 30, 'color', [0.66 0.76 1])
        %     lineHandle1 = line([1 2], [mean(barData2) mean(barData1)], 'linewidth', 3,        'color', [0.66 0.76 1])
        %     plot([1 2], [mean(barData4) mean(barData3)], '.', 'markersize', 30, 'color', [1 0.66 0.76])
        %     lineHandle2 = line([1 2], [mean(barData4) mean(barData3)], 'linewidth', 3,        'color', [1 0.66 0.76])
        %     xlim([0.5 2.5])
        % %     ylim([-2.5 0])
        % %     ylim([-2 1])
        %     ylim([-2 2])
        %     legend([lineHandle1 lineHandle2], {'TCT' 'TAU'})
        %     set(gca, 'xtick', [1 2], 'xticklabel', {'T0' 'T2'})
        %     set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)
        %     ylabel('t-scored rPDC')

    
    % Plot results
    pop_groupSIFT_viewResults_2x2ConditionSubtraction;
    set(gcf, 'Name', '(1.Top Left - 2.Top Right) - (3.Bottom Left - 4.Bottom Right); pop_groupSIFT_viewResults_2x2ConditionSubtraction')
    allAxesHandles = findobj(gcf, 'type', 'axes');
    tagList = [{get(allAxesHandles(1), 'Tag')}; {get(allAxesHandles(2), 'Tag')}; {get(allAxesHandles(3), 'Tag')}; {get(allAxesHandles(4), 'Tag')}];
    condition1AxesIdx = find(strcmp(tagList, 'condition1Axes'));
    condition2AxesIdx = find(strcmp(tagList, 'condition2Axes'));
    condition3AxesIdx = find(strcmp(tagList, 'condition3Axes'));
    condition4AxesIdx = find(strcmp(tagList, 'condition4Axes'));
    condition1AxesHandle = allAxesHandles(condition1AxesIdx);
    condition2AxesHandle = allAxesHandles(condition2AxesIdx);
    condition3AxesHandle = allAxesHandles(condition3AxesIdx);
    condition4AxesHandle = allAxesHandles(condition4AxesIdx);
    
    % Plot condition 1.
    axes(condition1AxesHandle);
    colormap(jet(256))
    colorScale = [-max(abs([tStatsTimeFreq1(:); tStatsTimeFreq2(:)])) max(abs([tStatsTimeFreq1(:); tStatsTimeFreq2(:)]))];
    imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), tStatsTimeFreq1, colorScale); axis xy
    set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickFreq); % Edit YTick. Decimal only < 10. 02/25/2019 Makoto.
    if min(connectivityData.latencies)<0 & max(connectivityData.latencies)>0
        zeroLatencyIdx = find(connectivityData.latencies>0,1);
        hold on
        line([0 0], ylim, 'color', [0 0 0], 'linewidth', 2)
    end
    set(get(gca, 'YLabel'), 'string', 'Frequency (Hz)')
    set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
    title('Condition 1: A in (A-B)-(C-D)', 'interpreter', 'none', 'fontsize', 10)
    
    % Plot condition 2.
    axes(condition2AxesHandle);
    imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), tStatsTimeFreq2, colorScale); axis xy
    set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickFreq); % Edit YTick. Decimal only < 10. 02/25/2019 Makoto.
    if min(connectivityData.latencies)<0 & max(connectivityData.latencies)>0
        zeroLatencyIdx = find(connectivityData.latencies>0,1);
        hold on
        line([0 0], ylim, 'color', [0 0 0], 'linewidth', 2)
    end
    set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
    title('Condition 2: B in (A-B)-(C-D)', 'interpreter', 'none', 'fontsize', 10)
    
    % Plot condition 3.
    axes(condition3AxesHandle);
    imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), tStatsTimeFreq3, colorScale); axis xy
    set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickFreq); % Edit YTick. Decimal only < 10. 02/25/2019 Makoto.
    if min(connectivityData.latencies)<0 & max(connectivityData.latencies)>0
        zeroLatencyIdx = find(connectivityData.latencies>0,1);
        hold on
        line([0 0], ylim, 'color', [0 0 0], 'linewidth', 2)
    end
    set(get(gca, 'XLabel'), 'string', 'Latency (s)')
    set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
    title('Condition 3: C in (A-B)-(C-D)', 'interpreter', 'none', 'fontsize', 10)
    
    % Plot condition 4.
    axes(condition4AxesHandle);
    imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), tStatsTimeFreq4, colorScale); axis xy
    set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickFreq); % Edit YTick. Decimal only < 10. 02/25/2019 Makoto.
    if min(connectivityData.latencies)<0 & max(connectivityData.latencies)>0
        zeroLatencyIdx = find(connectivityData.latencies>0,1);
        hold on
        line([0 0], ylim, 'color', [0 0 0], 'linewidth', 2)
    end
    set(get(gca, 'XLabel'), 'string', 'Latency (s)')
    set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
    title('Condition 4: D in (A-B)-(C-D)', 'interpreter', 'none', 'fontsize', 10)    
    
    
    
    % Show colorbar
    originalPosition = get(gca, 'position');
    colorbarHandle = colorbar;
    set(get(colorbarHandle, 'Title'), 'String', '  t-score')
    set(condition4AxesHandle, 'position', originalPosition)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show extra time-frequency plots if subtraction between two conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isfield(connectivityData, 'tStatistics_beforeSubtraction1') &...
       isfield(connectivityData, 'tStatistics_beforeSubtraction2')
    
    % Check if previous figure is open and is so, delete it
    figureHandleList = findall(0, 'type', 'figure');
    for n = 1:length(figureHandleList)
        tagList{n} = get(figureHandleList(n), 'Tag');
    end
    if any(strcmp(tagList, 'viewResultsTwoConditionSubtractionFigure'))
        subFigureIdx = find(strcmp(tagList, 'viewResultsTwoConditionSubtractionFigure'));
        close(figureHandleList(subFigureIdx));
    end
    
    % Obtain time-frequency t-statistics for Condition 1 and Condition 2.
    tStatsTimeFreq1 = squeeze(connectivityData.tStatistics_beforeSubtraction1(ponterIdxToTstatistics, :, :));
    tStatsTimeFreq2 = squeeze(connectivityData.tStatistics_beforeSubtraction2(ponterIdxToTstatistics, :, :));
    
    % For the case of continuous data analysis (06/28/2020)
    if size(tStatsTimeFreq1,1)==1 & size(tStatsTimeFreq2,1)==1
        tStatsTimeFreq1 = tStatsTimeFreq1';
        tStatsTimeFreq2 = tStatsTimeFreq2';
    end
    
    % Mask it with frequency and time.
    tStatsTimeFreq1 = tStatsTimeFreq1.*userSpecifiedMask;
    tStatsTimeFreq2 = tStatsTimeFreq2.*userSpecifiedMask;

    % Plot results
    pop_groupSIFT_viewResults_twoConditionSubtraction;
    set(gcf, 'Name', 'Condition 1 (Left) - Condition 2 (Right); pop_groupSIFT_viewResults_twoConditionSubtraction')
    allAxesHandles = findobj(gcf, 'type', 'axes');
    tagList = [{get(allAxesHandles(1), 'Tag')}; {get(allAxesHandles(2), 'Tag')}];
    condition1AxesIdx = find(strcmp(tagList, 'condition1Axes'));
    condition2AxesIdx = find(strcmp(tagList, 'condition2Axes'));
    condition1AxesHandle = allAxesHandles(condition1AxesIdx);
    condition2AxesHandle = allAxesHandles(condition2AxesIdx);

    axes(condition1AxesHandle);
    if size(tStatsTimeFreq1,2)==1 & size(tStatsTimeFreq2,2)==1
        plot(tStatsTimeFreq1, 'k', 'linewidth', 2)
        yLimitValues = [min([tStatsTimeFreq1;tStatsTimeFreq2]) max([tStatsTimeFreq1;tStatsTimeFreq2])];
        yRange = yLimitValues(2)-yLimitValues(1);
        yLimitValues(1) = yLimitValues(1)-yRange*0.05;
        yLimitValues(2) = yLimitValues(2)+yRange*0.05;
        xlim([1 length(tStatsTimeFreq1)])
        ylim(yLimitValues)
        set(gca, 'XTick', [idx1 idx2 idx3 idx4 idx5 idx6 idx7], 'XTickLabel', [2 4 8 13 20 30 50])
        set(get(gca, 'XLabel'), 'string', 'Frequency (Hz)')
        set(get(gca, 'YLabel'), 'string', 't-statistics') 
        set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
        title('Condition 1 (A in A-B)', 'interpreter', 'none', 'fontsize', 10)
        
        colorScale = [];
    else
        colormap(jet(256))
        colorScale = [-max(abs([tStatsTimeFreq1(:); tStatsTimeFreq2(:)])) max(abs([tStatsTimeFreq1(:); tStatsTimeFreq2(:)]))];
        imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), tStatsTimeFreq1, colorScale); axis xy
        set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickFreq); % Edit YTick. Decimal only < 10. 02/25/2019 Makoto.
        if min(connectivityData.latencies)<0 & max(connectivityData.latencies)>0
            zeroLatencyIdx = find(connectivityData.latencies>0,1);
            hold on
            line([0 0], ylim, 'color', [0 0 0], 'linewidth', 2)
        end
        set(get(gca, 'XLabel'), 'string', 'Latency (s)')
        set(get(gca, 'YLabel'), 'string', 'Frequency (Hz)')
        set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
        title('Condition 1 (A in A-B)', 'interpreter', 'none', 'fontsize', 10)
    end
    
    axes(condition2AxesHandle);
    if size(tStatsTimeFreq1,2)==1 & size(tStatsTimeFreq2,2)==1
        plot(freqs, tStatsTimeFreq2, 'k', 'linewidth', 2)
        xlim([1 length(tStatsTimeFreq2)])
        ylim(yLimitValues)
        set(gca, 'XTick', [idx1 idx2 idx3 idx4 idx5 idx6 idx7], 'XTickLabel', [2 4 8 13 20 30 50])
        set(get(gca, 'XLabel'), 'string', 'Frequency (Hz)')
        set(get(gca, 'YLabel'), 'string', 't-statistics') 
        set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
        title('Condition 2 (B in A-B)', 'interpreter', 'none', 'fontsize', 10)
    else
        imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), tStatsTimeFreq2, colorScale); axis xy
        set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickFreq); % Edit YTick. Decimal only < 10. 02/25/2019 Makoto.
        if min(connectivityData.latencies)<0 & max(connectivityData.latencies)>0
            zeroLatencyIdx = find(connectivityData.latencies>0,1);
            hold on
            line([0 0], ylim, 'color', [0 0 0], 'linewidth', 2)
        end
        set(get(gca, 'XLabel'), 'string', 'Latency (s)')
        set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
        title('Condition 2 (B in A-B)', 'interpreter', 'none', 'fontsize', 10)
    
        % Show colorbar
        originalPosition = get(gca, 'position');
        colorbarHandle = colorbar;
        set(get(colorbarHandle, 'Title'), 'String', '  t-score')
        set(condition2AxesHandle, 'position', originalPosition)
    end
end

% Store the results
timeFreqData = struct();
timeFreqData.fromIdx    = fromIdx;
timeFreqData.toIdx      = toIdx;
timeFreqData.colorScale = colorScale;
timeFreqData.unmasked   = tStatsTimeFreq;
timeFreqData.masked     = maskedStats;
timeFreqData.plotTitle  = timeFreqTitle;
timeFreqData.ponterIdxToTstatistics = ponterIdxToTstatistics;
set(handles.timeFrequencyAxes, 'UserData', timeFreqData);

% Update handles
guidata(hObject, handles); % This line should come in the end ??



% --- Executes on button press in applyMaskCheckbox.
function applyMaskCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to applyMaskCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% extract data
connectivityData = get(handles.connectivityAxes,  'UserData');
timeFreqData     = get(handles.timeFrequencyAxes, 'UserData');

% plot results
axes(handles.timeFrequencyAxes);

if length(connectivityData.latencies) == 1 % For the case of continuous data analysis.
    
    switch get(handles.applyMaskCheckbox, 'Value')
        case 0
            plot(timeFreqData.unmasked, 'k', 'linewidth', 2);
        case 1
            plotData = timeFreqData.masked;
            plotData(plotData==0) = NaN;
            plot(plotData, 'r', 'linewidth', 2);
    end
    xlim([1 length(timeFreqData.unmasked)])
    plotRange = abs(max(timeFreqData.unmasked)-min(timeFreqData.unmasked));
    ylim([min(timeFreqData.unmasked)-plotRange*0.05 max(timeFreqData.unmasked)+plotRange*0.05])
    
    % Obtain edge information.
    connectivityData = get(handles.connectivityAxes, 'UserData');
    timeFreqTitle = ['From ' connectivityData.roiLabels{timeFreqData.fromIdx} ' to ' connectivityData.roiLabels{timeFreqData.toIdx}];
    freqs = connectivityData.frequencies;
    
    % Plot axis labels
    set(get(gca, 'XLabel'), 'string', 'Frequency (Hz)')
    [~,idx1] = min(abs(freqs-2));
    [~,idx2] = min(abs(freqs-4));
    [~,idx3] = min(abs(freqs-8));
    [~,idx4] = min(abs(freqs-13));
    [~,idx5] = min(abs(freqs-20));
    [~,idx6] = min(abs(freqs-30));
    [~,idx7] = min(abs(freqs-50));
    set(gca, 'XTick', [idx1 idx2 idx3 idx4 idx5 idx6 idx7], 'XTickLabel', [2 4 8 13 20 30 50])
    set(get(gca, 'YLabel'), 'string', 't-Statistics')
    set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
    
    % Plot title
    title(timeFreqTitle, 'interpreter', 'none', 'fontsize', 10)
    
else % For time-frequency plots.
    
    colorScale = [-max(abs(timeFreqData.unmasked(:))) max(abs(timeFreqData.unmasked(:)))];
    switch get(handles.applyMaskCheckbox, 'Value')
        case 0
            imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), timeFreqData.unmasked, colorScale); axis xy
        case 1
            imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), timeFreqData.masked, colorScale);   axis xy
            %         % Plot contour
            %         imagesc(connectivityData.latencies, 1:length(connectivityData.frequencies), timeFreqData.unmasked, colorScale); axis xy
            %         [~,contHandle] = contour(connectivityData.latencies, 1:length(connectivityData.frequencies), logical(timeFreqData.masked), 1, 'k');
            %         set(contHandle, 'LineWidth', 2)
    end

    % Edit YTick. One decimal display when < 10 Hz. 12/13/2019 Makoto.
    % Worked on this again. 02/25/2020 Makoto.
    YTickInterval = round(length(connectivityData.frequencies)/7);
    YTickIdx      = [1:YTickInterval:length(connectivityData.frequencies)];
    YTickLabel_original            = connectivityData.frequencies(YTickIdx);
    YTickLabel_rounded_zeroDecimal = round(YTickLabel_original);
    YTickLabel_rounded_oneDecimal  = (round(YTickLabel_original*10))/10;
    largerThanTenIdx = find(YTickLabel_original>10);
    YTickFreq = YTickLabel_rounded_oneDecimal;
    YTickFreq(largerThanTenIdx) = YTickLabel_rounded_zeroDecimal(largerThanTenIdx);
    set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickFreq);
    
    % % edit YTick
    % YTickInterval = round(length(connectivityData.frequencies)/10);
    % YTickIdx      = YTickInterval:YTickInterval:length(connectivityData.frequencies);
    % YTickLabel    = round(connectivityData.frequencies(YTickIdx)*10)/10;
    % set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickLabel);
    
    % add vertical bar for latency zero
    if min(connectivityData.latencies)<0 & max(connectivityData.latencies)>0
        zeroLatencyIdx = find(connectivityData.latencies>0,1);
        hold on
        line([0 0], ylim, 'color', [0 0 0], 'linewidth', 2)
    end
    
    % Plot axis labels
    set(get(gca, 'XLabel'), 'string', 'Latency (s)')
    set(get(gca, 'YLabel'), 'string', 'Frequency (Hz)')
    set(findall(gca, '-property', 'fontsize'), 'fontsize', 10)
    
    % Plot title
    timeFreqTitle = timeFreqData.plotTitle;
    title(timeFreqTitle, 'interpreter', 'none', 'fontsize', 10)
end

% If subFigure is present, make it on top of others (05/31/2020 updated by Makoto)
figureHandleList = findall(0, 'type', 'figure');
for n = 1:length(figureHandleList)
    tagList{n} = get(figureHandleList(n), 'Tag');
end
if any(find(strcmp(tagList, 'viewResultsTwoConditionSubtractionFigure')));
    subFigureHandleIdx = strcmp(tagList, 'viewResultsTwoConditionSubtractionFigure');
    uistack(figureHandleList(subFigureHandleIdx), 'top')
elseif any(find(strcmp(tagList, 'viewResults2x2ConditionSubtractionFigure')))
    subFigureHandleIdx = strcmp(tagList, 'viewResults2x2ConditionSubtractionFigure');
    uistack(figureHandleList(subFigureHandleIdx), 'top')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply the mask for extra two plots for the case of subtraction. (02/25/2020 Makoto) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(connectivityData.latencies) ~= 1 % If continuous data, skip this part.
    
    axesHandleList = findall(0, 'type', 'axes');
    axesTitleList = [];
    for n = 1:length(axesHandleList)
        axesTitleList{n,1} = get(get(axesHandleList(n), 'Title'), 'String');
    end
    
    if any(strcmp(axesTitleList, 'Condition 1 (A in A-B)')) & any(strcmp(axesTitleList, 'Condition 2 (B in A-B)'))
        
        switch get(handles.applyMaskCheckbox, 'Value')
            case 0
                contourHandleList = findall(0, 'Type', 'contour');
                for contourHandleIdx = 1:length(contourHandleList)
                    delete(contourHandleList(contourHandleIdx));
                end
                
            case 1
                currentAxesIdx    = find(strcmp(axesTitleList, 'Condition 1 (A in A-B)'));
                currentAxesHandle = axesHandleList(currentAxesIdx);
                axes(currentAxesHandle)
                hold on
                significanceMask = logical(timeFreqData.masked);
                contour(currentAxesHandle, connectivityData.latencies, 1:length(connectivityData.frequencies), significanceMask, 2, 'k-')
                
                currentAxesIdx    = find(strcmp(axesTitleList, 'Condition 2 (B in A-B)'));
                currentAxesHandle = axesHandleList(currentAxesIdx);
                axes(currentAxesHandle)
                hold on
                significanceMask = logical(timeFreqData.masked);
                contour(currentAxesHandle, connectivityData.latencies, 1:length(connectivityData.frequencies), significanceMask, 2, 'k-')
        end
    end
    
    if any(     strcmp(axesTitleList, 'Condition 1: A in (A-B)-(C-D)')) & ...
            any(strcmp(axesTitleList, 'Condition 2: B in (A-B)-(C-D)')) & ...
            any(strcmp(axesTitleList, 'Condition 3: C in (A-B)-(C-D)')) & ...
            any(strcmp(axesTitleList, 'Condition 4: D in (A-B)-(C-D)'))
        
        switch get(handles.applyMaskCheckbox, 'Value')
            case 0
                contourHandleList = findall(0, 'Type', 'contour');
                for contourHandleIdx = 1:length(contourHandleList)
                    delete(contourHandleList(contourHandleIdx));
                end
            case 1
                for plotIdx = 1:4
                    switch plotIdx
                        case 1
                            currentAxesTag = 'Condition 1: A in (A-B)-(C-D)';
                        case 2
                            currentAxesTag = 'Condition 2: B in (A-B)-(C-D)';
                        case 3
                            currentAxesTag = 'Condition 3: C in (A-B)-(C-D)';
                        case 4
                            currentAxesTag = 'Condition 4: D in (A-B)-(C-D)';
                    end
                    
                    currentAxesIdx    = find(strcmp(axesTitleList, currentAxesTag));
                    currentAxesHandle = axesHandleList(currentAxesIdx);
                    axes(currentAxesHandle)
                    hold on
                    significanceMask = logical(timeFreqData.masked);
                    contour(currentAxesHandle, connectivityData.latencies, 1:length(connectivityData.frequencies), significanceMask, 2, 'k-')
                end
        end
    end
end

% pass data to the next plot
set(handles.timeFrequencyAxes, 'UserData', timeFreqData);
guidata(hObject, handles); % This line should come to the end ??



function freqEdit_Callback(hObject, eventdata, handles)
% hObject    handle to freqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqEdit as text
%        str2double(get(hObject,'String')) returns contents of freqEdit as a double


% --- Executes during object creation, after setting all properties.
function freqEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latencyEdit_Callback(hObject, eventdata, handles)
% hObject    handle to latencyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latencyEdit as text
%        str2double(get(hObject,'String')) returns contents of latencyEdit as a double


% --- Executes during object creation, after setting all properties.
function latencyEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latencyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in absPosiNegaTscorePopupmenu.
function absPosiNegaTscorePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to absPosiNegaTscorePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns absPosiNegaTscorePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from absPosiNegaTscorePopupmenu

% Refresh the screen
plotConnectivityButton_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function absPosiNegaTscorePopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to absPosiNegaTscorePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: absposinegatscorepopupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in moviePushbutton.
function moviePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to moviePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain *.set to save current connectivity data.
[loadFile, loadFolder] = uigetfile('*.set', 'Select any one of the .set files included');
if ~any(loadFolder)
    disp('Cancelled.')
    return
end

% Obtain the path to save data.
[saveFile, saveFolder] = uiputfile('', 'Select the path and enter the file name for saving');
if ~any(saveFolder)
    disp('Cancelled.')
    return
end

% Obtain data.
connectivityData       = get(handles.connectivityAxes, 'UserData');
tStatistics            = connectivityData.tStatistics;
finallySelectedEdgeIdx = connectivityData.finallySelectedEdgeIdx;
frequencies            = connectivityData.frequencies;
latencies              = connectivityData.latencies;
pValues                = connectivityData.pValues;
significantClusterMask = connectivityData.significantClusterMask;
symmetricRoiCentroids  = connectivityData.symmetricRoiCentroids;
roiLabels              = connectivityData.roiLabels;
connectivityType       = connectivityData.connectivityType;

% Populate connectivityMatrix with the masked data.
maskedTstatistics  = tStatistics.*significantClusterMask;
[toIdx, fromIdx]   = ind2sub([length(connectivityData.roiLabels) length(connectivityData.roiLabels)], finallySelectedEdgeIdx);
connectivityMatrix = zeros(length(connectivityData.roiLabels), length(connectivityData.roiLabels), size(tStatistics,2), size(tStatistics,3));
for n = 1:length(toIdx)
    connectivityMatrix(toIdx(n),fromIdx(n),:,:) = squeeze(maskedTstatistics(n,:,:));
end

% Load the donor dataset.
EEG = pop_loadset('filename', loadFile, 'filepath', loadFolder);
EEG.data   = single(zeros(size(EEG.data)));
EEG.icaact = single(zeros(size(EEG.icaact)));
disp('Time series data of the donor .set file are zeroed out to mark it as a dummy.')

% Swap dipole information.
EEG.dipfit.model(1, size(roiLabels,1)).rv = 0;
symmetricRoiCentroidsCell = num2cell(symmetricRoiCentroids, 2);
[EEG.dipfit.model.posxyz] = symmetricRoiCentroidsCell{:};
[EEG.dipfit.model.momxyz] = deal([0 0 0]);
[EEG.dipfit.model.rv]     = deal(0);

% Populate Conn and save group-mean data to .set file
EEG.CAT.curComponentNames = roiLabels;
EEG.CAT.curComps = 1:size(roiLabels,1);
EEG.CAT.nbchan = size(roiLabels,1);
if     strcmp(connectivityType, 'rPDC')
    EEG.CAT.Conn.RPDC = connectivityMatrix;
elseif strcmp(connectivityType, 'dDTF08')
    EEG.CAT.Conn.dDTF08 = connectivityMatrix;
else
    error('Who said feed me weired things?')
end
EEG.CAT.Conn.freqs = frequencies;
pop_saveset(EEG, 'filename', [saveFile '.set'], 'filepath', saveFolder);
eeglab redraw

% Display process end
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''6. View Results and Export for Movie'' finished. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))






% --- Executes when selected object is changed in directionPanel.
function directionPanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in directionPanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% Refresh the screen
plotConnectivityButton_Callback(hObject, eventdata, handles)



% --- Executes on button press in setBackgroundColorWhiteCheckbox.
function setBackgroundColorWhiteCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to setBackgroundColorWhiteCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of setBackgroundColorWhiteCheckbox

% Change background color to white
if     get(hObject,'Value') == 0
    set(gcf, 'color', [0.66 0.76 1.00])
elseif get(hObject,'Value') == 1
    set(gcf, 'color', [1.00 1.00 1.00])
end

% Apply the same if two time-freq plot window is present
if any(strcmp(get(findall(0, 'type', 'figure'), 'tag'), 'viewResultsTwoConditionSubtractionFigure'))
    twoPlotIdx = find(strcmp(get(findall(0, 'type', 'figure'), 'tag'), 'viewResultsTwoConditionSubtractionFigure'));
    figureHandles = findall(0, 'type', 'figure');
    
    if     get(hObject,'Value') == 0
        set(figureHandles(twoPlotIdx), 'color', [0.66 0.76 1.00])
    elseif get(hObject,'Value') == 1
        set(figureHandles(twoPlotIdx), 'color', [1.00 1.00 1.00])
    end
end


% --- Executes on button press in outputIndividualDataPushbutton.
function outputIndividualDataPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to outputIndividualDataPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain the current Mask.
edgeIdx          = handles.connectivityAxes.UserData.finallySelectedEdgeIdx;
significanceMask = handles.connectivityAxes.UserData.significantClusterMask; % This reflects time and freq limits. 12/13/2019 Makoto.
roiLabels        = handles.connectivityAxes.UserData.roiLabels;
latency          = handles.connectivityAxes.UserData.latencies;
freqs            = handles.connectivityAxes.UserData.frequencies;

% Load the _allSubjStack.mat
[FILENAME, PATHNAME, FILTERINDEX] = uigetfile('*_allSubjStack.mat', 'Select *_allSubjStack.mat');
disp('Loading data...')
load([PATHNAME FILENAME]);

% Decode edgeIdx to from & to anatomical ROIs.
[toIdx, fromIdx] = ind2sub(size(allConnectivityStack,1),edgeIdx);
    % [toIdx, fromIdx]   = ind2sub([length(connectivityData.roiLabels) length(connectivityData.roiLabels)], finallySelectedEdgeIdx);

% Obtain individual data.
individualSubjectConnectivity = struct();
for loopIdx = 1:length(edgeIdx)
    currentEdgeTF = squeeze(allConnectivityStack(toIdx(loopIdx), fromIdx(loopIdx), :,:,:));
    currentMask   = squeeze(significanceMask(loopIdx,:,:));
    [indexedBlobs, numBlobs] = bwlabeln(currentMask);
    for blobIdx = 1:numBlobs
        currentMask = indexedBlobs == blobIdx;
        
        % Transpose the mask if the data is continuous.
        if size(currentMask,1) == 1
            currentMask = currentMask';
        end
        
        % Obtain cluster metrics.
        timeFreqCentroid = regionprops(currentMask, 'Centroid');
        timeFreqCentroid = round(timeFreqCentroid.Centroid);
        clusterSize      = sum(currentMask(:));
        
        % Obtain blob-mean connectivity measures.
        if size(currentMask,2) == 1 % Already transposed.
            currentConn = squeeze(sum(bsxfun(@times, currentEdgeTF, currentMask)))/clusterSize;
        else
            currentConn = squeeze(sum(sum(bsxfun(@times, currentEdgeTF, currentMask))))/clusterSize;
        end
                
        % Obtain subject names.
        nonZeroIdx  = find(currentConn);
        currentSubj = fileNameList(nonZeroIdx);
        
        individualSubjectConnectivity.graphEdge(loopIdx).allBlobs           = indexedBlobs;
        individualSubjectConnectivity.graphEdge(loopIdx).blob(blobIdx).from = roiLabels{fromIdx(loopIdx)};
        individualSubjectConnectivity.graphEdge(loopIdx).blob(blobIdx).to   = roiLabels{toIdx(  loopIdx)};
        individualSubjectConnectivity.graphEdge(loopIdx).blob(blobIdx).centroidTimeS  = latency(timeFreqCentroid(1));
        individualSubjectConnectivity.graphEdge(loopIdx).blob(blobIdx).centroidFreqHz = freqs(timeFreqCentroid(2));
        individualSubjectConnectivity.graphEdge(loopIdx).blob(blobIdx).clusterSize   = clusterSize;
        individualSubjectConnectivity.graphEdge(loopIdx).blob(blobIdx).blobMeanConn  = currentConn(nonZeroIdx);
        individualSubjectConnectivity.graphEdge(loopIdx).blob(blobIdx).subjNames     = currentSubj;
    end
end

% Store the output to Matlab's 'base' workspace.
assignin('base', 'individualSubjectConnectivity', individualSubjectConnectivity);
disp('''individualSubjectConnectivity'' successfully generated as Matlab variable.')

% Preliminary solution to export to excel.
dataSheet = cell(1,1);
dataSheet{1,1} = 'from';
dataSheet{2,1} = 'to';
dataSheet{3,1} = 'centroidTimeS';
dataSheet{4,1} = 'centroidFreqHz';
dataSheet{5,1} = 'clusterSize';

columnIdx = 1; % The first column contains labels hence needs to be skipped (01/06/2020 Makoto).
for graphEdgeIdx = 1:length(individualSubjectConnectivity.graphEdge)
    for blobIdx = 1:length(individualSubjectConnectivity.graphEdge(graphEdgeIdx).blob)
        
        columnIdx = columnIdx+1;
        
        dataSheet{1,columnIdx} = individualSubjectConnectivity.graphEdge(graphEdgeIdx).blob(blobIdx).from;
        dataSheet{2,columnIdx} = individualSubjectConnectivity.graphEdge(graphEdgeIdx).blob(blobIdx).to;
        dataSheet{3,columnIdx} = individualSubjectConnectivity.graphEdge(graphEdgeIdx).blob(blobIdx).centroidTimeS;
        dataSheet{4,columnIdx} = individualSubjectConnectivity.graphEdge(graphEdgeIdx).blob(blobIdx).centroidFreqHz;
        dataSheet{5,columnIdx} = individualSubjectConnectivity.graphEdge(graphEdgeIdx).blob(blobIdx).clusterSize;
        dataLength = length(individualSubjectConnectivity.graphEdge(graphEdgeIdx).blob(blobIdx).blobMeanConn);
        dataSheet(6:6+dataLength-1,columnIdx)   = num2cell(individualSubjectConnectivity.graphEdge(graphEdgeIdx).blob(blobIdx).blobMeanConn);
        
        columnIdx = columnIdx+1;
        dataSheet(6:6+dataLength-1,columnIdx) = individualSubjectConnectivity.graphEdge(graphEdgeIdx).blob(blobIdx).subjNames;
    end
end
%assignin('base', 'dataSheet', dataSheet);

% Convert cell to table.
resultTable = cell2table(dataSheet);

% Write the table to the cwd.
writetable(resultTable, 'individualSubjectConnectivity.xlsx', 'FileType','spreadsheet', 'WriteVariableNames', false)
disp(['''individualSubjectConnectivity.xlsx'' successfully saved at ' pwd])

disp('Done.')


% --- Executes on button press in mccForGraphEdgesCheckbox.
function mccForGraphEdgesCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to mccForGraphEdgesCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.mccForGraphEdgesCheckbox, 'Value')==0
    set(handles.useGfwerCheckbox, 'Enable', 'off', 'Value', 0)
else
    set(handles.useGfwerCheckbox, 'Enable', 'on')
end

if get(handles.mccForGraphEdgesCheckbox, 'Value')==0
    msgbox(sprintf('Caution: Multiple comparison correction for graph edges is now disabled (that for time-frequency is still valid). Use it AT YOUR OWN RISK for, for example, hypothesis-driven ROI test.'))
end
% Hint: get(hObject,'Value') returns toggle state of mccForGraphEdgesCheckbox



% --- Executes on button press in useGfwerCheckbox.
function useGfwerCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to useGfwerCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useGfwerCheckbox
