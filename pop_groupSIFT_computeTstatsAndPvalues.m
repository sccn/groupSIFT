% pop_groupSIFT_computeTstatsAndPvalues(varargin)
% 
% History
% 08/03/2020 Makoto. Continuous data check now bypasses same latency test. 
% 06/26/2020 Makoto. empty baseline supported. filesep added.
% 06/22/2020 Makoto. -1 removed from pathName = fullPath(1:end-(length('_dipolePairDensity.mat')-1));
% 05/27/2020 Makoto. fileNameList = []; added for 1x2 and 2x2.
% 05/25/2020 Makoto. Updatad. 2x2 design supported.
% 02/25/2020 Makoto. Updated. '-v7.3' supported for saving > 2GB data.
% 12/20/2019 Makoto. Updated. Copyright updated. Minor graphics fixed.
% 03/06/2019 Makoto. Start display fixed. Modified single-condition baseline test supported.

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



function varargout = pop_groupSIFT_computeTstatsAndPvalues(varargin)
% POP_GROUPSIFT_COMPUTETSTATSANDPVALUES MATLAB code for pop_groupSIFT_computeTstatsAndPvalues.fig
%      POP_GROUPSIFT_COMPUTETSTATSANDPVALUES, by itself, creates a new POP_GROUPSIFT_COMPUTETSTATSANDPVALUES or raises the existing
%      singleton*.
%
%      H = POP_GROUPSIFT_COMPUTETSTATSANDPVALUES returns the handle to a new POP_GROUPSIFT_COMPUTETSTATSANDPVALUES or the handle to
%      the existing singleton*.
%
%      POP_GROUPSIFT_COMPUTETSTATSANDPVALUES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_GROUPSIFT_COMPUTETSTATSANDPVALUES.M with the given input arguments.
%
%      POP_GROUPSIFT_COMPUTETSTATSANDPVALUES('Property','Value',...) creates a new POP_GROUPSIFT_COMPUTETSTATSANDPVALUES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_groupSIFT_computeTstatsAndPvalues_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_groupSIFT_computeTstatsAndPvalues_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_groupSIFT_computeTstatsAndPvalues

% Last Modified by GUIDE v2.5 15-Aug-2020 10:08:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_groupSIFT_computeTstatsAndPvalues_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_groupSIFT_computeTstatsAndPvalues_OutputFcn, ...
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


% --- Executes just before pop_groupSIFT_computeTstatsAndPvalues is made visible.
function pop_groupSIFT_computeTstatsAndPvalues_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_groupSIFT_computeTstatsAndPvalues (see VARARGIN)

try
    condition1Path = dir('*_allSubjStack.mat');
    set(handles.condition1PathEdit, 'String', [pwd filesep condition1Path.name]);
end

% % Gray out the Condition 2 Editbox
% set(handles.condition2PathEdit, 'Enable', 'off')
% 
% % Hide some items
% set(handles.selectNewFolderButton,     'Visible', 'off')
% set(handles.selectNewFolderEdit,       'Visible', 'off')
% set(handles.twoConditionText, 'Visible', 'off')

% Choose default command line output for pop_groupSIFT_computeTstatsAndPvalues
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_groupSIFT_computeTstatsAndPvalues wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = pop_groupSIFT_computeTstatsAndPvalues_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in condition1Button.
function condition1Button_Callback(hObject, eventdata, handles)
% hObject    handle to condition1Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain *_allSubjStack.mat
[loadFile, workingFolder] = uigetfile('*_allSubjStack.mat');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end

% set *_allSubjStack.mat to Edit Box 1
set(handles.condition1PathEdit, 'string', [workingFolder loadFile]);

% Choose default command line output for pop_groupSIFT_computeTstatsAndPvalues
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in condition2Button.
function condition2Button_Callback(hObject, eventdata, handles)
% hObject    handle to condition2Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain *_allSubjStack.mat
[loadFile, workingFolder] = uigetfile('*_allSubjStack.mat');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end

% set *_allSubjStack.mat to Edit Box 1
set(handles.condition2PathEdit, 'string', [workingFolder loadFile]);

% Choose default command line output for pop_groupSIFT_computeTstatsAndPvalues
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



function condition1PathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to condition1PathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of condition1PathEdit as text
%        str2double(get(hObject,'String')) returns contents of condition1PathEdit as a double


% --- Executes during object creation, after setting all properties.
function condition1PathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condition1PathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function condition2PathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to condition2PathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of condition2PathEdit as text
%        str2double(get(hObject,'String')) returns contents of condition2PathEdit as a double


% --- Executes during object creation, after setting all properties.
function condition2PathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condition2PathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function baselineEdit_Callback(hObject, eventdata, handles)
% hObject    handle to baselineEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baselineEdit as text
%        str2double(get(hObject,'String')) returns contents of baselineEdit as a double


% --- Executes during object creation, after setting all properties.
function baselineEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baselineEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in selectNewFolderButton.
function selectNewFolderButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectNewFolderButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% obtain working folder path
workingFolder = uigetdir;

% store the result to the Editbox
set(handles.selectNewFolderEdit, 'String', workingFolder)

function selectNewFolderEdit_Callback(hObject, eventdata, handles)
% hObject    handle to selectNewFolderEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of selectNewFolderEdit as text
%        str2double(get(hObject,'String')) returns contents of selectNewFolderEdit as a double

% --- Executes during object creation, after setting all properties.
function selectNewFolderEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectNewFolderEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function pValEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pValEdit as text
%        str2double(get(hObject,'String')) returns contents of pValEdit as a double



% --- Executes during object creation, after setting all properties.
function pValEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pValEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iterationEdit_Callback(hObject, eventdata, handles)
% hObject    handle to iterationEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterationEdit as text
%        str2double(get(hObject,'String')) returns contents of iterationEdit as a double



% --- Executes during object creation, after setting all properties.
function iterationEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iterationEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% 06/26/2020 Updated. Now takes 2-D input for resting-state data. Makoto.
function output = excludeMissingValue(input)

if length(size(input)) == 2
    
    missingValueIdx = find(sum(abs(input),1)==0);
    dataPresentIdx  = setdiff(1:size(input,2), missingValueIdx);
    output          = input(:,dataPresentIdx);
    
elseif length(size(input)) == 3
    
    % Exclude subjects with missing values.
    missingValueIdx = find(squeeze(sum(sum(abs(input),1),2))==0);
    dataPresentIdx  = setdiff(1:size(input,3), missingValueIdx);
    output          = input(:,:,dataPresentIdx);
end


% --- Executes on button press in condition3Button.
function condition3Button_Callback(hObject, eventdata, handles)
% hObject    handle to condition3Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain *_allSubjStack.mat
[loadFile, workingFolder] = uigetfile('*_allSubjStack.mat');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end

% set *_allSubjStack.mat to Edit Box 1
set(handles.condition3PathEdit, 'string', [workingFolder loadFile]);

% Choose default command line output for pop_groupSIFT_computeTstatsAndPvalues
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in condition4Button.
function condition4Button_Callback(hObject, eventdata, handles)
% hObject    handle to condition4Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain *_allSubjStack.mat
[loadFile, workingFolder] = uigetfile('*_allSubjStack.mat');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end

% set *_allSubjStack.mat to Edit Box 1
set(handles.condition4PathEdit, 'string', [workingFolder loadFile]);

% Choose default command line output for pop_groupSIFT_computeTstatsAndPvalues
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


function condition3PathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to condition3PathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of condition3PathEdit as text
%        str2double(get(hObject,'String')) returns contents of condition3PathEdit as a double


% --- Executes during object creation, after setting all properties.
function condition3PathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condition3PathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function condition4PathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to condition4PathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of condition4PathEdit as text
%        str2double(get(hObject,'String')) returns contents of condition4PathEdit as a double


% --- Executes during object creation, after setting all properties.
function condition4PathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condition4PathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function numConditionPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numConditionPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in numConditionPopupmenu.
function numConditionPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to numConditionPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(handles.numConditionPopupmenu, 'Value')
    case 1
        set(handles.conditionDisplayText,'Visible', 'off')
        set(handles.condition1Button,    'Visible', 'off')
        set(handles.condition1PathEdit,  'Visible', 'off')
        set(handles.condition2Button,    'Visible', 'off')
        set(handles.condition2PathEdit,  'Visible', 'off')
        set(handles.condition3Button,    'Visible', 'off')
        set(handles.condition3PathEdit,  'Visible', 'off')
        set(handles.condition4Button,    'Visible', 'off')
        set(handles.condition4PathEdit,  'Visible', 'off')
        set(handles.selectNewFolderButton,'Visible', 'off')
        set(handles.selectNewFolderEdit,  'Visible', 'off')
    case 2
        set(handles.conditionDisplayText,'Visible', 'on', 'String', 'Single condition; t-test against mean baseline values.')
        set(handles.condition1Button,    'Visible', 'on')
        set(handles.condition1PathEdit,  'Visible', 'on')
        set(handles.condition2Button,    'Visible', 'off')
        set(handles.condition2PathEdit,  'Visible', 'off')
        set(handles.condition3Button,    'Visible', 'off')
        set(handles.condition3PathEdit,  'Visible', 'off')
        set(handles.condition4Button,    'Visible', 'off')
        set(handles.condition4PathEdit,  'Visible', 'off')
        set(handles.selectNewFolderButton,'Visible', 'off')
        set(handles.selectNewFolderEdit,  'Visible', 'off')
    case 3
        set(handles.conditionDisplayText,'Visible', 'on', 'String', 'Subtraction 1-2; repeated or non-repeated.')
        set(handles.condition1Button,    'Visible', 'on')
        set(handles.condition1PathEdit,  'Visible', 'on')
        set(handles.condition2Button,    'Visible', 'on')
        set(handles.condition2PathEdit,  'Visible', 'on')
        set(handles.condition3Button,    'Visible', 'off')
        set(handles.condition3PathEdit,  'Visible', 'off')
        set(handles.condition4Button,    'Visible', 'off')
        set(handles.condition4PathEdit,  'Visible', 'off')
        set(handles.selectNewFolderButton,'Visible', 'on')
        set(handles.selectNewFolderEdit,  'Visible', 'on')
    case 4
        set(handles.conditionDisplayText,'Visible', 'on', 'String', 'Interaction (1-2)-(3-4); repeated, non-repeated, or mixed.')
        set(handles.condition1Button,    'Visible', 'on')
        set(handles.condition1PathEdit,  'Visible', 'on')
        set(handles.condition2Button,    'Visible', 'on')
        set(handles.condition2PathEdit,  'Visible', 'on')
        set(handles.condition3Button,    'Visible', 'on')
        set(handles.condition3PathEdit,  'Visible', 'on')
        set(handles.condition4Button,    'Visible', 'on')
        set(handles.condition4PathEdit,  'Visible', 'on')
        set(handles.selectNewFolderButton,'Visible', 'on')
        set(handles.selectNewFolderEdit,  'Visible', 'on')
end
    
    
        

% --- Executes on button press in startButton.
function startButton_Callback(hObject, eventdata, handles)
% hObject    handle to startButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Display process start
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''4.Compute t-scores and p-values'' started. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))

% Set number of iterations for permutation test
numIterations = str2num(get(handles.iterationEdit, 'String'));

switch get(handles.numConditionPopupmenu, 'Value')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Single-condition process. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2    
        % Load allSubjStack.mat and dipoleDensity.mat
        fullPath       = get(handles.condition1PathEdit, 'string');
        fullPathParsed = parsetxt(fullPath, filesep);
        fileName       = fullPathParsed{end};
        pathName       = fullPath(1:end-(length(fileName))); 
        prefixName     = fileName(1:end-length('_allSubjStack.mat'));
        disp(sprintf('Loading %s...', fullPath))
        load(fullPath);
        disp(sprintf('Loaded successfully.\n'))
        disp(sprintf('Loading %s...', [pathName prefixName '_dipolePairDensity.mat']))
        load([pathName prefixName '_dipolePairDensity.mat']); 
        disp(sprintf('Loaded successfully.\n'))
        
        % Find baseline index.
        userInputBaselinePeriod = str2num(get(handles.baselineEdit, 'String'));
        baselineIdx = find(latencies>userInputBaselinePeriod(1) & latencies<userInputBaselinePeriod(2));
        
        % Prepare a list to store processing time.
        processTimeList = zeros(length(finallySelectedEdgeIdx),1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Single-condition t-test against baseline %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % start the loop
        tStatistics        = zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack,3), size(allConnectivityStack,4));
        pValues            = zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack,3), size(allConnectivityStack,4));
        surroMassOfCluster = zeros(numIterations*2, length(finallySelectedEdgeIdx));
        clusterMask        = zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack,3), size(allConnectivityStack,4));
        disp(sprintf('\n### Starting the main loop. ###\n'))
        for edgeIdxIdx = 1:length(finallySelectedEdgeIdx)
            
            % Start measuring time
            tic;
            
            % Convert edgeIdx to toIdx and fromIdx
            [toIdx, fromIdx] = ind2sub(length(roiLabels), finallySelectedEdgeIdx(edgeIdxIdx));
            
            % Extract current time-frequency-subject matrix
            tmpConnectivity = excludeMissingValue(squeeze(allConnectivityStack(toIdx, fromIdx, :, :, :)));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Compute cluster-level multiple-comparison correction (not corrected here yet) %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            baseSubtractedConnectivity = bsxfun(@minus, tmpConnectivity, mean(tmpConnectivity(:,baselineIdx,:),2));
            
            %         % Original idea.
            %         fullSizeBaseline = baseSubtractedConnectivity(:,baselineIdx,:);
            %         fullSizeBaseline = fullSizeBaseline(:,:);
            %         fullSizeBaseline = repmat(fullSizeBaseline, [1 1 size(baseSubtractedConnectivity,2)]);
            %         fullSizeBaseline = permute(fullSizeBaseline, [1 3 2]);
            %
            %         % Perform cluster-level statistics
            %         [edgeBlobMask, edgeTScore, edgePValue, edgeSurroMassOfCluster] = ...
            %             clusterLevelPermutationTest(baseSubtractedConnectivity,... % Post-stimulus datapoints
            %             fullSizeBaseline,...           % Pre-stimulus datapoints (baseline)
            %             2,...                          % This indicates two-sample t-test with forced to be equal sample size
            %             str2num(get(handles.pValEdit,      'String')),... % uncorr. p-val threhold for preselection (i.e. selecting pixels--this determines the cluster size)
            %             numIterations); % Number of iterations.
            
            % Alternative idea: to prepare a representative of the baseline
            % distribution so that one-sample t-test is available. 03/06/2019 Makoto.
            %
            % Re-visiting this idea, just taking mean across all baseline
            % time points makes the result close to zero, which does not
            % reflect variance of the baseline period. The proposed method
            % is to preserve within- and across-subject variance of baseline
            % value in the third dimension with resolution of number of ICs
            % available. An interesting idea. 05/22/2020 Makoto.
            baselineData        = baseSubtractedConnectivity(:,baselineIdx,:);
            baseline2D_sorted   = sort(baselineData(:,:),2);
            baseline3D          = reshape(baseline2D_sorted, size(baselineData));
            baseline3D_mean     = mean(baseline3D, 2);
            fullSizeBaseline    = repmat(baseline3D_mean, [1 size(baseSubtractedConnectivity,2) 1]);
            
            % Perform cluster-level statistics
            [edgeBlobMask, edgeTScore, edgePValue, edgeSurroMassOfCluster] = ...
                clusterLevelPermutationTest(baseSubtractedConnectivity,... % Post-stimulus datapoints
                fullSizeBaseline,...           % Pre-stimulus datapoints (baseline)
                1,...                          % This indicates paired t-test.
                str2num(get(handles.pValEdit, 'String')),... % uncorr. p-val threhold for preselection (i.e. selecting pixels--this determines the cluster size)
                numIterations); % Number of iterations.
            
            % Store the results
            clusterMask(edgeIdxIdx, :, :)    = edgeBlobMask;
            tStatistics(edgeIdxIdx, :, :)    = edgeTScore;
            pValues(    edgeIdxIdx, :, :)    = edgePValue;
            surroMassOfCluster(:,edgeIdxIdx) = edgeSurroMassOfCluster(:);

            % Display time elapsed.
            timeElapsed = toc;
            processTimeList = processLoopTime(timeElapsed, edgeIdxIdx, processTimeList, finallySelectedEdgeIdx);
        end
        
        % Save results
        dimensionLabels = dimensionLabels(1,1:4);
        fileNameParsed  = parsetxt(fileName, '_');
        filePath        = get(handles.condition1PathEdit, 'String');
        savePathParsed  = parsetxt(filePath, filesep);
        savePath        = filePath(1:strfind(filePath, savePathParsed{end})-2);
        save([savePath filesep fileNameParsed{1} '_tStatistics'],...
            'tStatistics',...
            'pValues', 'finallySelectedEdgeIdx', 'connectivityType', ...
            'latencies', 'frequencies', 'dimensionLabels', 'fileNameList', 'baselineIdx', ...
            'clusterMask', 'surroMassOfCluster', '-v7.3');
     
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Two-sample t-test (either paired or unpaired) across conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
       
        % Load allSubjStack.mat and dipoleDensity.mat as data1.
        fullPath       = get(handles.condition1PathEdit, 'string');
        fullPathParsed = parsetxt(fullPath, filesep);
        fileName       = fullPathParsed{end};
        pathName       = fullPath(1:end-(length(fileName))); 
        prefixName     = fileName(1:end-length('_allSubjStack.mat'));
        disp(sprintf('Loading %s...', fullPath))
        load(fullPath);
        disp(sprintf('Loaded successfully.\n'))
        disp(sprintf('Loading %s...', [pathName prefixName '_dipolePairDensity.mat']))
        load([pathName prefixName '_dipolePairDensity.mat']); 
        disp(sprintf('Loaded successfully.\n'))
        allConnectivityStack1      = allConnectivityStack;
        connectivityType1          = connectivityType;
        dipolePairDensity1         = dipolePairDensity;
        dipoleProbabilityInRegion1 = dipoleProbabilityInRegion;
        fileNameList1              = fileNameList;
        frequencies1               = frequencies;
        workingFolder1             = pathName;
        latencies1                 = latencies;
        fileName1                  = fileName;
        edgeIdx1                   = finallySelectedEdgeIdx;
        preselectedRoiIdx1         = preselectedRoiIdx;
        
        % Load allSubjStack.mat and dipoleDensity.mat as data2.
        fullPath       = get(handles.condition2PathEdit, 'string');
        fullPathParsed = parsetxt(fullPath, filesep);
        fileName       = fullPathParsed{end};
        pathName       = fullPath(1:end-(length(fileName))); 
        prefixName     = fileName(1:end-length('_allSubjStack.mat'));
        disp(sprintf('Loading %s...', fullPath))
        load(fullPath);
        disp(sprintf('Loaded successfully.\n'))
        disp(sprintf('Loading %s...', [pathName prefixName '_dipolePairDensity.mat']))
        load([pathName prefixName '_dipolePairDensity.mat']); 
        disp(sprintf('Loaded successfully.\n'))
        allConnectivityStack2      = allConnectivityStack;
        connectivityType2          = connectivityType;
        dipolePairDensity2         = dipolePairDensity;
        dipoleProbabilityInRegion2 = dipoleProbabilityInRegion;
        fileNameList2              = fileNameList;
        frequencies2               = frequencies;
        workingFolder2             = pathName;
        latencies2                 = latencies;
        fileName2                  = fileName;
        edgeIdx2                   = finallySelectedEdgeIdx;
        preselectedRoiIdx2         = preselectedRoiIdx;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Check the overlap of the two finallySelectedEdgeIdx %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Take the overlap of the two connectivity matrices
        finallySelectedEdgeIdx = intersect(edgeIdx1, edgeIdx2);
            
        if length(edgeIdx1) == sum(ismember(edgeIdx1, edgeIdx2))
            repeatedMeasureFlag = 1;
            disp('Will perform paired t-test.')
        else
            repeatedMeasureFlag = 0;
            disp('Will perform two-sample t-test.')

            % Ask user to proceed or not with the current edges
            qstring = sprintf('%.0f edges overlapped between the two conditions. Continue?', length(finallySelectedEdgeIdx));
            userInput = questdlg(qstring, 'Confirmation');
            if ~strcmp(userInput,'Yes')
                disp('Aborted.')
                return
            end
            
            % Apply the intersectMask to allConnectivityStack to mask out edges that are unique to either of condition.
            intersectMask = zeros(size(allConnectivityStack1,1)*size(allConnectivityStack1,2),1);
            intersectMask(finallySelectedEdgeIdx) = 1;
            intersectMask = reshape(intersectMask, [size(allConnectivityStack1,1) size(allConnectivityStack1,2)]);
            allConnectivityStack1 = bsxfun(@times, allConnectivityStack1, intersectMask);
            allConnectivityStack2 = bsxfun(@times, allConnectivityStack2, intersectMask);
            dipolePairDensity1 = bsxfun(@times, dipolePairDensity1, intersectMask);
            dipolePairDensity2 = bsxfun(@times, dipolePairDensity2, intersectMask);
        end
        
        % Create the list of elapsed time.
        processTimeList = zeros(length(finallySelectedEdgeIdx),1);
        
        % Check consistency in connectivity algorithm
        connectivityTest = strcmp(connectivityType1, connectivityType2);
        if connectivityTest == 1
            disp('Check1 OK: Connectivity algorithm consistent.')
        else
            error('Check1 NG: Connectivity algorithm inconsistent.')
        end
        
        % check consistency in frequencies
        frequencyTest = ~(frequencies1 == frequencies2);
        if ~any(frequencyTest)
            disp('Check2 OK: Frequency range consistent.')
        else
            error('Check2 NG: Frequency range inconsistent.')
        end
        
        % check consistency in latency.
        latencyTest = ~(latencies1 == latencies2);
        if ~any(latencyTest)
            disp('Check3 OK: Latencies consistent.')
        else
            if size(allConnectivityStack1,4)==1
                disp('Check3 OK: Continuous data detected.')
            else
                error('Check3 NG: Latencies inconsistent.')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Unify the common parameters. %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        connectivityType  = connectivityType1;
        frequencies       = frequencies1;
        latencies         = latencies1;
        fileNameList      = [];
        fileNameList{1,1} = fileNameList1;
        fileNameList{1,2} = fileNameList2;
        dipoleProbabilityInRegion = [dipoleProbabilityInRegion1; dipoleProbabilityInRegion2];
        preselectedRoiIdx = intersect(preselectedRoiIdx1, preselectedRoiIdx2);
        clear connectivityType1 connectivityType2 frequencies1 frequencies2 latencies1 latencies2 fileNameList1 fileNameList2
        
        % Find baseline index.
        userInputBaselinePeriod = str2num(get(handles.baselineEdit, 'String'));
        if isempty(userInputBaselinePeriod)
            baselineIdx = [];
            disp(sprintf('\n\n Empty baseline specified: will perform continuous data analysis.\n\n'))
            
            % If continuous data, compute a baseline value.
            allConnCombined = cat(5, allConnectivityStack1, allConnectivityStack2);
            allConnCombined = permute(allConnCombined, [3, 5, 1, 2, 4]);
            allConnCombined = allConnCombined(:,:);
            nonZeroIdx      = find((sum(allConnCombined)~=0));
            baselineForContinuous = mean(allConnCombined(:,nonZeroIdx),2);
            clear allConnCombined nonZeroIdx
        else
            baselineIdx = find(latencies>userInputBaselinePeriod(1) & latencies<userInputBaselinePeriod(2));
        end
        
        % Start the loop.
        tStatistics = single(zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4)));
        tStatistics_beforeSubtraction1 = single(zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4)));
        tStatistics_beforeSubtraction2 = single(zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4)));
        pValues            = zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4));
        surroMassOfCluster = zeros(numIterations*2, length(finallySelectedEdgeIdx));
        clusterMask        = zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4));
        
        disp(sprintf('\n### Starting the main loop. ###\n'))
        for edgeIdxIdx = 1:length(finallySelectedEdgeIdx)
            
            % Start measuring time
            tic
            
            % Convert edgeIdx to toIdx and fromIdx
            [toIdx, fromIdx] = ind2sub(length(roiLabels), finallySelectedEdgeIdx(edgeIdxIdx));
            
            % Extract current time-frequency-subject matrix
            input1 = squeeze(allConnectivityStack1(toIdx, fromIdx, :, :, :));
            input2 = squeeze(allConnectivityStack2(toIdx, fromIdx, :, :, :));
            tmpConnectivity1 = excludeMissingValue(input1);
            tmpConnectivity2 = excludeMissingValue(input2);
            clear input1 input2
            
            if isempty(baselineIdx) % This is for resting-state analysis. 06/26/2020 Makoto.

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Compute Condition 1 and 2 t-scores against baseline %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                meanSubtracted1 = bsxfun(@minus, tmpConnectivity1, baselineForContinuous);
                [H1,P1,CI1,STATS1] = ttest(meanSubtracted1');
                tStatistics_beforeSubtraction1(edgeIdxIdx, :) = STATS1.tstat;
                
                meanSubtracted2 = bsxfun(@minus, tmpConnectivity2, baselineForContinuous);
                [H2,P2,CI2,STATS2] = ttest(meanSubtracted2');
                tStatistics_beforeSubtraction2(edgeIdxIdx, :) = STATS2.tstat;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Compute (Condition 1) - (Condition 2) t-scores %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Insert singleton dimension that represents time.
                connectivityToTest1 = permute(tmpConnectivity1, [1 3 2]);
                connectivityToTest2 = permute(tmpConnectivity2, [1 3 2]);
                
                [edgeBlobMask, edgeTScore, edgePValue, edgeSurroMassOfCluster] = ...
                    clusterLevelPermutationTest(connectivityToTest1,...
                                                connectivityToTest2,...
                                                repeatedMeasureFlag,...  % Repeated measures flag
                                                str2num(get(handles.pValEdit, 'String')),... % uncorr. p-val threhold for preselection (i.e. selecting pixels--this determines the cluster size)
                                                numIterations); % Number of iterations.
                
            else
                % Subtract mean baseline value
                baseSubtractedConnectivity1 = bsxfun(@minus, tmpConnectivity1, mean(tmpConnectivity1(:,baselineIdx,:),2));
                baseSubtractedConnectivity2 = bsxfun(@minus, tmpConnectivity2, mean(tmpConnectivity2(:,baselineIdx,:),2));
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Compute Condition 1 and 2 t-scores against baseline %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tStats = testAgainstBaseline(baseSubtractedConnectivity1, baselineIdx);
                tStatistics_beforeSubtraction1(edgeIdxIdx, :, :) = tStats;
                tStats = testAgainstBaseline(baseSubtractedConnectivity2, baselineIdx);
                tStatistics_beforeSubtraction2(edgeIdxIdx, :, :) = tStats;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Compute (Condition 1) - (Condition 2) t-scores %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [edgeBlobMask, edgeTScore, edgePValue, edgeSurroMassOfCluster] = ...
                    clusterLevelPermutationTest(baseSubtractedConnectivity1,...
                    baseSubtractedConnectivity2,...
                    repeatedMeasureFlag,...  % Repeated measures flag
                    str2num(get(handles.pValEdit, 'String')),... % uncorr. p-val threhold for preselection (i.e. selecting pixels--this determines the cluster size)
                    numIterations); % Number of iterations.
            end
            
            % Store the results
            clusterMask(edgeIdxIdx, :, :)    = edgeBlobMask;
            tStatistics(edgeIdxIdx, :, :)    = edgeTScore;
            pValues(    edgeIdxIdx, :, :)    = edgePValue;
            surroMassOfCluster(:,edgeIdxIdx) = edgeSurroMassOfCluster(:);
            
            % Display time elapsed.
            timeElapsed = toc;
            processTimeList = processLoopTime(timeElapsed, edgeIdxIdx, processTimeList, finallySelectedEdgeIdx);
        end
        
        % Save results.
        dimensionLabels = dimensionLabels(1,1:4);
        fileNameParsed1 = parsetxt(fileName1, '_');
        fileNameParsed2 = parsetxt(fileName2, '_');
        saveFileName = [fileNameParsed1{1} '_minus_' fileNameParsed2{1}];
        savePath = get(handles.selectNewFolderEdit, 'String');
        save([savePath filesep saveFileName '_tStatistics'],...
            'tStatistics', 'tStatistics_beforeSubtraction1', 'tStatistics_beforeSubtraction2',...
            'pValues', 'finallySelectedEdgeIdx', 'connectivityType', ...
            'latencies', 'frequencies', 'dimensionLabels', 'fileNameList', 'baselineIdx',...
            'clusterMask', 'surroMassOfCluster', '-v7.3');
        
        % Save the merged '_dipolePairDensity' for this new folder
        save([savePath filesep saveFileName '_dipolePairDensity'],...
            'roiLabels', 'symmetricRoiCentroids',...
            'fileNameList', 'dipolePairDensityFromIcSquareToRoiSquare', 'linearizedMeasure',...
            'preselectedRoiIdx', 'dipoleProbabilityInRegion', '-v7.3');
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interaction either paired, unpaired, or mixed. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
        
        % Load allSubjStack.mat and dipoleDensity.mat as data1.
        fullPath       = get(handles.condition1PathEdit, 'string');
        fullPathParsed = parsetxt(fullPath, filesep);
        fileName       = fullPathParsed{end};
        pathName       = fullPath(1:end-(length(fileName))); 
        prefixName     = fileName(1:end-length('_allSubjStack.mat'));
        disp(sprintf('Loading %s...', fullPath))
        load(fullPath);
        disp(sprintf('Loaded successfully.\n'))
        disp(sprintf('Loading %s...', [pathName prefixName '_dipolePairDensity.mat']))
        load([pathName prefixName '_dipolePairDensity.mat']); 
        disp(sprintf('Loaded successfully.\n'))
        allConnectivityStack1      = allConnectivityStack;
        connectivityType1          = connectivityType;
        dipolePairDensity1         = dipolePairDensity;
        dipoleProbabilityInRegion1 = dipoleProbabilityInRegion;
        fileNameList1              = fileNameList;
        frequencies1               = frequencies;
        workingFolder1             = pathName;
        latencies1                 = latencies;
        fileName1                  = fileName;
        edgeIdx1                   = finallySelectedEdgeIdx;
        preselectedRoiIdx1         = preselectedRoiIdx;
        
        % Load allSubjStack.mat and dipoleDensity.mat as data2.
        fullPath       = get(handles.condition2PathEdit, 'string');
        fullPathParsed = parsetxt(fullPath, filesep);
        fileName       = fullPathParsed{end};
        pathName       = fullPath(1:end-(length(fileName))); 
        prefixName     = fileName(1:end-length('_allSubjStack.mat'));
        disp(sprintf('Loading %s...', fullPath))
        load(fullPath);
        disp(sprintf('Loaded successfully.\n'))
        disp(sprintf('Loading %s...', [pathName prefixName '_dipolePairDensity.mat']))
        load([pathName prefixName '_dipolePairDensity.mat']); 
        disp(sprintf('Loaded successfully.\n'))
        allConnectivityStack2      = allConnectivityStack;
        connectivityType2          = connectivityType;
        dipolePairDensity2         = dipolePairDensity;
        dipoleProbabilityInRegion2 = dipoleProbabilityInRegion;
        fileNameList2              = fileNameList;
        frequencies2               = frequencies;
        workingFolder2             = pathName;
        latencies2                 = latencies;
        fileName2                  = fileName;
        edgeIdx2                   = finallySelectedEdgeIdx;
        preselectedRoiIdx2         = preselectedRoiIdx;
        
        % Load allSubjStack.mat and dipoleDensity.mat as data3.
        fullPath       = get(handles.condition3PathEdit, 'string');
        fullPathParsed = parsetxt(fullPath, filesep);
        fileName       = fullPathParsed{end};
        pathName       = fullPath(1:end-(length(fileName))); 
        prefixName     = fileName(1:end-length('_allSubjStack.mat'));
        disp(sprintf('Loading %s...', fullPath))
        load(fullPath);
        disp(sprintf('Loaded successfully.\n'))
        disp(sprintf('Loading %s...', [pathName prefixName '_dipolePairDensity.mat']))
        load([pathName prefixName '_dipolePairDensity.mat']); 
        disp(sprintf('Loaded successfully.\n'))
        allConnectivityStack3      = allConnectivityStack;
        connectivityType3          = connectivityType;
        dipolePairDensity3         = dipolePairDensity;
        dipoleProbabilityInRegion3 = dipoleProbabilityInRegion;
        fileNameList3              = fileNameList;
        frequencies3               = frequencies;
        workingFolder3             = pathName;
        latencies3                 = latencies;
        fileName3                  = fileName;
        edgeIdx3                   = finallySelectedEdgeIdx;
        preselectedRoiIdx3         = preselectedRoiIdx;
        
        % Load allSubjStack.mat and dipoleDensity.mat as data4.
        fullPath       = get(handles.condition4PathEdit, 'string');
        fullPathParsed = parsetxt(fullPath, filesep);
        fileName       = fullPathParsed{end};
        pathName       = fullPath(1:end-(length(fileName))); 
        prefixName     = fileName(1:end-length('_allSubjStack.mat'));
        disp(sprintf('Loading %s...', fullPath))
        load(fullPath);
        disp(sprintf('Loaded successfully.\n'))
        disp(sprintf('Loading %s...', [pathName prefixName '_dipolePairDensity.mat']))
        load([pathName prefixName '_dipolePairDensity.mat']); 
        disp(sprintf('Loaded successfully.\n'))
        allConnectivityStack4      = allConnectivityStack;
        connectivityType4          = connectivityType;
        dipolePairDensity4         = dipolePairDensity;
        dipoleProbabilityInRegion4 = dipoleProbabilityInRegion;
        fileNameList4              = fileNameList;
        frequencies4               = frequencies;
        workingFolder4             = pathName;
        latencies4                 = latencies;
        fileName4                  = fileName;
        edgeIdx4                   = finallySelectedEdgeIdx;
        preselectedRoiIdx4         = preselectedRoiIdx;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Check the overlap of the two finallySelectedEdgeIdx %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Take the overlap of the two connectivity matrices
        finallySelectedEdgeIdx = intersect(intersect(edgeIdx1, edgeIdx2), intersect(edgeIdx3, edgeIdx4));
        
        % All repeated measures
        if length(finallySelectedEdgeIdx) == sum(ismember(edgeIdx1, edgeIdx2)) & sum(ismember(edgeIdx1, edgeIdx3)) == sum(ismember(edgeIdx2, edgeIdx4))
            repeatedMeasureFlag = 1; % Repeated measures ANOVA.
            disp('Will perform fixed-effect test.')
            
        elseif length(edgeIdx1) == sum(ismember(edgeIdx1, edgeIdx2)) & length(edgeIdx3) == sum(ismember(edgeIdx3, edgeIdx4))
            repeatedMeasureFlag = 2; % Mixed design ANOVA.
            disp('Will perform mixed-effect test; within-subject paires are (1,2) and (3,4).')
        
        else
            repeatedMeasureFlag = 3; % Non-repeated measures ANOVA.
            disp('Will perform random-effect test.')
        end
        
        % Ask user whether to proceed with the given edges.
        qstring = sprintf('%.0f edges overlapped across conditions. Continue?', length(finallySelectedEdgeIdx));
        userInput = questdlg(qstring, 'Confirmation');
        if ~strcmp(userInput,'Yes')
            disp('Aborted.')
            return
        end
        
        % Apply the intersectMask to allConnectivityStack to mask out edges that are unique to either of condition.
        intersectMask = zeros(size(allConnectivityStack1,1)*size(allConnectivityStack1,2),1);
        intersectMask(finallySelectedEdgeIdx) = 1;
        intersectMask = reshape(intersectMask, [size(allConnectivityStack1,1) size(allConnectivityStack1,2)]);
        allConnectivityStack1 = bsxfun(@times, allConnectivityStack1, intersectMask);
        allConnectivityStack2 = bsxfun(@times, allConnectivityStack2, intersectMask);
        dipolePairDensity1 = bsxfun(@times, dipolePairDensity1, intersectMask);
        dipolePairDensity2 = bsxfun(@times, dipolePairDensity2, intersectMask);
        dipolePairDensity3 = bsxfun(@times, dipolePairDensity3, intersectMask);
        dipolePairDensity4 = bsxfun(@times, dipolePairDensity4, intersectMask);

        % Create the list of elapsed time.
        processTimeList = zeros(length(finallySelectedEdgeIdx),1);
        
        % Check consistency in connectivity algorithm
        connectivityTest = strcmp(connectivityType1, connectivityType2);
        if connectivityTest == 1
            disp('Check1 OK: Connectivity algorithm consistent.')
        else
            error('Check1 NG: Connectivity algorithm inconsistent.')
        end
        
        % check consistency in frequencies
        frequencyTest = ~(frequencies1 == frequencies2) | ~(frequencies3 == frequencies4) | ~(frequencies1 == frequencies3);
        if ~any(frequencyTest)
            disp('Check2 OK: Frequency range consistent.')
        else
            error('Check2 NG: Frequency range inconsistent.')
        end
        
        % check consistency in latency
        latencyTest = ~(latencies1 == latencies2) | ~(latencies3 == latencies4) | ~(latencies1 == latencies3);
        if ~any(latencyTest)
            disp('Check3 OK: Latencies consistent.')
        else
            error('Check3 NG: Latencies inconsistent.')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Unify the common parameters. %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        connectivityType  = connectivityType1;
        frequencies       = frequencies1;
        latencies         = latencies1;
        fileNameList      = [];
        fileNameList{1,1} = fileNameList1;
        fileNameList{1,2} = fileNameList2;
        fileNameList{1,3} = fileNameList3;
        fileNameList{1,4} = fileNameList4;        
        dipoleProbabilityInRegion = [dipoleProbabilityInRegion1; dipoleProbabilityInRegion2; dipoleProbabilityInRegion3; dipoleProbabilityInRegion4];
        preselectedRoiIdx = intersect(preselectedRoiIdx1, preselectedRoiIdx2);
        clear connectivityType1 connectivityType2 connectivityType3 connectivityType4 frequencies1 frequencies2 frequencies3 frequencies4
        clear latencies1 latencies2 latencies3 latencies4 fileNameList1 fileNameList2 fileNameList3 fileNameList4
        
        % Find baseline index.
        userInputBaselinePeriod = str2num(get(handles.baselineEdit, 'String'));
        baselineIdx = find(latencies>userInputBaselinePeriod(1) & latencies<userInputBaselinePeriod(2));
        
        % Start the loop.
        tStatistics = single(zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4)));
        tStatistics_beforeSubtraction1 = single(zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4)));
        tStatistics_beforeSubtraction2 = single(zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4)));
        tStatistics_beforeSubtraction3 = single(zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4)));
        tStatistics_beforeSubtraction4 = single(zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4)));
        pValues            = zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4));
        surroMassOfCluster = zeros(numIterations*2, length(finallySelectedEdgeIdx));
        clusterMask        = zeros(length(finallySelectedEdgeIdx), size(allConnectivityStack1,3), size(allConnectivityStack1,4));
        
        disp(sprintf('\n### Starting the main loop. ###\n'))
        for edgeIdxIdx = 1:length(finallySelectedEdgeIdx)
            
            % Start measuring time
            tic
            
            % Convert edgeIdx to toIdx and fromIdx
            [toIdx, fromIdx] = ind2sub(length(roiLabels), finallySelectedEdgeIdx(edgeIdxIdx));
            
            % Extract current time-frequency-subject matrix
            tmpConnectivity1 = excludeMissingValue(squeeze(allConnectivityStack1(toIdx, fromIdx, :, :, :)));
            tmpConnectivity2 = excludeMissingValue(squeeze(allConnectivityStack2(toIdx, fromIdx, :, :, :)));
            tmpConnectivity3 = excludeMissingValue(squeeze(allConnectivityStack3(toIdx, fromIdx, :, :, :)));
            tmpConnectivity4 = excludeMissingValue(squeeze(allConnectivityStack4(toIdx, fromIdx, :, :, :)));
                        
            % Subtract mean baseline value
            baseSubtractedConnectivity1 = bsxfun(@minus, tmpConnectivity1, mean(tmpConnectivity1(:,baselineIdx,:),2));
            baseSubtractedConnectivity2 = bsxfun(@minus, tmpConnectivity2, mean(tmpConnectivity2(:,baselineIdx,:),2));
            baseSubtractedConnectivity3 = bsxfun(@minus, tmpConnectivity3, mean(tmpConnectivity3(:,baselineIdx,:),2));
            baseSubtractedConnectivity4 = bsxfun(@minus, tmpConnectivity4, mean(tmpConnectivity4(:,baselineIdx,:),2));
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Compute Condition 1, 2, 3, 4 t-scores against baseline %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tStats = testAgainstBaseline(baseSubtractedConnectivity1, baselineIdx);
            tStatistics_beforeSubtraction1(edgeIdxIdx, :, :) = tStats;
            tStats = testAgainstBaseline(baseSubtractedConnectivity2, baselineIdx);
            tStatistics_beforeSubtraction2(edgeIdxIdx, :, :) = tStats;
            tStats = testAgainstBaseline(baseSubtractedConnectivity3, baselineIdx);
            tStatistics_beforeSubtraction3(edgeIdxIdx, :, :) = tStats;
            tStats = testAgainstBaseline(baseSubtractedConnectivity4, baselineIdx);
            tStatistics_beforeSubtraction4(edgeIdxIdx, :, :) = tStats;            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Compute (Condition 1) - (Condition 2) t-scores %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [edgeBlobMask, edgeTScore, edgePValue, edgeSurroMassOfCluster] = clusterLevelPermutationTest2x2(...
                baseSubtractedConnectivity1,...
                baseSubtractedConnectivity2,...
                baseSubtractedConnectivity3,...
                baseSubtractedConnectivity4,...
                repeatedMeasureFlag,...  % Repeated measures flag
                str2num(get(handles.pValEdit, 'String')),... % uncorr. p-val threhold for preselection (i.e. selecting pixels--this determines the cluster size)
                numIterations); % Number of iterations.
            
            % Store the results
            clusterMask(edgeIdxIdx, :, :)    = edgeBlobMask;
            tStatistics(edgeIdxIdx, :, :)    = edgeTScore;
            pValues(    edgeIdxIdx, :, :)    = edgePValue;
            surroMassOfCluster(:,edgeIdxIdx) = edgeSurroMassOfCluster(:);
            
            % Display time elapsed.
            timeElapsed = toc;
            processTimeList = processLoopTime(timeElapsed, edgeIdxIdx, processTimeList, finallySelectedEdgeIdx);
        end
        
        % Save results.
        dimensionLabels = dimensionLabels(1,1:4);
        saveFileName = 'doubleSubtraction';
        savePath = get(handles.selectNewFolderEdit, 'String');
        save([savePath filesep saveFileName '_tStatistics'],...
            'tStatistics',...
            'tStatistics_beforeSubtraction1',...
            'tStatistics_beforeSubtraction2',...
            'tStatistics_beforeSubtraction3',...
            'tStatistics_beforeSubtraction4',...
            'pValues', 'finallySelectedEdgeIdx', 'connectivityType', ...
            'latencies', 'frequencies', 'dimensionLabels', 'fileNameList', 'baselineIdx',...
            'clusterMask', 'surroMassOfCluster', '-v7.3');
        
        % Save the merged '_dipolePairDensity' for this new folder
        save([savePath filesep saveFileName '_dipolePairDensity'],...
            'roiLabels', 'symmetricRoiCentroids',...
            'fileNameList', 'dipolePairDensityFromIcSquareToRoiSquare', 'linearizedMeasure',...
            'preselectedRoiIdx', 'dipoleProbabilityInRegion', '-v7.3');
end

cd(savePath)

% Display process end
disp(newline)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''4.Compute t-scores and p-values'' finished. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))





% 05/20/2020 Makoto. Isolating this loop-time prcoessing part.
function processTimeList = processLoopTime(timeElapsed, edgeIdxIdx, processTimeList, finallySelectedEdgeIdx)

% Remove < 1 sec results (generated by non-significant results)
if timeElapsed > 1
    processTimeList(edgeIdxIdx) = timeElapsed;
end
nonzeroProcessTimes = nonzeros(processTimeList);

% Skip if nonzeroProcessTimes is empty (i.e. so far, results are all non-significant)
if isempty(nonzeroProcessTimes)
    disp(sprintf('%.0f/%.0f graph edges done.', edgeIdxIdx, length(finallySelectedEdgeIdx)));
    return
end

% Use robust measure
if length(nonzeroProcessTimes) >= 4
    robustIdx = find(nonzeroProcessTimes >= prctile(nonzeroProcessTimes,25) & nonzeroProcessTimes <= prctile(nonzeroProcessTimes,75));
    meanValueSoFar = mean(nonzeroProcessTimes(robustIdx));
else
    meanValueSoFar = mean(nonzeroProcessTimes);
end

% Display time elapsed.
timeRemaining = (length(processTimeList)-edgeIdxIdx)*meanValueSoFar;
if     timeRemaining >= 3600
    disp(sprintf('%.0f/%.0f graph edges done, maximum %.1f hours remaining.',   edgeIdxIdx, length(finallySelectedEdgeIdx), timeRemaining/3600));
elseif timeRemaining >  60
    disp(sprintf('%.0f/%.0f graph edges done, maximum %.1f minutes remaining.', edgeIdxIdx, length(finallySelectedEdgeIdx), timeRemaining/60));
else
    disp(sprintf('%.0f/%.0f graph edges done, maximum %.0f seconds remaining.', edgeIdxIdx, length(finallySelectedEdgeIdx), timeRemaining));
end





function tStats = testAgainstBaseline(input, baselineIdx) % input must be [freq x time x subjs]

% Build the data structure that represents baseline.
baselineData      = input(:,baselineIdx,:);
baseline2D_sorted = sort(baselineData(:,:),2);
baseline3D        = reshape(baseline2D_sorted, size(baselineData));
baseline3D_mean   = mean(baseline3D, 2);
fullSizeBaseline  = repmat(baseline3D_mean, [1 size(input,2) 1]);

% Perform paired t-test against baseline data.
inputA = (reshape(input, [size(input,1)*size(input,2) size(input,3)]))';
inputB = (reshape(fullSizeBaseline, [size(fullSizeBaseline,1)*size(fullSizeBaseline,2) size(fullSizeBaseline,3)]))';
[~,~,~,STATS] = ttest(inputA-inputB);   % paired t-test.
%[~,~,~,STATS] = ttest2(inputA, inputB); % two-sample t-test.
tStats = reshape(STATS.tstat, [size(input,1) size(input,2)]);            
