% pop_groupSIFT_runSiftBatch(varargin)
%
% History
% 06/24/2020 Makoto. Supporting resting-state analysis.
% 05/10/2020 Makoto. 'rectwin' selected for 300-ms sliding window analysis.
% 02/20/2020 Makoto. Used.
% 02/05/2020 Makoto. Used.
% 12/09/2019 Makoto. Updated.
% 01/30/2018 Makoto. Modified.

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

function varargout = pop_groupSIFT_runSiftBatch(varargin)
% POP_GROUPSIFT_RUNSIFTBATCH MATLAB code for pop_groupSIFT_runSiftBatch.fig
%      POP_GROUPSIFT_RUNSIFTBATCH, by itself, creates a new POP_GROUPSIFT_RUNSIFTBATCH or raises the existing
%      singleton*.
%
%      H = POP_GROUPSIFT_RUNSIFTBATCH returns the handle to a new POP_GROUPSIFT_RUNSIFTBATCH or the handle to
%      the existing singleton*.
%
%      POP_GROUPSIFT_RUNSIFTBATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_GROUPSIFT_RUNSIFTBATCH.M with the given input arguments.
%
%      POP_GROUPSIFT_RUNSIFTBATCH('Property','Value',...) creates a new POP_GROUPSIFT_RUNSIFTBATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_groupSIFT_runSiftBatch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_groupSIFT_runSiftBatch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_groupSIFT_runSiftBatch

% Last Modified by GUIDE v2.5 24-Jun-2020 17:14:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_groupSIFT_runSiftBatch_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_groupSIFT_runSiftBatch_OutputFcn, ...
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


% --- Executes just before pop_groupSIFT_runSiftBatch is made visible.
function pop_groupSIFT_runSiftBatch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_groupSIFT_runSiftBatch (see VARARGIN)

% Choose default command line output for pop_groupSIFT_runSiftBatch
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_groupSIFT_runSiftBatch wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pop_groupSIFT_runSiftBatch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function windowLengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to windowLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowLengthEdit as text
%        str2num(get(hObject,'String')) returns contents of windowLengthEdit as a double


% --- Executes during object creation, after setting all properties.
function windowLengthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function windowStepSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to windowStepSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowStepSizeEdit as text
%        str2num(get(hObject,'String')) returns contents of windowStepSizeEdit as a double


% --- Executes during object creation, after setting all properties.
function windowStepSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowStepSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freqRangeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to freqRangeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqRangeEdit as text
%        str2num(get(hObject,'String')) returns contents of freqRangeEdit as a double


% --- Executes during object creation, after setting all properties.
function freqRangeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqRangeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freqBinEdit_Callback(hObject, eventdata, handles)
% hObject    handle to freqBinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqBinEdit as text
%        str2num(get(hObject,'String')) returns contents of freqBinEdit as a double


% --- Executes during object creation, after setting all properties.
function freqBinEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqBinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in windowingPopupmenu.
function windowingPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to windowingPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns windowingPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from windowingPopupmenu


% --- Executes during object creation, after setting all properties.
function windowingPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowingPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in optimumOrderAlgoPopupmenu.
function optimumOrderAlgoPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to optimumOrderAlgoPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns optimumOrderAlgoPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from optimumOrderAlgoPopupmenu


% --- Executes during object creation, after setting all properties.
function optimumOrderAlgoPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to optimumOrderAlgoPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in singleWindowCheckbox.
function singleWindowCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to singleWindowCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(handles.singleWindowCheckbox, 'Value')
    case 0
        set(handles.windowLengthEdit,   'Enable', 'on')
        set(handles.windowStepSizeEdit, 'Enable', 'on')
    case 1
        set(handles.windowLengthEdit,   'Enable', 'off')
        set(handles.windowStepSizeEdit, 'Enable', 'off')
end



% --- Executes on button press in selectFilesButton.
function selectFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain multiple .set files
[allFiles, workingFolder] = uigetfile('*.set', 'MultiSelect', 'on');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end

% Convert char to cell if n = 1.
if ischar(allFiles)
    subjName = allFiles;
    clear allFiles
    allFiles{1,1} = subjName;    
end

% Display process start
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''1.Run SIFT batch'' started. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))

% Move to the working folder
cd(workingFolder)

% Start the batch
for setIdx = 1:length(allFiles)
    tic
    loadName = allFiles{setIdx};
    dataName = loadName(1:end-4);
    
    %% STEP1: load data
    EEG = pop_loadset('filename', loadName, 'filepath', workingFolder);

    %% STEP 2: Define key Processing Parameters
    Components        = 1:size(EEG.icaweights);       % these are the components/channels to which we'll fit our multivariate model
    
    switch get(handles.singleWindowCheckbox, 'Value')
        case 0
            WindowLengthSec   = str2num(get(handles.windowLengthEdit, 'String'));   % sliding window length in seconds
            WindowStepSizeSec = str2num(get(handles.windowStepSizeEdit, 'String')); % sliding window step size in seconds
        case 1
            WindowLengthSec   = EEG.xmax;
            WindowStepSizeSec = 1;
    end
    NewSamplingRate   = [];                           % new sampling rate (if downsampling)
    EpochTimeRange    = [EEG.xmin EEG.xmax];          % this is the time range (in seconds) to analyze (relative to event at t=0)
    GUI_MODE          = 'nogui';                      % whether or not to show the Graphical User Interfaces. Can be 'nogui' or anything else (to show the gui)
    VERBOSITY_LEVEL   = 1;                            % Verbosity Level (0=no/minimal output, 2=graphical output)
    freqEdges         = str2num(get(handles.freqRangeEdit, 'String'));
    
    % This is a modified version of log scale by Nattapong and Makoto (03/04/2017)
    deviationFromLog = 5;
    freqs            = logspace(log10(freqEdges(1)+deviationFromLog), log10(freqEdges(2)+deviationFromLog), str2num(get(handles.freqBinEdit, 'String')))-deviationFromLog;
    %freqs = 10:50;
    
    %% STEP 3: Pre-process the data
    disp('===================================')
    disp('PRE-PROCESSING DATA');
    
    % select time range
    EEG = pop_select( EEG,'time',EpochTimeRange );

    % convert list of components to cell array of strings
    ComponentNames = strtrim(cellstr(num2str(Components')));
    
    % apply the command to pre-process the data
    EEG = pop_pre_prepData(EEG,GUI_MODE, ...
        'VerbosityLevel',VERBOSITY_LEVEL,   ...
        'SignalType',{'Components'},  ...
        'VariableNames',ComponentNames,   ...
        'NormalizeData',    ...
        {'verb' 0       ...
        'method' {'time' 'ensemble'}},   ...
         'Detrend',  ...
        {'verb' VERBOSITY_LEVEL ...
        'method' {'linear'}},  ...
        'resetConfigs',true,    ...
        'badsegments',[],       ...
        'newtrials',[],         ...
        'equalizetrials',false);
    
    disp('===================================')
    
%         'Detrend',  ...
%         {'verb' VERBOSITY_LEVEL ...
%         'method' {'linear'} ...
%         'piecewise' ...
%         {'seglength' 1.00   ...
%         'stepsize' 0.25} ...
%         'plot' false},  ...

    %% STEP 4: Identify the optimal model order
    disp('===================================')
    disp('MODEL ORDER IDENTIFICATION');
    
    % Here we compute various model order selection criteria for varying model
    % orders (e.g. 1 to 30) and visualize the results
    
    % Obtain user input for windowing type.
    switch get(handles.windowingPopupmenu, 'Value')
        case 1
            windowingType = 'rectwin';
        case 2
            windowingType = 'hamming';
        case 3
            windowingType = 'blackmanharris';
    end
    
    % compute model order selection criteria...
    EEG = pop_est_selModelOrder(EEG,GUI_MODE, ...
        'modelingApproach',         ...
        {'Segmentation VAR'     ...
        'algorithm' {'Vieira-Morf'} ...
        'winStartIdx' []    ...
        'winlen'  WindowLengthSec    ...
        'winstep' WindowStepSizeSec  ...
        'taperfcn' windowingType  ... % 'rectwin' 'hamming' 'blackmanharris'
        'epochTimeLims' []      ...
        'prctWinToSample' 100   ...
        'normalize' {'method' {'time' 'ensemble'}} ...
        'detrend' {'method' 'linear'} ...
        'verb' VERBOSITY_LEVEL},      ...
        'morderRange',[1 20] ,  ...
        'downdate',true,        ...
        'runPll',[],            ...
        'icselector',{'sbc' 'aic' 'fpe' 'hq'},  ...
        'winStartIdx',[],       ...
        'epochTimeLims',[],     ...
        'prctWinToSample',100,   ...
        'plot', [], ...
        'verb', VERBOSITY_LEVEL);
    
    % To plot the results, use this:
    %     handles = vis_plotOrderCriteria(EEG.CAT.IC,{'conditions' []    ...
    %         'icselector' {'sbc','aic','fpe','hq'}  ...
    %         'minimizer' {'min'} ...
    %         'prclim' 90});
    
    % If you want to save this figure you can uncomment the following lines:
    %
    % for i=1:length(handles)
    %     saveas(handles(i),sprintf('orderResults%d.fig',i));
    % end
    % close(handles);
    
    % Finally, we can automatically select the model order which minimizes one
    % of the criteria (or you can set this manually based on above figure)
    
    % Obtain user input for windowing type.
    switch get(handles.optimumOrderAlgoPopupmenu, 'Value')
        case 1
            ModelOrder = ceil(mean(EEG.CAT.IC.hq.popt));
        case 2
            ModelOrder = ceil(mean(EEG.CAT.IC.hq.pelbow));
    end
    
    % As an alternative to using the minimum of the selection criteria over
    % model order, you can find the "elbow" in the plot of model order versus
    % selection criterion value. This is useful in cases where the selection
    % criterion does not have a clear minimum. For example, the lines below
    % plot and select the elbow location (averaged across windows) for the AIC
    % criterion
    %
    % vis_plotOrderCriteria(EEG(1).CAT.IC,{},{},'elbow');
    % ModelOrder = ceil(mean(EEG(1).CAT.IC.aic.pelbow));
    
    disp('===================================')
    
    %% STEP 5: Fit the VAR model

    disp('===================================')
    disp('MODEL FITTING');

    % Here we can check that our selected parameters make sense
    fprintf('===================================================\n');
    fprintf('MVAR PARAMETER SUMMARY FOR CONDITION: %s\n',EEG.condition);
    fprintf('===================================================\n');
    est_dispMVARParamCheck(EEG,struct('morder', ModelOrder', 'winlen', WindowLengthSec, 'winstep', WindowStepSizeSec,'verb', VERBOSITY_LEVEL));

    % Once we have identified our optimal model order, we can fit our VAR model.

    % Fit a model using the options specifed for model order selection (STEP 4)
    EEG = pop_est_fitMVAR(EEG,GUI_MODE, ...
                EEG.CAT.configs.est_selModelOrder.modelingApproach, ...
                'ModelOrder',ModelOrder);

    % Note that EEG.CAT.MODEL now contains the model structure with
    % coefficients (in MODEL.AR), prediction errors (MODEL.PE) and other
    % self-evident information

    % Alternately, we can fit the VAR parameters using a Kalman filter (see
    % doc est_fitMVARKalman for more info on arguments)
    %
    % EEG.CAT.MODEL = est_fitMVARKalman(EEG,0,'updatecoeff',0.0005,'updatemode',2,'morder',ModelOrder,'verb',2,'downsampleFactor',50);

    disp('===================================')

    %% STEP 6: Validate the fitted model

    disp('===================================')
    disp('MODEL VALIDATION');

    % Here we assess the quality of the fit of our model w.r.t. the data. This
    % step can be slow.

    % We can obtain statistics for residual whiteness, percent consistency, and
    % model stability ...
    [EEG] = pop_est_validateMVAR(EEG,GUI_MODE,...
                                'checkWhiteness', ...
                                    {'alpha' 0.05 ...
                                     'statcorrection' 'none' ...
                                     'numAcfLags' 50         ...
                                     'whitenessCriteria' {'Ljung-Box' 'ACF' 'Box-Pierce' 'Li-McLeod'} ...
                                     'winStartIdx' [] ...
                                     'prctWinToSample' 100  ...
                                     'verb' 0}, ...
                                 'checkResidualVariance',...
                                    {'alpha' 0.05 ...
                                     'statcorrection' 'none' ...
                                     'numAcfLags' 50    ...
                                     'whitenessCriteria' {}  ...
                                     'winStartIdx' []        ...
                                     'prctWinToSample' 100   ...
                                     'verb' 0}, ...
                                 'checkConsistency',    ...
                                    {'winStartIdx' []   ...
                                     'prctWinToSample' 100 ...
                                     'Nr' []                ...
                                     'donorm' 0         ...
                                     'nlags' []         ...
                                     'verb' 0}, ...
                                 'checkStability',  ...
                                    {'winStartIdx' []   ...
                                     'prctWinToSample' 100 ...
                                     'verb' 0},     ...
                                 'prctWinToSample',100,  ...
                                 'winStartIdx',[],      ...
                                 'verb',VERBOSITY_LEVEL,...
                                 'plot',false);

    %     % ... and then plot the results
    %     handles = [];
    %     for k=1:length(EEG)
    %         handles(k) = vis_plotModelValidation(EEG(k).CAT.VALIDATION.whitestats, ...
    %                                              EEG(k).CAT.VALIDATION.PC,         ...
    %                                              EEG(k).CAT.VALIDATION.stability);
    %     end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The code above has a bug. Use below instead. 01/27/2016 Makoto
    %         tmpWhitestats     = EEG.CAT.VALIDATION.whitestats;
    %         tmpPCstats        = EEG.CAT.VALIDATION.PCstats;
    %         tmpStabilitystats = EEG.CAT.VALIDATION.stabilitystats;
    %         EEG.CAT.VALIDATION.whitestats     = {};
    %         EEG.CAT.VALIDATION.PCstats        = {};
    %         EEG.CAT.VALIDATION.stabilitystats = {};
    %         EEG.CAT.VALIDATION.whitestats{1}     = tmpWhitestats;
    %         EEG.CAT.VALIDATION.PCstats{1}        = tmpPCstats;
    %         EEG.CAT.VALIDATION.stabilitystats{1} = tmpStabilitystats;
    %         
    %          vis_plotModelValidation(EEG.CAT.VALIDATION.whitestats, ...
    %                                  EEG.CAT.VALIDATION.PCstats,         ...
    %                                  EEG.CAT.VALIDATION.stabilitystats);

                         
    % If you want to save this figure you can uncomment the following lines:
    %
    % for i=1:length(handles)
    %     saveas(handles(i),sprintf('validationResults%d.fig',i));
    % end
    % close(handles);


    % To automatically determine whether our model accurately fits the data you
    % can write a few lines as follows (replace 'acf' with desired statistic):
    %
    % if ~all(EEG(1).CAT.VALIDATION.whitestats.acf.w)
    %     msgbox('Residuals are not completely white!');
    % end

    disp('===================================')


    %% STEP 7: Compute Connectivity

    disp('===================================')
    disp('CONNECTIVITY ESTIMATION');

    % Next we will compute various dynamical quantities, including connectivity,
    % from the fitted VAR model. We can compute these for a range of
    % frequencies (here 1-40 Hz). See 'doc est_mvarConnectivity' for a complete
    % list of available connectivity and spectral estimators.

    EEG = pop_est_mvarConnectivity(EEG,GUI_MODE, ...
                'connmethods',{'dDTF08' 'RPDC'}, ...
                'absvalsq',true,           ...
                'spectraldecibels',true,   ...
                'freqs', freqs,        ...
                'verb',VERBOSITY_LEVEL);

    %% STEP 8: Save data
    pop_saveset(EEG, 'filename', dataName, 'filepath', workingFolder);
    
    %% STEP 9: Report time lapse
    timeLapse = toc;
    disp(sprintf('%2.0d/%2.0d subjects done (%0.1d sec lapsed for this one)', setIdx, length(allFiles), round(timeLapse)));
end

% Display process end
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''1.Run SIFT batch'' finished. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))
