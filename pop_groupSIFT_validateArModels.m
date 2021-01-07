% pop_groupSIFT_validateArModels(varargin)
%
% History:
% 06/24/2020 Makoto. Compatible with single subject and/or single-window data for resting state.
% 12/09/2019 Makoto. Result graphics updated.

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

function varargout = pop_groupSIFT_validateArModels(varargin)
% POP_GROUPSIFT_VALIDATEARMODELS MATLAB code for pop_groupSIFT_validateArModels.fig
%      POP_GROUPSIFT_VALIDATEARMODELS, by itself, creates a new POP_GROUPSIFT_VALIDATEARMODELS or raises the existing
%      singleton*.
%
%      H = POP_GROUPSIFT_VALIDATEARMODELS returns the handle to a new POP_GROUPSIFT_VALIDATEARMODELS or the handle to
%      the existing singleton*.
%
%      POP_GROUPSIFT_VALIDATEARMODELS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_GROUPSIFT_VALIDATEARMODELS.M with the given input arguments.
%
%      POP_GROUPSIFT_VALIDATEARMODELS('Property','Value',...) creates a new POP_GROUPSIFT_VALIDATEARMODELS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_groupSIFT_validateArModels_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_groupSIFT_validateArModels_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_groupSIFT_validateArModels

% Last Modified by GUIDE v2.5 18-Jun-2016 00:21:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_groupSIFT_validateArModels_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_groupSIFT_validateArModels_OutputFcn, ...
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


% --- Executes just before pop_groupSIFT_validateArModels is made visible.
function pop_groupSIFT_validateArModels_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_groupSIFT_validateArModels (see VARARGIN)

% Choose default command line output for pop_groupSIFT_validateArModels
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_groupSIFT_validateArModels wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pop_groupSIFT_validateArModels_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in validateButton.
function validateButton_Callback(hObject, eventdata, handles)
% hObject    handle to validateButton (see GCBO)
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
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''2.Validate AR models'' started. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))

% Move to the working folder
cd(workingFolder)
disp(sprintf('=========== Loading one set to obtain EEG structure ==========='))
EEG = pop_loadset('filename', allFiles{1}, 'filepath', workingFolder, 'loadmode', 'info');
dataLength = size(EEG.CAT.VALIDATION.whitestats.ljungbox.pval,2);

% Separate load log (Thanks Clement Lee!)
disp(sprintf('\n\n\n'))
disp(sprintf('=========== Loading files below for validation ==========='))

% Load all .set data located under the working folder
validationResults = zeros(dataLength, 11, length(allFiles)); 
svdResults = cell(length(allFiles),1);
for setIdx = 1:length(allFiles)
    loadName = allFiles{setIdx};
    EEG = pop_loadset('filename', loadName, 'filepath', workingFolder, 'loadmode', 'info');
    validationResults(:,1,setIdx)  = repmat(EEG.trials, [size(validationResults,1),1]);
    validationResults(:,2,setIdx)  = repmat(size(EEG.icawinv,2), [size(validationResults,1),1]);
    validationResults(:,3,setIdx)  = ones(size(validationResults,1),1);
    %validationResults(:,4,n)  = EEG.CAT.IC.hq.popt';
    validationResults(:,4,setIdx)  = EEG.CAT.MODEL.morder;
    validationResults(:,5,setIdx)  = size(EEG.icaweights,1)^2*EEG.CAT.MODEL.morder / (EEG.srate*EEG.CAT.MODEL.winlen*EEG.trials); % from est_checkMVARParams.m
    validationResults(:,6,setIdx)  = EEG.CAT.VALIDATION.whitestats.ljungbox.pval';
    validationResults(:,7,setIdx)  = EEG.CAT.VALIDATION.whitestats.acf.pval';
    validationResults(:,8,setIdx)  = EEG.CAT.VALIDATION.whitestats.boxpierce.pval';
    validationResults(:,9,setIdx)  = EEG.CAT.VALIDATION.whitestats.limcleod.pval';
    validationResults(:,10,setIdx) = EEG.CAT.VALIDATION.PCstats.PC';
    validationResults(:,11,setIdx) = max(real(EEG.CAT.VALIDATION.stabilitystats.lambda),[],2)';
    % EEG.CAT.MODEL.AR{1,x} contains IC x (IC x model order) that represent
    % AR coefficient of each time lag.
    
    svdResults{setIdx,1} = svd(EEG.icawinv);
end

% Take average across windows.
meanResults = squeeze(mean(validationResults, 1));


% stdResults  = squeeze(std(validationResults, 0, 1));
switch size(validationResults,1) % If single-window analysis.
    case 1
        steResults = zeros(size(meanResults));
    case ~1
        steResults = squeeze(std(validationResults, 0, 1)/sqrt(size(validationResults,1)));
end

% If single-subject, transpose data.
if size(meanResults,1) == 1
    meanResults = meanResults';
    steResults  = steResults';
end

conditionMean = mean(meanResults,2);
conditionStd  = std(meanResults, 0, 2);
% conditionStd  = std(meanResults, 0, 2)/sqrt(size(meanResults,2));

% Generate bar graphs
barGraphColors = colormap(cool(11));
figHandle = figure;
set(figHandle, 'color', [0.93 0.96 1])
for setIdx = 1:11
    figure(figHandle);
    subplot(4,3,setIdx);
    
    if setIdx == 3
        concatenateValue = [];
        concatenateLabel = [];
        medianSingularValueList = zeros(length(svdResults),1);
        for m = 1:length(svdResults)
            tmpValue = svdResults{m};
            tmpLabel = repmat(m, [length(tmpValue) 1]);
            concatenateValue = [concatenateValue; tmpValue];
            concatenateLabel = [concatenateLabel; tmpLabel];
            medianSingularValueList(m,:) = median(tmpValue);
        end
        svGrandMean = mean(concatenateValue);
        svStd         = std(concatenateValue);
        concatenateValue(end+1) = NaN;
        concatenateLabel(end+1) = m+1;
        %         concatenateValue(end+1:end+length(singularValueList)) = medianSingularValueList;
        %         concatenateLabel(end+1:end+length(singularValueList)) = m+2;
        concatenateValue(end+1:end+length(concatenateValue)-1) = concatenateValue(1:end-1);
        concatenateLabel(end+1:length(concatenateValue)) = m+2;
        
        currentPosition = get(gca, 'OuterPosition');
        boxplot(concatenateValue, num2str(concatenateLabel), 'boxstyle', 'filled', 'colors', barGraphColors(setIdx,:));
        set(gca, 'OuterPosition', currentPosition);
        if length(unique(concatenateLabel)) > 10
            cellLabels{end} = 'Median';
            set(gca, 'XTick', 1:length(cellLabels), 'XTickLabel', cellLabels)
        end
    else
        barHandle = bar([meanResults(setIdx,:) 0 conditionMean(setIdx)]');
        barLabels = get(barHandle, 'XData');
        cellLabels = num2cell(barLabels);
        cellLabels{end-1} = '';
        cellLabels{end} = 'Mean';
  
        if length(cellLabels) > 10
            reducedIdx = 0:5:length(cellLabels)-2;
            reducedIdx = reducedIdx(2:end);
            deleteIdx  = setdiff(1:length(cellLabels)-2, reducedIdx);
            [cellLabels{deleteIdx}] = deal('');
        end
        
        set(gca, 'XTick', 1:size(meanResults,2)+2, 'XtickLabel', cellLabels')
        set(barHandle, 'facecolor', barGraphColors(setIdx,:))
        hold on
        errorbarHandle = errorbar([meanResults(setIdx,:) 0 conditionMean(setIdx)], [steResults(setIdx,:) 0 conditionStd(setIdx)]);
        set(errorbarHandle, 'linestyle', 'none', 'Color', [0 0 0])
        
        if setIdx == 10
            set(get(gca, 'xlabel'), 'string', 'Subjects')
        end
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Arrange axis labels %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    xlim([0 size(meanResults,2)+3])
    if     setIdx==1
        title(sprintf('Num trials. M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        ylim([0 max(meanResults(setIdx,:))*1.2])
    elseif setIdx==2
        title(sprintf('Num ICs. M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        ylim([0 max(meanResults(setIdx,:))*1.2])
    elseif setIdx==3
        title(sprintf('Eigval of icawinv (high?). M %.2f (SD %.2f)', svGrandMean, svStd));
    elseif setIdx==4
        title(sprintf('Model order. M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        ylim([0 max(meanResults(setIdx,:))*1.2])
    elseif setIdx==5
        title(sprintf('Param/Datap (<0.1?) M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        ylim([0 max(meanResults(setIdx,:))*1.2])
        line([0 length(barLabels)+1], [0.1  0.1],  'color', [1 0 0])
    elseif setIdx==6
        title(sprintf('Ljung-Box pval (>0.05?) M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        line([0 length(barLabels)+1], [0.05 0.05], 'color', [1 0 0])
    elseif setIdx==7
        title(sprintf('ACF pval (>0.95?) M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        ylim([0 1])
        line([0 length(barLabels)+1], [0.95 0.95], 'color', [1 0 0])
    elseif setIdx==8
        title(sprintf('Box-Pierce pval (>0.05?) M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        line([0 length(barLabels)+1], [0.05 0.05], 'color', [1 0 0])
    elseif setIdx==9
        title(sprintf('Li-McLeod pval (>0.05?) M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        line([0 length(barLabels)+1], [0.05 0.05], 'color', [1 0 0])
    elseif setIdx==10
        title(sprintf('Percent consist (high?) M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        ylim([0 max(meanResults(setIdx,:))*1.2]);
    elseif setIdx==11
        title(sprintf('Stability idx (<0?) M %.2f (SD %.2f)', conditionMean(setIdx), conditionStd(setIdx)));
        ylim([min(meanResults(setIdx,:))*1.2 0]);
        line([0 length(barLabels)+1], [0 0], 'color', [1 0 0])
    end
    
    line([length(barLabels)-1 length(barLabels)-1], ylim, 'color', [0 0 0], 'linestyle', '--')
    
end
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)

% Display process end
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''2.Validate AR models'' finished. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))
