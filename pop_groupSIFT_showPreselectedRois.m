% pop_groupSIFT_showPreselectedRois(varargin)
%
% History
% 12/20/2019 Makoto. Updated. The copyright description is updated.

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

function varargout = pop_groupSIFT_showPreselectedRois(varargin)
% POP_GROUPSIFT_SHOWPRESELECTEDROIS MATLAB code for pop_groupSIFT_showPreselectedRois.fig
%      POP_GROUPSIFT_SHOWPRESELECTEDROIS, by itself, creates a new POP_GROUPSIFT_SHOWPRESELECTEDROIS or raises the existing
%      singleton*.
%
%      H = POP_GROUPSIFT_SHOWPRESELECTEDROIS returns the handle to a new POP_GROUPSIFT_SHOWPRESELECTEDROIS or the handle to
%      the existing singleton*.
%
%      POP_GROUPSIFT_SHOWPRESELECTEDROIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_GROUPSIFT_SHOWPRESELECTEDROIS.M with the given input arguments.
%
%      POP_GROUPSIFT_SHOWPRESELECTEDROIS('Property','Value',...) creates a new POP_GROUPSIFT_SHOWPRESELECTEDROIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_groupSIFT_showPreselectedRois_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_groupSIFT_showPreselectedRois_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_groupSIFT_showPreselectedRois

% Last Modified by GUIDE v2.5 15-Mar-2017 18:20:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_groupSIFT_showPreselectedRois_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_groupSIFT_showPreselectedRois_OutputFcn, ...
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


% --- Executes just before pop_groupSIFT_showPreselectedRois is made visible.
function pop_groupSIFT_showPreselectedRois_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_groupSIFT_showPreselectedRois (see VARARGIN)

% Choose default command line output for pop_groupSIFT_showPreselectedRois
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_groupSIFT_showPreselectedRois wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pop_groupSIFT_showPreselectedRois_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in runBrainBlobBrowserPushbutton.
function runBrainBlobBrowserPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runBrainBlobBrowserPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load *_dipolePairDensity.mat
[loadFile, workingFolder] = uigetfile('*_dipolePairDensity.mat');
if ~any(workingFolder)
    disp('Cancelled.')
    return
end
load([workingFolder loadFile])

% Display process start
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''5.Show preselected ROIs'' started. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make a mask to determine ROIs that receives user-specified percent of unique subjects. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
roiDipoleDensiyPct = 100*sum(dipoleProbabilityInRegion(:,preselectedRoiIdx))/sum(dipoleProbabilityInRegion(:));



%%%%%%%%%%%%%%%%%%%%%%
%%% Report resutls %%%
%%%%%%%%%%%%%%%%%%%%%%
[sortedRoiValuesPct, sortIdxPct] = sort(roiDipoleDensiyPct, 'descend');
sortedRoiLabels                  = roiLabels(preselectedRoiIdx(sortIdxPct));
sumDipoleDensity                 = sum(sortedRoiValuesPct);
reportCell = cell(length(sortedRoiLabels)+1,2);
reportCell(1,1) = {'-ROI labels-'};
reportCell(1,2) = {'-Dipole Density (%)-'};
reportCell(2:end,1) = sortedRoiLabels;
for n = 1:length(sortedRoiValuesPct')
    reportCell{n+1,2} = sprintf('%.2f', sortedRoiValuesPct(n));
end
disp(reportCell)
preselectionLog = sprintf('\n\n%.0f/%.0f anatomical ROIs passed the selection.\nTotal of %.1f%% dipole density is accounted for.',...
                            length(preselectedRoiIdx), length(roiLabels), sumDipoleDensity);
disp(preselectionLog)



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visualize the ROIs. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:length(sortedRoiLabels)
    tmpRoi = pr.regionOfInterestFromAnatomy(pr.headGrid, sortedRoiLabels{n});
    if n == 1
        roiSum = double(tmpRoi.membershipCube)*sortedRoiValuesPct(n);
    else
        roiSum = roiSum + double(tmpRoi.membershipCube)*sortedRoiValuesPct(n);
    end
end

% Adjust roiSum position by 1 slide lower.
roiSumAdjusted = cat(3, roiSum(:,:,2:end), zeros(size(roiSum,1), size(roiSum,2)));
[Xq, Yq, Zq]       = meshgrid(23/91:23/91:23, 27/109:27/109:27, 23/91:23/91:23); % Resize roiSum to 109 x 91 x 91, which is the size brainBlobBrowser takes.
interpolatedRoiSum = permute(interp3(roiSumAdjusted, Xq, Yq, Zq), [2 1 3]);
nanMask            = isnan(interpolatedRoiSum);
interpolatedRoiSum(nanMask) = 0;

% Show dipole density using brainBlobBrowser.
brainBlobBrowserCustom('data', interpolatedRoiSum, 'roiDipoleDensityReport', reportCell)

% Display process end.
disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%% ''5.Show preselected ROIs'' finished. %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))