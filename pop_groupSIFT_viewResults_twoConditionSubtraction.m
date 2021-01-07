% pop_groupSIFT_viewResults_twoConditionSubtraction(varargin)
% 
% History
% 12/20/2019 Makoto. Updated. Copyright updated. Minor graphics fixed.

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

function varargout = pop_groupSIFT_viewResults_twoConditionSubtraction(varargin)
% POP_GROUPSIFT_VIEWRESULTS_TWOCONDITIONSUBTRACTION MATLAB code for pop_groupSIFT_viewResults_twoConditionSubtraction.fig
%      POP_GROUPSIFT_VIEWRESULTS_TWOCONDITIONSUBTRACTION, by itself, creates a new POP_GROUPSIFT_VIEWRESULTS_TWOCONDITIONSUBTRACTION or raises the existing
%      singleton*.
%
%      H = POP_GROUPSIFT_VIEWRESULTS_TWOCONDITIONSUBTRACTION returns the handle to a new POP_GROUPSIFT_VIEWRESULTS_TWOCONDITIONSUBTRACTION or the handle to
%      the existing singleton*.
%
%      POP_GROUPSIFT_VIEWRESULTS_TWOCONDITIONSUBTRACTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_GROUPSIFT_VIEWRESULTS_TWOCONDITIONSUBTRACTION.M with the given input arguments.
%
%      POP_GROUPSIFT_VIEWRESULTS_TWOCONDITIONSUBTRACTION('Property','Value',...) creates a new POP_GROUPSIFT_VIEWRESULTS_TWOCONDITIONSUBTRACTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_groupSIFT_viewResults_twoConditionSubtraction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_groupSIFT_viewResults_twoConditionSubtraction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_groupSIFT_viewResults_twoConditionSubtraction

% Last Modified by GUIDE v2.5 01-Jun-2016 11:54:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_groupSIFT_viewResults_twoConditionSubtraction_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_groupSIFT_viewResults_twoConditionSubtraction_OutputFcn, ...
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


% --- Executes just before pop_groupSIFT_viewResults_twoConditionSubtraction is made visible.
function pop_groupSIFT_viewResults_twoConditionSubtraction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_groupSIFT_viewResults_twoConditionSubtraction (see VARARGIN)

% Choose default command line output for pop_groupSIFT_viewResults_twoConditionSubtraction
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_groupSIFT_viewResults_twoConditionSubtraction wait for user response (see UIRESUME)
% uiwait(handles.viewResultsTwoConditionSubtractionFigure);


% --- Outputs from this function are returned to the command line.
function varargout = pop_groupSIFT_viewResults_twoConditionSubtraction_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
