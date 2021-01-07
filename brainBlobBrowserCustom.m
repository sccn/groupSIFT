function varargout = brainBlobBrowserCustom(varargin)
% BRAINBLOBBROWSERCUSTOM MATLAB code for brainBlobBrowserCustom.fig
%
%      Example:
%      brainBlobBrowserCustom('data', rand(91,109,91), 'color', 'red')
%
%      Input:
%      'data'  - input data. Size must be [91 109 91]
%      'color' - color scheme in plot. ['red'|'blue'|'green'|'magenta'|'blue-green'|'purple-green']
%      'mri'   - mri data
%
%      BRAINBLOBBROWSERCUSTOM, by itself, creates a new BRAINBLOBBROWSERCUSTOM or raises the existing
%      singleton*.
%
%      H = BRAINBLOBBROWSERCUSTOM returns the handle to a new BRAINBLOBBROWSERCUSTOM or the handle to
%      the existing singleton*.
%
%      BRAINBLOBBROWSERCUSTOM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAINBLOBBROWSERCUSTOM.M with the given input arguments.
%
%      BRAINBLOBBROWSERCUSTOM('Property','Value',...) creates a new BRAINBLOBBROWSERCUSTOM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before brainBlobBrowserCustom_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to brainBlobBrowserCustom_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help brainBlobBrowserCustom

% Last Modified by GUIDE v2.5 15-Mar-2017 10:42:59

% Begin initialization code - DO NOT EDIT

% History
% 03/06/2017 Makoto. Modified.


gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @brainBlobBrowserCustom_OpeningFcn, ...
                   'gui_OutputFcn',  @brainBlobBrowserCustom_OutputFcn, ...
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


% --- Executes just before brainBlobBrowserCustom is made visible.
function brainBlobBrowserCustom_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to brainBlobBrowserCustom (see VARARGIN)

% this uses Christian's function
userInput = hlp_varargin2struct(varargin);

% Show report on the brainBlobBrowserCustom GUI
reportText = userInput.roiDipoleDensityReport;
set(handles.roiLabelText,      'String', reportText(:,1))
set(handles.dipoleDensityText, 'String', reportText(:,2))

% add spm8 to the path
% addpath /data/projects/makoto/Tools/spm8

% check input
if ~isfield(userInput, 'data')
    error('No input data provided.')
%     disp(sprintf('\n\nNo input data.'));
%     disp('SupraMarginal_L will be shown for demo. ')
%     disp('The input data must be 91x109x91 of 2x2x2 voxels, ranges -88:92, -128:90, -74:108')
%     % load AAL-segmented brain from nii (ComeOnJohnAshburner)
%     tmpAAL = spm_vol('/data/projects/makoto/Tools/spm8/toolbox/wfu_pickatlas/MNI_atlas_templates/aal_MNI_V4.nii');
%     Vols = zeros(numel(tmpAAL),1);
%     for j=1:numel(tmpAAL),
%         tot = 0;
%         for i=1:tmpAAL(1).dim(3),
%                 img     = spm_slice_vol(tmpAAL(j),spm_matrix([0 0 i]),tmpAAL(j).dim(1:2),0);
%                 pileimg(:,:,i) = img; 
%         end;
%     end
%     handles.cubeAAL = pileimg;
%     handles.roiAAL = handles.cubeAAL==63; % SupraMarginal_L
%     handles.inputData = smooth3(handles.roiAAL, 'gaussian', [7 7 7]);
else
    if isequal(size(userInput.data),[91 109 91])
        handles.inputData = userInput.data;
    else
        error('Invalid input data.')
    end
end

if ~isfield(userInput, 'mri')
    dipfitdefs;
    load('-mat', template_models(1).mrifile); % load mri variable
else
    if isequal(size(userInput.mri),[91 109 91])
        mri = userInput.mri;
    else
        error('Invalid input data.')
    end
end

% read MNI template header
% hdr = spm_read_hdr('/data/projects/makoto/Tools/spm8/canonical/single_subj_T1.nii');
handles.dimension = [91 109 91];
% handles.voxelSize = hdr.dime.pixdim(2:4);
% handles.origin    = [0 0 0];% if any(handles.origin); error('Origin is not [0 0 0].'); end
% http://imaging.mrc-cbu.cam.ac.uk/imaging/FormatAnalyze
% Note that if the Origin is set to 0 0 0, then SPM routines will assume that the origin is in fact the central voxel of the image.
% handles.dimensionInMillimeter = [-(handles.dimension.*handles.voxelSize)/2; (handles.dimension.*handles.voxelSize)/2];
handles.currentPointer = [46 64 37];

% Note that origine in voxel space is [46 64 37]
set(handles.axialXSlider,'Value', 46)
set(handles.axialXSlider,'Min', 1)
set(handles.axialXSlider,'Max', handles.dimension(1))

set(handles.axialYSlider,'Value', 64)
set(handles.axialYSlider,'Min', 1)
set(handles.axialYSlider,'Max', handles.dimension(2))

set(handles.sagittalYSlider,'Value', 64)
set(handles.sagittalYSlider,'Min', 1)
set(handles.sagittalYSlider,'Max', handles.dimension(2))

set(handles.coronalZSlider,'Value', 37)
set(handles.coronalZSlider,'Min', 1)
set(handles.coronalZSlider,'Max', handles.dimension(3))


% load SPM single-subject T1 from nii (ComeOnJohnAshburner)
% tmpT1 = spm_vol('/data/projects/makoto/Tools/spm8/canonical/single_subj_T1.nii');
% V.mat   - a 4x4 affine transformation matrix mapping from
%           voxel coordinates to real world coordinates.
handles.affinMat = mri.transform; % don't forget to add the last '1' after the voxel coordinate!
% handles.affinMat = tmpT1.mat; % don't forget to add the last '1' after the voxel coordinate!
handles.affinMat(1,1) = handles.affinMat(1,1)*-1; % non-radiological convention
handles.affinMat(1,4) = handles.affinMat(1,4)*-1; % non-radiological convention
Vols = zeros(numel(mri),1);
for j=1:numel(mri),
    tot = 0;
    for i=1:mri(1).dim(3),
%         img     = spm_slice_vol(tmpT1(j),spm_matrix([0 0 i]),tmpT1(j).dim(1:2),0);
        img = mri.anatomy(:,:,i);
        pileimg(:,:,i) = img;
    end;
end
handles.cubeT1 = pileimg;

% define color map
if ~isfield(userInput,'color')
    handles.cmap = colormap(jet(195-74));
    handles.cmap(1,:) = 0;
    disp(sprintf('\n\nNo color specified.'));
    disp('Color scheme is jet.')
else
    handles.cmap = colormap(hot(200));
    handles.cmap = handles.cmap(75:195,:);
    if     strcmp(userInput.color, 'red')
        disp('Color scheme is red.')
    elseif strcmp(userInput.color, 'blue')
        handles.cmap = handles.cmap(:, [3 2 1]); 
        disp('Color scheme is blue.')
    elseif strcmp(userInput.color, 'green')
        handles.cmap = handles.cmap(:, [2 1 3]);
        disp('Color scheme is blue.')
    elseif strcmp(userInput.color, 'magenta')
        handles.cmap = handles.cmap(:, [1 3 2]);
        disp('Color shceme is magenta')
    elseif strcmp(userInput.color, 'blue-green')
        handles.cmap = handles.cmap(:, [3 1 2]);
        disp('Color shceme is blue-green')
    elseif strcmp(userInput.color, 'purple-blue')
        handles.cmap = handles.cmap(:, [2 3 1]);
        disp('Color shceme is purple-blue')
    else
        error('Unsupported color.')
    end
end

% Choose default command line output for brainBlobBrowserCustom
handles.output = hObject;

% update data
handles = prepare_dens(handles);
handles.alpha = .4;

% update handles
guidata(hObject, handles);

% draw maps
drawAxial(   hObject, eventdata, handles)
drawSagittal(hObject, eventdata, handles)
drawCoronal( hObject, eventdata, handles)

% Set color bar.
axes(handles.colorbarAxes)
colorbarMax = max(userInput.data(:));
colorbarData = [colorbarMax:-colorbarMax/127:0]';
imagesc(colorbarData)
yTicks = [0:length(colorbarData)/8:length(colorbarData)];
yTicks(1) = 1;
yTickLabels = round([colorbarMax:-colorbarMax/(length(yTicks)-1):0]*100)/100;
set(handles.colorbarAxes, 'XTick', [], 'YAxisLocation', 'right', 'YTick', yTicks, 'YTickLabel', yTickLabels)

% set(gcf, 'tool', 'figure')

% Change the figure title on the top bar.
set(gcf, 'name', 'BrainiBlobBrowser--std_dipoleDensity()')

% draw colorbar
% colorbarHandle = colorbar(handles.sagittalAxes, 'EastOutside');
% set(colorbarHandle, 'Ticks', [])
% warning('off','MATLAB:colorbar:DeprecatedV6Argument');

% colormap(handles.cmap); % probably...
% set(handles.axes4, 'YTick', [1 200], 'YTickLabel', {'0' handles.maxdens})
% set(handles.axes4, 'YTick', [1 200], 'YTickLabel', {'1' num2str(userInput.normCoeff)})

% UIWAIT makes brainBlobBrowserCustom wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = brainBlobBrowserCustom_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function [handles maxdens] = prepare_dens(handles)

handles.maxdens = max(handles.inputData(:));
% ncolors = size(handles.cmap,1);
% handles.inputData(~handles.inputData) = nan;
% 
% handles.inputData    = round((handles.inputData)/(handles.maxdens)*(ncolors-1))+1; % project desnity image into the color space: [1:ncolors]
% handles.inputData( find(handles.inputData > ncolors) ) = ncolors;
% handles.inputData( find(handles.inputData < 1))        = 1; % added by Makoto
% % newprob3d = zeros(size(handles.inputData,1), size(handles.inputData,2), size(handles.inputData,3), 3);
% 
% % outOfBrainMask = find(isnan(handles.inputData)); % place NaNs in a mask, NaNs are assumed for points outside the brain
% % handles.inputData(outOfBrainMask) = 1;
handles.inputData = handles.inputData./max(handles.inputData(:));
% handles.inputData(isnan(handles.inputData)) = 0;


% tmp = handles.cmap(handles.inputData,1); newprob3d(:,:,:,1) = reshape(tmp, size(handles.inputData));
% tmp = handles.cmap(handles.inputData,2); newprob3d(:,:,:,2) = reshape(tmp, size(handles.inputData));
% tmp = handles.cmap(handles.inputData,3); newprob3d(:,:,:,3) = reshape(tmp, size(handles.inputData));


% --- Executes on slider movement.
function axialYSlider_Callback(hObject, eventdata, handles)
% hObject    handle to axialYSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currentPointer(2) = round(get(hObject,'Value'));
set(handles.sagittalYSlider, 'Value', get(handles.axialYSlider,'Value'));
% update handles
guidata(hObject, handles);
drawCoronal( hObject, eventdata, handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function axialYSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axialYSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function axialXSlider_Callback(hObject, eventdata, handles)
% hObject    handle to axialXSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currentPointer(1) = round(get(hObject,'Value'));
% update handles
guidata(hObject, handles);
drawSagittal( hObject, eventdata, handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function axialXSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axialXSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function coronalZSlider_Callback(hObject, eventdata, handles)
% hObject    handle to coronalZSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currentPointer(3) = round(get(hObject,'Value'));
% update handles
guidata(hObject, handles);
drawAxial( hObject, eventdata, handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function coronalZSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coronalZSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function drawAxial(hObject, eventdata, handles)

% Select the background brain MRI slice.
currentZ = handles.currentPointer(3);
tmpT1Slice      = rot90(handles.cubeT1(:,:,currentZ)); % rotate to frontal-up
tmpT1SliceColor = repmat(tmpT1Slice, [1 1 3]);
axes(handles.axialAxes);

% Overlay the density blob.
inputDataSlice = squeeze(rot90(handles.inputData(:,:,currentZ)));
inputDataSlice_colorIdx = round(inputDataSlice*(195-75));
inputDataSlice_colorIdx(inputDataSlice_colorIdx==0)=1;
inputDataSliceR = reshape(handles.cmap(inputDataSlice_colorIdx,1), size(inputDataSlice_colorIdx));
inputDataSliceG = reshape(handles.cmap(inputDataSlice_colorIdx,2), size(inputDataSlice_colorIdx));
inputDataSliceB = reshape(handles.cmap(inputDataSlice_colorIdx,3), size(inputDataSlice_colorIdx));
inputDataSliceColor = cat(3, inputDataSliceR, inputDataSliceG, inputDataSliceB);
image(handles.alpha*inputDataSliceColor + (1 - handles.alpha)*tmpT1SliceColor);
set(gca,'XTickLabel', '', 'YTickLabel', '')
currentPointerInRealWorld = handles.affinMat*[handles.currentPointer';1];
set(handles.coordinateZText, 'String', num2str(currentPointerInRealWorld(3)), 'FontSize', 16);

drawCrossHair(hObject, eventdata, handles)



function drawSagittal(hObject, eventdata, handles)

% Select the background brain MRI slice.
currentX = handles.currentPointer(1);
tmpT1Slice      = rot90(squeeze(handles.cubeT1(currentX,:,:))); % rotate to frontal-up
tmpT1SliceColor = repmat(tmpT1Slice, [1 1 3]);
axes(handles.sagittalAxes);

% Overlay the density blob.
inputDataSlice = rot90(squeeze(handles.inputData(currentX,:,:)));
inputDataSlice_colorIdx = round(inputDataSlice*(195-75));
inputDataSlice_colorIdx(inputDataSlice_colorIdx==0)=1;
inputDataSliceR = reshape(handles.cmap(inputDataSlice_colorIdx,1), size(inputDataSlice_colorIdx));
inputDataSliceG = reshape(handles.cmap(inputDataSlice_colorIdx,2), size(inputDataSlice_colorIdx));
inputDataSliceB = reshape(handles.cmap(inputDataSlice_colorIdx,3), size(inputDataSlice_colorIdx));
inputDataSliceColor = cat(3, inputDataSliceR, inputDataSliceG, inputDataSliceB);
inputDataSliceColor = cat(3, inputDataSliceR, inputDataSliceG, inputDataSliceB);
image(handles.alpha*inputDataSliceColor + (1 - handles.alpha)*tmpT1SliceColor);
set(gca,'XTickLabel', '', 'YTickLabel', '')
currentPointerInRealWorld = handles.affinMat*[handles.currentPointer';1];
set(handles.coordinateXText, 'String', num2str(currentPointerInRealWorld(1)*-1), 'FontSize', 16);
axes(handles.sagittalAxes);
% originalPosition = get(gca, 'position');
% colorbarHandle = colorbar;
% set(colorbarHandle, 'Box', 'off', 'YTick', [])
% set(gca, 'position', originalPosition)

drawCrossHair(hObject, eventdata, handles)



function drawCoronal(hObject, eventdata, handles)

% Select the background brain MRI slice.
currentY = handles.currentPointer(2);
tmpT1Slice      = rot90(squeeze(handles.cubeT1(:,currentY,:))); % rotate to frontal-up
tmpT1SliceColor = repmat(tmpT1Slice, [1 1 3]);
axes(handles.coronalAxes);

% Overlay the density blob.
inputDataSlice = rot90(squeeze(handles.inputData(:,currentY,:)));
inputDataSlice_colorIdx = round(inputDataSlice*(195-75));
inputDataSlice_colorIdx(inputDataSlice_colorIdx==0)=1;
inputDataSliceR = reshape(handles.cmap(inputDataSlice_colorIdx,1), size(inputDataSlice_colorIdx));
inputDataSliceG = reshape(handles.cmap(inputDataSlice_colorIdx,2), size(inputDataSlice_colorIdx));
inputDataSliceB = reshape(handles.cmap(inputDataSlice_colorIdx,3), size(inputDataSlice_colorIdx));
inputDataSliceColor = cat(3, inputDataSliceR, inputDataSliceG, inputDataSliceB);
inputDataSliceColor = cat(3, inputDataSliceR, inputDataSliceG, inputDataSliceB);
image(handles.alpha*inputDataSliceColor + (1 - handles.alpha)*tmpT1SliceColor);
set(gca,'XTickLabel', '', 'YTickLabel', '')
currentPointerInRealWorld = handles.affinMat*[handles.currentPointer';1];
set(handles.coordinateYText, 'String', num2str(currentPointerInRealWorld(2)), 'FontSize', 16);

drawCrossHair(hObject, eventdata, handles)



% --- Executes on slider movement.
function sagittalYSlider_Callback(hObject, eventdata, handles)
% hObject    handle to sagittalYSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currentPointer(2) = round(get(hObject,'Value'));
set(handles.axialYSlider, 'Value', get(handles.sagittalYSlider,'Value'));
% update handles
guidata(hObject, handles);
drawCoronal( hObject, eventdata, handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sagittalYSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sagittalYSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function drawCrossHair(hObject, eventdata, handles)

lineHandles = findobj(gcf, 'Type', 'line');
delete(lineHandles);
    
axes(handles.coronalAxes)
hold on
line([handles.currentPointer(1) handles.currentPointer(1)], [0 handles.dimension(3)], 'color', [1 1 1]);
line([0 handles.dimension(1)], [handles.dimension(3)-handles.currentPointer(3) handles.dimension(3)-handles.currentPointer(3)], 'color', [1 1 1]);
hold off
 
axes(handles.axialAxes)
hold on
line([handles.currentPointer(1) handles.currentPointer(1)], [0 handles.dimension(2)], 'color', [1 1 1]);
line([0 handles.dimension(1)], [handles.dimension(2)-handles.currentPointer(2) handles.dimension(2)-handles.currentPointer(2)], 'color', [1 1 1]);
hold off
 
axes(handles.sagittalAxes)
hold on
line([handles.currentPointer(2) handles.currentPointer(2)], [0 handles.dimension(3)], 'color', [1 1 1]);
line([0 handles.dimension(2)], [handles.dimension(3)-handles.currentPointer(3) handles.dimension(3)-handles.currentPointer(3)], 'color', [1 1 1]);
hold off
