function varargout = fmri_mip_color_gui(varargin)
% 
% GUI for "fmri_mip_color.m"
% 
%__________________________________________________________________________
% Copyright (C) 2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2014-Jul-02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fmri_mip_color_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @fmri_mip_color_gui_OutputFcn, ...
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

function fmri_mip_color_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for fmri_mip_color_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function varargout = fmri_mip_color_gui_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

%==========================================================================

function stats_Padj_radiobutton1_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.stats_Padj_radiobutton2, 'Value', 0)
else
    set(hObject, 'Value', 1)
end
%--------------------------------------------------------------------------
function stats_Padj_radiobutton2_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.stats_Padj_radiobutton1, 'Value', 0)
else
    set(hObject, 'Value', 1)
end
%--------------------------------------------------------------------------
function spmmat_pushbutton_Callback(hObject, eventdata, handles)
prev = get(hObject, 'UserData');
if isempty(prev)
    prev = {''};
end
[t, sts] = spm_select([1 Inf], '^SPM\.mat$', 'Select SPM.mat', prev);
if sts
    set(hObject, 'UserData', cellstr(t))
end
%--------------------------------------------------------------------------
function run_pushbutton_Callback(hObject, eventdata, handles)
if get(handles.stats_Padj_radiobutton1, 'Value')
    adj = 'none';
else
    adj = 'FWE';
end
fmri_mip_color(get(handles.spmmat_pushbutton, 'UserData'), ...
    str2num(get(handles.spmmat_con_edit, 'String')), ...
    str2double(get(handles.stats_P_edit, 'String')), ...
    adj, ...
    str2double(get(handles.stats_ET_edit, 'String')), ...
    str2num(get(handles.spmmat_con_color_edit, 'String')))
