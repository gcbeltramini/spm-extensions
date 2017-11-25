function varargout = connext_gui(varargin)
% 
% Graphical user interface to plot the signal of a chosen voxel after some
% processing options. It can also perform functional connectivity analysis.
% 
%__________________________________________________________________________
% Copyright (C) 2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2014-Jul-15


% Begin initialization code - DO NOT EDIT ---------------------------------
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @connext_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @connext_gui_OutputFcn, ...
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
% End initialization code - DO NOT EDIT -----------------------------------


function connext_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% Executes just before connext_gui is made visible

% Choose default command line output for connext_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

function varargout = connext_gui_OutputFcn(hObject, eventdata, handles) 
% Outputs from this function are returned to the command line
varargout{1} = handles.output;


%==========================================================================
%                             BASIC INPUT
%==========================================================================

% COREGISTERED IMAGE
%==========================================================================
function input_coreg_img_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    status = 'on';
    if ~isempty(get(handles.images_coreg_button, 'UserData'))
        set(handles.run_button, 'Enable', 'on')
    else
        set(handles.run_button, 'Enable', 'off')
    end
else
    status = 'off';
    if ~isempty(get(handles.images_fmri_button, 'UserData'))
        set(handles.run_button, 'Enable', 'on')
    end
end
set(handles.images_coreg_button, 'Enable', status)
set(handles.images_coreg_txt   , 'Enable', status)


% TR
%==========================================================================
function input_TR_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'TR', 'real', 0)


% REGRESSORS
%==========================================================================

function input_regressors_checkbox_Callback(hObject, eventdata, handles)
aux_enablehndl(hObject, handles.input_regressors_button)


function input_regressors_button_Callback(hObject, eventdata, handles)
tmp = get(hObject, 'UserData');
if isempty(tmp)
    tmp = '';
end
[tmp, sts] = spm_select([1 Inf], 'any', ...
    'Regressor file(s): time along the rows, regressors along the columns', cellstr(tmp), '', '.*txt');
set(hObject, 'UserData', tmp)
if isempty(tmp) % no file was chosen
    set(handles.input_regressors_checkbox, 'Value', 0)
    input_regressors_checkbox_Callback(handles.input_regressors_checkbox, eventdata, handles)
end


% DETREND
%==========================================================================
%function input_detrend_checkbox_Callback(hObject, eventdata, handles)



% FREQUENCY FILTER
%==========================================================================

function input_freqfilt_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.input_freqfilt_panel, 'BorderType', 'beveledout')
    set(handles.input_freqfilt_menu              , 'Enable', 'on')
    set(handles.input_freqfilt_lowcutoff_checkbox, 'Enable', 'on')
    set(handles.input_freqfilt_lowcutoff_edit    , 'Enable', 'on')
    input_freqfilt_menu_Callback(handles.input_freqfilt_menu, eventdata, handles)
else
    set(handles.input_freqfilt_panel, 'BorderType', 'beveledin')
    set(handles.input_freqfilt_menu               , 'Enable', 'off')
    set(handles.input_freqfilt_lowcutoff_checkbox , 'Enable', 'off')
    set(handles.input_freqfilt_lowcutoff_edit     , 'Enable', 'off')
    set(handles.input_freqfilt_highcutoff_checkbox, 'Enable', 'off')
    set(handles.input_freqfilt_highcutoff_edit    , 'Enable', 'off')
    set(handles.input_freqfilt_order_checkbox     , 'Enable', 'off')
    set(handles.input_freqfilt_order_edit         , 'Enable', 'off')
    set(handles.input_freqfilt_dir_checkbox       , 'Enable', 'off')
    set(handles.input_freqfilt_dir_edit           , 'Enable', 'off')
end


function input_freqfilt_menu_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')==1 % high-pass filter using DCT
    set(handles.input_freqfilt_highcutoff_checkbox, 'Enable', 'off')
    set(handles.input_freqfilt_highcutoff_edit    , 'Enable', 'off')
    set(handles.input_freqfilt_order_checkbox     , 'Enable', 'off')
    set(handles.input_freqfilt_order_edit         , 'Enable', 'off')
    set(handles.input_freqfilt_dir_checkbox       , 'Enable', 'off')
    set(handles.input_freqfilt_dir_edit           , 'Enable', 'off')
else
    set(handles.input_freqfilt_highcutoff_checkbox, 'Enable', 'on')
    set(handles.input_freqfilt_highcutoff_edit    , 'Enable', 'on')
    if get(hObject, 'Value')==2 % rectangular band-pass filter
        set(handles.input_freqfilt_order_checkbox, 'Enable', 'off')
        set(handles.input_freqfilt_order_edit    , 'Enable', 'off')
        set(handles.input_freqfilt_dir_checkbox  , 'Enable', 'off')
        set(handles.input_freqfilt_dir_edit      , 'Enable', 'off')
    else
        if any(cellfun(@isempty, {which('filtfilt'), which('firls'), which('butter')}))
            warndlg('You must have the Signal Processing Toolbox to use this option.', ...
                'Freq. filter: Invalid option', 'modal')
            set(hObject, 'Value', 1)
            input_freqfilt_menu_Callback(hObject, eventdata, handles)
            return
        end
        set(handles.input_freqfilt_order_checkbox, 'Enable', 'on')
        if get(handles.input_freqfilt_order_checkbox, 'Value')
            set(handles.input_freqfilt_order_edit, 'Enable', 'on')
        end
        set(handles.input_freqfilt_dir_checkbox, 'Enable', 'on')
        if get(handles.input_freqfilt_dir_checkbox, 'Value')
            set(handles.input_freqfilt_dir_edit, 'Enable', 'on')
        end
    end
end


function input_freqfilt_lowcutoff_checkbox_Callback(hObject, eventdata, handles)
aux_enablehndl(hObject, handles.input_freqfilt_lowcutoff_edit)


function input_freqfilt_highcutoff_checkbox_Callback(hObject, eventdata, handles)
aux_enablehndl(hObject, handles.input_freqfilt_highcutoff_edit)


function input_freqfilt_lowcutoff_edit_Callback(hObject, eventdata, handles)
tmp = str2double(get(hObject, 'String'));
if isnan(tmp) || isinf(tmp) || ~isreal(tmp) || tmp<0
    warndlg('Insert a non-negative real number', ...
        'Low cutoff: Invalid input', 'modal')
    set(hObject, 'String', num2str(get(hObject, 'UserData')))
elseif tmp>=get(handles.input_freqfilt_highcutoff_edit, 'UserData')
    warndlg('The low cutoff must be smaller than the high cutoff', ...
        'Low cutoff: Invalid input', 'modal')
    set(hObject, 'String', num2str(get(hObject, 'UserData')))
else
    set(hObject, 'UserData', tmp)
end
if get(hObject, 'UserData')==0
    set(handles.input_freqfilt_lowcutoff_checkbox, 'Value', 0)
    set(handles.input_freqfilt_lowcutoff_edit    , 'Enable', 'off')
end


function input_freqfilt_highcutoff_edit_Callback(hObject, eventdata, handles)
tmp = str2double(get(hObject, 'String'));
if isnan(tmp) || ~isreal(tmp) || tmp<=0
    warndlg('Insert a positive real number', ...
        'High cutoff: Invalid input', 'modal')
    set(hObject, 'String', num2str(get(hObject, 'UserData')))
elseif tmp<=get(handles.input_freqfilt_highcutoff_edit, 'UserData')
    warndlg('The high cutoff must be greater than the low cutoff', ...
        'High cutoff: Invalid input', 'modal')
    set(hObject, 'String', num2str(get(hObject, 'UserData')))
else
    set(hObject, 'UserData', tmp)
end
if isinf(tmp)
    set(handles.input_freqfilt_highcutoff_checkbox, 'Value', 0)
    set(handles.input_freqfilt_highcutoff_edit    , 'Enable', 'off')
end


function input_freqfilt_order_checkbox_Callback(hObject, eventdata, handles)
aux_enablehndl(hObject, handles.input_freqfilt_order_edit)


function input_freqfilt_order_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'Filter order', 'integer', 1)


function input_freqfilt_dir_checkbox_Callback(hObject, eventdata, handles)
aux_enablehndl(hObject, handles.input_freqfilt_dir_edit)


% MASK
%==========================================================================

function input_mask_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.input_mask_panel, 'BorderType', 'beveledout')
    status = 'on';
else
    set(handles.input_mask_panel, 'BorderType', 'beveledin')
    status = 'off';
end
set(handles.input_mask_button     , 'Enable', status)
set(handles.input_mask_thresh_txt , 'Enable', status)
set(handles.input_mask_thresh_edit, 'Enable', status)


function input_mask_button_Callback(hObject, eventdata, handles)
tmp = get(hObject, 'UserData');
if isempty(tmp)
    tmp = fullfile(fileparts(which('spm')), 'apriori', 'brainmask.nii');
    if exist(tmp, 'file')~=2 % file was not found
        tmp = '';
    end
end
[tmp, sts] = spm_select(1, 'image', 'Select a MASK file', cellstr(tmp));
set(hObject, 'UserData', tmp)
if isempty(tmp) % no file was chosen
    set(handles.input_mask_checkbox, 'Value', 0)
    input_mask_checkbox_Callback(handles.input_mask_checkbox, eventdata, handles)
end


function input_mask_thresh_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'Threshold', 'real', 0)


% MOVING AVERAGE
%==========================================================================
function input_movavg_checkbox_Callback(hObject, eventdata, handles)
aux_enablehndl(hObject, handles.input_movavg_window_txt)
aux_enablehndl(hObject, handles.input_movavg_window_edit)

function input_movavg_window_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'Mov. avg. (window size)', 'integer', 1)


% AVERAGE ACROSS CONDITIONS
%==========================================================================

function input_avgblock_checkbox_Callback(hObject, eventdata, handles)
aux_enablehndl(hObject, handles.input_avgblock_button)


function input_avgblock_button_Callback(hObject, eventdata, handles)
TR = get(handles.input_TR_edit, 'UserData');
if TR>0 % TR is known
    txt = 'seconds';
else
    txt = 'scans';
end

% Get the beginning and duration of the current block
tmp = get(hObject, 'UserData');
if isempty(tmp) % no previous answer
    tmp = {sprintf('0\n0'), '0'};
else % the user has already made choices
    % Beginning and duration are recorded in units of scans for convenience
    % when changing the TR
    if TR>0 % TR is known
        tmp = {num2str((tmp(:,1)-1)*TR), num2str((tmp(1,2)-tmp(1,1))*TR)};
        % doubts about using (tmp(1,2)-tmp(1,1)+1)*TR or (tmp(1,2)-tmp(1,1))*TR
        % scans 1-2 => duration = 1 or 2 scans?
    else
        tmp = {num2str(tmp(:,1)), num2str(tmp(1,2)-tmp(1,1))};
        % doubts about using tmp(1,2)-tmp(1,1) or tmp(1,2)-tmp(1,1)+1
    end
end
answer = inputdlg({sprintf('Beginning (units: %s):\nScan 1 <--> t=0 s', txt), ...
    sprintf('Duration (units: %s):', txt)}, ...
    'Blocks information', [10 30; 1 30], tmp, 'on');

if isempty(answer) % clicked the Cancel button or closed the window
    return
    % Disable the option
    %set(handles.input_avgblock_checkbox, 'Value', 0)
    %set(hObject, 'Enable', 'off')
else
    
    [sts, ons] = aux_onset_dur(str2num(answer{1}), str2double(answer{2}), TR);
    if ~sts
        return
    end
    
    if ~isequal(floor(ons), ceil(ons)) || any(ons(:, 1)==0)
        if TR>0
            msg = ['When dividing the events onset and ending by the TR, ' ...
                'the selected blocks must refer to positive integer indices.'];
        else
            msg = 'The selected blocks must refer to positive integer indices.';
        end
        waitfor(warndlg(msg, 'Warning', 'modal'))
        set(handles.input_avgblock_checkbox, 'Value', 0)
        input_avgblock_checkbox_Callback(handles.input_avgblock_checkbox, eventdata, handles)
        return
    end
        
    set(hObject, 'UserData', ons)
end


% SCALE THE SIGNAL
%==========================================================================

function input_scale_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.input_scale_panel, 'BorderType', 'beveledout')
    status = 'on';
else
    set(handles.input_scale_panel, 'BorderType', 'beveledin')
    status = 'off';
end
set(handles.input_scale_Z_button   , 'Enable', status)
set(handles.input_scale_rng_button , 'Enable', status)
set(handles.input_scale_mean_button, 'Enable', status)
if get(hObject, 'Value') && ~get(handles.input_scale_rng_button, 'Value')
    status = 'off';
end
set(handles.input_scale_rng_min_edit, 'Enable', status)
set(handles.input_scale_rng_max_edit, 'Enable', status)
set(handles.input_scale_rng_min_txt , 'Enable', status)
set(handles.input_scale_rng_max_txt , 'Enable', status)


function input_scale_Z_button_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.input_scale_rng_button  , 'Value', 0)
    set(handles.input_scale_mean_button , 'Value', 0)
    set(handles.input_scale_rng_min_edit, 'Enable', 'off')
    set(handles.input_scale_rng_max_edit, 'Enable', 'off')
    set(handles.input_scale_rng_min_txt , 'Enable', 'off')
    set(handles.input_scale_rng_max_txt , 'Enable', 'off')
else
    set(hObject, 'Value', 1)
end
    

function input_scale_rng_button_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.input_scale_Z_button    , 'Value', 0)
    set(handles.input_scale_rng_min_edit, 'Enable', 'on')
    set(handles.input_scale_rng_max_edit, 'Enable', 'on')
    set(handles.input_scale_rng_min_txt , 'Enable', 'on')
    set(handles.input_scale_rng_max_txt , 'Enable', 'on')
    set(handles.input_scale_mean_button , 'Value', 0)
else
    set(hObject, 'Value', 1)
end


function input_scale_mean_button_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.input_scale_Z_button    , 'Value', 0)
    set(handles.input_scale_rng_button  , 'Value', 0)
    set(handles.input_scale_rng_min_edit, 'Enable', 'off')
    set(handles.input_scale_rng_max_edit, 'Enable', 'off')
    set(handles.input_scale_rng_min_txt , 'Enable', 'off')
    set(handles.input_scale_rng_max_txt , 'Enable', 'off')
else
    set(hObject, 'Value', 1)
end


function input_scale_rng_min_edit_Callback(hObject, eventdata, handles)
tmp = str2double(get(hObject, 'String'));
if isnan(tmp) || isinf(tmp) || ~isreal(tmp)
    warndlg('Insert a real number', ...
        'Scale: Invalid input', 'modal')
    set(hObject, 'String', num2str(get(hObject, 'UserData')))
elseif tmp>=get(handles.input_scale_rng_max_edit, 'UserData')
    warndlg('The minimum value must be smaller than the maximum value', ...
        'Scale: Invalid input', 'modal')
    set(hObject, 'String', num2str(get(hObject, 'UserData')))
else
    set(hObject, 'UserData', tmp)
end


function input_scale_rng_max_edit_Callback(hObject, eventdata, handles)
tmp = str2double(get(hObject, 'String'));
if isnan(tmp) || isinf(tmp) || ~isreal(tmp)
    warndlg('Insert a real number', ...
        'Scale: Invalid input', 'modal')
    set(hObject, 'String', num2str(get(hObject, 'UserData')))
elseif tmp<=get(handles.input_scale_rng_min_edit, 'UserData')
    warndlg('The maximum value must be greater than the minimum value', ...
        'Scale: Invalid input', 'modal')
    set(hObject, 'String', num2str(get(hObject, 'UserData')))
else
    set(hObject, 'UserData', tmp)
end


%==========================================================================
%                            SIGNAL PLOT
%==========================================================================

% The plot checkbox was removed because via GUI an image will always be
% displayed.
% 
% function plot_checkbox_Callback(hObject, eventdata, handles)
% if get(hObject, 'Value')
%     status = 'on';
%     plot_roi_checkbox_Callback(handles.plot_roi_checkbox, eventdata, handles)
%     plot_events_checkbox_Callback(handles.plot_events_checkbox, eventdata, handles)
%     plot_hist_checkbox_Callback(handles.plot_hist_checkbox, eventdata, handles)
% else
%     status = 'off';
%     set(handles.plot_roi_x_txt              , 'Enable', 'off')
%     set(handles.plot_roi_y_txt              , 'Enable', 'off')
%     set(handles.plot_roi_z_txt              , 'Enable', 'off')
%     set(handles.plot_roi_x_edit             , 'Enable', 'off')
%     set(handles.plot_roi_y_edit             , 'Enable', 'off')
%     set(handles.plot_roi_z_edit             , 'Enable', 'off')
%     set(handles.plot_events_ncond_txt       , 'Enable', 'off')
%     set(handles.plot_events_ncond_edit      , 'Enable', 'off')
%     set(handles.plot_events_cond_info_button, 'Enable', 'off')
%     set(handles.plot_hist_nbins_button      , 'Enable', 'off')
%     set(handles.plot_hist_nbins_edit        , 'Enable', 'off')
%     set(handles.plot_hist_binsize_button    , 'Enable', 'off')
%     set(handles.plot_hist_binsize_edit      , 'Enable', 'off')
% end
% set(handles.plot_roi_checkbox   , 'Enable', status)
% set(handles.plot_events_checkbox, 'Enable', status)
% set(handles.plot_hist_checkbox  , 'Enable', status)


% RECTANGULAR ROI
%==========================================================================

function plot_roi_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.plot_roi_panel, 'BorderType', 'beveledout')
    status = 'on';
else
    set(handles.plot_roi_panel, 'BorderType', 'beveledin')
    status = 'off';
end
set(handles.plot_roi_x_txt , 'Enable', status)
set(handles.plot_roi_y_txt , 'Enable', status)
set(handles.plot_roi_z_txt , 'Enable', status)
set(handles.plot_roi_x_edit, 'Enable', status)
set(handles.plot_roi_y_edit, 'Enable', status)
set(handles.plot_roi_z_edit, 'Enable', status)

function plot_roi_x_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'ROI', 'integer', 0)

function plot_roi_y_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'ROI', 'integer', 0)

function plot_roi_z_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'ROI', 'integer', 0)


% PLOT EVENTS
%==========================================================================

function plot_events_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.plot_events_panel, 'BorderType', 'beveledout')
    status = 'on';
else
    set(handles.plot_events_panel, 'BorderType', 'beveledin')
    status = 'off';
end
set(handles.plot_events_ncond_txt       , 'Enable', status)
set(handles.plot_events_ncond_edit      , 'Enable', status)
set(handles.plot_events_cond_info_button, 'Enable', status)


function plot_events_ncond_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'Plot events', 'integer', 1)


function plot_events_cond_info_button_Callback(hObject, eventdata, handles)
ncond = get(handles.plot_events_ncond_edit, 'UserData'); % number of conditions
tmp   = cell(ncond, 1);
cond  = struct('name', tmp, 'time', tmp, 'color', tmp);

prev_ans = get(hObject, 'UserData');
if length(prev_ans)<ncond
    prev_ans = [];
elseif length(prev_ans)>ncond
    prev_ans = prev_ans(1:ncond);
end

% Onset and duration are recorded in units of scans for convenience when
% changing the TR
TR = get(handles.input_TR_edit, 'UserData');
if TR>0 % TR is known
    for cc=1:length(prev_ans)
        prev_ans(cc).time = (prev_ans(cc).time-1)*TR;
    end
else
   % do nothing
end

cc = 0; % counter for the conditions
while cc<ncond % loop for the conditions
    
    cc = cc + 1;
    
    if TR>0 % TR is known
        txt = 'seconds';
    else
        txt = 'scans';
    end
    
    % Get condition name, onset and duration
    if isempty(prev_ans) % no previous answer
        tmp = {sprintf('Condition %d', cc), sprintf('0\n0'), '0'};
        col = [.8 .8 .8];
    else % the user has already made choices
        tmp = prev_ans(cc);
        col = tmp.color;
        tmp = {tmp.name, num2str(tmp.time(:,1)), num2str(tmp.time(:,2)-tmp.time(:,1))};
    end
    answer = inputdlg({'Cond. name:', ...
        sprintf('Onset (units: %s):\nScan 1 <--> t=0 s', txt), ...
        sprintf('Duration (units: %s):', txt)}, ...
        sprintf('Cond. %d', cc), [1 10 10], tmp, 'on');
    
    if isempty(answer) % clicked the Cancel button or closed the window
        cond(cc:ncond) = []; % ignore this and further conditions
        ncond = cc - 1;
        if ncond==0
            return
        else
            break
        end
    end
    
    % Condition name
    cond(cc).name = answer{1};
    
    % Onset and end of the condition
    [sts, tmp] = aux_onset_dur(str2num(answer{2}), str2num(answer{3}), TR);
    if ~sts
        cc = cc - 1;
        continue
    end
    cond(cc).time = tmp;
    
    % Color for the rectangles
    cond(cc).color = uisetcolor(col, sprintf('Rectangle color %d', cc));
    
end

set(handles.plot_events_ncond_edit, 'String'  , num2str(ncond))
set(handles.plot_events_ncond_edit, 'UserData', ncond)
set(hObject, 'UserData', cond)


% HISTOGRAM
%==========================================================================

function plot_hist_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.plot_hist_panel, 'BorderType', 'beveledout')
    status = 'on';
else
    set(handles.plot_hist_panel, 'BorderType', 'beveledin')
    status = 'off';
end
set(handles.plot_hist_nbins_button  , 'Enable', status)
set(handles.plot_hist_binsize_button, 'Enable', status)
set(handles.plot_hist_nbins_edit    , 'Enable', status)
set(handles.plot_hist_binsize_edit  , 'Enable', status)
if get(hObject, 'Value')
    if get(handles.plot_hist_nbins_button, 'Value')==1
        set(handles.plot_hist_binsize_edit, 'Enable', 'off')
    else
        set(handles.plot_hist_nbins_edit, 'Enable', 'off')
    end
end


function plot_hist_nbins_button_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.plot_hist_nbins_edit    , 'Enable', 'on')
    set(handles.plot_hist_binsize_button, 'Value', 0)
    set(handles.plot_hist_binsize_edit  , 'Enable', 'off')
else
    set(hObject, 'Value', 1)
end


function plot_hist_binsize_button_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.plot_hist_nbins_button, 'Value', 0)
    set(handles.plot_hist_nbins_edit  , 'Enable', 'off')
    set(handles.plot_hist_binsize_edit, 'Enable', 'on')
else
    set(hObject, 'Value', 1)
end


function plot_hist_nbins_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'Histogram', 'integer', 1)


function plot_hist_binsize_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'Histogram', 'real', 1)


%==========================================================================
%                       FUNCTIONAL CONNECTIVITY
%==========================================================================

function fc_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    status = 'on';
    
    fc_seed_pos_menu_Callback(handles.fc_seed_pos_menu, eventdata, handles)
    fc_seed_roi_checkbox_Callback(handles.fc_seed_roi_checkbox, eventdata, handles)
    
    % Histogram
    set(handles.plot_hist_checkbox      , 'Enable', 'off')
    set(handles.plot_hist_nbins_button  , 'Enable', 'off')
    set(handles.plot_hist_binsize_button, 'Enable', 'off')
    set(handles.plot_hist_nbins_edit    , 'Enable', 'off')
    set(handles.plot_hist_binsize_edit  , 'Enable', 'off')
    
else
    status = 'off';
    set(handles.fc_seed_pos_edit      , 'Enable', 'off')
    set(handles.fc_seed_pos_units_txt , 'Enable', 'off')
    set(handles.fc_seed_pos_units_menu, 'Enable', 'off')
    set(handles.fc_seed_roi_x_txt     , 'Enable', 'off')
    set(handles.fc_seed_roi_x_edit    , 'Enable', 'off')
    set(handles.fc_seed_roi_y_txt     , 'Enable', 'off')
    set(handles.fc_seed_roi_y_edit    , 'Enable', 'off')
    set(handles.fc_seed_roi_z_txt     , 'Enable', 'off')
    set(handles.fc_seed_roi_z_edit    , 'Enable', 'off')
    
    % Histogram
    set(handles.plot_hist_checkbox, 'Enable', 'on')
    plot_hist_checkbox_Callback(handles.plot_hist_checkbox, eventdata, handles)
    
end
set(handles.fc_text             , 'Enable', status)
set(handles.fc_seed_pos_menu    , 'Enable', status)
set(handles.fc_seed_roi_checkbox, 'Enable', status)


% INITIAL SEED POSITION
%==========================================================================

function fc_seed_pos_menu_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')==1 % only ctrl+click
    sts_pos_edit = 'off';
    sts_seed_roi = 'on';
elseif get(hObject, 'Value')==2 % initial seed position
    sts_pos_edit = 'on';
    sts_seed_roi = 'on';
elseif get(hObject, 'Value')==3 % time series
    sts_pos_edit = 'off';
    prev = get(hObject, 'UserData');
    if isempty(prev)
        prev = '';
    end
    [tmp, sts] = spm_select(1, 'any', 'Select the time series file', ...
        cellstr(prev), '', '.*txt');
    if ~sts && isempty(prev) % user closed the window and there's no previous file
        set(hObject, 'Value', 1)
        fc_seed_pos_menu_Callback(handles.fc_seed_pos_menu, eventdata, handles)
        return
    end
    if ~isempty(tmp)
        set(hObject, 'UserData', tmp)
    end
    sts_seed_roi = 'off';
end
set(handles.fc_seed_pos_edit      , 'Enable', sts_pos_edit)
set(handles.fc_seed_pos_units_txt , 'Enable', sts_pos_edit)
set(handles.fc_seed_pos_units_menu, 'Enable', sts_pos_edit)
set(handles.fc_seed_roi_checkbox  , 'Enable', sts_seed_roi)
if strcmp(sts_seed_roi, 'on')
    fc_seed_roi_checkbox_Callback(handles.fc_seed_roi_checkbox, eventdata, handles)
else
    set(handles.fc_seed_roi_x_txt , 'Enable', 'off')
    set(handles.fc_seed_roi_y_txt , 'Enable', 'off')
    set(handles.fc_seed_roi_z_txt , 'Enable', 'off')
    set(handles.fc_seed_roi_x_edit, 'Enable', 'off')
    set(handles.fc_seed_roi_y_edit, 'Enable', 'off')
    set(handles.fc_seed_roi_z_edit, 'Enable', 'off')
end


function fc_seed_pos_edit_Callback(hObject, eventdata, handles)
tmp = str2num(get(hObject, 'String'));
ok = 1;
if length(tmp)~=3 || any(isnan(tmp)) || any(isinf(tmp)) || any(~isreal(tmp))
    txt = 'Insert 3 real numbers';
    ok  = 0;
elseif get(handles.fc_seed_pos_units_menu, 'Value')==2 && ...
        (any(ceil(tmp)~=floor(tmp)) || any(tmp<=0))
    txt = 'Insert 3 positive integers';
    ok  = 0;
    set(hObject, 'String', num2str([1 1 1]))
    set(hObject, 'UserData', [1 1 1])
else
    set(hObject, 'UserData', tmp)
end
if ~ok
    warndlg(txt, ...
        'Seed pos.: Invalid input', 'modal')
    tmp = num2str(get(hObject, 'UserData')); % get the previous value
    tmp = regexprep(tmp, '\s+', ' '); % remove repeated spaces
    set(hObject, 'String', tmp)
end


function fc_seed_pos_units_menu_Callback(hObject, eventdata, handles)
switch get(hObject, 'Value')
    case 1 % mm
        % do nothing
    case 2 % voxel
        fc_seed_pos_edit_Callback(handles.fc_seed_pos_edit, eventdata, handles)
end
fc_seed_pos_edit_Callback(handles.fc_seed_pos_edit, eventdata, handles)


% SEED AS A RECTANGULAR ROI
%==========================================================================

function fc_seed_roi_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.fc_seed_roi_panel, 'BorderType', 'beveledout')
    status = 'on';
else
    set(handles.fc_seed_roi_panel, 'BorderType', 'beveledin')
    status = 'off';
end
set(handles.fc_seed_roi_x_txt , 'Enable', status)
set(handles.fc_seed_roi_y_txt , 'Enable', status)
set(handles.fc_seed_roi_z_txt , 'Enable', status)
set(handles.fc_seed_roi_x_edit, 'Enable', status)
set(handles.fc_seed_roi_y_edit, 'Enable', status)
set(handles.fc_seed_roi_z_edit, 'Enable', status)


function fc_seed_roi_x_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'ROI (seed)', 'integer', 0)


function fc_seed_roi_y_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'ROI (seed)', 'integer', 0)


function fc_seed_roi_z_edit_Callback(hObject, eventdata, handles)
aux_chcknmbr(hObject, 'ROI (seed)', 'integer', 0)


%==========================================================================
%                                IMAGES
%==========================================================================

function images_fmri_button_Callback(hObject, eventdata, handles)
tmp = get(hObject, 'UserData');
if isempty(tmp)
    tmp = '';
end
[tmp, sts] = spm_select(1, 'image', ...
    'Select the FUNCTIONAL image', cellstr(tmp));
if sts % user did not close the window
    set(handles.images_fmri_button, 'UserData', tmp)
    set(handles.images_fmri_txt   , 'String'  , '1 file selected')
    if ~get(handles.input_coreg_img_checkbox, 'Value') || ...
            ~isempty(get(handles.images_coreg_button, 'UserData'))
        set(handles.run_button, 'Enable' ,'on')
    end
end


function images_coreg_button_Callback(hObject, eventdata, handles)
tmp = get(hObject, 'UserData');
if isempty(tmp)
    tmp = '';
end
[tmp, sts] = spm_select(1, 'image', ...
    'Select the COREGISTERED image (e.g., T1W)', cellstr(tmp));
if sts % user did not close the window
    set(handles.images_coreg_button, 'UserData', tmp)
    set(handles.images_coreg_txt   , 'String'  , '1 file selected')
    if ~isempty(get(handles.images_fmri_button, 'UserData'))
        set(handles.run_button, 'Enable' ,'on')
    end
end


%==========================================================================
%                              SAVE & RUN
%==========================================================================

function save_button_Callback(hObject, eventdata, handles)
opt = read_opt(handles);
waitfor(warndlg({'The variable "opt" will be saved.',...
    '',...
    ['By setting "opt.plot=0" later, you can disable the image display ' ...
    'done by the GUI. This can be useful when analysing many images via script.']}, ...
    'Attention', 'modal'))
uisave('opt', 'connext_params')

function run_button_Callback(hObject, eventdata, handles)

opt = read_opt(handles);

set(handles.main_window, 'Pointer', 'watch')
set(handles.run_txt    , 'Visible', 'on')
connext(opt) % Run program
set(handles.run_txt    , 'Visible', 'off')
set(handles.main_window, 'Pointer', 'arrow')


%==========================================================================
% Auxiliary functions
%==========================================================================

function aux_chcknmbr(hObject, msg, type, zero)
% Check if the number in the "String" field of hObject is a non-negative
% number of type "type". If the input is valid, the "UserData" field will
% be updated; otherwise, "String" will receive the previous value.
% 
% hObject: handle of the object of style "edit"
% msg    : identifier to indicate where the error occurred
% type   : type of number ('integer' or 'real')
% zero   : 0 to include zero, 1 for only positive numbers
% 
tmp = str2double(get(hObject, 'String'));
if zero % zero is not allowed
    zero   = tmp<=0;
    txt{1} = 'positive';
else
    zero   = tmp<0;
    txt{1} = 'non-negative';
end
switch type
    case 'integer'
        type   = ceil(tmp)~=floor(tmp);
        txt{2} = 'integer';
    case 'real'
        type   = 0;
        txt{2} = 'real number';
end
if isnan(tmp) || isinf(tmp) || ~isreal(tmp) || type || zero
    warndlg(['Insert a ' txt{1} ' ' txt{2}], ...
        [msg ': Invalid input'], 'modal')
    set(hObject, 'String', num2str(get(hObject, 'UserData')))
else
    set(hObject, 'UserData', tmp)
end


function aux_enablehndl(handle_get, handle_chng)
% If the checkbox identified by the handle "handle_get" is checked, the
% object identified by "handle_chng" is enabled; otherwise, it is disabled.
if get(handle_get, 'Value')
    status = 'on';
else
    status = 'off';
end
set(handle_chng, 'Enable', status)


function [sts, ons] = aux_onset_dur(ons, dur, TR)
% Returns sts=0 when an error occurred and sts=1 otherwise.
% "ons" contains the onsets in column 1 and the end of the events in col. 2

sts = 1;

% Onset
if min(size(ons))~=1 || any(isnan(ons)) || any(isinf(ons)) || any(ons<0)
    waitfor(warndlg('Invalid onset matrix.', 'Warning', 'modal'))
    sts = 0;
    return
end
ons = ons(:);

% Duration
if min(size(dur))~=1 || any(isnan(dur)) || any(isinf(dur)) || any(dur<=0)
    waitfor(warndlg('Invalid duration matrix.', 'Warning', 'modal'))
    sts = 0;
    return
end
dur = dur(:);
if size(dur, 1)==1 % same duration for all conditions
    dur = repmat(dur, size(ons,1), 1);
elseif size(dur, 1)~=size(ons, 1)
    waitfor(warndlg('The number of duration values must be 1 or the same number of onsets.', ...
        'Warning', 'modal'))
    sts = 0;
    return
end

% doubts about the duration: scans 1-2 => duration = 1 or 2 scans?
% if TR>0
%     dur = dur - TR;
% else
%     dur = dur - 1;
% end

ons = [ons ons+dur];

if TR>0 % the elements were given in time units
    ons = ons/TR + 1; % scan 1 <--> t = 0 s
end


function opt = read_opt(handles)
% Read the options from the GUI

% Input
%------

% Coregistered image
opt.coreg.use = get(handles.input_coreg_img_checkbox, 'Value');
opt.coreg.img = '';

% TR
opt.TR = get(handles.input_TR_edit, 'UserData');

% Regressors
opt.reg.use = get(handles.input_regressors_checkbox, 'Value');
if opt.reg.use && isempty(get(handles.input_regressors_button, 'UserData'))
    input_regressors_button_Callback(handles.input_regressors_button, [], handles)
end
opt.reg.file = get(handles.input_regressors_button, 'UserData');

% Detrend
opt.detrend = get(handles.input_detrend_checkbox, 'Value');

% Frequency filter
opt.filt.use  = get(handles.input_freqfilt_checkbox, 'Value');
opt.filt.type = get(handles.input_freqfilt_menu    , 'Value');
switch opt.filt.type
    case 1
        opt.filt.type = 'dct';
    case 2
        opt.filt.type = 'rectbandpass';
    case 3
        opt.filt.type = 'but';
    case 4
        opt.filt.type = 'fir';
    case 5
        opt.filt.type = 'firls';
    otherwise
        opt.filt.type = [];
end
opt.filt.cutoff = [0 Inf];
if get(handles.input_freqfilt_lowcutoff_checkbox, 'Value')
    opt.filt.cutoff(1) = get(handles.input_freqfilt_lowcutoff_edit, 'UserData');
end
if get(handles.input_freqfilt_highcutoff_checkbox, 'Value')
    opt.filt.cutoff(2) = get(handles.input_freqfilt_highcutoff_edit, 'UserData');
end
if get(handles.input_freqfilt_order_checkbox, 'Value')
    opt.filt.order = get(handles.input_freqfilt_order_edit, 'UserData');
else
    opt.filt.order = [];
end
if get(handles.input_freqfilt_dir_checkbox, 'Value')
    opt.filt.dir = get(handles.input_freqfilt_dir_edit, 'Value');
    switch opt.filt.dir
        case 1
            opt.filt.dir = 'onepass';
        case 2
            opt.filt.dir = 'onepass-reverse';
        case 3
            opt.filt.dir = 'twopass';
        case 4
            opt.filt.dir = 'twopass-reverse';
        case 5
            opt.filt.dir = 'twopass-average';
        otherwise
            opt.filt.dir = [];
    end
else
    opt.filt.dir = [];
end

% Scale
opt.scale.use = get(handles.input_scale_checkbox, 'Value');
if get(handles.input_scale_Z_button, 'Value') % Z-score
    opt.scale.type = 1;
elseif get(handles.input_scale_rng_button, 'Value') % specified range
    opt.scale.type = 2;
    opt.scale.rng = [get(handles.input_scale_rng_min_edit, 'UserData') ; ...
        get(handles.input_scale_rng_max_edit, 'UserData')];
elseif get(handles.input_scale_mean_button, 'Value') % variation around the mean
    opt.scale.type = 3;
end

% Mask
opt.mask.use = get(handles.input_mask_checkbox, 'Value');
if opt.mask.use && isempty(get(handles.input_mask_button, 'UserData'))
    input_mask_button_Callback(handles.input_mask_button, [], handles)
end
opt.mask.file   = get(handles.input_mask_button     , 'UserData');
opt.mask.thresh = get(handles.input_mask_thresh_edit, 'UserData');

% Average across conditions
opt.avgblock.use = get(handles.input_avgblock_checkbox, 'Value');
opt.avgblock.ons = get(handles.input_avgblock_button  , 'UserData'); % units of scans
if opt.avgblock.use && isempty(opt.avgblock.ons)
    input_avgblock_button_Callback(handles.input_avgblock_button, [], handles)
    opt.avgblock.ons = get(handles.input_avgblock_button, 'UserData');
end

% Moving average
opt.movavg.use      = get(handles.input_movavg_checkbox, 'Value');
opt.movavg.windowsz = get(handles.input_movavg_window_edit, 'UserData');
if opt.movavg.windowsz==1
    opt.movavg.use = 0;
end


% Signal plot
%------------
opt.plot = 1; %get(handles.plot_checkbox, 'Value');

% ROI
opt.ROI.use = get(handles.plot_roi_checkbox, 'Value');
opt.ROI.x   = get(handles.plot_roi_x_edit  , 'UserData');
opt.ROI.y   = get(handles.plot_roi_y_edit  , 'UserData');
opt.ROI.z   = get(handles.plot_roi_z_edit  , 'UserData');

% Plot events
opt.plot_events.use  = get(handles.plot_events_checkbox, 'Value');
opt.plot_events.cond = get(handles.plot_events_cond_info_button, 'UserData');
if opt.plot_events.use && length(opt.plot_events.cond)~=str2double(get(handles.plot_events_ncond_edit, 'String'))
    plot_events_cond_info_button_Callback(handles.plot_events_cond_info_button, [], handles)
    opt.plot_events.cond = get(handles.plot_events_cond_info_button, 'UserData');
end
if opt.TR>0 % change to unit of time
    for cc=1:length(opt.plot_events.cond)
        opt.plot_events.cond(cc).time = (opt.plot_events.cond(cc).time-1)*opt.TR;
    end
end

% Histogram
opt.hist.use        = get(handles.plot_hist_checkbox      , 'Value');
opt.hist.nbins.use  = get(handles.plot_hist_nbins_button  , 'Value');
opt.hist.nbins.bins = get(handles.plot_hist_nbins_edit    , 'UserData');
opt.hist.binsz.use  = get(handles.plot_hist_binsize_button, 'Value');
opt.hist.binsz.size = get(handles.plot_hist_binsize_edit  , 'UserData');


% Functional connectivity
%------------------------
opt.funccon.use      = get(handles.fc_checkbox         , 'Value');
opt.funccon.seed.use = get(handles.fc_seed_pos_menu    , 'Value') - 1;
opt.funccon.seed.pos = [get(handles.fc_seed_pos_edit   , 'UserData') get(handles.fc_seed_pos_units_menu, 'Value')-1];
opt.funccon.roi.use  = get(handles.fc_seed_roi_checkbox, 'Value');
opt.funccon.roi.x    = get(handles.fc_seed_roi_x_edit  , 'UserData');
opt.funccon.roi.y    = get(handles.fc_seed_roi_y_edit  , 'UserData');
opt.funccon.roi.z    = get(handles.fc_seed_roi_z_edit  , 'UserData');
if opt.funccon.seed.use==2
    opt.funccon.seed.series = load(get(handles.fc_seed_pos_menu, 'UserData'));
else
    opt.funccon.seed.series = [];
end


% Check if there is anything to do
%=================================
if ~opt.reg.use && ~opt.filt.use && ~opt.scale.use && ~opt.mask.use && ...
        ~opt.avgblock.use && ~opt.movavg.use && ~opt.plot && ~opt.funccon.use
    warndlg('Nothing to do. Choose at least one option.', 'Warning', 'modal')
    opt = [];
    return
end


% Save
%=====
opt.save.use   = 0; % don't save the resulting images
opt.save.fname = '';


% Images
%=======
opt.img       = get(handles.images_fmri_button , 'UserData');
opt.coreg.img = get(handles.images_coreg_button, 'UserData');


% Order the fields
%=================
opt = orderfields(opt, {'img', 'coreg', 'TR', 'reg', 'detrend', 'filt', ...
    'movavg', 'avgblock', 'scale', 'mask', 'plot', 'ROI', 'plot_events', ...
    'hist', 'funccon', 'save'});
