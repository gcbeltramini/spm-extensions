function varargout = safe(varargin)
% 
% SAFE opens the GUI for the fMRI analysis
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2013-Aug-22


% TEMPLATE
% --------
% function XXX_Callback(hObject, eventdata, handles)
% hObject    handle to XXX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% Hints:
% -----
% Togglebutton, Radiobutton, Checkbox: get(hObject,'Value') returns toggle state of XXX
% Edit: get(hObject,'String') returns contents of XXX as text
%       str2double(get(hObject,'String')) returns contents of XXX as a double
% Slider: get(hObject,'Value') returns position of slider
%         get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% Listbox, Popupmenu: contents = cellstr(get(hObject,'String')) returns XXX contents as cell array
%                     contents{get(hObject,'Value')} returns selected item
%                     from XXX


%==========================================================================
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @safe_OpeningFcn, ...
                   'gui_OutputFcn',  @safe_OutputFcn, ...
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

end
% End initialization code - DO NOT EDIT

%--------------------------------------------------------------------------

function safe_OpeningFcn(hObject, eventdata, handles, varargin)
% Executes just before "safe" is made visible
% hObject: handle to figure
    
    % Load the default values
    %------------------------
    handles = safe_read_GUI(handles,'default');
    
    handles.output = handles.opt_safe; % default fields
    tmp = sprintf('Positive numbers are in mm, and negative numbers\nare the multiplicity factor relative to the original size.');
    set(handles.preproc_norm_factor_struc_edit,'TooltipString',tmp)
    set(handles.preproc_norm_factor_func_edit ,'TooltipString',tmp)
    set(handles.preproc_smoothfwhm_edit       ,'TooltipString',tmp)
    % this cannot be set in the FIG
    guidata(hObject, handles);
    
    % Change window size to fit the screen
    %-------------------------------------
    oldunits = get(0,'Units'); set(0,'Units','pixels')
    res = get(0,'ScreenSize'); res = res(3:4);
    oldfigunits = get(hObject,'Units'); set(hObject,'Units','pixels')
    pos = get(hObject,'Position'); pos = pos(1:2);
    set(hObject,'Position',[pos res(1)*790/1280 res(2)*560/1024])
    set(0,'Units',oldunits), set(hObject,'Units',oldfigunits) % reset units
    
    % Display text in the command window
    %-----------------------------------
    
    
    %uiwait(handles.safe_main_window); % wait for user response
end

%--------------------------------------------------------------------------

function varargout = safe_OutputFcn(hObject, eventdata, handles)
    varargout{1} = handles.opt_safe; %[];
    %delete(hObject); % must be in state from uiwait, otherwise figure will
                     % be immediately closed
end
%==========================================================================



%==========================================================================
%                         AUXILIARY FUNCTIONS
%==========================================================================

%--------------------------------------------------------------------------
function handles = safe_read_GUI(handles,option)
% Create the necessary structure for the fMRI analysis
    
    % BASIC INPUT
    %======================================================================
    if strcmp(option,'default')
        % These options cannot be changed via GUI:
        handles.opt_safe.basic.parent_dir = ''; % Parent directory for the subject folder
        handles.opt_safe.basic.fmri_path  = ''; % Folders inside "parent_dir"\"sbj_folder" (same for all subjects)
    end
    
    handles.opt_safe.basic.xlsfile = get(handles.basic_input_xlsfile_txt,'UserData'); % xlsfile{FF}: XLS file FF
    if isempty(handles.opt_safe.basic.xlsfile) % when opening the GUI
        handles.opt_safe.basic.xlsfile = {''};
    end
    
    handles.opt_safe.basic.sbj_folder = get(handles.basic_input_datadir_txt,'UserData');
    % sbj_folder{FF}{WW}: XLS file FF, worksheet tab WW
    if isempty(handles.opt_safe.basic.sbj_folder)
        handles.opt_safe.basic.sbj_folder = {{''}};
    end
    
    if strcmp(option,'default')
        handles.opt_safe.basic.worksheet = {{''}};
        set(handles.basic_input_worksheet_txt,'UserData',handles.opt_safe.basic.worksheet)
        % this cannot be set in the FIG
    else
        handles.opt_safe.basic.worksheet = get(handles.basic_input_worksheet_txt,'UserData');
    end
    % worksheet{FF}{WW}: XLS file FF, worksheet tab WW (same size as "opt_safe.basic.sbj_folder")
    
    handles.opt_safe.basic.eeg_fmri_exp = get(handles.iseegfmri,'Value');
    
    
    % PREPROCESSING
    %======================================================================
    preproc.run           = get(handles.preproc_chckbx                ,'Value');
    preproc.stc           = get(handles.preproc_slice_timing          ,'Value');
    preproc.norm.run      = get(handles.preproc_normalize             ,'Value');
    preproc.norm.struc_vx = get(handles.preproc_norm_factor_struc_edit,'UserData');
    preproc.norm.func_vx  = get(handles.preproc_norm_factor_func_edit ,'UserData');
    preproc.smooth.run    = get(handles.preproc_smooth                ,'Value');
    preproc.smooth.fwhm   = get(handles.preproc_smoothfwhm_edit       ,'UserData');
    handles.opt_safe.preproc = preproc;
    
    
    % DESIGN & ESTIMATE
    %======================================================================
    stat_des.run        = get(handles.design_chckbx         ,'Value');
    stat_des.pref_func  = get(handles.design_funcprefix_edit,'String');
    stat_des.pref_struc = get(handles.res_strucprefix_edit  ,'String');
    stat_des.tag_struc  = get(handles.res_struc_extra_edit  ,'String');
    stat_des.use_rp     = get(handles.design_useRP          ,'Value');
    stat_des.hpf        = get(handles.design_hpf_edit       ,'UserData');
    switch get(handles.design_hrf_menu,'Value')
        case 1
            stat_des.bf = 'can_hrf';
        case 2
            stat_des.bf = 'gamma';
        case 3
            stat_des.bf = 'fourier';
        case 4
            stat_des.bf = 'fourier_han';
        case 5
            stat_des.bf = 'fir';
    end
    stat_des.bf_length = get(handles.design_hrf_length_edit,'UserData');
    stat_des.bf_order  = get(handles.design_hrf_nmbrbf_edit,'UserData');
    stat_des.delay     = get(handles.design_hrf_delay_edit ,'UserData');
    handles.opt_safe.stat_des = stat_des;
    
    
    % CONTRASTS
    %======================================================================
    if get(handles.res_chckbx,'Value')
        handles.opt_safe.stat_con = get(handles.contr_chckbx,'Value');
    else
        handles.opt_safe.stat_con = 0;
    end
    
    
    % RESULTS
    %======================================================================
    
    handles.opt_safe.stat_res.run = get(handles.res_chckbx,'Value');
    
    % P-value
    %--------
    switch get(handles.res_FWE,'Value')
        case 0
            stat_res.thresh.adj = 'none';
        case 1
            stat_res.thresh.adj = 'FWE';
    end
    stat_res.thresh.P  = get(handles.res_pvalue_edit    , 'UserData');
    stat_res.thresh.ET = get(handles.res_vxl_thresh_edit, 'UserData');
    handles.opt_safe.stat_res.thresh = stat_res.thresh;
    
    % Table and PostScript file
    %--------------------------
    if get(handles.res_chckbx,'Value')
        handles.opt_safe.stat_res.ps_save = get(handles.res_PS,'Value');
    else
        handles.opt_safe.stat_res.ps_save = 0;
    end
    
    switch option
        case 'default'
            
            % Save image with slices
            %-----------------------
            if get(handles.res_chckbx,'Value')
                stat_res.slices.save = get(handles.res_slice_view,'Value');
            else
                stat_res.slices.save = 0;
            end
            stat_res.slices.img_ext = 'png';
            stat_res.slices.step    = 0;
            stat_res.slices.number  = 36;
            stat_res.slices.slices  = -50:4:85;
            stat_res.slices.orient  = {'axial'};
            handles.opt_safe.stat_res.slices = stat_res.slices;
            
            % Group of images
            %----------------
            if get(handles.res_chckbx,'Value')
                stat_res.group.save = get(handles.res_union,'Value');
            else
                stat_res.group.save = 0;
            end
            stat_res.group.scale.adj   = 1;
            stat_res.group.scale.ftest = 'max';
            stat_res.group.scale.ttest = 'max';
            stat_res.group.number      = 3;
            handles.opt_safe.stat_res.group = stat_res.group;
            
            % Sequence of images
            %-------------------
            if get(handles.res_chckbx,'Value')
                stat_res.event_seq.save = get(handles.res_seq,'Value');
            else
                stat_res.event_seq.save = 0;
            end
            stat_res.event_seq.scale.ftest   = 'max';
            stat_res.event_seq.scale.posbold = 'max';
            stat_res.event_seq.scale.negbold = 'max';
            handles.opt_safe.stat_res.event_seq = stat_res.event_seq;
            
        case 'read'
            % The other options were set already
            if get(handles.res_chckbx,'Value')
                handles.opt_safe.stat_res.slices.save    = get(handles.res_slice_view,'Value');
                handles.opt_safe.stat_res.group.save     = get(handles.res_union     ,'Value');
                handles.opt_safe.stat_res.event_seq.save = get(handles.res_seq       ,'Value');
            else
                handles.opt_safe.stat_res.slices.save    = 0;
                handles.opt_safe.stat_res.group.save     = 0;
                handles.opt_safe.stat_res.event_seq.save = 0;
            end
    end
    
    % Auxiliary files
    %----------------
    if strcmp(option,'default')
        handles.opt_safe.stat_res.save_thr    = 0;
        handles.opt_safe.stat_res.save_slover = 0;
    end
    
    % Mask
    %-----
    if strcmp(option,'default')
        stat_res.mask.use        = 0;
        stat_res.mask.imgs_text  = {'mwc*'};
        stat_res.mask.incl_excl  = 'inclusive';
        stat_res.mask.bin        = 1;
        stat_res.mask.bin_thresh = 0.2;
        handles.opt_safe.stat_res.mask = stat_res.mask;
    end
    
    % Clean up
    %---------
    handles.opt_safe.stat_res.cleanup = get(handles.res_cleanup,'Value');
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = safe_load(handles,opt_safe)
% Check parameters while loading

ok = 1;
try
    
    % BASIC INPUT
    %======================================================================
    %opt_safe.basic.parent_dir  --  cannot change via GUI => just check
    %opt_safe.basic.fmri_path   --  cannot change via GUI => just check
    tmp = opt_safe.basic.parent_dir;
    if ~ischar(tmp)
        ok   = 0;
        warn = 'opt_safe.basic.parent_dir';
    end
    tmp = opt_safe.basic.fmri_path;
    if ~ischar(tmp)
        ok   = 0;
        warn = 'opt_safe.basic.fmri_path';
    end
    
    tmp = opt_safe.basic.xlsfile;
    nmbr_xls_xlsfile = 0;
    enable_work = 'off';
    if ~iscellstr(tmp)
        ok = 0;
    else
        nmbr_xls_xlsfile = length(tmp);
        if nmbr_xls_xlsfile>1
            tmp_text    = sprintf('%d files selected',nmbr_xls_xlsfile);
            enable_work = 'on';
        elseif nmbr_xls_xlsfile==1
            tmp_text    = '1 file selected';
            enable_work = 'on';
        else %nmbr_xls_xlsfile=0
            tmp_text = '(no file selected)';
            tmp      = {''}; % guarantee that it is a 1x1 cell array
        end
        set(handles.basic_input_xlsfile_txt,'UserData',tmp)
        %basic_input_fit_text(1,handles.basic_input_xlsfile_txt,tmp_text)
        set(handles.basic_input_xlsfile_txt,'String',tmp_text)
    end
    
    tmp = opt_safe.basic.worksheet;
    set(handles.basic_input_worksheet_button,'Enable',enable_work)
    set(handles.basic_input_worksheet_txt   ,'Enable',enable_work)
    enable_dir = enable_work;
    if ~iscell(tmp)
        ok = 0;
    else
        num_xls_wkrsht = length(tmp);
        if num_xls_wkrsht~=nmbr_xls_xlsfile % invalid input
            warndlg('Basic input: Inconsistent number of XLS files listed in worksheet names.',...
                'Invalid input','modal')
            num_xls_wkrsht = nmbr_xls_xlsfile;
            tmp = cell(1,nmbr_xls_xlsfile);
            for ii=1:nmbr_xls_xlsfile
                tmp{ii} = {''};
            end
        end
        if num_xls_wkrsht>=1
            if num_xls_wkrsht==1
                tmp_text = '(1 XLS file)';
            else
                tmp_text = sprintf('(%d XLS files)',num_xls_wkrsht);
            end
            cnt = 0;
            for ii=1:num_xls_wkrsht
                if iscellstr(tmp{ii})
                    cnt = cnt + sum(~cellfun(@isempty,tmp{ii})); % count non-empty cells
                else
                    ok = 0;
                    %enable_dir = 'off';
                    break
                end
            end
            if cnt>1
                tmp_text = sprintf('%d worksheets selected %s',cnt,tmp_text);
                %enable_dir = 'on';
            elseif cnt==1
                tmp_text = sprintf('1 worksheet selected %s',tmp_text);
                %enable_dir = 'on';
            else %cnt=0
                tmp_text = '(no worksheet selected)';
            end
            set(handles.basic_input_worksheet_txt,'UserData',tmp)
            %basic_input_fit_text(2,handles.basic_input_worksheet_txt,tmp_text)
            set(handles.basic_input_worksheet_txt,'String',tmp_text)
        else % {}
            ok  = 0;
        end
    end
    
    tmp = opt_safe.basic.sbj_folder; % similar to opt_safe.basic.worksheet
    set(handles.basic_input_datadir_button,'Enable',enable_dir)
    set(handles.basic_input_datadir_txt   ,'Enable',enable_dir)
    if ~iscell(tmp)
        ok = 0;
    else
        nmbr_xls_dir = length(tmp);
        if nmbr_xls_dir~=nmbr_xls_xlsfile % invalid input
            warndlg('Basic input: Inconsistent number of XLS files listed in data directories.',...
                'Invalid input','modal')
            nmbr_xls_dir = nmbr_xls_xlsfile;
            tmp = cell(1,nmbr_xls_xlsfile);
            for ii=1:nmbr_xls_xlsfile
                tmp{ii} = {''};
            end
        end
        if nmbr_xls_dir>=1
            if nmbr_xls_dir==1
                tmp_text = '(1 XLS file)';
            else
                tmp_text = sprintf('(%d XLS files)',nmbr_xls_dir);
            end
            cnt = 0;
            for ii=1:nmbr_xls_dir
                if iscellstr(tmp{ii})
                    cnt = cnt + sum(~cellfun(@isempty,tmp{ii})); % count non-empty cells
                else
                    ok = 0;
                    break
                end
            end
            if cnt>1
                tmp_text = sprintf('%d directories selected %s',cnt,tmp_text);
            elseif cnt==1
                tmp_text = sprintf('1 directory selected %s',tmp_text);
            else %cnt=0
                tmp_text = '(no directory selected)';
            end
            set(handles.basic_input_datadir_txt,'UserData',tmp)
            %basic_input_fit_text(3,handles.basic_input_datadir_txt,tmp_text)
            set(handles.basic_input_datadir_txt,'String',tmp_text)
        else % {}
            ok  = 0;
        end
    end
    
    if any(opt_safe.basic.eeg_fmri_exp==[0 1])
        set(handles.iseegfmri,'Value',opt_safe.basic.eeg_fmri_exp)
    else
        ok = 0;
    end
    if ~ok && exist('warn','var')~=1
        warn = 'Basic input';
    end
    
    
    % PREPROCESSING
    %======================================================================
    if any(opt_safe.preproc.run==[0 1]) && ...
            any(opt_safe.preproc.stc==[0 1]) && ...
            any(opt_safe.preproc.norm.run==[0 1]) && ...
            any(opt_safe.preproc.smooth.run==[0 1])
        set(handles.preproc_chckbx      ,'Value',opt_safe.preproc.run)
        set(handles.preproc_slice_timing,'Value',opt_safe.preproc.stc)
        set(handles.preproc_normalize   ,'Value',opt_safe.preproc.norm.run)
        set(handles.preproc_smooth      ,'Value',opt_safe.preproc.smooth.run)
    else
        ok = 0;
    end
    
    set(handles.preproc_norm_factor_struc_edit,'String',num2str(opt_safe.preproc.norm.struc_vx))
    set(handles.preproc_norm_factor_func_edit ,'String',num2str(opt_safe.preproc.norm.func_vx ))
    set(handles.preproc_smoothfwhm_edit       ,'String',num2str(opt_safe.preproc.smooth.fwhm  ))
    
    ok_tmp = preproc_norm_factor_struc_edit_Callback(handles.preproc_norm_factor_struc_edit,[],[]);
    if ~ok_tmp
        ok = 0;
    end
    ok_tmp = preproc_norm_factor_func_edit_Callback(handles.preproc_norm_factor_func_edit ,[],[]);
    if ~ok_tmp
        ok = 0;
    end
    ok_tmp = preproc_smoothfwhm_edit_Callback(handles.preproc_smoothfwhm_edit,[],[]);
    if ~ok_tmp
        ok = 0;
    end
    
    preproc_chckbx_Callback(handles.preproc_chckbx,[],handles)
    if ~ok && exist('warn','var')~=1
        warn = 'Preprocessing';
    end
    
    
    % DESIGN & ESTIMATE
    %======================================================================
    if any(opt_safe.stat_des.run==[0 1])
        set(handles.design_chckbx,'Value',opt_safe.stat_des.run)
        design_chckbx_Callback(handles.design_chckbx,[],handles)
    else
        ok = 0;
    end
    
    set(handles.design_funcprefix_edit,'String',opt_safe.stat_des.pref_func )
    set(handles.res_strucprefix_edit  ,'String',opt_safe.stat_des.pref_struc)
    set(handles.res_struc_extra_edit  ,'String',opt_safe.stat_des.tag_struc )
    
    if any(opt_safe.stat_des.use_rp==[0 1])
        set(handles.design_useRP,'Value',opt_safe.stat_des.use_rp)
    else
        ok = 0;
    end
    
    set(handles.design_hpf_edit,'String',opt_safe.stat_des.hpf)
    ok_tmp = design_hpf_edit_Callback(handles.design_hpf_edit,[],[]);
    if ~ok_tmp
        ok = 0;
    end
    
    switch opt_safe.stat_des.bf
        case 'can_hrf'
            set(handles.design_hrf_menu,'Value',1)
        case 'gamma'
            set(handles.design_hrf_menu,'Value',2)
        case 'fourier'
            set(handles.design_hrf_menu,'Value',3)
        case 'fourier_han'
            set(handles.design_hrf_menu,'Value',4)
        case 'fir'
            set(handles.design_hrf_menu,'Value',5)
        otherwise
            ok = 0;
    end
    design_hrf_menu_Callback(handles.design_hrf_menu,[],handles)
    
    set(handles.design_hrf_length_edit,'String',num2str(opt_safe.stat_des.bf_length))
    ok_tmp = design_hrf_length_edit_Callback(handles.design_hrf_length_edit,[],[]);
    if ~ok_tmp
        ok = 0;
    end
    
    set(handles.design_hrf_nmbrbf_edit,'String',num2str(opt_safe.stat_des.bf_order))
    design_hrf_nmbrbf_edit_Callback(handles.design_hrf_nmbrbf_edit,[],handles)
    
    set(handles.design_hrf_delay_edit,'String',num2str(opt_safe.stat_des.delay))
    ok_tmp = design_hrf_delay_edit_Callback(handles.design_hrf_delay_edit,[],handles);
    if ~ok_tmp
        ok = 0;
    else
        set(handles.design_hrf_delay_edit, 'String',...
            regexprep(get(handles.design_hrf_delay_edit, 'String'), '\s+', ' ')) % remove repeated spaces
        % Change tooltip string
        gui_chg_tooltip(handles.design_hrf_delay_edit, ...
            handles.design_hrf_delay_txt, 'Add an additional delay to the HRF (can be a set of numbers)')
    end
    
    if ~ok && exist('warn','var')~=1
        warn = 'Design';
    end
    
    
    % CONTRASTS
    %======================================================================
    if any(opt_safe.stat_con==[0 1])
        set(handles.contr_chckbx,'Value',opt_safe.stat_con)
    elseif exist('warn','var')~=1
        warn = 'Contrasts';
        ok   = 0;
    end
    
    
    % RESULTS
    %======================================================================
    
    % Results checkbox
    %-----------------
    if any(opt_safe.stat_res.run==[0 1])
        set(handles.res_chckbx,'Value',opt_safe.stat_res.run)
        res_chckbx_Callback(handles.res_chckbx,[],handles)
    else
        ok = 0;
    end
    
    % P-value
    %--------
    switch opt_safe.stat_res.thresh.adj
        case 'none'
            set(handles.res_FWE,'Value',0)
        case 'FWE'
            set(handles.res_FWE,'Value',1)
        otherwise
            ok = 0;
    end
    
    set(handles.res_pvalue_edit, 'String', opt_safe.stat_res.thresh.P)
    ok_tmp = res_pvalue_edit_Callback(handles.res_pvalue_edit,[],[]);
    if ~ok_tmp
        ok = 0;
    end
    
    set(handles.res_vxl_thresh_edit, 'String', opt_safe.stat_res.thresh.ET)
    ok_tmp = res_vxl_thresh_edit_Callback(handles.res_vxl_thresh_edit,[],[]);
    if ~ok_tmp
        ok = 0;
    end
    
    % Table and PostScript file
    %--------------------------
    if any(opt_safe.stat_res.ps_save==[0 1])
        set(handles.res_PS,'Value',opt_safe.stat_res.ps_save)
    else
        ok   = 0;
        warn = 'Save PS file';
    end
    
    % Auxiliary files
    %----------------
    %opt_safe.stat_res.save_slover -- cannot change via GUI => just check
    %opt_safe.stat_res.save_thr    -- cannot change via GUI => just check
    if ~any(opt_safe.stat_res.save_slover==[0 1])
        ok   = 0;
        warn = 'opt_safe.stat_res.save_slover';
    end
    if ~any(opt_safe.stat_res.save_thr==[0 1])
        ok   = 0;
        warn = 'opt_safe.stat_res.save_thr';
    end
    
    % Save image with slices
    %-----------------------
    if any(opt_safe.stat_res.slices.save==[0 1])
        set(handles.res_slice_view,'Value',opt_safe.stat_res.slices.save)
    else
        ok = 0;
    end
    %opt_safe.stat_res.slices.img_ext -- cannot change via GUI => just check
    tmp = opt_safe.stat_res.slices.img_ext;
    if ~ischar(tmp) || size(tmp,1)~=1 || size(tmp,2)>4
        ok   = 0;
        warn = 'opt_safe.stat_res.slices.img_ext';
    end
    tmp = opt_safe.stat_res.slices.step;
    if ~isnumeric(tmp) || length(tmp)~=1 || tmp<0
        % must be a positive number (but step=0 => ignore it)
        ok = 0;
    end
    tmp = opt_safe.stat_res.slices.number;
    if ~isnumeric(tmp) || length(tmp)~=1 || floor(tmp)~=ceil(tmp) || tmp<0
        % must be a positive integer number (but number=0 => ignore it)
        ok = 0;
    end
    tmp = opt_safe.stat_res.slices.slices;
    if ~isnumeric(tmp) || isempty(tmp) || min(size(tmp))~=1
        % must be an integer one-dimensional matrix
        ok = 0;
    end
    tmp = opt_safe.stat_res.slices.orient;
    if ~iscell(tmp) || isempty(tmp)
        ok = 0;
    else
        tmp = sort(tmp);
        while ~isempty(tmp)
            if strcmp(tmp{1},'axial') || ...
                    strcmp(tmp{1},'coronal') || ...
                    strcmp(tmp{1},'sagittal')
                tmp(1) = [];
            else
                ok = 0;
            end
        end
    end
    if ~ok && exist('warn','var')~=1
        warn = 'Slice view';
    end
    
    % Group of images
    %----------------
    if any(opt_safe.stat_res.group.save==[0 1])
        set(handles.res_union,'Value',opt_safe.stat_res.group.save)
    else
        ok = 0;
    end
    tmp = opt_safe.stat_res.group.number;
    if ~isnumeric(tmp) || length(tmp)~=1 || floor(tmp)~=ceil(tmp) || tmp<1
        % must be a positive integer number
        ok = 0;
    end
    if ~any(opt_safe.stat_res.group.scale.adj==[0 1])
        ok = 0;
    end
    for ii=1:4
        switch ii
            case 1
                tmp = opt_safe.stat_res.group.scale.ftest;
            case 2
                tmp = opt_safe.stat_res.event_seq.scale.ftest;
            case 3
                tmp = opt_safe.stat_res.event_seq.scale.posbold;
            case 4
                tmp = opt_safe.stat_res.event_seq.scale.negbold;
        end
        if ~strcmp(tmp,'max') && ( ...
                min(size(tmp))~=1 || max(size(tmp))~=2 || ...
                ~isnumeric(tmp) || any(tmp<0) || tmp(1)>tmp(2) )
            % must be a 1x2 or 2x1 numeric matrix, with non-negative
            % numbers, and the 1st element smaller than the 2nd one
            ok = 0;
        end
    end
    tmp = opt_safe.stat_res.group.scale.ttest;
    if ~strcmp(tmp,'max') && (...
            min(size(tmp))~=1 || max(size(tmp))~=4 || ~isnumeric(tmp) || ...
            any(tmp<0) || tmp(1)>tmp(2) || tmp(3)>tmp(4) )
        % must be a 1x4 or 4x1 numeric matrix, with non-negative numbers,
        % the first element smaller than the second one, and the third
        % element smaller than the fourth one
        ok = 0;
    end
    if ~ok && exist('warn','var')~=1
        warn = 'Union or Sequence of delays';
    end
    
    % Sequence of images
    %-------------------
    if any(opt_safe.stat_res.event_seq.save==[0 1])
        set(handles.res_seq,'Value',opt_safe.stat_res.event_seq.save)
    else
        ok = 0;
    end
    if ~ok && exist('warn','var')~=1
        warn = 'Sequence of delays';
    end
    % the statistical thresholds were tested above (with "Group ...")
    
    
    % Results checkbox
    %-----------------
    res_chckbx_Callback(handles.res_chckbx,[],handles)
    
    % Mask
    %-----
    %opt_safe.stat_res.mask.use        -- cannot change via GUI => just check
    %opt_safe.stat_res.mask.imgs_text  -- cannot change via GUI => just check
    %opt_safe.stat_res.mask.incl_excl  -- cannot change via GUI => just check
    %opt_safe.stat_res.mask.bin        -- cannot change via GUI => just check
    %opt_safe.stat_res.mask.bin_thresh -- cannot change via GUI => just check
    if ~any(opt_safe.stat_res.mask.use==[0 1])
        ok = 0;
    end
    if ~iscellstr(opt_safe.stat_res.mask.imgs_text)
        ok = 0;
    end
    tmp = opt_safe.stat_res.mask.incl_excl;
    if ~strcmp(tmp,'inclusive') && ~strcmp(tmp,'exclusive')
        ok = 0;
    end
    if ~any(opt_safe.stat_res.mask.bin==[0 1])
        ok = 0;
    end
    tmp = opt_safe.stat_res.mask.bin_thresh;
    if ~isnumeric(tmp) || length(tmp)~=1 || tmp<=0
        % must be a positive number
        ok = 0;
    end
    if ~ok && exist('warn','var')~=1
        warn = 'Results (mask)';
    end
    
    if ~any(opt_safe.stat_res.cleanup==[0 1])
        ok = 0;
    else
        set(handles.res_cleanup,'Value',opt_safe.stat_res.cleanup)
    end
    
    if ~ok && exist('warn','var')~=1
        warn = 'Results';
    end
    
    safe_load_finish(handles,opt_safe)
    
    if ~ok
        warndlg(sprintf('%s: Problem with input data.',warn),'Loading...',...
            'modal')
    end
catch
    ok = 0;
    safe_load_finish(handles,opt_safe)
    warndlg('Invalid parameter structure.','Loading...','modal')
end

end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function safe_load_finish(handles,opt_safe)
% Enable/Disable "Basic input" panel
%basic_input_panel(handles)

% Check number of data directories and worksheet names & Results options
ok = check_param_part(opt_safe,handles);

% Enable/Disable "Options" (for Results) button
res_opt_button_enable(handles)

end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% function basic_input_panel(handles)
% 
% This function is commented because the "Basic input" panel will stay
% enabled.
% 
% % Enable/Disable "Basic input" panel
%     if get(handles.preproc_chckbx,'Value') || ...
%             get(handles.design_chckbx,'Value') || ...
%             get(handles.res_chckbx,'Value')
%         set(handles.basic_input_panel         ,'BorderType','beveledout')
%         set(handles.basic_input_xlsfile_button,'Enable','on')
%         set(handles.basic_input_xlsfile_txt   ,'Enable','on')
%         if ~isempty(get(handles.basic_input_xlsfile_txt,'UserData'))
%             set(handles.basic_input_datadir_button,'Enable','on')
%             set(handles.basic_input_datadir_txt   ,'Enable','on')
%             if ~isempty(get(handles.basic_input_datadir_txt,'UserData'))
%                 set(handles.basic_input_worksheet_button ,'Enable','on')
%                 set(handles.basic_input_worksheet_txt,'Enable','on')
%             end
%         end
%         set(handles.iseegfmri,'Enable','on')
%     else
%         set(handles.basic_input_panel           ,'BorderType','beveledin')
%         set(handles.basic_input_xlsfile_button  ,'Enable','off')
%         set(handles.basic_input_xlsfile_txt     ,'Enable','off')
%         set(handles.basic_input_datadir_button  ,'Enable','off')
%         set(handles.basic_input_datadir_txt     ,'Enable','off')
%         set(handles.basic_input_worksheet_button,'Enable','off')
%         set(handles.basic_input_worksheet_txt   ,'Enable','off')
%         set(handles.iseegfmri,'Enable','off')
%     end
% 
% end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% function basic_input_fit_text(N,txt_handle,txt)
% 
% This function is commented because the text should always fit the space,
% and the window is resizable.
% 
% % Fit the text TXT in 'String' field of handle TXT_HANDLE in the "Basic
% % input" panel. N is 1, 2 or 3, dependending on the vertical position of
% % the text.
% 
% % This function is useful when the names of the XLS file, worksheet and
% % directory are written in the GUI. The units below are in points for the
% % fonts and in characteres for the position and size.
% 
% % Load defaults
% %--------------
% font = 10; % 'FontUnits' = points
% if N==1 % 'Units' = characters
%     pos = [28.6 8.5 43.2 2];
% elseif N==2
%     pos = [28.6 5.75 43.2 2];
% else % N=3
%     pos = [28.6 3 43.2 2];
% end
% 
% % Change text size and position
% %------------------------------
% if length(txt)>36
%     font_chg = -2;
%     pos_chg  = -0.1;
% else
%     font_chg = 0;
%     pos_chg  = 0;
% end
% 
% % Adjust text and write
% %----------------------
% set(txt_handle,'FontSize',font+font_chg)
% set(txt_handle,'Position',pos+[0 pos_chg 0 0])
% set(txt_handle,'String',txt)
% 
% end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = check_param_part(opt_safe,handles)
% Check some parameters

ok = 1;
try
    
    % Check the consistency in the basic input
    %-----------------------------------------
    
    % 1) Number of XLS files listed in data directory
    tmp = length(opt_safe.basic.sbj_folder);
    if tmp~=length(opt_safe.basic.xlsfile) && ~isempty(opt_safe.basic.xlsfile)
        ok = 0;
        warndlg('Basic input: Inconsistent number of XLS files listed in data directories.',...
            'Invalid input','modal')
        return
    end
    
    % 2) Number of XLS files listed in data dir. and in worksheet names
    if tmp~=length(opt_safe.basic.worksheet)
        ok = 0;
        warndlg('Basic input: Inconsistent number of XLS files listed in worksheet names and in data directories.',...
            'Invalid input','modal')
        return
    end
    
    % 3) Number of data dirs. and worksheet names for each XLS file
    for ii=1:tmp
        if ~iscellstr(opt_safe.basic.sbj_folder{ii}) || ...
                ~iscellstr(opt_safe.basic.worksheet{ii}) || ...
                length(opt_safe.basic.sbj_folder{ii})~=length(opt_safe.basic.worksheet{ii}) || ...
                ~all( cellfun(@isempty,opt_safe.basic.sbj_folder{ii})==...
                cellfun(@isempty,opt_safe.basic.worksheet{ii}) )
            ok = 0;
            break
        end
    end
    if ~ok
        warndlg('Basic input: The number of data directories and worksheet names must be the same for each XLS file.',...
            'Invalid input','modal')
        return
    end
    
    % Check if any result option was selected
    %----------------------------------------
    if get(handles.res_chckbx,'Value')
        if ~opt_safe.stat_con                     && ...
                ~opt_safe.stat_res.ps_save        && ...
                ~opt_safe.stat_res.slices.save    && ...
                ~opt_safe.stat_res.group.save     && ...
                ~opt_safe.stat_res.event_seq.save && ...
                ~opt_safe.stat_res.cleanup
            warndlg('Results: You must select at least one output option.',...
                'Invalid input','modal')
            ok = 0;
        elseif (opt_safe.stat_res.group.save || ...
                opt_safe.stat_res.event_seq.save) && ...
                length(opt_safe.stat_des.delay)==1
            warndlg('Results: You must select at least two delays to make a union or sequence of delays.',...
                'Invalid input','modal')
            ok = 0;
%             if  ~opt_safe.stat_res.ps_save && ...
%                     ~opt_safe.stat_res.slices.save
%                 ok = 0;
%             end
        end
    end
    
catch
    ok = 0;
end

end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function files_prefix(handles)
    if get(handles.preproc_normalize,'Value')
        set(handles.design_funcprefix_edit, 'String', 'swa')
        set(handles.res_strucprefix_edit  , 'String', 'wm')
    else
        set(handles.design_funcprefix_edit, 'String', 'sar')
        set(handles.res_strucprefix_edit  , 'String', '')
    end
    tmp = get(handles.design_funcprefix_edit, 'String');
    if ~get(handles.preproc_slice_timing,'Value')
        tmp = strrep(tmp,'a','');
    end
    if ~get(handles.preproc_smooth,'Value')
        tmp = strrep(tmp,'s','');
    end
    set(handles.design_funcprefix_edit, 'String', tmp)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% function struc_pref_enable(handles)
% % Enable/Disable "Structural file prefix" in Design
%     if ~get(handles.preproc_chckbx,'Value') && ...
%             ( get(handles.design_chckbx,'Value') || ...
%             ( get(handles.res_chckbx,'Value') && ...
%             ( get(handles.res_PS,'Value') || ...
%             get(handles.res_slice_view,'Value') || ...
%             get(handles.res_union,'Value') || ...
%             get(handles.res_seq,'Value') || ...
%         state = 'on';
%     else
%         state = 'off';
%     end
%     set(handles.res_strucprefix_txt ,'Enable',state)
%     set(handles.res_strucprefix_edit,'Enable',state)
%     set(handles.res_struc_extra_txt ,'Enable',state)
%     set(handles.res_struc_extra_edit,'Enable',state)
% end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function hrf_delays_enable(handles)
% Enable/Disable "HRF delays" in Design
    if get(handles.design_chckbx,'Value') || get(handles.res_chckbx,'Value')
        state = 'on';
    else
        state = 'off';
    end
    set(handles.design_hrf_delay_edit, 'Enable', state)
    set(handles.design_hrf_delay_txt , 'Enable', state)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function res_opt_button_enable(handles)
% Enable/Disable "Options" (for Results) button
    if get(handles.res_slice_view,'Value') || ...
            get(handles.res_union,'Value') || ...
            get(handles.res_seq,'Value')
        state = 'on';
    else
        state = 'off';
    end
    set(handles.res_opt_button,'Enable',state)
end
%--------------------------------------------------------------------------


%==========================================================================
%                             BASIC INPUT
%==========================================================================

%--------------------------------------------------------------------------
function basic_input_xlsfile_button_Callback(hObject, eventdata, handles)
    
    % Choose file(s)
    %---------------
%     [xlsfile,fpath] = uigetfile({'*.xls;*.xlsx','Excel files (*.xls,*.xlsx)'},...
%         'Select the XLS file',pwd,'MultiSelect','on');
    xlsfile = get(handles.basic_input_xlsfile_txt,'UserData');
    curr_number = length(xlsfile);
    xlsfile = uipickfiles(...
        'FilterSpec',pwd,...
        'REFilter'  ,'\.xls$|\.xlsx$',...
        'Prompt'    ,'Select the XLS file',...
        'Append'    ,xlsfile,...
        'Output'    ,'cell');
    if isequal(xlsfile,0) % clicked Cancel
        return
    end
    xlsfile = xlsfile(:);
    
    % Count number of XLS files
    %--------------------------
    cnt = length(xlsfile);
    if cnt>1
        tmp_text = sprintf('%d files selected',cnt);
    elseif cnt==1
        tmp_text = '1 file selected';
    else %cnt=0
        tmp_text = '(no file selected)';
        xlsfile  = {''}; % guarantee that it is a 1x1 cell array
    end
    
    % Set data
    %---------
    %basic_input_fit_text(1,handles.basic_input_xlsfile_txt,tmp_text)
    set(handles.basic_input_xlsfile_txt,'String'  ,tmp_text)
    set(handles.basic_input_xlsfile_txt,'UserData',xlsfile)
    
    % Enable/Disable worksheet and directory selection
    %-------------------------------------------------
    if cnt>0
        state = 'on';
        if curr_number~=cnt % changed the number of XLS-files
            worksheet = cell(1,length(xlsfile));
            for xx=1:length(xlsfile)
                worksheet{xx} = {''};
            end
            set(handles.basic_input_worksheet_txt ,'String'  ,'(no worksheet selected)')
            set(handles.basic_input_worksheet_txt ,'UserData',worksheet)
            %set(handles.basic_input_datadir_button,'Enable'  ,'off')
            %set(handles.basic_input_datadir_txt   ,'Enable'  ,'off')
            set(handles.basic_input_datadir_txt   ,'String'  ,'(no directory selected)')
            set(handles.basic_input_datadir_txt   ,'UserData',worksheet)
        end
    else
        state = 'off';
    end
    set(handles.basic_input_worksheet_button,'Enable',state)
    set(handles.basic_input_worksheet_txt   ,'Enable',state)
    set(handles.basic_input_datadir_button  ,'Enable',state)
    set(handles.basic_input_datadir_txt     ,'Enable',state)
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function basic_input_worksheet_button_Callback(hObject, eventdata, handles)
	
    xlsfile = get(handles.basic_input_xlsfile_txt,'UserData');
    num_xls = length(xlsfile);
    
    % Get only file name from XLS file(s)
    %------------------------------------
    xlsfile_name = xlsfile;
    %[tmp,xlsfile_name] = cellfun(@fileparts,xlsfile,'UniformOutput',false); % file name without extension
    for xx=1:num_xls
        tmp = strfind(xlsfile{xx},filesep);
        tmp = tmp(end);
        xlsfile_name{xx} = xlsfile{xx}(tmp+1:end);
    end
    
    % Check for previous data
    %------------------------
    worksheet = get(handles.basic_input_worksheet_txt,'UserData');
    if isempty(worksheet) || ~iscell(worksheet) || length(worksheet)~=num_xls % wrong format
        worksheet = cell(1,num_xls);
        for xx=1:num_xls
            worksheet{xx} = '';
        end
    else
        for xx=1:num_xls
            if isempty(worksheet{xx}) || ~iscellstr(worksheet{xx}) % wrong format
                worksheet{xx} = '';
            else
                % When Option 1 below is used:
%                 tmp = [];
%                 for ww=1:length(worksheet{xx})
%                     tmp = [tmp strrep(worksheet{xx}{ww},' ','') ' '];
%                 end
%                 tmp(end) = [];
%                 worksheet{xx} = tmp;
            end
        end
    end
    
    % Apparent number of XLS files
    %-----------------------------
    num_xls_wkrsht = length(worksheet);
    if num_xls_wkrsht==1
        tmp_text = '(1 XLS file)';
    else
        tmp_text = sprintf('(%d XLS files)',num_xls_wkrsht);
    end
    
    % Get worksheet names
    %--------------------
    % Option 1) Input dialog: the user must insert the spreadhsheet names
    % manually
%     opt_dlg.Resize = 'on';
%     num_prompts = 5; % use figures with up to 5 prompts
%     num_figs = ceil(num_xls_wkrsht/num_prompts);
%     for xx=1:num_figs
%         tmp = ((xx-1)*num_prompts + 1) : (xx*num_prompts);
%         if xx==num_figs
%             tmp(tmp>num_xls_wkrsht) = [];
%         end
%         %worksheet_tmp{xx} = inputdlg(xlsfile_name(tmp),...
%         try
%             worksheet(tmp) = inputdlg(xlsfile_name(tmp),...
%                 'Worksheet names (separate with spaces)',1,worksheet(tmp),opt_dlg);
%         catch % clicked "Cancel" => returns: {}
%             return
%         end
%     end
    
    % Option 2) List dialog: MATLAB reads the spreadsheet names with XLSFINFO
    warning off MATLAB:xlsfinfo:ActiveX % disable warning
    cnt = 0;
    set(handles.safe_main_window,'Pointer','watch')
    drawnow; pause(.05)
    for xx=1:num_xls % loop for all XLS files
        [tmp,sheets] = xlsfinfo(xlsfile{xx});
        if isempty(tmp) % "xlsread" cannot read xlsfile{xx}
            error('XLSREAD is unable to read %s',xlsfile{xx})
            return
        end
        initval = [];
        for ww=1:length(worksheet{xx})
            initval = [initval find(strcmp(worksheet{xx}{ww},sheets))];
        end
        if isempty(initval)
            initval = 1;
        end
        [worksheet{xx},ok] = listdlg(...
            'ListString',sheets,...
            'SelectionMode','multiple',...
            'ListSize',[200 300],...
            'InitialValue',initval,...
            'Name','',...
            'PromptString',{'Select the spreadsheets for',sprintf('%s:',xlsfile_name{xx})});
        % default size ('ListSize'): [160 300]
        if ~ok % clicked on Cancel or closed the dialog box
            set(handles.safe_main_window,'Pointer','arrow')
            drawnow; pause(.05)
            return
        end
        cnt = cnt + size(worksheet{xx},2);
        worksheet{xx} = sheets(worksheet{xx});
        worksheet{xx} = worksheet{xx}(:);
    end
    set(handles.safe_main_window,'Pointer','arrow')
    drawnow; pause(.05)
    warning on MATLAB:xlsfinfo:ActiveX % restore warning state
    
    % Count number of worksheets
    %---------------------------
    % When Option 1 above is used:
%     cnt = 0;
%     for xx=1:num_xls_wkrsht
%         if ~isempty(worksheet{xx})
%             worksheet{xx} = textscan(worksheet{xx},'%s');
%             worksheet{xx} = worksheet{xx}{1}';
%         else
%             worksheet{xx} = {''};
%         end
%         cnt = cnt + sum(~cellfun(@isempty,worksheet{xx})); % count non-empty cells
%     end
    if cnt>1
        tmp_text = sprintf('%d worksheets selected %s',cnt,tmp_text);
    elseif cnt==1
        tmp_text = sprintf('1 worksheet selected %s',tmp_text);
    else %cnt=0
        tmp_text = '(no worksheet selected)';
    end
    
    % Set data
    %---------
    %basic_input_fit_text(2,handles.basic_input_worksheet_txt,tmp_text)
    set(handles.basic_input_worksheet_txt,'String',tmp_text)
    set(handles.basic_input_worksheet_txt,'UserData',worksheet)
    
    % Enable/Disable data directory selection
    %----------------------------------------
%     if cnt>0
%         state = 'on';
%     else
%         state = 'off';
%     end
%     set(handles.basic_input_datadir_button,'Enable',state)
%     set(handles.basic_input_datadir_txt   ,'Enable',state)
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function basic_input_datadir_button_Callback(hObject, eventdata, handles)
    
    % Get previous data
    %------------------
    xlsfile   = get(handles.basic_input_xlsfile_txt  ,'UserData');
    worksheet = get(handles.basic_input_worksheet_txt,'UserData');
    data_dir  = get(handles.basic_input_datadir_txt  ,'UserData');
    if isempty(data_dir) || ~iscell(data_dir) || length(data_dir)~=length(xlsfile) % wrong format
        data_dir = cell(1,length(xlsfile));
    end
    
    % Select directories
    %-------------------
    max_char = 50;
    for xx=1:length(xlsfile) % loop to select folders for each XLS file
        
        tmp = [];
        for ww=1:length(worksheet{xx})
            tmp = [tmp strrep(worksheet{xx}{ww},' ','') ', '];
        end
        if length(tmp)>2
            tmp(end-1:end) = [];
            worksheet{xx} = tmp;
            prompt = 'Add the directories for';
        else
            prompt = '';
        end
        
        tmp = xlsfile{xx};
        if length(tmp)>max_char
            tmp = ['...' tmp(length(tmp)-(max_char-4):length(tmp))];
        end
        if ~strcmp(prompt,'')
            prompt = sprintf('%s: %s (%s)',prompt,worksheet{xx},tmp);
        else
            prompt = sprintf('Add the directories for %s',tmp);
        end
        tmp = uipickfiles(...
            'FilterSpec',pwd,...
            'REFilter'  ,'$',...
            'Prompt'    ,prompt,...
            'Append'    ,data_dir{xx},...
            'Output'    ,'cell');
        if ~isempty(tmp) % clicked "Done", but left something selected; or clicked "Cancel" or closed the window (tmp=0)
            if ~isequal(tmp,0) % didn't click "Cancel"
                data_dir{xx} = tmp(:);
            elseif isempty(data_dir{xx})
                data_dir{xx} = {''};
            end
        else
            data_dir{xx} = {''};
        end
        
    end
    
    % Apparent number of XLS files
    %-----------------------------
    nmbr_xls_dir = length(data_dir);
    if nmbr_xls_dir==1
        tmp_text = '(1 XLS file)';
    else
        tmp_text = sprintf('(%d XLS files)',nmbr_xls_dir);
    end
    
    % Count number of directories
    %----------------------------
    cnt = 0;
    for xx=1:length(data_dir)
        if isempty(data_dir{xx})
            data_dir{xx} = {''};
        end
        cnt = cnt + sum(~cellfun(@isempty,data_dir{xx})); % count non-empty cells
    end
    if cnt>1
        tmp_text = sprintf('%d directories selected %s',cnt,tmp_text);
    elseif cnt==1
        tmp_text = sprintf('1 directory selected %s',tmp_text);
    else %cnt=0
        tmp_text = '(no directory selected)';
    end
    
    % Set data
    %---------
    %basic_input_fit_text(3,handles.basic_input_datadir_txt,tmp_text)
    set(handles.basic_input_datadir_txt,'String'  ,tmp_text)
    set(handles.basic_input_datadir_txt,'UserData',data_dir)
    
end
%--------------------------------------------------------------------------



%==========================================================================
%                             PREPROCESSING
%==========================================================================

%--------------------------------------------------------------------------
function preproc_chckbx_Callback(hObject, eventdata, handles)
    
    if get(hObject,'Value') % preproc = 1
        set(handles.preproc_panel,'BorderType','beveledout')
        state = 'on';
        if get(handles.preproc_normalize,'Value')
            set(handles.preproc_norm_factor_struc_txt ,'Enable','on')
            set(handles.preproc_norm_factor_struc_edit,'Enable','on')
            set(handles.preproc_norm_factor_func_txt  ,'Enable','on')
            set(handles.preproc_norm_factor_func_edit ,'Enable','on')
        end
        if get(handles.preproc_smooth,'Value')
            set(handles.preproc_smoothfwhm_edit,'Enable','on')
        end
        set(handles.design_funcprefix_txt  ,'Enable','off')
        set(handles.design_funcprefix_edit ,'Enable','off')
    else % preproc = 0
        set(handles.preproc_panel,'BorderType','beveledin')
        state = 'off';
        set(handles.preproc_norm_factor_struc_txt ,'Enable','off')
        set(handles.preproc_norm_factor_struc_edit,'Enable','off')
        set(handles.preproc_norm_factor_func_txt  ,'Enable','off')
        set(handles.preproc_norm_factor_func_edit ,'Enable','off')
        set(handles.preproc_smoothfwhm_edit       ,'Enable','off')
        if get(handles.design_chckbx,'Value')
            set(handles.design_funcprefix_txt  ,'Enable','on')
            set(handles.design_funcprefix_edit ,'Enable','on')
        end
    end
    set(handles.preproc_slice_timing,'Enable',state)
    set(handles.preproc_normalize   ,'Enable',state)
    set(handles.preproc_smooth      ,'Enable',state)
    
    %struc_pref_enable(handles)
    %basic_input_panel(handles)
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function preproc_slice_timing_Callback(hObject, eventdata, handles)
    files_prefix(handles)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function preproc_normalize_Callback(hObject, eventdata, handles)
    files_prefix(handles)
    if get(hObject,'Value')
        state = 'on';
    else
        state = 'off';
    end
    set(handles.preproc_norm_factor_struc_txt ,'Enable',state)
    set(handles.preproc_norm_factor_struc_edit,'Enable',state)
    set(handles.preproc_norm_factor_func_txt  ,'Enable',state)
    set(handles.preproc_norm_factor_func_edit ,'Enable',state)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = preproc_norm_factor_struc_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    ok  = 1;
    if length(tmp)~=3 || ~isreal(tmp) || any(tmp==0)
        warndlg('Normalization factor: Insert a 1x3 non-zero real numeric matrix',...
            'Preprocessing: Invalid input','modal')
        set(hObject,'String','[-1 -1 -1]')
        tmp = [-1 -1 -1];
        ok  = 0;
    end
    set(hObject,'UserData',tmp)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = preproc_norm_factor_func_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    ok  = 1;
    if length(tmp)~=3 || ~isreal(tmp) || any(tmp==0)
        warndlg('Normalization factor: Insert a 1x3 non-zero real numeric matrix',...
            'Preprocessing: Invalid input','modal')
        set(hObject,'String','[-1 -1 -1]')
        ok  = 0;
    end
    set(hObject,'UserData',tmp)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function preproc_smooth_Callback(hObject, eventdata, handles)
    files_prefix(handles)
    if get(hObject,'Value')
        state = 'on';
    else
        state = 'off';
    end
    set(handles.preproc_smoothfwhm_edit,'Enable',state)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = preproc_smoothfwhm_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    ok  = 1;
    if length(tmp)~=3 || ~isreal(tmp) || all(tmp==0)
        warndlg('FWHM: Insert a 1x3 real numeric matrix',...
            'Preprocessing: Invalid input','modal')
        set(hObject,'String','[-2 -2 -2]')
        ok  = 0;
    end
    set(hObject,'UserData',tmp)
end
%--------------------------------------------------------------------------



%==========================================================================
%                               DESIGN
%==========================================================================

%--------------------------------------------------------------------------
function design_chckbx_Callback(hObject, eventdata, handles)
    
    if get(hObject,'Value')
        state = 'on';
        set(handles.design_panel,'BorderType','beveledout')
        if ~get(handles.preproc_chckbx,'Value') % not preprocessing
            set(handles.design_funcprefix_txt ,'Enable','on')
            set(handles.design_funcprefix_edit,'Enable','on')
            %set(handles.res_strucprefix_txt   ,'Enable','on')
            %set(handles.res_strucprefix_edit  ,'Enable','on')
            %set(handles.res_struc_extra_txt   ,'Enable','on')
            %set(handles.res_struc_extra_edit  ,'Enable','on')
        end
        if get(handles.design_hrf_menu,'Value')>1 % not the canonical HRF
            set(handles.design_hrf_length_txt ,'Enable','on')
            set(handles.design_hrf_length_edit,'Enable','on')
        end
    else
        state = 'off';
        set(handles.design_panel,'BorderType','beveledin')
        set(handles.design_funcprefix_txt  ,'Enable','off')
        set(handles.design_funcprefix_edit ,'Enable','off')
        set(handles.design_hrf_length_txt  ,'Enable','off')
        set(handles.design_hrf_length_edit ,'Enable','off')
        %struc_pref_enable(handles)
    end
    set(handles.design_hpf_txt        ,'Enable',state)
    set(handles.design_hpf_edit       ,'Enable',state)
    set(handles.design_useRP          ,'Enable',state)
    set(handles.design_hrf_menu       ,'Enable',state)
    set(handles.design_hrf_nmbrbf_txt ,'Enable',state)
    set(handles.design_hrf_nmbrbf_edit,'Enable',state)
    hrf_delays_enable(handles)
    
    %basic_input_panel(handles)
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = design_hpf_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    ok  = 1;
    if length(tmp)~=1 || tmp<=0
        warndlg('High-pass filter: Insert a positive number',...
            'Design: Invalid input','modal')
        set(hObject,'String','128')
        tmp = 128;
        ok  = 0;
    end
    set(hObject,'UserData',tmp)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function design_hrf_menu_Callback(hObject, eventdata, handles)
    if get(hObject,'Value')==1
        state = 'off';
        if str2double(get(handles.design_hrf_nmbrbf_edit,'String'))>3
            set(handles.design_hrf_nmbrbf_edit,'String','3')
        end
    else
        state = 'on';
    end
    set(handles.design_hrf_length_txt ,'Enable',state)
    set(handles.design_hrf_length_edit,'Enable',state)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function design_hrf_nmbrbf_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    ok  = 1;
    if length(tmp)~=1 || tmp<=0 || ceil(tmp)~=floor(tmp)
        warndlg('Basis functions: Insert a positive integer',...
            'Design: Invalid input','modal')
        set(hObject,'String','1')
        tmp = 1;
        ok  = 0;
    end
    if get(handles.design_hrf_menu,'Value')==1 && tmp>3
        warndlg('The maximum number of basis functions is 3',...
            'Design: Invalid input','modal')
        set(hObject,'String','3')
        tmp = 3;
        ok  = 0;
    end
    set(hObject,'UserData',tmp)
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = design_hrf_length_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    ok  = 1;
    if length(tmp)~=1 || tmp<=0
        warndlg('HRF length: Insert a positive number',...
            'Design: Invalid input','modal')
        set(hObject,'String','32')
        tmp = 32;
        ok  = 0;
    end
    set(hObject,'UserData',tmp)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = design_hrf_delay_edit_Callback(hObject, eventdata, handles)
    
    tmp = str2num(get(hObject,'String'));
    ok  = 1;
    if length(tmp)<1
        warndlg('Delays: Insert a real numeric matrix',...
            'Design: Invalid input','modal')
        set(hObject,'String','0')
        tmp = 0;
        ok  = 0;
    end
    set(hObject,'UserData',tmp)
    
    % Change tooltip string
    gui_chg_tooltip(hObject, handles.design_hrf_delay_txt, ...
        'Add an additional delay to the HRF (can be a set of numbers)')
    
end
%--------------------------------------------------------------------------



%==========================================================================
%                               RESULTS
%==========================================================================

%--------------------------------------------------------------------------
function res_chckbx_Callback(hObject, eventdata, handles)
    
    if get(hObject,'Value')
        set(handles.res_panel,'BorderType','beveledout')
        state = 'on';
        res_opt_button_enable(handles)
    else
        set(handles.res_panel,'BorderType','beveledin')
        state = 'off';
        set(handles.res_opt_button,'Enable','off')
    end
    set(handles.contr_chckbx        ,'Enable',state)
    set(handles.res_pvalue_txt      ,'Enable',state)
    set(handles.res_pvalue_edit     ,'Enable',state)
    set(handles.res_FWE             ,'Enable',state)
    set(handles.res_vxl_thresh_txt  ,'Enable',state)
    set(handles.res_vxl_thresh_edit ,'Enable',state)
    set(handles.res_PS              ,'Enable',state)
    set(handles.res_slice_view      ,'Enable',state)
    set(handles.res_union           ,'Enable',state)
    set(handles.res_seq             ,'Enable',state)
    set(handles.res_cleanup         ,'Enable',state)
    set(handles.res_strucprefix_txt ,'Enable',state)
    set(handles.res_strucprefix_edit,'Enable',state)
    set(handles.res_struc_extra_txt ,'Enable',state)
    set(handles.res_struc_extra_edit,'Enable',state)
    hrf_delays_enable(handles)
    
    %struc_pref_enable(handles)
    %basic_input_panel(handles)
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = res_pvalue_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    ok  = 1;
    if length(tmp)~=1 || tmp<=0 || tmp>1
        warndlg('P-value: Insert a number between 0 and 1',...
            'Results: Invalid input','modal')
        set(hObject, 'String', '0.001')
        tmp = 0.001;
        ok  = 0;
    end
    set(hObject,'UserData',tmp)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = res_vxl_thresh_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    ok  = 1;
    if length(tmp)~=1 || tmp<0
        warndlg('Cluster threshold: Insert a non-negative integer',...
            'Results: Invalid input','modal')
        set(hObject, 'String', '0')
        tmp = 0;
        ok  = 0;
    end
    tmp = round(tmp);
    set(hObject, 'String', num2str(tmp))
    set(hObject, 'UserData', tmp)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function res_PS_Callback(hObject, eventdata, handles)
    %struc_pref_enable(handles)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function res_slice_view_Callback(hObject, eventdata, handles)
    res_opt_button_enable(handles)
    %struc_pref_enable(handles)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function res_union_Callback(hObject, eventdata, handles)
    res_opt_button_enable(handles)
    %struc_pref_enable(handles)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function res_seq_Callback(hObject, eventdata, handles)
    res_opt_button_enable(handles)
    %struc_pref_enable(handles)
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function res_opt_button_Callback(hObject, eventdata, handles)
    handles.opt_safe.stat_res.slices.save    = get(handles.res_slice_view,'Value');
    handles.opt_safe.stat_res.group.save     = get(handles.res_union     ,'Value');
    handles.opt_safe.stat_res.event_seq.save = get(handles.res_seq       ,'Value');
    [stat_res button_press] = safe_results_opts(handles.opt_safe.stat_res);
    if ~strcmp(button_press,'ok')
        return
    end
    handles.opt_safe.stat_res.slices    = stat_res.slices;
    handles.opt_safe.stat_res.group     = stat_res.group;
    handles.opt_safe.stat_res.event_seq = stat_res.event_seq;
    guidata(handles.safe_main_window, handles)
end
%--------------------------------------------------------------------------



%==========================================================================
%                                BUTTONS
%==========================================================================

%--------------------------------------------------------------------------
function load_button_Callback(hObject, eventdata, handles)

    % Get file
    %---------
    [filename, tmp] = uigetfile({'*.mat','MAT-files (*.mat)'},'Load parameters...');
    if isequal(tmp,0) % "Cancel"
        return
    end
    
    % Load parameters ("opt_safe" variable)
    %--------------------------------------
    load(fullfile(tmp,filename))
    
    % Check and load parameters to GUI
    %---------------------------------
    ok = safe_load(handles,opt_safe);
    if ok
        handles.opt_safe = opt_safe;
        guidata(handles.safe_main_window,handles)
    end
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function save_button_Callback(hObject, eventdata, handles)
    
    % Get the parameters
    %-------------------
    handles  = safe_read_GUI(handles,'read');
    opt_safe = handles.opt_safe;
    
    % Check some parameters
    %----------------------
    tmp = check_param_part(opt_safe,handles);
%     if ~tmp % if there is some problem with the parameters, stop
%         return
%     end
    
    % Save parameters
    %----------------
    waitfor(findobj(0,'Name','Invalid input',...
        'Tag','Msgbox_Invalid input',...
        'Type','figure'))
    tmp = clock;
    tmp(6) = round(tmp(6));
    tmp = sprintf('%d_%.2d_%.2d_%.2d_%.2d_%.2d',tmp);
    uisave('opt_safe',sprintf('SAfE_%s',tmp));
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function run_button_Callback(hObject, eventdata, handles)
    
    % Basic check
    %------------
    if ~get(handles.preproc_chckbx,'Value') && ...
            ~get(handles.design_chckbx,'Value') && ...
            ~get(handles.res_chckbx,'Value')
        warndlg('You must select at least one output option.',...
            'Invalid input','modal')
        return
    end
    
    % Get the parameters
    %-------------------
    handles  = safe_read_GUI(handles,'read');
    opt_safe = handles.opt_safe;
    
    % Check some parameters
    %----------------------
    tmp = check_param_part(opt_safe,handles);
    if ~tmp % if there is some problem with the parameters, stop
        return
    end
    
    % Check if anterior commissure is at the origin
    %----------------------------------------------
    question = 'Yes';
    if get(handles.preproc_chckbx,'Value') && get(handles.preproc_normalize,'Value')
        question = questdlg(...
            {'Has the anterior commissure been marked in the images?',...
            'You need to do this with "SPM-Display" before running the analysis.',...
            '',...
            'If you choose Yes, the analysis will continue.'},...
            'Attention!','Yes', 'No', 'Yes');
    end
    if ~strcmp(question,'Yes') % 'No' or empty (closed the window)
        msgbox('Please save the parameters and reorient the images with SPM-Display.','The analysis will not continue')
        return
    end
    
    % Run
    %----
    set(handles.safe_main_window,'Pointer','watch')
    drawnow; pause(.05)
    safe_preproc_stat(opt_safe)
    set(handles.safe_main_window,'Pointer','arrow')
    drawnow; pause(.05)
    %set(handles.safe_main_window,'WindowStyle','normal')
    
    %uiresume(handles.safe_main_window); % will close the main figure
    msgbox('The process has been completed!','SAfE: finish','modal')
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function safe_main_window_CloseRequestFcn(hObject, eventdata, handles)
% Executes when user attempts to close safe_main_window
    question = questdlg('Do you really want to quit?','Exit confirmation','Yes','No','No');
    if strcmp(question,'Yes')
        if isequal(get(hObject, 'waitstatus'), 'waiting')
            % The GUI is still in UIWAIT, use UIRESUME and return
            uiresume(hObject);
        else
            % The GUI is no longer waiting, so destroy it now.
            delete(hObject); % close the figure
        end
    end
end
%--------------------------------------------------------------------------
