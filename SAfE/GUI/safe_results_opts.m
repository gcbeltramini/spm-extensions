function varargout = safe_results_opts(varargin)
% 
% SAFE_RESULTS_OPTS opens a GUI with the options for the results
% 
%__________________________________________________________________________
% Copyright (C) 2012 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2012-Dec-27, 07:15 pm

%==========================================================================
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @safe_results_opts_OpeningFcn, ...
                   'gui_OutputFcn',  @safe_results_opts_OutputFcn, ...
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

function safe_results_opts_OpeningFcn(hObject, eventdata, handles, varargin)
% Executes just before safe_results_opts is made visible

dontOpen = false;

if length(varargin)==1
    handles.stat_res = varargin{1};
    safe_load(handles)
else
    dontOpen = true;
    disp('Input missing')
end

% Update handles structure
guidata(hObject, handles);

if ~dontOpen
    disable_panels(handles) % adjust what options are enabled
    uiwait(hObject);
    % WindowStyle property for the figure must be 'modal' to block all
    % other windows in Matlab
end

end

%--------------------------------------------------------------------------

function varargout = safe_results_opts_OutputFcn(hObject, eventdata, handles)
    varargout{1} = [];
    try
        varargout{2} = get(handles.OK_button,'UserData');
        if strcmp(varargout{2},'ok')
            varargout{1} = handles.stat_res;
        end
    catch
        varargout{2} = 'cancel';
    end
    delete(hObject)
end
%==========================================================================



%==========================================================================
%                           AUXILIARY FUNCTIONS
%==========================================================================

%--------------------------------------------------------------------------
function safe_load(handles)
% Load parameters into GUI
    
    stat_res = handles.stat_res;
    
    % 1) Slice view
    %--------------
    if isequal(stat_res.slices.step,0)
        if isequal(stat_res.slices.number,0) % coordinates of the slices
            set(handles.slice_view_slices_menu,'Value'   ,3)
            set(handles.slice_view_slices_edit,'String'  ,num2str(stat_res.slices.slices))
            set(handles.slice_view_slices_edit,'UserData',stat_res.slices.slices)
            set(handles.slice_view_slices_edit,'String',...
                regexprep(get(handles.slice_view_slices_edit,'String'),'\s+',' ')) % remove repeated spaces
            gui_chg_tooltip(handles.slice_view_slices_edit,[],'') % change tooltip string
        else % number of slices
            set(handles.slice_view_slices_menu,'Value'   ,2)
            set(handles.slice_view_slices_edit,'String'  ,num2str(stat_res.slices.number))
            set(handles.slice_view_slices_edit,'UserData',stat_res.slices.number)
        end
    else % step (mm)
        set(handles.slice_view_slices_menu,'Value'   ,1)
        set(handles.slice_view_slices_edit,'String'  ,num2str(stat_res.slices.step))
        set(handles.slice_view_slices_edit,'UserData',stat_res.slices.step)
    end
    
    %stat_res.slices.orient = sort(stat_res.slices.orient);
    set(handles.slice_view_slice_orient_axi,'Value',0)
    set(handles.slice_view_slice_orient_cor,'Value',0)
    set(handles.slice_view_slice_orient_sag,'Value',0)
    for oo=1:length(stat_res.slices.orient)
        switch stat_res.slices.orient{oo}
            case 'axial'
                set(handles.slice_view_slice_orient_axi,'Value',1)
            case 'coronal'
                set(handles.slice_view_slice_orient_cor,'Value',1)
            case 'sagittal'
                set(handles.slice_view_slice_orient_sag,'Value',1)
        end
    end
    
    % 2) Union of delays
    %-------------------
    set(handles.union_delays_per_group_edit,'String'  ,num2str(stat_res.group.number))
    set(handles.union_delays_per_group_edit,'UserData',stat_res.group.number)
    set(handles.union_adj_scale,'Value',stat_res.group.scale.adj)
    if stat_res.group.scale.adj
        set(handles.union_adj_scale_max    ,'Enable','on')
        set(handles.union_adj_scale_specify,'Enable','on')
    else
        set(handles.union_adj_scale_max      ,'Enable','off')
        set(handles.union_adj_scale_specify  ,'Enable','off')
        set(handles.union_adj_scale_F_txt    ,'Enable','off')
        set(handles.union_adj_scale_F_edit   ,'Enable','off')
        set(handles.union_adj_scale_Tpos_txt ,'Enable','off')
        set(handles.union_adj_scale_Tpos_edit,'Enable','off')
        set(handles.union_adj_scale_Tneg_txt ,'Enable','off')
        set(handles.union_adj_scale_Tneg_edit,'Enable','off')
    end
    if strcmp(stat_res.group.scale.ftest,'max') || ...
            strcmp(stat_res.group.scale.ttest,'max')
        set(handles.union_adj_scale_max    ,'Value',1)
        set(handles.union_adj_scale_specify,'Value',0)
    else
        if stat_res.group.scale.adj
            set(handles.union_adj_scale_F_txt    ,'Enable','on')
            set(handles.union_adj_scale_F_edit   ,'Enable','on')
            set(handles.union_adj_scale_Tpos_txt ,'Enable','on')
            set(handles.union_adj_scale_Tpos_edit,'Enable','on')
            set(handles.union_adj_scale_Tneg_txt ,'Enable','on')
            set(handles.union_adj_scale_Tneg_edit,'Enable','on')
        end
        set(handles.union_adj_scale_max      ,'Value',0)
        set(handles.union_adj_scale_specify  ,'Value',1)
        set(handles.union_adj_scale_F_edit   ,'String'  ,num2str(stat_res.group.scale.ftest(2)))
        set(handles.union_adj_scale_Tpos_edit,'String'  ,num2str(stat_res.group.scale.ttest(2)))
        set(handles.union_adj_scale_Tneg_edit,'String'  ,num2str(stat_res.group.scale.ttest(4)))
        set(handles.union_adj_scale_F_edit   ,'UserData',stat_res.group.scale.ftest(2))
        set(handles.union_adj_scale_Tpos_edit,'UserData',stat_res.group.scale.ttest(2))
        set(handles.union_adj_scale_Tneg_edit,'UserData',stat_res.group.scale.ttest(4))
    end
    
    % 3) Sequence of delays
    %----------------------
    if strcmp(stat_res.event_seq.scale.ftest,'max') || ...
            strcmp(stat_res.event_seq.scale.posbold,'max') || ...
            strcmp(stat_res.event_seq.scale.negbold,'max')
        set(handles.seq_adj_scale_max      ,'Value',1)
        set(handles.seq_adj_scale_specify  ,'Value',0)
        set(handles.seq_adj_scale_F_txt    ,'Enable','off')
        set(handles.seq_adj_scale_F_edit   ,'Enable','off')
        set(handles.seq_adj_scale_Tpos_txt ,'Enable','off')
        set(handles.seq_adj_scale_Tpos_edit,'Enable','off')
        set(handles.seq_adj_scale_Tneg_txt ,'Enable','off')
        set(handles.seq_adj_scale_Tneg_edit,'Enable','off')
    else
        set(handles.seq_adj_scale_max      ,'Value',0)
        set(handles.seq_adj_scale_specify  ,'Value',1)
        set(handles.seq_adj_scale_F_txt    ,'Enable','on')
        set(handles.seq_adj_scale_F_edit   ,'Enable','on')
        set(handles.seq_adj_scale_Tpos_txt ,'Enable','on')
        set(handles.seq_adj_scale_Tpos_edit,'Enable','on')
        set(handles.seq_adj_scale_Tneg_txt ,'Enable','on')
        set(handles.seq_adj_scale_Tneg_edit,'Enable','on')
        set(handles.seq_adj_scale_F_edit   ,'String',num2str(stat_res.event_seq.scale.ftest(2)))
        set(handles.seq_adj_scale_Tpos_edit,'String',num2str(stat_res.event_seq.scale.posbold(2)))
        set(handles.seq_adj_scale_Tneg_edit,'String',num2str(stat_res.event_seq.scale.negbold(2)))
        set(handles.seq_adj_scale_F_edit   ,'UserData',stat_res.event_seq.scale.ftest(2))
        set(handles.seq_adj_scale_Tpos_edit,'UserData',stat_res.event_seq.scale.posbold(2))
        set(handles.seq_adj_scale_Tneg_edit,'UserData',stat_res.event_seq.scale.negbold(2))
    end
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function disable_panels(handles)
    
    % 1) Slice view
    %--------------
    if ~handles.stat_res.slices.save && ~handles.stat_res.group.save && ...
            ~handles.stat_res.event_seq.save
        set(handles.slice_view_panel           ,'BorderType','beveledin')
        set(handles.slice_view_panel           ,'FontWeight','normal')
        set(handles.slice_view_slices_menu     ,'Enable','off')
        set(handles.slice_view_slices_edit     ,'Enable','off')
        set(handles.slice_view_slice_orient    ,'Enable','off')
        set(handles.slice_view_slice_orient_axi,'Enable','off')
        set(handles.slice_view_slice_orient_cor,'Enable','off')
        set(handles.slice_view_slice_orient_sag,'Enable','off')
    end
    
    % 2) Union of delays
    %-------------------
    if ~handles.stat_res.group.save
        set(handles.union_panel                ,'BorderType','beveledin')
        set(handles.union_panel                ,'FontWeight','normal')
        set(handles.union_delays_per_group_txt ,'Enable','off')
        set(handles.union_delays_per_group_edit,'Enable','off')
        set(handles.union_adj_scale            ,'Enable','off')
        set(handles.union_adj_scale_max        ,'Enable','off')
        set(handles.union_adj_scale_specify    ,'Enable','off')
        set(handles.union_adj_scale_F_txt      ,'Enable','off')
        set(handles.union_adj_scale_F_edit     ,'Enable','off')
        set(handles.union_adj_scale_Tpos_txt   ,'Enable','off')
        set(handles.union_adj_scale_Tpos_edit  ,'Enable','off')
        set(handles.union_adj_scale_Tneg_txt   ,'Enable','off')
        set(handles.union_adj_scale_Tneg_edit  ,'Enable','off')
    end
    
    % 3) Sequence of delays
    %----------------------
    if ~handles.stat_res.event_seq.save
        set(handles.seq_panel              ,'BorderType','beveledin')
        set(handles.seq_panel              ,'FontWeight','normal')
        set(handles.seq_adj_scale          ,'Enable','off')
        set(handles.seq_adj_scale_max      ,'Enable','off')
        set(handles.seq_adj_scale_specify  ,'Enable','off')
        set(handles.seq_adj_scale_F_txt    ,'Enable','off')
        set(handles.seq_adj_scale_F_edit   ,'Enable','off')
        set(handles.seq_adj_scale_Tpos_txt ,'Enable','off')
        set(handles.seq_adj_scale_Tpos_edit,'Enable','off')
        set(handles.seq_adj_scale_Tneg_txt ,'Enable','off')
        set(handles.seq_adj_scale_Tneg_edit,'Enable','off')
    end
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function [stat_res warn] = safe_read_GUI(handles)
% Read options
    
    zero_thresh = 10^(-10);
    
    warn = 0;
    
    % 1) Slice view
    %--------------
    stat_res.slices.save    = handles.stat_res.slices.save;
    stat_res.slices.img_ext = 'png';
    stat_res.slices.step    = 3;
    stat_res.slices.number  = 36;
    stat_res.slices.slices  = -50:4:85;
    switch get(handles.slice_view_slices_menu,'Value')
        case 1 % Step (mm)
            stat_res.slices.step = get(handles.slice_view_slices_edit,'UserData');
        case 2 % Number of slices
            stat_res.slices.step   = 0;
            stat_res.slices.number = get(handles.slice_view_slices_edit,'UserData');
        case 3 % Coordinates of the slices
            stat_res.slices.step   = 0;
            stat_res.slices.number = 0;
            stat_res.slices.slices = get(handles.slice_view_slices_edit,'UserData');
    end
    
    % Check slice orientation
    stat_res.slices.orient = [];
    if get(handles.slice_view_slice_orient_axi,'Value')
        stat_res.slices.orient = ['axial   '];
    end
    if get(handles.slice_view_slice_orient_cor,'Value')
        stat_res.slices.orient = [stat_res.slices.orient; 'coronal '];
    end
    if get(handles.slice_view_slice_orient_sag,'Value')
        stat_res.slices.orient = [stat_res.slices.orient; 'sagittal'];
    end
    if isempty(stat_res.slices.orient)
        warndlg('Slice view: Orientation is missing','Invalid input')
        warn = 1;
        stat_res.slices.orient = '';
    end
    stat_res.slices.orient = cellstr(stat_res.slices.orient);
    
    % 2) Union of delays
    %-------------------
    stat_res.group.save   = handles.stat_res.group.save;
    stat_res.group.number = get(handles.union_delays_per_group_edit,'UserData');
    if stat_res.group.save
        if get(handles.union_adj_scale,'Value') && get(handles.union_adj_scale_specify,'Value')
            % check if F & T are not empty
            ok = F_T_allowed(handles.union_adj_scale_F_edit);
            if ~ok, warn = 1; return, end
            ok = F_T_allowed(handles.union_adj_scale_Tpos_edit);
            if ~ok, warn = 1; return, end
            ok = F_T_allowed(handles.union_adj_scale_Tneg_edit);
            if ~ok, warn = 1; return, end
        end
    end
    stat_res.group.scale.adj = get(handles.union_adj_scale,'Value');
    if get(handles.union_adj_scale_max,'Value')
        stat_res.group.scale.ftest = 'max';
        stat_res.group.scale.ttest = 'max';
    else
        stat_res.group.scale.ftest = [zero_thresh get(handles.union_adj_scale_F_edit   ,'UserData')];
        stat_res.group.scale.ttest = [zero_thresh get(handles.union_adj_scale_Tpos_edit,'UserData') ...
            zero_thresh get(handles.union_adj_scale_Tneg_edit,'UserData')];
    end
    
    % 3) Sequence of delays
    %----------------------
    stat_res.event_seq.save = handles.stat_res.event_seq.save;
    if stat_res.event_seq.save
        if get(handles.seq_adj_scale_specify,'Value')
            % check if F & T are not empty
            ok = F_T_allowed(handles.seq_adj_scale_F_edit);
            if ~ok, warn = 1; return, end
            ok = F_T_allowed(handles.seq_adj_scale_Tpos_edit);
            if ~ok, warn = 1; return, end
            ok = F_T_allowed(handles.seq_adj_scale_Tneg_edit);
            if ~ok, warn = 1; return, end
        end
    end
    if get(handles.seq_adj_scale_max,'Value')
        stat_res.event_seq.scale.ftest   = 'max';
        stat_res.event_seq.scale.posbold = 'max';
        stat_res.event_seq.scale.negbold = 'max';
    else
        stat_res.event_seq.scale.ftest   = [zero_thresh get(handles.seq_adj_scale_F_edit   ,'UserData')];
        stat_res.event_seq.scale.posbold = [zero_thresh get(handles.seq_adj_scale_Tpos_edit,'UserData')];
        stat_res.event_seq.scale.negbold = [zero_thresh get(handles.seq_adj_scale_Tneg_edit,'UserData')];
    end
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function adjust_scale(handles,panel,object,state)
    tmp = sprintf('set(handles.%s_adj_scale_',panel);
    if ~isempty(strfind(get(object,'Tag'),'_specify')) % specify values
        tmp2 = [1 0];
    else % maximum value
        tmp2 = [0 1];
    end
    if get(object,'Value')
        eval(sprintf('%sspecify  ,''Value'' ,%d)'    ,tmp,tmp2(1)))
        eval(sprintf('%smax      ,''Value'' ,%d)'    ,tmp,tmp2(2)))
        eval(sprintf('%sF_txt    ,''Enable'',''%s'')',tmp,state))
        eval(sprintf('%sF_edit   ,''Enable'',''%s'')',tmp,state))
        eval(sprintf('%sTpos_txt ,''Enable'',''%s'')',tmp,state))
        eval(sprintf('%sTpos_edit,''Enable'',''%s'')',tmp,state))
        eval(sprintf('%sTneg_txt ,''Enable'',''%s'')',tmp,state))
        eval(sprintf('%sTneg_edit,''Enable'',''%s'')',tmp,state))
    else
        set(object,'Value',1) % cannot disable option by clicking on it
    end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function ok = F_T_allowed(hObject)
    ok = 1;
    tmp = get(hObject,'Tag');
    switch tmp(1:3)
        case 'uni'
            warntxt = 'Union of delays';
        case 'seq'
            warntxt = 'Sequence of delays';
    end
    tmp = str2num(get(hObject,'String'));
    if length(tmp)~=1 || tmp<=0
        warndlg(sprintf('%s: Insert a positive number',warntxt),'Invalid input')
        set(hObject,'String'  ,'')
        set(hObject,'UserData',[])
        ok = 0;
    else
        set(hObject,'UserData',tmp)
    end
end
%--------------------------------------------------------------------------



%==========================================================================
%                              SLICE VIEW
%==========================================================================

%--------------------------------------------------------------------------
function slice_view_slices_menu_Callback(hObject, eventdata, handles)
    switch get(hObject,'Value')
        case 1 % Step (mm)
            set(handles.slice_view_slices_edit,'String','3')
            set(handles.slice_view_slices_edit,'UserData',3)
        case 2 % Number of slices
            set(handles.slice_view_slices_edit,'String','36')
            set(handles.slice_view_slices_edit,'UserData',36)
        case 3 % Coordinates of the slices
            set(handles.slice_view_slices_edit,'String','-50:4:85')
            set(handles.slice_view_slices_edit,'UserData',-50:4:85)
    end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function slice_view_slices_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    switch get(handles.slice_view_slices_menu,'Value')
        case 1 % Step (mm)
            if length(tmp)~=1 || tmp<=0
                warndlg('Slice view: Insert a positive number','Invalid input')
                set(hObject,'String','3')
                set(hObject,'UserData',3)
                return
            end
        case 2 % Number of slices
            if length(tmp)~=1 || floor(tmp)~=ceil(tmp) || tmp<1
                warndlg('Slice view: Insert a positive integer number','Invalid input')
                set(hObject,'String','36')
                set(hObject,'UserData',36)
                return
            end
        case 3 % Coordinates of the slices
            if isempty(tmp) || any(floor(tmp)~=ceil(tmp))
                warndlg('Slice view: Insert an integer matrix','Invalid input')
                set(hObject,'String','-50:4:85')
                set(hObject,'UserData',-50:4:85)
                return
            end
            % Change tooltip string
            gui_chg_tooltip(hObject,[],'')
    end
    set(hObject,'UserData',tmp)
end
%--------------------------------------------------------------------------



%==========================================================================
%                            UNION OF DELAYS
%==========================================================================

%--------------------------------------------------------------------------
function union_delays_per_group_edit_Callback(hObject, eventdata, handles)
    tmp = str2num(get(hObject,'String'));
    if length(tmp)~=1 || floor(tmp)~=ceil(tmp) || tmp<1
        warndlg('Union of delays: Insert a positive integer number','Invalid input')
        set(hObject,'String','3')
        tmp = 3;
    end
    set(hObject,'UserData',tmp)
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function union_adj_scale_Callback(hObject, eventdata, handles)
    if get(hObject,'Value')
        set(handles.union_adj_scale_max    ,'Enable','on')
        set(handles.union_adj_scale_specify,'Enable','on')
        if get(handles.union_adj_scale_specify,'Value')
            set(handles.union_adj_scale_F_txt    ,'Enable','on')
            set(handles.union_adj_scale_F_edit   ,'Enable','on')
            set(handles.union_adj_scale_Tpos_txt ,'Enable','on')
            set(handles.union_adj_scale_Tpos_edit,'Enable','on')
            set(handles.union_adj_scale_Tneg_txt ,'Enable','on')
            set(handles.union_adj_scale_Tneg_edit,'Enable','on')
        end
    else
        set(handles.union_adj_scale_max      ,'Enable','off')
        set(handles.union_adj_scale_specify  ,'Enable','off')
        set(handles.union_adj_scale_F_txt    ,'Enable','off')
        set(handles.union_adj_scale_F_edit   ,'Enable','off')
        set(handles.union_adj_scale_Tpos_txt ,'Enable','off')
        set(handles.union_adj_scale_Tpos_edit,'Enable','off')
        set(handles.union_adj_scale_Tneg_txt ,'Enable','off')
        set(handles.union_adj_scale_Tneg_edit,'Enable','off')
    end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function union_adj_scale_max_Callback(hObject, eventdata, handles)
    adjust_scale(handles,'union',hObject,'off')
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function union_adj_scale_specify_Callback(hObject, eventdata, handles)
    adjust_scale(handles,'union',hObject,'on')
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function union_adj_scale_F_edit_Callback(hObject, eventdata, handles)
    F_T_allowed(hObject);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function union_adj_scale_Tpos_edit_Callback(hObject, eventdata, handles)
    F_T_allowed(hObject);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function union_adj_scale_Tneg_edit_Callback(hObject, eventdata, handles)
    F_T_allowed(hObject);
end
%--------------------------------------------------------------------------



%==========================================================================
%                          SEQUENCE OF DELAYS
%==========================================================================

%--------------------------------------------------------------------------
function seq_adj_scale_max_Callback(hObject, eventdata, handles)
    adjust_scale(handles,'seq',hObject,'off')
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function seq_adj_scale_specify_Callback(hObject, eventdata, handles)
    adjust_scale(handles,'seq',hObject,'on')
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function seq_adj_scale_F_edit_Callback(hObject, eventdata, handles)
    F_T_allowed(hObject);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function seq_adj_scale_Tpos_edit_Callback(hObject, eventdata, handles)
    F_T_allowed(hObject);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function seq_adj_scale_Tneg_edit_Callback(hObject, eventdata, handles)
    F_T_allowed(hObject);
end
%--------------------------------------------------------------------------


%==========================================================================
%                                BUTTONS
%==========================================================================

%--------------------------------------------------------------------------
function OK_button_Callback(hObject, eventdata, handles)
    set(handles.OK_button,'UserData','ok') % remember which button was pressed
    [stat_res warn] = safe_read_GUI(handles);
    if ~warn
        handles.stat_res = stat_res;
        guidata(handles.safe_res_opt_window,handles)
        uiresume(handles.safe_res_opt_window);
    end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function cancel_button_Callback(hObject, eventdata, handles)
%     question = questdlg('Do you really want to quit? The options will be lost',...
%         'Exit confirmation','Yes','No','No');
%     if strcmp(question,'Yes')
        set(handles.OK_button,'UserData','cancel') % remember which button was pressed
        uiresume(handles.safe_res_opt_window);
%     end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function safe_res_opt_window_CloseRequestFcn(hObject, eventdata, handles)
%     question = questdlg('Do you really want to quit? The options will be lost',...
%         'Exit confirmation','Yes','No','No');
%     if strcmp(question,'Yes')
        set(handles.OK_button,'UserData','close') % remember which button was pressed
        uiresume(hObject);
%     end
end
%--------------------------------------------------------------------------
