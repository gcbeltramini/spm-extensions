function [fighandle mov_pars] = fmri_plot_rp(fmri_files,cond_time,color)
% 
% Plot motion parameters from fMRI files
% 
% 
%USAGE
%-----
% [fighandle mov_pars] = fmri_plot_rp(fmri_files,cond_time,color)
% 
% 
%INPUT
%-----
% >> fmri_plot_rp
% >> fmri_plot_rp(fmri_files)
% >> fmri_plot_rp(fmri_files,cond_time,color)
% 
% - Empty: file(s) will be asked
% - FMRI_FILES: Nx1 cell array with the name of the unpreprocessed fMRI
%   files, where N is the number of sessions
% - COND_TIME: Cx1 cell array, where C is the number of conditions.
%   COND_TIME{cc} is an EVx2 matrix with the onsets in column 1 and the end
%   of the stimulus in column 2 (EV is the number of events)
%   N.B.: The units must in scans
% - COLOR: Cx3 double array with the color in RGB for condition C
%   - [1 1 1]/3 is the default gray
% 
% For all cases:
% - The order is important
% - Choose one scan per session
% 
% 
%OUTPUT
%------
% >> fmri_plot_rp(...)
% >> fighandle = fmri_plot_rp(...)
% >> [fighandle mov_pars] = fmri_plot_rp(...)
% 
% - FIGHANDLE: figure handle
% - MOV_PARS : Nx6 matrix with the motion parameters. Columns 1-3 with the
%   translational movement (in mm) in the x, y and z direction. Columns 4-6
%   with the rotational movement (in degrees) for pitch, roll and yaw.
% 
% 
% Based on function plot_parameters from spm_realign.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2012-May-22, 00:15 am
% 
% Changelog
%--------------------------------------------------------------------------
% 2012-Oct-10, 12:55 pm: added cond_time option
% 2013-Feb-23, 01:07 pm: different conditions are allowed and their color
%   can be different


% Input
%==========================================================================
if nargin==0
%     filt = ['^rp_*','.*\.txt$']; % rp file
%     file = spm_select([1 Inf],'any','Select realignment parameters files in the correct order',[],pwd, filt);
%     file = cellstr(file);
    file = spm_select([1 10],'image','Realigned files (order is important; 1 scan per session; files without prefix)',[],pwd);
    if isempty(file)
        fighandle = []; mov_pars = [];
        return
    end
    file = cellstr(file);
else
    if ~iscellstr(fmri_files)
        fmri_files = cellstr(fmri_files);
    end
    file = fmri_files(:);
    
    if nargin==2
        color = [1 1 1]/3;
    end
    
end


% Realignment parameters
%==========================================================================
% old: params = load(file);
params = [];
for ff=1:size(file,1) % loop for all chosen files
    
%     RP file (rp_...txt):
%     [fpath,fname,fext] = fileparts(file{ff,:});
%     fname = fullfile(fpath,[fname(4:end) '.nii']);
%     if exist(fname,'file')~=2
%         fname = fullfile(fpath,[fname(4:end) '.img']);
%     end
%     V = spm_vol(fname);
    
    tmp = file{ff,:};
    if nargin==0
        tmp = tmp(1:end-2); % exclude ",1"
    else
        if exist(tmp,'file')~=2
            fprintf('File %s not found\n',tmp)
            fighandle = []; mov_pars = [];
            return
        end
        if isempty(fileparts(tmp))
            tmp = fullfile(pwd,tmp);
        end
    end
    
    V = spm_vol(tmp);
    
    if ff==1
        mat0 = V(1).mat;
    end
    
    tmp = zeros(size(V,1),12);
    for vv=1:size(V,1)
        tmp(vv,:) = spm_imatrix(V(vv).mat/mat0);
    end
    
    params = [params ; tmp];
    
end


% Plot
%==========================================================================
fig = figure('Name','Motion parameters', 'Visible','on');
[fpath,fname,fext] = fileparts(file{1,:});

if nargout>0
    fighandle = fig;
end
if nargout>1
    mov_pars = params(:,1:6);
    mov_pars(:,4:6) = mov_pars(:,4:6)*180/pi;
end


% Translation
%--------------------------------------------------------------------------
subplot(2,1,1)
    plot(params(:,1:3),'LineWidth',2)
    if size(file,1)==1
        title(['Data from ' fname fext(1:4)],...
            'Interpreter','none',...
            'FontSize',14) % "fext" may be ".nii,1"
    else
        title(['Data from ' fname fext(1:4) ', etc.'],...
            'Interpreter','none',...
            'FontSize',14)
    end
    grid on, axis tight
    xlabel('Scan','FontSize',12)
    ylabel('Translation (mm)','FontSize',12)
    legend('x translation','y translation','z translation')
    
    % Plot events
    %------------
    if nargin>1
        plot_events(fig,cond_time,color)
    end


% Rotation
%--------------------------------------------------------------------------
subplot(2,1,2)
    plot(params(:,4:6)*180/pi,'LineWidth',2)
    grid on, axis tight
    xlabel('Scan','FontSize',12)
    ylabel('Rotation (degrees)','FontSize',12)
    legend('pitch','roll','yaw')
    
    % Plot events
    %------------
    if nargin>1
        plot_events(fig,cond_time,color)
    end


drawnow; pause(.1)