function fmri_mip_color(spmmat, con, P, adj, ET, color)
% 
% fmri_mip_color(spmmat, con, P, adj, ET, color)
% 
% Create figure with the maximum intensity projection (MIP) on a glass
% brain in three orthogonal planes using colors.
% 
% 
%USAGE
%-----
% fmri_mip_color: open GUI
% fmri_mip_color(spmmat)
% fmri_mip_color(spmmat, con)
% fmri_mip_color(spmmat, con, P)
% fmri_mip_color(spmmat, con, P, adj)
% fmri_mip_color(spmmat, con, P, adj, ET)
% fmri_mip_color(spmmat, con, P, adj, ET, color)
% 
% If any of the input variables is empty or missing, the default value will
% be used.
% 
% 
%INPUT
%-----
% - SPMMAT: full path for the file SPM.mat [default: SPM.mat in the current folder]
%     SPMMAT can be cell array of file paths, but only one figure will be
%     created.
% - CON   : contrast number (1, 2, 3, ...)       [default: 1]
%     CON can be a vector. Only one figure will be created with all the
%     contrasts applied to all SPM.mat files.
% - P     : P-value (0 < P < 1)                  [default: 0.001]
% - ADJ   : 'FWE' or 'none' (P-value adjustment) [default: 'none']
% - ET    : extent threshold in voxels (ET >= 0) [default: 0]
% - COLOR : 1 (red), 2 (blue), 3 (green)         [default: 1]
%     COLOR must have the same dimensions as CON or one value (all
%     conditions will have the same color).
% 
% 
%OUTPUT
%------
% - Figure with the MIP of the chosen contrast and color
%   - Figure tag  = 'fig_mip_color'
%   - Figure name = 'MIP: P=0.xxxx, adj=xxx, ET=x'
% 
% 
%EXAMPLES
%--------
% fmri_mip_color('D:\MyData\Sbj14\fMRI\SPM.mat', 2, .05, 'FWE', 0, 1)
% fmri_mip_color('D:\MyData\Sbj14\fMRI\SPM.mat', [2 4], .05, 'FWE', 10, [1 2])
% fmri_mip_color({'D:\MyData\Sbj01\fMRI\SPM.mat', 'D:\MyData\Sbj02\fMRI\SPM.mat'}, ...
%     [1 2], [], [], [], 3)
% 
%__________________________________________________________________________
% Copyright (C) 2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2014-Jul-01


max_time   = 20; % maximum tolerance (in seconds) to calculate the height threshold (see "spm_results_ui.m")
cmaplevels = 64; % number of levels of the colormap


% Input
%==========================================================================

% GUI
%----
if nargin==0
    fmri_mip_color_gui;
    return
end

% Color
%------
if nargin<6 || isempty(color)
    color = 1:3;
elseif ~isequal(intersect(color, [1 2 3]), color)
    error('"color" can be 1, 2 or 3 (or a vector with these numbers)')
end

% Extent threshold
%-----------------
if nargin<5 || isempty(ET)
    ET = 0;
elseif ET<0
    error('"ET" must be a non-negative number')
end
ET = round(ET);

% P-value adjustment
%-------------------
if nargin<4 || isempty(adj)
    adj = 'none';
elseif all(~strcmpi(adj, {'none', 'FWE'}))
    error('"adj" must be ''none'' or ''FWE''')
elseif strcmpi(adj, 'none')
    adj = 'none';
else
    adj = 'FWE';
end

% P-value
%--------
if nargin<3 || isempty(P)
    P = 0.001;
elseif P<=0
    error('"P" must be a positive number')
end

% Contrast
%---------
if nargin<2 || isempty(con)
    con = 1;
elseif any(floor(con)~=ceil(con)) || any(con<1)
    error('"con" must be a positive integer number (or a vector with up to 3 values)')
end
if length(con)<=length(color)
    color = color(1:length(con));
else
    color = repmat(color(1), size(con));
end

% SPM.mat file
%-------------
if isempty(spmmat)
    spmmat = {fullfile(pwd, 'SPM.mat')};
else
    spmmat = cellstr(spmmat);
    spmmat = spmmat(:);
    spmmat(strcmp(spmmat, 'SPM.mat')) = {fullfile(pwd, 'SPM.mat')};
end


% Figure visibility: if SPM is open, do nothing
%------------------
if isempty(spm_figure('FindWin', 'Menu')) % SPM is not open
    visible = 'off';
else
    visible = 'on';
end


% Create image
%==========================================================================

spm('defaults', 'FMRI'); % set the default global variables
Ncond = 0; % number of conditions
curr_dir = pwd;

for ss=1:length(spmmat) % loop for the SPM.mat files
    
    fprintf('\n%d) File:\n  %s\n', ss, spmmat{ss})
    
    % Load the SPM.mat file
    %----------------------
    if exist(spmmat{ss}, 'file')~=2
        fprintf('  SPM.mat not found\n')
        continue
    end
    load(spmmat{ss})
    
    con(con>length(SPM.xCon)) = [];
    if isempty(con)
        fprintf('  The maximum number allowed for "con" is %d\n', length(SPM.xCon))
        continue
    end
    
    % Create variable xSPM
    %---------------------
    xSPM.swd   = fileparts(spmmat{ss});  % SPM.mat folder
    xSPM.u     = P;                      % Threshold
    xSPM.Im    = [];                     % Masking
    xSPM.Ex    = [];                     % exclusive or inclusive
    xSPM.k     = ET;                     % Extent (voxels)
    xSPM.units = {'mm','mm','mm'};       % Data type: Volumetric (2D/3D)
    
    
    for cc=1:length(con) % loop for the contrasts
        
        xSPM.title     = SPM.xCon(con(cc)).name; % Results title
        xSPM.Ic        = con(cc);                % Contrast(s)
        xSPM.thresDesc = adj; % Threshold type: 'none', 'FWE' (must be inside the for-loop)
        
        % Test time to calculate "height threshold" in spm_results_ui (from spm_uc_RF.m)
        %------------------------------------------------------------------
        stop = test_time(SPM, xSPM, max_time);
        if stop
            fprintf('%.1f s have elapsed. Condition skipped.\n', max_time)
            continue
        end
        
        fprintf('\n  %d.%d) Condition: %s\n', ss, cc, SPM.xCon(con(cc)).name)
        
        Ncond = Ncond + 1;
        if Ncond==1
            fig_title = xSPM.title; % figure title
        end
        
        fprintf('       Color = ')
        switch color(cc)
            case 1
                fprintf('red\n')
            case 2
                fprintf('blue\n')
            case 3
                fprintf('green\n')
        end
        
        [hReg, xSPM, SPM] = spm_results_ui('Setup', xSPM); % show results
        if strcmp(visible, 'off')
            set(spm_figure('FindWin','Graphics'), 'Visible', visible) % big window
            set(spm_figure('FindWin','Interactive'), 'Visible', visible) % small window
        end
        
        % Create MIP image
        %-----------------
        [mip, d, grid_mask] = spm_mip_mod(xSPM.Z, xSPM.XYZmm, SPM.xVol.M);
        if exist('img', 'var')~=1
            img = zeros([size(d), 3]);
        end
        img(:, :, color(cc)) = max(img(:, :, color(cc)), 1-d);
        
    end
    
end

cd(curr_dir)

if exist('img', 'var')~=1 % "img" was not created
    return
end

clear SPM xSPM % clean-up


% Create figure
%==========================================================================

figure('Name', sprintf('MIP: P=%.4f, adj=%s, ET=%d', P, adj, ET), ...
    'Color', [1 1 1], ...
    'Colormap', repmat(linspace(0, 1, cmaplevels)', 1, 3), ...
    'NumberTitle', 'off', ...
    'Tag', 'fig_mip_color')

% Colored image
%--------------
img(:, :, [2 3]) = img(:, :, [3 2]); % 1: red, 2: blue, 3: green, but RGB = [1 2 3]
% This order was chosen because most of the time red and blue will be used
BG = img(:,:,1)==0 & img(:,:,2)==0 & img(:,:,3)==0; % only activation map = 0
image(img + repmat(grid_mask.*BG,[1,1,3]))

% Axis and title
%---------------
axis tight; axis off;
if Ncond>1
    fig_title = sprintf('%s, ... (%d conditions)', fig_title, Ncond);
end
title(fig_title)


% Other possibilities for the images
%==========================================================================

% Only grid and mask (gray scale)
%-------------------
%image(repmat(grid_mask, [1,1,3]));
%image(grid_mask*cmaplevels); % better solution

% Only the functional map (one map)
%------------------------
%image(repmat(d, [1,1,3])); % gray scale
%image(d*cmaplevels); % gray scale (better solution)
%BG = ones(size(d)); BG(d~=1) = 0; % background (1: background; 0: brain)
%image(cat(3, d, BG, BG)); % red
%image(cat(3, BG, d, BG)); % green
%image(cat(3, BG, BG, d)); % blue

% Grid, mask and functional map (one map)
%------------------------------
%image(repmat(mip, [1,1,3])) % gray scale
%image(mip*cmaplevels) % gray scale (better solution)

% Binary images (one map)
%--------------
%BG = mip==d; % background(1: background and activation map; 0: grid and mask)
%image(BG*cmaplevels) % 0: black; 1: white
%image(cmaplevels*BG.*mip)
% With colors:
% BG = mip==grid_mask; % background(1: background, lines and mask; 0: activation map)
% if color==1 % red
%     image(cat(3, mip, BG.*mip, BG.*mip))
% elseif color==2 % blue
%     image(cat(3, BG.*mip, BG.*mip, mip))
% elseif color==3 % green
%     image(cat(3, BG.*mip, mip, BG.*mip))
% end


% Close SPM windows
%------------------
if strcmp(visible, 'off')
    tmp = spm_figure('FindWin', 'Graphics'); % findobj('Tag', 'Graphics')
    close(tmp);
    tmp = spm_figure('FindWin', 'Interactive'); % findobj('Tag', 'Interactive')
    close(tmp);
end


% =========================================================================
%                           AUXILIARY FUNCTIONS
%==========================================================================

function stop = test_time(SPM, xSPM, max_time)
% When calculating the height threshold (at spm_uc_RF.m), it can take a
% long time because of problems in the design. Since I couldn't find a
% better way, I copied the commands from spm_uc_RF and inserted a time
% counter to check the elapsed time in the problematic loop.

% Based on spm_getSPM.m and spm_uc_RF.m

a    = 0.05;
stop = 0;

Ic   = xSPM.Ic;
xCon = SPM.xCon;
STAT = xCon(Ic(1)).STAT;

R = SPM.xVol.R;

xX = SPM.xX;
df = [xCon(Ic(1)).eidf xX.erdf];

try
    n = xSPM.n;
catch
    n = 1;
end

u  = spm_u((a/max(R))^(1/n),df,STAT);
du = 1e-6;

% Approximate estimate using E{m}
%--------------------------------
d = 1;
tic
while abs(d) > 1e-6
    if toc>max_time, stop = 1; break, end
    [P P p] = spm_P_RF(1,0,u,df,STAT,R,n);
    [P P q] = spm_P_RF(1,0,u + du,df,STAT,R,n);
    d       = (a - p)/((q - p)/du);
    u       = u + d;
end
