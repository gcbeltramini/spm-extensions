function [mn_mx, slvr_obj] = safe_slover_create(opt_slvr_create)
% 
% SAFE_SLOVER_CREATE creates a slover object.
% 
% 
%USAGE
%-----
% [mn_mx, slvr_obj] = safe_slover_create(opt_slvr_create)
% 
% 
%INPUT
%-----
% - OPT_SLVR_CREATE: structure with the following fields:
%   - IMGS : Nx1 cell array containing the full path of the images (e.g.,
%   structural, positive and negative BOLD)
%   - TYPE : Nx1 cell array of the image type corresponding to IMGS. The
%   available options are: 'Structural', 'Truecolour', 'Blobs', 'Negative
%   blobs', 'Contours'. One of the cells must be 'Structural'.
%   - CMAP : Nx1 cell array of strings (structural image is necessarily
%   grayscale, so leave its corresponding cell empty). Possibilities are:
%   1) Filename of *.mat or *.lut file
%   a) *.lut files from ...\spm8\@slover\private:
%      blackbdy.lut, cardiac.lut, flow.lut    , ge_color.lut, gold.lut    ,
%      hotiron.lut ,  nih.lut   , nih_fire.lut, nih_ice.lut , rainramp.lut,
%      spectrum.lut, x_hot.lut  , x_rain.lut
%   2) Colormaps from Matlab:
%      'jet','hsv','hot','cool','spring','summer','autumn','winter','gray',
%      'bone','copper','pink','lines','colorcube','flag','prism','white'
%   3) Colour name:
%      'red','green','blue','cyan','magenta','yellow','black','white'
%   4) Matrix name
%      e.g., pure_red = [0 0 0; .2 0 0; .4 0 0; .6 0 0; .8 0 0; 1 0 0]
%   Default:
%     - 'Truecolour'    : 'flow.lut'
%     - 'Blobs'         : 'hot'
%     - 'Negative blobs': 'winter'
%     - 'Contours'      : 'white'
%   - ORIENT: slice orientation: 'axial', 'coronal', 'sagittal'
%   - SLC_STEP: spacing between slices (mm) (any positive number)
%   - SLC_NMBR: number of slices (to use this option, SLC_STEP must be 0)
%   - SLC_SLICES: coordinates of the slices (to use this option, SLC_STEP
%     and SLC_NMBR must be 0)
% 
% 
%OUTPUT
%------
% - MN_MX   : Nx2 matrix, with the minimum and maximum values of the N
%             images
% - SLVR_OBJ: slice overlay (slover) object
% 
% 
% Based on slover.m (slice overlay), pr_basic_ui.m
% 
% 
% See also SAFE_SLOVER_SAVE
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme C. Beltramini - 2013-Mar-16

% obj.img(X).range:
% Image value range for colormap: [0 149.177]
% (1st value in the colormap -- last value in the colormap)
% Example: Colors symmetrical around 0 (half the colormap is for values
% greater than 0, and half for less than 0); image range = [-4.81 7.54]
% => range = [-7.54 7.54]


zero_thresh = 10^(-10);


% Read options
%==========================================================================
imgs       = opt_slvr_create.imgs(:);
itype      = opt_slvr_create.type;
cmap       = opt_slvr_create.cmap;
orient     = opt_slvr_create.orient;
slc_step   = opt_slvr_create.slc_step;
slc_nmbr   = opt_slvr_create.slc_nmbr;
slc_slices = opt_slvr_create.slc_slices;


% Create slover object
%==========================================================================
obj      = slover;
cscale   = [];
obj.cbar = [];
mn_mx    = zeros(size(imgs,1),2);
for ii = 1:size(imgs,1)
    
    obj.img(ii).vol = spm_vol(imgs{ii});
    [mn_mx(ii,2) mn_mx(ii,1)] = slover('volmaxmin', obj.img(ii).vol); % max & min
    
    if strcmp(itype{ii},'Structural')
        
        obj.img(ii).type  = 'truecolour';
        obj.img(ii).cmap  = gray(64); % "64" prevents creating a new figure
        obj.img(ii).range = [mn_mx(ii,1) mn_mx(ii,2)];
        cscale = [cscale ii];
        
    else
        
        % No suprathresold voxels
        %------------------------------------------------------------------
        if mn_mx(ii,2)==0
            mn_mx(ii,2) = 1;
        end
        % obj.img(ii).range = [n1 n2]
        % - if n1>n2: obj.img(ii).outofrange = {1 0} (error)
        % - if n2>n1: obj.img(ii).outofrange = {0 size(obj.img(ii).cmap,1)}
        % - if n1=0 : image does not show the background
        % => choose n1>0 and n2>n1 => n1=zero_thresh & n2=1
        
        switch itype{ii}
            case 'Truecolour'
                obj.img(ii).type = 'truecolour';
                %dcmap            = 'flow.lut';
                drange           = [mn_mx(ii,1) mn_mx(ii,2)];
                cscale           = [cscale ii];
                obj.cbar         = [obj.cbar ii];
            case 'Blobs'
                obj.img(ii).type = 'split';
                %dcmap            = 'hot';
                %drange           = [0 mx];
                drange           = [zero_thresh mn_mx(ii,2)]; % [eps ...]
                obj.img(ii).prop = 1;
                obj.cbar         = [obj.cbar ii];
            case 'Negative blobs'
                obj.img(ii).type = 'split';
                %dcmap            = 'winter';
                %drange           = [0 mn];
                drange           = [zero_thresh mn_mx(ii,2)]; % [eps ...]
                obj.img(ii).prop = 1;
                obj.cbar         = [obj.cbar ii];
            case 'Contours'
                obj.img(ii).type = 'contour';
                %dcmap            = 'white';
                drange           = [mn_mx(ii,1) mn_mx(ii,2)];
                obj.img(ii).prop = 1;
        end
        
        obj.img(ii).cmap  = slover('getcmap',cmap{ii}); % Nx3 colormap matrix
        %close(gcf)
        obj.img(ii).range = drange;
        
    end
end

ncmaps = length(cscale);
if ncmaps == 1 % one structural or truecolour image
    obj.img(cscale).prop = 1;
else
    remcol = 1;
    for ii = 1:ncmaps
        ino = cscale(ii);
        obj.img(ino).prop = remcol/(ncmaps-ii+1);
        remcol = remcol - obj.img(ino).prop;
    end
end

obj.transform = orient;


% Slices for display (mm)
%--------------------------------------------------------------------------
obj    = fill_defaults(obj);
%close(gcf)
slices = obj.slices;
if slc_step~=0
    obj.slices = (slices(1):slc_step:slices(end));
elseif slc_nmbr~=0
    obj.slices = round(linspace(slices(1),slices(end),slc_nmbr));
    obj.slices(~diff(obj.slices)) = []; % eliminate the repeated slices because of "round"
else
    obj.slices = slc_slices;
end


% Close the empty figures (slover('getcmap') may create an empty figure)
%--------------------------------------------------------------------------
figs = get(0, 'Children');
for ff=1:numel(figs)
   if isempty(get(figs(ff),'Children'))
       close(figs(ff));
   end
end


% Output
%--------------------------------------------------------------------------
%mn_mx
slvr_obj = obj;