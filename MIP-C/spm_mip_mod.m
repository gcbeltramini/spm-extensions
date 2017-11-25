function [mip, d, grid_mask] = spm_mip_mod(Z, XYZ, M)
% 
% [mip, d, grid_mask] = spm_mip_mod(Z, XYZ, M)
% 
% SPM maximum intensity projection
% 
%USAGE
%-----
% [mip, d, grid_mask] = spm_mip(Z, XYZ, M);
% 
% 
%INPUT
%-----
% Z  : vector point list of SPM values for MIP
% XYZ: matrix of coordinates of points (mip coordinates)
% M  : voxels - > mip matrix or size of voxels (mm)
% 
% 
%OUTPUT
%------
% MIP      : maximum intensity projection
% D        : only the functional map
% GRID_MASK: only grid and brain mask
%__________________________________________________________________________
% 
% spm_mip creates a maximum intensity projection of a point list of voxel
% values (Z) and their location (XYZ) in three orthogonal views of the
% brain.  It is assumed voxel locations conform to the space defined in the
% atlas of Talairach and Tournoux (1988); unless the third dimension is
% time.
% 
% This routine loads a mip outline from MIP.mat. This is an image with
% contours and grids defining the space of Talairach & Tournoux (1988).
% mip95 corresponds to the Talairach atlas, mip96 to the MNI templates.
% The outline and grid are superimposed at intensity 0.4.
% 
% A customised mip outline can be used instead of the default.
% 
% A default colormap of 64 levels is assumed. The pointlist image is
% scaled to fit in the interval [0.25,1]*64 for display. Flat images
% are scaled to 1*64.
% 
% If M is not specified, it is assumed the XYZ locations are
% in Talairach mm.
% 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston et al.
% $Id: spm_mip.m 4860 2012-08-24 13:23:55Z volkmar $
% 
% Changelog
%--------------------------------------------------------------------------
% 2014-Jul-01, 08:11 pm: Modified by Guilherme Coco Beltramini to use with
%   the function "fmri_mip_color.m"


Grid = 0.4; % grid scaling


% Transpose locations if necessary
%==========================================================================
if size(XYZ,1) ~= 3, XYZ = XYZ';         end
if size(Z,1)   ~= 1, Z   = Z';           end
if size(M,1)   == 1, M   = speye(4,4)*M; end


% Scale & offset point list values to fit in [0.25,1]
%==========================================================================
Z    = Z - min(Z);
mx   = max(Z);
Scal = 8;
if isempty(mx),
    Z = [];
elseif isfinite(mx) && (numel(Z) ~= 1),
    Z = (1 + Scal*Z/mx)/(Scal + 1);
else
    Z = ones(1,length(Z));
end


% Load various grids, DXYZ, CXYZ, scale (see spm_mip_ui and spm_project)
%==========================================================================
%load(char(spm_get_defaults('stats.results.mipmat')));
load('MIP.mat');
clear mip95 mip96


% Grid and mask
%==========================================================================
grid_mask = 4*grid_all + mask_all;
grid_mask = grid_mask/max(grid_mask(:)); % sparse matrix
clear grid_all grid_time grid_trans mask_all mask_trans


% Maximum intensity projection (MIP)
%==========================================================================
c   = [0 0 0 ;
    0 0 1 ;
    0 1 0 ;
    0 1 1 ;
    1 0 0 ;
    1 0 1 ;
    1 1 0 ;
    1 1 1 ] - 0.5;
c   = c * M(1:3,1:3);
dim = [(max(c) - min(c)) size(grid_mask)];
if exist('DXYZ', 'var')==1 && exist('CXYZ', 'var')==1
    d = spm_project(Z, round(XYZ), dim, DXYZ, CXYZ);
else
    d = spm_project(Z, round(XYZ), dim);
end
mip = max(d, Grid*grid_mask);
clear Z XYZ dim DXYZ CXYZ


% MIP, grid and the combination of both
%==========================================================================
grid_mask = rot90(1-Grid*grid_mask); % grid and mask
d         = rot90(1-d);              % functional map
mip       = rot90(1-mip);            % grid, mask and functional map
