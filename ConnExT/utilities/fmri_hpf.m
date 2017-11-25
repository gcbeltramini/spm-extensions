function yf = fmri_hpf(y,TR,hp_period)
% 
% High-pass filter using DCT (discrete cosine transform)
% 
% 
%USAGE
%-----
% yf = fmri_hpf(y,TR,hp_period)
% 
% 
%INPUT
%-----
% - Y        : NxS matrix of S series with N time points
% - TR       : repetition time (in seconds)
% - HP_PERIOD: cut-off period (in seconds)
% 
% 
%OUTPUT
%------
% - YF: filtered data (same size as Y)
% 
% 
% See also SPM_FILTER, SPM_DCTMTX

% Guilherme Coco Beltramini - 2014-Mar-29, 06:08 pm


K(1).RT     = TR;          % [s]
K(1).row    = 1:size(y,1); % filter all data (select all rows)
K(1).HParam = hp_period;   % [s]
K  = spm_filter(K);        % create filter structure (set K.X0)
yf = spm_filter(K,y);      % filter data