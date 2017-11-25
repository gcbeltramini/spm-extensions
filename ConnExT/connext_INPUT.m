%                         =============
%                            CONNEXT
%                         =============
% 
% Input variables of ConnExT. For more details, see the help text of
% "connext".
% 
% Edit this file with your data. "connext_INPUT_template.m" contains some
% dummy values.


% Images
%==========================================================================
opt.img = 'SbjXX_EPI.nii'; % functional image file name
opt.coreg.use = 0; % 1 or 0
opt.coreg.img = 'SbjXX_T1W.nii'; % coregistered image file name

% Basic input
%==========================================================================

% Repetition time
%----------------
opt.TR = 2; % TR (in seconds). Use 0 to ignore TR (units: scans)

% Regressors
%-----------
opt.reg.use  = 0; % 1 or 0
opt.reg.file = 'covariates_times_series.txt'; % regressors file names (character array or cell array of strings)
%  For example, files containing: realignment parameters; white matter and
%  CSF signals; expected hemodynamic response

% Detrend
%--------
opt.detrend = 0; % 1 or 0

% Frequency filter
%-----------------
opt.filt.use  = 0; % 1 or 0
opt.filt.type = 1; % 1 (high-pass filter using DCT), 2 (rectangular band-pass
%  filter), 3 (Butterworth IIR filter), 4 (FIR filter), 5 (FIR filter
%  using least-squares error minimization)
opt.filt.cutoff = [0 Inf]; % [low_cutoff high_cutoff]
%  High-pass filter: low_cutoff=0
%  Low-pass filter : high_cutoff=Inf
opt.filt.order = []; % positive integer or [] (use default)
opt.filt.dir   = ''; % 'onepass', 'onepass-reverse', 'twopass',
%  'twopass-reverse', 'twopass-average', '' (use default)

% Moving average
%---------------
opt.movavg.use      = 0; % 1 or 0
opt.movavg.windowsz = 2; % window size (integer greater than 1)

% Average across blocks
%----------------------
opt.avgblock.use = 0; % 1 or 0
opt.avgblock.ons = [0 10]; % Nblocks x 2 matrix with the first and last
%  elements to take the average of the signal (units: scans)

% Scale the signal
%-----------------
opt.scale.use  = 0; % 1 or 0. If opt.funccon.use=1, opt.scale.use is ignored.
opt.scale.type = 1; % 1 (Z-score), 2 (specified range), 3 (percent change from the mean)
opt.scale.rng  = [-1 1]; % [minimum maximum] (when opt.scale.type=2)

% Mask
%-----
opt.mask.use    = 0;   % 1 or 0
opt.mask.file   = 'brainmask.nii'; % mask file name
opt.mask.thresh = 0.3; % threshold value to binarize the mask


% Signal plot
%==========================================================================

opt.plot = 0; % 1 (display the images [3 orthogonal planes and time series]) or 0

% Rectangular ROI
%----------------
opt.plot.ROI.use = 1; % 1 or 0
opt.plot.ROI.x   = 0; % no. of voxels for increasing x
opt.plot.ROI.y   = 0; % no. of voxels for increasing y
opt.plot.ROI.z   = 0; % no. of voxels for increasing z

% Plot events
%------------
opt.plot.plot_events.use = 0; % 1 or 0
opt.plot.plot_events.cond(cond1).name = 'Name of condition 1';
opt.plot.plot_events.cond(cond1).time = [10 30; 40 60; 70 90];
% no. events x 2 matrix: onsets in col. 1 and end of the stimulus in col. 2
%  The units must be in seconds if TR>0, or in scans if TR=0.
opt.plot.plot_events.cond(cond1).color = [.8 .8 .8]; % 1x3 matrix with color in RGB scale
% opt.plot.plot_events.cond(condN): analogous


% Histogram
%----------
opt.plot.hist.use        = 0; % 1 or 0. If opt.funccon.use=1, opt.hist.use is ignored.
opt.plot.hist.binsz.use  = 0; % 1 or 0
opt.plot.hist.binsz.size = 50; % bin size (in signal units)
opt.plot.hist.nbins.use  = 1; % 1 or 0. If opt.hist.binsz.use=1, this option is ignored.
opt.plot.hist.nbins.bins = 10;


% Functional connectivity
%==========================================================================

opt.funccon.use = 0; % 1 or 0. If opt.funccon.use=1, opt.scale.use and opt.hist.use are ignored.

% Type of seed
%-------------
opt.funccon.seed.use = 1; % 0 (only ctrl+click to select the seed position),
%  1 (use initial seed position), 2 (time series)

% Seed position (opt.funccon.seed.use = 1)
%--------------
opt.funccon.seed.pos = [0 -54 27 0]; % seed position
% If funccon.seed.pos(4)=0, the coordinates (first three elements) are in
% milimeters; if funccon.seed.pos(4)=1, the coordinates are in voxels of
% the functional scan.

% Seed as a rectangular ROI (opt.funccon.seed.use = 0 or 1)
%--------------------------
opt.funccon.roi.use = 0;
opt.funccon.roi.x   = 0;
opt.funccon.roi.y   = 0;
opt.funccon.roi.z   = 0;

% Time series (opt.funccon.seed.use = 2)
%------------
opt.funccon.seed.series = []; % Nscans x 1 matrix with the seed time series


% Save
%==========================================================================
opt.save.use   = 0; % 1 or 0
opt.save.fname = 'output_file.nii'; % output file name


% Run
%==========================================================================
connext(opt)