function connext(opt)
% 
% Plot dynamically the MR signal and the histogram. Perform functional
% connectivity analysis.
% 
% 
%USAGE
%-----
% connext: open GUI
% connext(opt): run analysis
% 
% 
%INPUT
%-----
% - OPT can come from the GUI called by "connext_gui.m" or can
%   be created explicitly, with the following fields:
% 
% IMAGES
% - img: functional image file name
% 
% - coreg.use: 1 (use coregistered image [e.g., T1W]), 0 (don't use). Only
%   works when opt.plot=1 (see below).
% - coreg.img: coregistered image file name
% 
% - TR: repetition time (in seconds). Use 0 to ignore this option.
% 
% BASIC INPUT (the processing is performed in this order)
% - reg.use : 1 (use regressors), 0 (don't use)
% - reg.file: regressors file names (character array or cell array of strings)
%   For example, files containing: realignment parameters; white matter and
%   CSF signals; expected hemodynamic response of a particular condition
% 
% - detrend: 1 (apply detrend to the data), 0 (don't use)
% 
% - filt.use   : 1 (use frequency filter), 0 (don't use)
% - filt.type  : 1 (high-pass filter using DCT), 2 (rectangular band-pass
%   filter), 3 (Butterworth IIR filter), 4 (FIR filter), 5 (FIR filter
%   using least-squares error minimization)
% - filt.cutoff: [low_cutoff high_cutoff] (for a high-pass filter, use
%   low_cutoff=0; for a low-pass filter, use high_cutoff=Inf)
% - filt.order : filter order (positive integer or [] [use default])
% - filt.dir   : 'onepass', 'onepass-reverse', 'twopass',
%   'twopass-reverse', 'twopass-average', '' (use default)
% 
% - movavg.use     : 1 (perform moving average), 0 (don't perform)
% - movavg.windowsz: window size (integer greater than 1)
% 
% - avgblock.use: 1 (average across blocks), 0 (don't average)
% - avgblock.ons: Nblocks x 2 matrix with the first and last elements to
%   take the average of the signal (units: scans)
% 
% - scale.use : 1 (scale the signal), 0 (use raw signal). If
%   opt.funccon.use=1, opt.scale.use is ignored.
% - scale.type: 1 (Z-score), 2 (specified range), 3 (percent change from
%   the mean)
% - scale.rng : [minimum maximum] (only used when opt.scale.type=2)
% 
% - mask.use   : 1 (use mask), 0 (don't use)
% - mask.file  : mask file name
% - mask.thresh: threshold value to binarize the mask
% 
% SIGNAL PLOT
% - plot: 1 (plot signal), 0 (don't plot)
% 
%   - ROI.use: 1 (plot signal from rectangular ROI), 0 (don't use ROI)
%   - ROI.x  : no. of voxels for increasing x
%   - ROI.y  : no. of voxels for increasing y
%   - ROI.z  : no. of voxels for increasing z
% 
%   - plot_events.use : 1 (show rectangles indicating events on the MR
%     signal plot), 0 (don't show)
%   - plot_events.cond: Ncond structure array with the fields:
%     - name : condition name (string)
%     - time : Nevents x 2 matrix with the onsets in column 1 and the end
%       of the stimulus in column 2. The units must be in seconds if TR>0,
%       or in scans if TR=0. Scan 1 corresponds to instant 0, scan 2 to
%       instant TR, and so on.
%     - color: 1x3 matrix with color in RGB scale
% 
%   - hist.use       : 1 (plot histogram of the plotted signal), 0 (don't
%     plot). If opt.funccon.use=1, opt.hist.use is ignored.
%   - hist.binsz.use : 1 (use bin size option), 0 (don't use)
%   - hist.binsz.size: bin size (in signal units)
%   - hist.nbins.use : 1 (use number of bins option), 0 (don't use). If
%     opt.hist.binsz.use=1, this option is ignored.
%   - hist.nbins.bins: number of bins
% 
% FUNCTIONAL CONNECTIVITY
% - funccon.use        : 1 (do functional connectivity analysis), 0 (don't
%   do). If opt.funccon.use=1, opt.scale.use and opt.hist.use are
%   ignored.
% - funccon.seed.use   : 0 (only ctrl+click to select the seed position), 1
%   (use initial seed position), 2 (time series)
% - funccon.seed.pos   : 1x4 matrix with the seed position. If
%   funccon.seed.pos(4)=0, the coordinates (first three elements) are in
%   milimeters; if funccon.seed.pos(4)=1, the coordinates are in voxels of
%   the functional scan.
% - funccon.roi.use    : 1 (use seed as rectangular ROI), 0 (don't use).
%   This option is ignored when opt.funccon.seed.use=2.
% - funccon.roi.x      : no. of voxels for increasing x
% - funccon.roi.y      : no. of voxels for increasing y
% - funccon.roi.z      : no. of voxels for increasing z
% - funccon.seed.series: Nscans x 1 matrix with the time series of a seed
%   region. This option only works when opt.funccon.seed.use=2.
% 
% SAVE OUTPUT IMAGE
% - save.use  : 1 (save the processed image or the functional connectivity
%   image), 0 (don't save).
% - save.fname: output file name
% 
% 
% IMPORTANT: all voxel coordinates are rounded to the closest integer
% 
% 
%OUTPUT
%------
% - The signal of all voxels will be processed according to the chosen
% options.
% - If opt.plot=1 (it always is when running the GUI):
%   - Two windows:
%     1) SPM Display window
%        - If opt.funccon.use=1, blobs of the correlation coefficient are
%          overlaid on the functional or coregistered image.
%     2) Figure with the signal plot, and the histogram (if it was chosen)
%        or the seed voxel time series for functional connectivity (this
%        option overrides the histogram option). All the voxel coordinates
%        are in mm.
% - If opt.save.use=1: new file with name opt.save.fname. If the extension
%   or path are not given, they will be obtained from the functional file.
%   - If opt.funccon.use=0, it will contain the processed signal from
%     opt.img.
%   - If opt.funccon.use=1, it will contain the correlation coefficients
%     (note that a seed position or a time series must be given, that is,
%     opt.funccon.seed.use=1 or 2).
% 
%__________________________________________________________________________
% Copyright (C) 2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2014-Apr-27


if nargin==0
    connext_gui;
    return
end

%==========================================================================
% Initialize
%==========================================================================

thresh_minmax = 10^(-5);  % threshold to consider the signal as having the
% same value, when scaling the data for a specified range
thresh_std    = 10^(-10); % threshold to consider the standard deviation = 0
opt.show      = opt.plot; % show image when plotting the signal
if opt.show
    global st
end


if nargin<1
    disp('Invalid number of input parameters')
    return
end

% Check part of the input
%------------------------
if any(~isfield(opt, {'img', 'coreg', 'reg', 'filt', 'scale', 'mask', ...
        'avgblock', 'movavg', 'plot', 'ROI', 'plot_events', 'hist', 'funccon', 'save'}))
    % The script will try to read the TR from the image header
    disp('Invalid structure as input. Check the function''s help for more information.')
    return
end

% Check if there is anything to do
%---------------------------------
if ~opt.plot
    if ~opt.save.use
        disp('Nothing to do. You must either plot and save the results.')
        return
    elseif ~opt.reg.use && ~opt.filt.use && ~opt.scale.use && ~opt.mask.use && ...
            ~opt.avgblock.use && ~opt.movavg.use && ~opt.funccon.use
        disp('Nothing to do. Choose at least one processing option.')
        return
    end
    opt.coreg.use = 0;
end


%==========================================================================
% SPM display
%==========================================================================
if opt.show % image is shown
    
    % Initialize
    %-----------
    %spm_image('Reset');
    spm_figure('GetWin', 'Graphics');
    spm_figure('Clear', 'Graphics');
    spm_orthviews('Reset'); % variable "st" is updated
    set(st.fig, 'Pointer', 'watch')
    
    % Show images
    %------------
    if ~opt.coreg.use % only the functional image was given
        opt.coreg.use = 0;
        tmp = spm_vol(opt.img);
        spm_image('Init', tmp);
        nii = st.vols{1}.private; % NIfTI object
    else % 2 images
        tmp = spm_vol(opt.coreg.img);
        spm_image('Init', tmp);
        nii = nifti(opt.img); % read NIfTI header
    end
    set(st.fig, 'Name', 'Signal plot and functional connectivity')
    
    % Hide the figure 
    %----------------
    % To revent the user from clicking on the figure before it loads completely
    set(st.fig, 'Visible', 'off')
    drawnow
    
    % Change the ButtonDownFcn property
    %----------------------------------
    % To change the behaviour when clicking on the images
    tmp = findobj(get(st.fig, 'Children'), 'Type', 'axes'); % find axes objects
    if length(tmp)~=3
        error('Axes object could not be found inside the "Graphics" window')
    end
    set(tmp, 'ButtonDownFcn', @repos_start_new);
    
    % Change the Callback property
    %-----------------------------
    curr_callbackmp = get(st.mp, 'Callback');
    curr_callbackvp = get(st.vp, 'Callback');
    set(st.mp, 'Callback', {@callbackmp}) % coordinate in mm
    set(st.vp, 'Callback', {@callbackvp}) % coordinate in vx
    tmp = findobj(st.fig, 'TooltipString', 'move crosshairs to origin', 'Type', 'uicontrol');
    if ~isempty(tmp) % button was found
        curr_callbackorig = get(tmp(1), 'Callback');
        set(tmp(1), 'Callback', {@callbackorig}) % origin button
    end
    
    % Position for the three slices (function "bbox" in "spm_orthviews"):
    %------------------------------
    % - Transverse: [offx offy s*Dims(1) s*Dims(2)]
    % - Coronal   : [offx offy+s*Dims(2)+sky s*Dims(1) s*Dims(3)]
    % - Sagittal  : if st.mode=0: [offx+s*Dims(1)+skx offy s*Dims(3) s*Dims(2)]
    %               if st.mode=1: [offx+s*Dims(1)+skx offy+s*Dims(2)+sky s*Dims(2) s*Dims(3)]
    
    tmp = get(st.vols{1}.ax{1}.ax, 'Units');
    set(st.vols{1}.ax{1}.ax, 'Units', 'pixels')
    disp_axi = get(st.vols{1}.ax{1}.ax, 'Position');
    set(st.vols{1}.ax{1}.ax, 'Units', tmp)
    
    %tmp = get(st.vols{1}.ax{2}.ax, 'Units');
    %set(st.vols{1}.ax{2}.ax, 'Units', 'pixels')
    %disp_cor = get(st.vols{1}.ax{2}.ax, 'Position');
    %set(st.vols{1}.ax{2}.ax, 'Units', tmp)
    
    tmp = get(st.vols{1}.ax{3}.ax, 'Units');
    set(st.vols{1}.ax{3}.ax, 'Units', 'pixels')
    disp_sag = get(st.vols{1}.ax{3}.ax, 'Position');
    set(st.vols{1}.ax{3}.ax, 'Units', tmp)
    
    % Create save button
    %-------------------
    if ~opt.funccon.use
        tmp = 'Save the processed image';
    else
        tmp = 'Save the correlation coefficients';
    end
    tmp = uicontrol(...
        'Callback'     , @save_img, ...
        'String'       , 'Save image', ...
        'Style'        , 'pushbutton', ...
        'TooltipString', tmp, ...
        'Parent'       , st.fig, ...
        'Tag'          , 'disp_save_button');
    bb_WS = spm('WinScale');
    sz_button = [100*bb_WS(3) 30*bb_WS(4)];
    set(tmp, 'Units', 'pixels', ...
        'Position', [disp_sag(1)+disp_sag(3)-sz_button(1) ...
        disp_axi(2)+disp_axi(4)-sz_button(2) ...
        sz_button])
    set(tmp, 'Units', 'normalized')
    % The size of the button is based on the size of the texts from the
    % display window (in "spm_image.m", search for "uicontrol")
    
else
    nii = nifti(opt.img); % read NIfTI header
end

if length(nii.dat.dim)<4 || nii.dat.dim(4)==1
    waitfor(warndlg('3D image. Please select a 4D image.', 'Warning', 'modal'))
    if opt.show
        delete(st.fig) % close figure
    end
    return
end


%==========================================================================
% Read the NIfTI file
%==========================================================================
fprintf('Reading data........................ ')
dims   = nii.dat.dim(1:3);    % spatial dimensions
Nscans = nii.dat.dim(4);      % number of scans
data   = nii.dat(:, :, :, :); % read NIfTI data
data   = reshape(data, prod(dims), Nscans)'; % time along the rows and voxels along the columns
Mm2v   = inv(nii.mat);        % for the converstion from mm to voxel coordinates
Mv2m   = nii.mat;             % for the converstion from voxel to mm coordinates
if isfield(opt, 'TR')
    TR = opt.TR;
else
    TR = nii.hdr.pixdim(5);
end
fprintf('Done!\n')


%==========================================================================
% Basic input
%==========================================================================

% Mask
%-----
if opt.mask.use
    
    fprintf('Applying mask....................... ')
    
    Vmask = spm_vol(opt.mask.file);
%     mask = nifti(opt.mask.file);
%     if isequal([mask.dat.dim(1:3) mask.hdr.pixdim(2:4)], [nii.dat.dim(1:3) nii.hdr.pixdim(2:4)])
    if ~isequal([Vmask.dim(1:3) Vmask.private.hdr.pixdim(2:4)], [nii.dat.dim(1:3) nii.hdr.pixdim(2:4)])
        % functional image and mask have different dimensions or voxel size
        
        % From spm_imcalc.m via spm_imcalc_ui.m:
        mask = zeros(nii.dat.dim(1:3));
        for pp = 1:nii.dat.dim(3) % loop over planes
            tmp = spm_matrix([0 0 -pp 0 0 0 1 1 1]);
            tmp = inv(tmp*inv(nii.mat)*Vmask.mat);
            tmp = spm_slice_vol(Vmask, tmp, nii.dat.dim(1:2), [1,NaN]);
            if prod(nii.dat.dim(1:2)) ~= numel(tmp) % incompatible image
                opt.mask.use = 0;
                break
            end
            mask(:,:,pp) = reshape(tmp, nii.dat.dim(1:2));
        end
%     
%     tmp = mri_same_res(opt.img, opt.mask.file);
%     if ~isempty(tmp) % new mask file was created with different dimensions
%         opt.mask.file = tmp;
%         %delete(tmp) % remove created file
%     end

    else
%         mask = nifti(opt.mask.file);
%         mask = mask.dat(:,:,:);
        mask = spm_read_vols(Vmask); % (same result but ~2x slower)
    end
    
    % Apply mask
    if opt.mask.use
        mask = mask>opt.mask.thresh; % binarize
        data = reshape(data', dims(1), dims(2), dims(3), Nscans); % restore data shape
        data = bsxfun(@times, data, mask);
        data = reshape(data, prod(dims), Nscans)'; % reshape again
        %if opt.funccon.use % run this when the mask is applied in the end
        %    dataFC = bsxfun(@times, dataFC, mask);
        %end
    end
    
    clear Vmask mask
    fprintf('Done!\n')
    
end


% Regressors & Detrend
%---------------------
% The regressors can also be removed after both the data and the regressors
% have been filtered ("Simult" approach described in Hallquist et al.,
% NeuroImage, 82:208-225, 2013), but the results are very similar when the
% regressors are removed before filtering ("RegBp"). Do not perform the
% correction in the other direction ("BpReg": filtering before removing
% covariates), because frequencies outside the bandpass may reappear. More
% about this can be read in the above mentioned article.
% 
% Another issue is the order between temporal filtering and detrending. The
% filter cannot remove the linear trend very well, because linear trend is
% not a single slow frequency component (it will "leak" and generate
% frequencies outside the passband if done after the filtering). So detrend
% first and then do temporal filtering.

if opt.reg.use
    opt.reg.file = cellstr(opt.reg.file);
else
    opt.reg.file = '';
end
if opt.detrend
    opt.reg.use = 1; % detrending is done through GLM (same result as detrend(data))
    nuisance    = (1:Nscans)'/Nscans;
else
    nuisance    = [];
end
if opt.reg.use
    
    fprintf('Removing regressors................. ')
    
    try
        
        % Load files
        for ff=1:length(opt.reg.file)
            nuisance = [nuisance load(opt.reg.file{ff})];
        end
        
        % Remove the mean
        nuisance = bsxfun(@minus, nuisance, sum(nuisance, 1)./Nscans);
        % //or// nuisance = spm_detrend(nuisance, 0);
        % //or// nuisance = detrend(nuisance, 'constant');
        
        % Get the Z-score of the nuisance factors (to put everything in the same scale)
        tmp = sum(abs(nuisance).^2, 1)./(Nscans-1); % abs guarantees a real result
        tmp = sqrt(tmp);                            % standard deviation along time (row vector)
        tmp(tmp<=thresh_std) = 1;                   % eliminate points where std=0
        nuisance = bsxfun(@rdivide, nuisance, tmp); % Z-score
        
        % Add a column of 1's
        nuisance = [nuisance ones(Nscans, 1)];
        
        % Adapted from "regress.m" (similar to "y_regress_ss.m" from DPARSF):
        % 
        % Use the rank-revealing QR to remove dependent columns of "tmp"
        [Q, R, perm] = qr(nuisance, 0);
        Ncol = size(nuisance, 2);        
        p = sum(abs(diag(R)) > max(size(nuisance, 1), Ncol)*eps(R(1)) );
        if p < Ncol
            %warning('stats:regress:RankDefDesignMat', ...
            %    'X is rank deficient to within machine precision.');
            R    = R(1:p, 1:p);
            Q    = Q(:, 1:p);
            perm = perm(1:p);
        end
        
        % Compute the LS coefficients, filling in zeros in elements
        % corresponding to rows of "nuisance" that were thrown out
        b = zeros(Ncol, size(data,2));
        b(perm, :) = R \ (Q'*data);
        
        % Calculating directly gives similar results:
        %b2 = nuisance\data;
        % In one test, the maximum difference between b and b2 is of order
        % 10^(-9). And the difference between the corrected data signal
        % using b and b2 is of order 10^(-12).
        
        % Use a linear model to remove unwanted signal variation
        mu   = sum(data, 1)./Nscans;    % mean along time (row vector)
        data = data - nuisance*b;
        data = bsxfun(@plus, data, mu); % return the mean (it's useful for
        % displaying results because it maintains the color intensities)
        
        clear nuisance Q R perm Ncol p b mu
        
    catch
        waitfor(warndlg(['Regressors will not be used. Probably the size ' ...
            'of the matrices in the files do not match.'], 'Warning', 'modal'))
        opt.reg.use = 0;
    end
    
    fprintf('Done!\n')
    
end


% Filtering
%----------
if opt.filt.use && TR>0
    
    fprintf('Frequency filter.................... ')
    
    if strcmp(opt.filt.type, 'dct') % High-pass filter using DCT (discrete cosine transform)
        % The error distribution becomes normal
        if opt.filt.cutoff(1)>0
            data = fmri_hpf(data, TR, 1/opt.filt.cutoff(1));
            %nuisance = fmri_hpf(nuisance, TR, 1/opt.filt.cutoff(1));
        end
    elseif strcmp(opt.filt.type, 'rectbandpass') % Rectangular band-pass filter
        if isinf(opt.filt.cutoff(2)) % ignore this frequency
            opt.filt.cutoff(2) = 0;
        end
        data = bandpass_ideal_rect_wdw(data, TR, opt.filt.cutoff(1), opt.filt.cutoff(2), 1);
        %nuisance = bandpass_ideal_rect_wdw(nuisance, TR, opt.filt.cutoff(1), opt.filt.cutoff(2), 1);
    else % IIR or FIR filter
        tmp = bandpass_filt(data, 1/TR, opt.filt.cutoff, opt.filt.order, opt.filt.type, opt.filt.dir, 1);
        if ~isempty(tmp) % everything went fine
            data = tmp;
            %nuisance = bandpass_filt(nuisance, 1/TR, opt.filt.cutoff, opt.filt.order, opt.filt.type, opt.filt.dir, 1);
        else
            waitfor(warndlg(['Unable to filter the data. Try a higher ' ...
                'cutoff frequency or a different type/order of filter.'], ...
                'Frequency filter', 'modal'))
        end
    end
    
    fprintf('Done!\n')
    
end
% The regressors can also be removed after both the data and the regressors
% have been filtered ("RegBp" approach described in Hallquist et al.,
% NeuroImage, 82:208-225, 2013). By default, this script removes the
% regressors first.


% Moving average
%---------------
if opt.movavg.use
    fprintf('Moving average...................... ')
    data = mov_avg(data, opt.movavg.windowsz);
    fprintf('Done!\n')
end


% Average across events
%----------------------
if opt.avgblock.use && ~isempty(opt.avgblock.ons)
    fprintf('Average across blocks............... ')
    data   = block_average(data, opt.avgblock.ons, 0);
    Nscans = size(data, 1); % update the number of scans
    fprintf('Done!\n')
end


% Data scaling
%-------------
if opt.funccon.use % use Z-score to calculate the correlation coefficient
    opt.scale.use  = 1;
    opt.scale.type = 1;
end
if opt.scale.use && opt.scale.type>0
    
    fprintf('Data scaling........................ ')
    
    switch opt.scale.type
        case 1 % Z-score
            if opt.funccon.use
                % Back-up data to plot the raw signal
                tmp = data;
            end
            mu   = sum(data, 1)./Nscans;             % mean along time (row vector)
            data = bsxfun(@minus, data, mu);         % remove mean
            S    = sum(abs(data).^2, 1)./(Nscans-1); % abs guarantees a real result
            S    = sqrt(S);                          % standard deviation along time (row vector)
            S(S<=thresh_std) = 1;                    % eliminate points where S=0
            data = bsxfun(@rdivide, data, S);        % Z-score
            clear mu S % clean-up memory
            if ~opt.funccon.use
                ylabeltxt = 'MR signal (Z-score)';
            else % Z-score is also used to calculate the correlation coefficient
                ylabeltxt = 'MR signal (a.u.)';
                dataFC    = data/sqrt(Nscans-1); % used to calculate the correlation coefficient
                % using dataFC takes a lot of memory, but speeds up the
                % calculation of the correlation coefficient
                data      = tmp; % raw data
            end
        case 2 % Specified range
            rng = opt.scale.rng;
            minvox = min(data,[],1); % minimum value for each voxel
            maxvox = max(data,[],1); % maximum value for each voxel
            data = bsxfun(@minus, data, minvox);
            data = bsxfun(@rdivide, data, (maxvox-minvox))*(rng(2)-rng(1)) + rng(1);
            tmp = abs(minvox-maxvox)<=thresh_minmax; % where the signal is constant
            data(:, tmp) = rng(1);
            ylabeltxt = 'MR signal (scaled)';
            clear minvox maxvox
        case 3 % Percent change from the mean
            mu   = sum(data, 1)/Nscans;            % mean along time (row vector)
            data = bsxfun(@minus, data, mu);       % remove mean
            data = 100*bsxfun(@rdivide, data, mu); % divide by the mean
            ylabeltxt = 'MR signal (% change from the mean)';
    end
    
    fprintf('Done!\n')
    
else % no scaling will be applied => raw data
    ylabeltxt = 'MR signal (a.u.)';
end


%==========================================================================
% Signal plot
%==========================================================================

if opt.plot
    
    fprintf('Preparing to plot................... ')
    
    % Choose rectangular ROI
    %-----------------------
    if opt.ROI.use % plot rectangle
        plot_roi_sz = [opt.ROI.x opt.ROI.y opt.ROI.z]';
    else % do not plot rectangle
        plot_roi_sz = [0 0 0]';
    end
    
    % Plot events
    %------------
    if opt.plot_events.use % plot events
        cond  = opt.plot_events.cond;
        ncond = length(cond);
    else % do not plot events
        cond  = [];
        ncond = 0;
    end
    
    % Histogram
    %----------
    plot_hist_binsz = opt.hist.binsz.use;
    if plot_hist_binsz % histogram bin size was given
        plot_hist_opt = opt.hist.binsz.size;
    elseif opt.hist.nbins.use % number of bins was given
        plot_hist_opt = opt.hist.nbins.bins;
    else % just in case...
        plot_hist_binsz = 0;
    end
    
end


%==========================================================================
% Figure
%==========================================================================
% The figure must be created before the functional connectivity analysis

if opt.plot
    
    % Create figure
    %--------------
    tmp = findobj('Tag', 'SPM_MR_signal'); % find other figures from this script
    for tt=1:length(tmp) % and close all of them
        delete(tmp(tt));
    end
    fig = figure(...
        'Name'       , [mfilename ': MR signal'], ...
        'NumberTitle', 'off', ...
        'Tag'        , 'SPM_MR_signal',...
        'Visible'    , 'off');
    %'CloseRequestFcn', @my_closereq, ...
    
    % Figure axes
    %------------
    fontsz = 12;
    
    if opt.hist.use || opt.funccon.use
        subplot(2,1,1)
    end
    
    % x-axis
    if TR
        t   = (0:(Nscans-1))*TR;
        tmp = xlabel('Time (s)');
    else % no TR could be found in the header
        t   = 1:Nscans;
        tmp = xlabel('Scan');
    end
    set(tmp, ...
        'FontSize'        , fontsz, ...
        'FontWeight'      , 'bold', ...
        'HandleVisibility', 'off')
    
    ylabel(ylabeltxt, ...
        'FontSize'        , fontsz, ...
        'FontWeight'      , 'bold', ...
        'HandleVisibility', 'off')
    
    set(gca, 'FontSize', fontsz, 'NextPlot', 'replacechildren')
    
    % Histogram
    %----------
    if opt.hist.use
        subplot(2,1,2)
        xlabel(ylabeltxt, ...
            'FontSize'        , fontsz, ...
            'FontWeight'      , 'bold', ...
            'HandleVisibility', 'off')
        ylabel('Count', ...
            'FontSize'        , fontsz, ...
            'FontWeight'      , 'bold', ...
            'HandleVisibility', 'off')
    end
    
    % Functional connectivity
    %------------------------
    if opt.funccon.use
        
        subplot(2,1,1)
        xlabel('')
        
        subplot(2,1,2)
        if TR
            tmp = xlabel('Time (s)');
        else
            tmp = xlabel('Scan');
        end
        set(tmp, ...
            'FontSize'        , fontsz, ...
            'FontWeight'      , 'bold', ...
            'HandleVisibility', 'off')
        ylabel('MR signal (Z-score)', ...
            'FontSize'        , fontsz, ...
            'FontWeight'      , 'bold', ...
            'HandleVisibility', 'off')
        
    end
    
    set(gca, 'FontSize', fontsz, 'NextPlot', 'replacechildren')
    
end


%==========================================================================
% Functional connectivity
%==========================================================================

if opt.funccon.use % functional connectivity was chosen
    
    if opt.show % image is shown
        
        % Create slider
        %--------------
        hslider = uicontrol(...
            'Callback'     , @funccon_slider, ...
            'Max'          , 1, 'Min', 0, ...
            'Style'        , 'slider', ...
            'SliderStep'   , [.01 .1], ...
            'TooltipString', 'Set the threshold for |r|', ...
            'Value'        , 0.3, ...
            'Parent'       , st.fig, ...
            'Tag'          , 'fc_slider');
        set(hslider, 'Units', 'pixels', ...
            'Position', [disp_sag(1) disp_axi(2) disp_axi(3)/12 disp_axi(4)])
        
        % Text for the slider
        %--------------------
        hslidertxt = uicontrol(...
            'BackgroundColor'    , [1 1 1], ...
            'HorizontalAlignment', 'left', ...
            'String'             , sprintf('|r| > %.4f', get(hslider, 'Value')), ...
            'Style'              , 'text', ...
            'Parent'             , st.fig, ...
            'Tag'                , 'fc_slider_txt');
        tmp = get(hslider, 'Position');
        set(hslidertxt, 'Units', 'pixels', ...
            'Position', [tmp(1)+tmp(3)*1.5 tmp(2)+tmp(4)*.5 100*bb_WS(3) 16*bb_WS(4)])
        set(hslider, 'Units', 'normalized')
        set(hslidertxt, 'Units', 'normalized')
        
        if opt.hist.use
            opt.hist.use = 0; % unable to plot the histogram
        end
        
    end
    
    % Seed as rectangular ROI
    %------------------------
    if opt.funccon.roi.use
        seed_roi = [opt.funccon.roi.x opt.funccon.roi.y opt.funccon.roi.z]';
    else
        seed_roi = [0 0 0]';
    end
    
    switch opt.funccon.seed.use
        case 0 % Click to select the image
        case 1 % Initial seed
            seed_pos = opt.funccon.seed.pos(1:3);
            % opt.funccon.seed.pos(4) = 0 (mm coord.) or 1 (voxel coord.)
            if opt.funccon.seed.pos(4)==0 % mm
                seed_pos = Mm2v(1:3, 1:3)*seed_pos(:) + Mm2v(1:3, 4);
                % Alternative: seed_pos = Mm2v*[seed_pos(:); 1]; seed_pos = seed_pos(1:3);
                seed_pos = round(seed_pos); % round to the closest voxel
            end
            dataFC = reshape(dataFC', dims(1), dims(2), dims(3), Nscans); % restore data shape
            funccon_seed; % create "corrFC" and show blobs
        case 2 % Given time series
            dataFC = reshape(dataFC', dims(1), dims(2), dims(3), Nscans); % restore data shape
            funccon_seed; % create "corrFC" and show blobs
        otherwise
            error('Unknown option for the seed.')
    end
    
    seed_pos = [];
    
end


%==========================================================================
% Finish
%==========================================================================

% Restore data shape
%-------------------
if length(size(data))==2
    data = reshape(data', dims(1), dims(2), dims(3), Nscans);
end
if opt.funccon.use
    if length(size(dataFC))==2
        dataFC = reshape(dataFC', dims(1), dims(2), dims(3), Nscans);
    end
end

if opt.show % image is shown
    if opt.funccon.seed.use==0
    % //or// if exist('corrFC', 'var')~=1
        corrFC = []; % initialize the variable
    end
    pause(1) % it takes some time to load everything, so let's wait
    set(st.fig, 'Visible', 'on') % make figure visible again
    set(fig   , 'Visible', 'on')
    set(st.fig, 'Pointer', 'arrow')
    drawnow
    fprintf('Done!\n') % assuming opt.show = opt.plot
end


%==========================================================================
% Save results
%==========================================================================
if opt.save.use
    
    fprintf('Saving the results.................. ')
    
    % File name
    %----------
    [fpath fname fext]    = fileparts(opt.img);
    [fpath2 fname2 fext2] = fileparts(opt.save.fname);
    if isempty(fpath2) % if path is not given, save the new file in the same folder
        fpath2 = fpath;
    end
    if isempty(fext2) % if extension is not found, save the new file with the same extension
        fext2 = fext(1:4); % sometimes fext='.nii,1'
    end
    
    save_img(fpath2, [fname2 fext2(1:4)]);
    fprintf('Done!\n')
    
end


% Clean up
clear bb_WS sz_button tmp

fprintf('\n') % skip one line at the end


%==========================================================================
% Auxiliary functions
%==========================================================================
function repos_start_new(varargin)
    % When moving the crosshair in the images
    % From "spm_orthviews.m"
    % don't use right mouse button to start reposition
    try
        if ~strcmpi(get(gcbf,'SelectionType'), 'alt')
            % left-click: 'normal'; double click: 'open';
            % Shift + left-click: 'extend'
            set(gcbf, 'windowbuttonmotionfcn', @repos_move, 'windowbuttonupfcn', @repos_end);
            spm_orthviews('reposition'); % move the crosshair
            if opt.plot
                plot_signal;
            end
        elseif opt.funccon.use % Ctrl + left-click or right-click
            spm_orthviews('reposition'); % move the crosshair
            funccon_seed;
            plot_signal;
        end
    catch % sometimes moving the crosshair results in errors
        return
    end
end
%__________________________________________________________________________
%__________________________________________________________________________
function repos_move(varargin)
    % From "spm_orthviews.m"
    spm_orthviews('reposition');
    if opt.plot
        plot_signal;
    end
end
%__________________________________________________________________________
%__________________________________________________________________________
function repos_end(varargin)
    % From "spm_orthviews.m"
    set(gcbf, 'windowbuttonmotionfcn', '', 'windowbuttonupfcn', '');
end
%__________________________________________________________________________
%__________________________________________________________________________
function callbackmp(varargin)
    % Coordinate in mm
    eval(curr_callbackmp)
    if opt.plot
        plot_signal
    end
end
%__________________________________________________________________________
%__________________________________________________________________________
function callbackvp(varargin)
    % Coordinate in voxel
    eval(curr_callbackvp)
    if opt.plot
        plot_signal
    end
end
%__________________________________________________________________________
%__________________________________________________________________________
function callbackorig(varargin)
    % Button that takes to the origin
    eval(curr_callbackorig)
    if opt.plot
        plot_signal
    end
end
%__________________________________________________________________________
%__________________________________________________________________________
function plot_signal
    % Get voxel position and plot MR signal
    
    % Get position in voxel coordinates
    %----------------------------------
    if ~opt.coreg.use % only the functional image
        
        pos = round(spm_orthviews('Pos', 1)); % get(st.vp, 'String')
        pos(pos<1) = 1; % minimum voxel position is 1
        
    else % opt.coreg.use = 1
        
        % Get position in mm coordinates of the coregistered image
        %pos = spm_orthviews('Pos') = get(st.mp, 'String') = st.centre
        
        % If the images are coregistered, the coordinates in mm should be
        % the same for both images
        
        % For the coregistered image:
        % Mm2v = inv(st.vols{1}.premul*st.vols{1}.mat);
        % (from function "pos" in "spm_orthviews.m")
        % (premul = eye(4), in "specify_image" from "spm_orthviews")
        % Alternatives:
        % 1) nii = st.vols{1}.private;
        % 2) nii.mat = spm_get_space(st.vols{1}.fname)
        % 3) V = spm_vol(st.vols{1}.fname);
        %    nii.mat = V(1).mat
        
        % For the functional image:
        %Mm2v = inv(nii.mat); % this is done when the image is loaded
        
        % mm to voxel coordinate:
        pos = Mm2v(1:3,1:3)*st.centre(:) + Mm2v(1:3,4);
        % Alternative: pos = Mm2v*[pos(:) ; 1]; pos = pos(1:3);
        pos = round(pos); % round to the closest voxel
        
        % voxel to mm coordinate:
        % pos = Mv2m*[vox(:) ; 1]; pos = pos(1:3);
        
    end
    
    % Get the voxels from the rectangular ROI
    %----------------------------------------
    vx = cell(3, 1);
    for dd=1:3 % loop for the 3 coordinates
        vx{dd} = pos(dd) : (pos(dd)+plot_roi_sz(dd)); % maximum voxel position is the image dimension
        vx{dd}(vx{dd}>dims(dd)) = [];
    end
    
    % Get the signal
    %---------------
    data_roi = data(vx{1}, vx{2}, vx{3}, :); % signal from ROI
    if numel(data_roi)>Nscans % more than 1 voxel
        data_roi = squeeze(sum( sum( sum(data_roi, 1), 2), 3)) / prod(cellfun(@length, vx)); % average signal
        % prod(cellfun(@length, vx)) = size(data_roi, 1) * size(data_roi, 2) * size(data_roi, 3)
        %txt = sprintf('Voxels: [%d %d %d] - [%d %d %d]', ...
        %    vx{1}(1), vx{2}(1), vx{3}(1), vx{1}(end), vx{2}(end), vx{3}(end)); % voxel coord.
        txt = Mv2m*[vx{1}(1) vx{1}(end); vx{2}(1) vx{2}(end); vx{3}(1) vx{3}(end); 1 1];
        txt = txt(1:3, :);
        txt = sprintf('Voxels: [%.1f %.1f %.1f] - [%.1f %.1f %.1f]', ...
            txt(1, 1), txt(2, 1), txt(3, 1), txt(1, 2), txt(2, 2), txt(3, 2)); % voxel coord.
        if opt.funccon.use && ~isempty(corrFC) % write also the correlation coefficient
            tmp = squeeze(corrFC(vx{1}, vx{2}, vx{3}));
            txt = sprintf('%s (r = %.2f - %.2f)', txt, ...
                min(min(min(tmp, [], 1), [], 2), [], 3), ...
                max(max(max(tmp, [], 1), [], 2), [], 3));
        end
    else % only 1 voxel
        data_roi = squeeze(data_roi);
        %txt = sprintf('Voxel: [%d %d %d]', vx{1}, vx{2}, vx{3}); % voxel coord.
        txt = Mv2m*[[vx{:}]' ; 1]; txt = txt(1:3);
        txt = sprintf('Voxel: [%.1f %.1f %.1f]', txt(1), txt(2), txt(3)); % mm coord.
        if opt.funccon.use && ~isempty(corrFC) % write also the correlation coefficient
            txt = sprintf('%s (r = %.4f)', txt, corrFC(vx{1}, vx{2}, vx{3}));
        end
    end
    
    % Plot the signal
    %----------------
    figure(fig)
    if opt.hist.use || opt.funccon.use
        subplot(2,1,1)
    end
    plot(t, data_roi)
    xlim([t(1) t(end)])
    title(txt, ...
        'FontSize'  , fontsz, ...
        'FontWeight', 'bold')
    grid on
    
    % Plot events
    %------------
    if ncond>0
        
        tmp = cell(ncond, 1);
        [tmp{:}] = deal(cond.time);
        
        figure(fig) % to guarantee that the events will be on the correct subplot
        if opt.hist.use
            subplot(2,1,1)
        end
        
        plot_events(fig, tmp, reshape([cond.color], 3, ncond)')
        
        if ncond>1 % for 1 condition, the legend is unnecessary
            tmp = cell(ncond, 1);
            [tmp{:}] = deal(cond.name);
            figure(fig)
            legend(cat(1, 'MR signal', tmp), 'FontSize', 8)
        end
        
    end
    
    % Histogram
    %----------
    if opt.hist.use
        figure(fig)
        subplot(2,1,2)
        try
            if ~plot_hist_binsz
                if ~isempty(which('histfit')) % statistics toolbox is required
                    histfit(data_roi, plot_hist_opt, 'normal')
                else
                    hist(data_roi, plot_hist_opt)
                end
            else
                tmp = min(data_roi):plot_hist_opt:max(data_roi);
                if length(tmp)>1
                    hist(data_roi, tmp)
                else % constant value
                    hist(data_roi, 5)
                end
            end
        catch
            waitfor(warndlg('Unable to plot histogram. Check bin size.', 'Warning', 'modal'))
        end
        grid on
    end
    
end
%__________________________________________________________________________
%__________________________________________________________________________
function funccon_seed
    % Get seed position, plot seed time series and add blobs to the image
    
    if opt.funccon.seed.use~=2 % time series was not given
        
        % Get position in voxel coordinates
        %----------------------------------
        if opt.show && isempty(seed_pos) % image is shown
            if ~opt.coreg.use % only the functional image
                seed_pos = round(spm_orthviews('Pos', 1)); % get(st.vp, 'String')
                seed_pos(seed_pos<1) = 1; % minimum voxel position is 1
            else % coregistered image
                % mm to voxel coordinate:
                seed_pos = Mm2v(1:3,1:3)*st.centre(:) + Mm2v(1:3,4);
                % Alternative: seed_pos = Mm2v*[seed_pos(:) ; 1]; seed_pos = seed_pos(1:3);
                seed_pos = round(seed_pos); % round to closest voxel
            end
        end
        opt.funccon.seed.pos = seed_pos; % update seed position (in voxel coordinates)
        
        seed_vx = cell(3, 1);
        for dd=1:3 % loop for the 3 coordinates
            seed_vx{dd} = seed_pos(dd) : (seed_pos(dd)+seed_roi(dd)); % maximum voxel position is the image dimension
            seed_vx{dd}(seed_vx{dd}>dims(dd)) = [];
        end
        seed_pos = []; % reset seed position
        
        % Get the signal
        %---------------
        data_seed = dataFC(seed_vx{1}, seed_vx{2}, seed_vx{3}, :); % signal from ROI
        if numel(data_seed)>Nscans % more than 1 voxel
            data_seed = squeeze(sum( sum( sum(data_seed, 1), 2), 3)) / prod(cellfun(@length, seed_vx)); % average signal
            % prod(cellfun(@length, seed_vx)) = size(data_seed, 1) * size(data_seed, 2) * size(data_seed, 3)
            %txt = sprintf('Seed voxels: [%d %d %d] - [%d %d %d]', ...
            %    seed_vx{1}(1), seed_vx{2}(1), seed_vx{3}(1), seed_vx{1}(end), seed_vx{2}(end), seed_vx{3}(end)); % voxel coord.
            txt = Mv2m*[seed_vx{1}(1) seed_vx{1}(end); seed_vx{2}(1) seed_vx{2}(end); seed_vx{3}(1) seed_vx{3}(end); 1 1];
            txt = txt(1:3, :);
            txt = sprintf('Seed voxels: [%.1f %.1f %.1f] - [%.1f %.1f %.1f]', ...
                txt(1, 1), txt(2, 1), txt(3, 1), txt(1, 2), txt(2, 2), txt(3, 2)); % voxel coord.
        else % only 1 voxel
            data_seed = squeeze(data_seed);
            %txt = sprintf('Seed voxel: [%d %d %d]', seed_vx{1}, seed_vx{2}, seed_vx{3});
            txt = Mv2m*[[seed_vx{:}]' ; 1]; txt = txt(1:3);
            txt = sprintf('Voxel: [%.1f %.1f %.1f]', txt(1), txt(2), txt(3)); % mm coord.
        end
        
    else % time series was given
        
        if ~isvector(opt.funccon.seed.series)
            waitfor(warndlg('The seed time series must be a vector.', 'Warning', 'modal'))
            return
        end
        data_seed = opt.funccon.seed.series(:);
        if size(data_seed, 1)~=Nscans
            waitfor(warndlg(['The number of elements in the seed time series ' ...
                'does not match the number of scans.'], 'Warning', 'modal'))
            return
        end
        
        data_seed = data_seed - sum(data_seed, 1)./Nscans; % remove mean
        tmp       = sum(abs(data_seed).^2, 1);             % abs guarantees a real result
        tmp       = sqrt(tmp);
        tmp(tmp<=thresh_std) = 1;                          % eliminate points where tmp=0
        data_seed = data_seed / tmp;
        txt = 'Seed time series';
        
    end
    
    % Seed time series
    %-----------------
    if opt.plot
        figure(fig)
        subplot(2,1,2)
        plot(t, data_seed*sqrt(Nscans-1)) % multiply by sqrt(Nscans-1) to get the Z-score
        xlim([t(1) t(end)])
        title(txt, ...
            'FontSize'  , fontsz, ...
            'FontWeight', 'bold')
        grid on
    end
    
    % Correlation coefficient
    %------------------------
    dataFC = reshape(dataFC, prod(dims), Nscans)'; % time along the rows and voxels along the columns
    corrFC = sum(bsxfun(@times, dataFC, data_seed), 1);
    
    % Restore data shape
    %-------------------
    dataFC = reshape(dataFC', dims(1), dims(2), dims(3), Nscans);
    corrFC = reshape(corrFC', dims(1), dims(2), dims(3));
    
    % Add blobs to the image
    %-----------------------
    if opt.show % image is shown
        funccon_addblobs
    end
    
end
%__________________________________________________________________________
%__________________________________________________________________________
function funccon_addblobs(varargin)
	% Add blobs to the image [from spm_image('addblobs') and
	% spm_orthviews('AddColouredBlobs', ...)]
    
    % Indices where the correlation coefficients are above the threshold
    %-------------------------------------------------------------------
    %tmp = find(abs(corrFC)>=get(hslider, 'Value'));
    tmp = {find(corrFC>=get(hslider, 'Value')) find(corrFC<=-get(hslider, 'Value'))}; % positive and negative correlations
    if isempty(tmp{1}) && isempty(tmp{2}) % no suprathreshold voxel was found
        return
    end
    
    % Remove all blobs
    %-----------------
    spm_orthviews('RemoveBlobs', 1);
    %set(st.blobber, 'String', 'Add Blobs', 'Callback', 'spm_image(''AddBlobs'');');
    %spm_orthviews('RemoveContext', 1); 
    %spm_orthviews('Redraw');
    
    color = [1 0 0 ; 0 0 1]; % red for positive and blue for negative correlations
    for ii=1:2 % positive and negative correlations
        
        if isempty(tmp{ii}) % no suprathreshold voxel
            continue
        end
        
        % Get the equivalent subscript values
        %------------------------------------
        sub = zeros(3, length(tmp{ii}));
        [sub(1,:) sub(2,:) sub(3,:)] = ind2sub([dims(1) dims(2) dims(3)], tmp{ii});
        
        % Correlation coefficients above the threshold
        %---------------------------------------------
        FCthresh = corrFC(tmp{ii});
        
        % Add coloured blobs
        %-------------------
        spm_orthviews('AddColouredBlobs', 1, sub, FCthresh, Mv2m, color(ii, :));
        
    end
    
    set(st.blobber, 'String', 'Remove Blobs', 'Callback', 'spm_image(''RemoveBlobs'');'); % change button behaviour
    spm_orthviews('AddContext', 1); % add context menu (right-button click on the image) and the blobs
    spm_orthviews('Redraw');
    
end
%__________________________________________________________________________
%__________________________________________________________________________
function funccon_slider(varargin)
    set(hslidertxt, 'String', ['|r| > ' sprintf('%.4f', get(hslider, 'Value'))])
    if isempty(corrFC) % seed was not selected, but the slider was moved
        return
    end
    funccon_addblobs
end
%__________________________________________________________________________
%__________________________________________________________________________
function save_img(varargin)
    % Save NIfTI image
    
    % File name
    %----------
    if ishandle(varargin{1}) % button was clicked
        tmp = findobj(st.fig, 'Tag', 'disp_save_button');
        if isempty(tmp) % button was not found
            tmp = 'Save image';
        else
            tmp = get(tmp, 'TooltipString');
        end
        [fname fpath] = uiputfile({'*.nii', 'NIfTI files (*.nii)';
            '*.img', 'NIfTI files with separated header (*.img)'}, tmp); % get file name
        if isequal(fname, 0) || isequal(fpath, 0) % user clicked cancel or closed the window
            return
        end
    else % via script call
        fpath = varargin{1};
        fname = varargin{2};
    end
    
    [tmp fname_orig fext] = fileparts(opt.img);
    fname_orig = [fname_orig fext(1:4)];
    
    % NIfTI structure
    %----------------
    nii_new = nii; % new NIfTI structure
    nii_new.dat.fname = fullfile(fpath, fname);
    if ~opt.funccon.use % save the processed data
        nii_new.dat.dim(4)    = size(data, 4);
        nii_new.dat.dtype     = 'FLOAT32-LE'; % if integers are used, the values must be rounded
        nii_new.dat.scl_slope = 1;
        nii_new.dat.scl_inter = 0;
        txt = fname_orig;
        if length(txt)>60 % too long (the limit length for the description is 80 characters)
            txt = [txt(1:20) '...' txt((end-19):end)]; % => length = 43
        end
        nii_new.descrip      = sprintf('Data from %s', txt);
        nii_new.dat(:,:,:,:) = data;
    elseif opt.funccon.use % save the functional connectivity data
        if isempty(corrFC) % correlation coefficient was not calculated
            waitfor(warndlg(['The correlation coefficient has not been calculated yet. ' ...
                'Please choose a seed region first.'], 'Warning', 'modal'))
            return
        end
        nii_new.dat.dim       = nii.dat.dim(1:3);
        nii_new.dat.dtype     = 'FLOAT32-LE';
        nii_new.dat.scl_slope = 1;
        nii_new.dat.scl_inter = 0;
        nii_new.descrip       = sprintf('Corr. coeff. with seed on [%d %d %d]', opt.funccon.seed.pos);
        nii_new.dat(:,:,:)    = corrFC; % file is created
    end
    create(nii_new) % the header is adjusted
    
end
%__________________________________________________________________________
%__________________________________________________________________________
% function my_closereq(varargin) % 2 input arguments are necessary
%     % Display a question dialog box upon closing the figure
%     selection = questdlg('Close this figure?', ...
%         '', ...
%         'Yes', 'No', 'No');
%     switch selection
%         case 'Yes'
%             delete(fig)
%         case 'No'
%             return
%     end
% end

end
