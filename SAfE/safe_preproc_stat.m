function safe_preproc_stat(opt_safe)
% 
% Preprocessing and statistical analysis for various subjects and delays.
% 
% 
%==========================================================================
%                                 USAGE
%==========================================================================
% safe_preproc_stat(opt_safe)
% 
% 
%==========================================================================
%                                 INPUT
%==========================================================================
% 
% All the variables must be in small caps
% 
% OPT_SAFE: structure with the following fields:
% 
% 1) Basic input:
% BASIC: structure with the following fields:
% - PARENT_DIR: parent directory for all subjects
% - X = 1, 2, ... XLS file number identification
%   - XLSFILE{X}: full path for the XLS files where the conditions are
%   - WORKSHEET{X}: list of worksheet names from the cell XLSFILE{X} (cell
%     array)
%   - SBJ_FOLDER{X}: subject folder name inside PARENT_DIR for each cell in
%     WORKSHEET{X} (cell array)
% - FMRI_PATH: path where the fMRI files are inside the subject folders
% - EEG_FMRI_EXP: 1 (EEG-fMRI experiment) or 0 (regular fMRI experiment)
% 
% Therefore, the fMRI files for subject S in the XLS file X are in:
% PARENT_DIR\SBJ_FOLDER{X}{S}\FMRI_PATH
% 
% Subject's data is in XLSFILE{X}, WORKSHEET{X}{S}.
% 
% 
% 2) Preprocessing
%    -------------
% PREPROC: structure with the following fields:
% - RUN : 1 (preprocess) or 0 (do not preprocess)
% - STC : 1 (do slice timing) or 0 (do not do slice timing)
% - NORM: structure with the following fields:
%   - RUN: 1 (normalize) or 0 (do not normalize)
%   If RUN=1:
%   - STRUC_VX: voxel size of the normalized structural image (any
%     positive or negative number)
%     If negative numbers are used in some dimension, to get the voxel size
%     of the normalized structural image, their opposite sign, that is, a
%     positive number, is multiplied by the voxel size of the structural
%     image in that dimension. Zero is not allowed.
%   - FUNC_VX : voxel size of the normalized functional image. Analogous.
% - SMOOTH: structure with the following fields:
%   - RUN: 1 (smooth) or 0 (do not smooth)
%   If RUN=1:
%   - FWHM: FWHM of the Gaussian kernel. Analogous. Zero is allowed.
% 
% 
% 3) Statistical analysis: Specify and estimate model
%    ------------------------------------------------
% STAT_DES: structure with the following fields:
% - RUN       : 1 (specify and estimate model) 0 (do not specify)
% - PREF_FUNC : prefix for the preprocessed functional files
% - PREF_STRUC: prefix for the structural file
%   If OPT_SAFE.STC=1 & OPT_SAFE.NORM.RUN=1: PREF_FUNC='swa' & PREF_STRUC='wm'
%   If OPT_SAFE.STC=1 & OPT_SAFE.NORM.RUN=0: PREF_FUNC='sar' & PREF_STRUC='m'
%   If OPT_SAFE.STC=0 & OPT_SAFE.NORM.RUN=1: PREF_FUNC='sw'  & PREF_STRUC='wm'
%   If OPT_SAFE.STC=0 & OPT_SAFE.NORM.RUN=0: PREF_FUNC='sr'  & PREF_STRUC='m'
% - TAG_STRUC: text to identify the structural file; useful when PREF_STRUC
%   is empty or if there are more files with the prefix PREF_STRUC (e.g.,
%   '', 'VBM', 'T1W3D')
% - USE_RP: 1 (use realignment parameters) or 0 (do not use)
% - HPF   : high-pass filter [s]
% - BF    : basis function ('can_hrf', 'fourier', 'fourier_han', 'gamma' or
%           'fir')
%   If BF = 'fourier', 'fourier_han', 'gamma', 'fir':
%   - BF_LENGTH: post-stimulus window length [s]
% - BF_ORDER: number of basis functions (maximum 3 for BF='can_hrf')
% - DELAY   : time in seconds from the stimulus onset to the beginning of
%   the hemodynamic response function. DELAY is a vector of any real number
%   (positive, negative or zero)
% 
% 
% 4) Statistical analysis: Contrasts
%    -------------------------------
% - STAT_CON: 1 (create contrasts) or 0 (do not create)
% 
% 
% 5) Statistical analysis: Results
%    -----------------------------
%   (the variable DELAY is also used in the Results)
% 
% STAT_RES: structure with the following fields:
% 
% - RUN: 1 (perform any of the operations below) or 0 (do not perform)
% 
% - THRESH.ADJ: 'FWE' or 'none' (P-value adjustment)
% - THRESH.P  : P-value (0 <= P <= 1)
% - THRESH.ET : extend threshold in voxels (ET >= 0)
% 
% - PS_SAVE: 1 (create PS file) or 0 (do not create)
% 
% - SAVE_THR: 1 (create thresholded image) or 0 (do not create)
% It is an intermediate step when SAVE_SLOVER=1 or
% SLICES.SAVE=1 or GROUP.SAVE=1 or EVENT_SEQ.SAVE=1. So in those
% cases it does not matter whether SAVE_THR=1 or 0.
% 
% - SAVE_SLOVER: 1 (create slover object) or 0 (do not create)
%   - If SAVE_THR=1 or SLICES.SAVE=1, the slover object will be created
%   anyway (but with different fields).
%   - If EVENT_SEQ.SAVE=1 or GROUP.SAVE=1, SAVE_SLOVER will be set to 1.
% 
% - SLICES.SAVE: 1 (create image with slices) or 0 (do not create)
% If SLICES.SAVE=1 or EVENT_SEQ.SAVE=1 or GROUP.SAVE=1:
%   - SLICES.IMG_EXT: image extension ('bmp', 'fig', 'jpg', 'png', 'tif')
% If SLICES.SAVE=1 or GROUP.SAVE=1:
%   - SLICES.STEP  : spacing between slices (mm) (any positive number)
%   - SLICES.NUMBER: number of slices (to use this option, STEP must be 0)
%   - SLICES.SLICES: coordinates of the slices (to use this option, STEP
%   and NUMBER must be 0)
%   - SLICES.ORIENT: cell array of string, containig at least one of the
%     following options: 'axial', 'coronal', 'sagittal'
% 
% - GROUP.SAVE: 1 (group images from different delays) or 0 (do not group)
% If GROUP.SAVE=1:
%   - Options from SLICES will be applied here.
%   - GROUP.NUMBER   : number of delays per group
%   - GROUP.SCALE.ADJ: 1 (adjust scale) or 0 (do not adjust)
%   If GROUP.SCALE.ADJ=1:
%     - GROUP.SCALE.FTEST: 'max' (the first value will be 10^(-10) and the
%     last value will be the maximum value from all delays); or 1x2 matrix 
%     containing the first and last value in the colormap (the first value
%     should be a very small positive number)
%     - GROUP.SCALE.TTEST: 'max' (the first value will be 10^(-10) and the
%     last value will be the maximum value from all delays); or 1x4 matrix
%     containing the first and last value in the colormap for positive
%     (SCALE.TTEST([1 2])) and negative (SCALE.TTEST([3 4])) BOLD
%     (SCALE.TTEST(1) and SCALE.TTEST(3) should be a very small positive
%     number)
% 
% - EVENT_SEQ.SAVE: 1 (create sequence of images) or 0 (do not create)
%   It is similar to SLICES.SAVE, but the images are organized in the same
%   folder, with the names in a sequential order, and with the scale
%   adjusted.
% If EVENT_SEQ.SAVE=1:
%   - Options from SLICES will be applied here.
%   - EVENT_SEQ.SCALE.FTEST  : same as in GROUP.SCALE.FTEST
%   - EVENT_SEQ.SCALE.POSBOLD: same as in GROUP.SCALE.FTEST
%   - EVENT_SEQ.SCALE.NEGBOLD: same as in GROUP.SCALE.FTEST
% 
% - MASK.USE: 1 (apply masking) or 0 (do not apply) -- DISABLED
%   If MASK.use=1:
%   - MASK.IMGS_TEXT: cell array with the text to find the mask images
%   (e.g., {'mwc*'} or {'mwc1*','mwc2*'})
%   - MASK.INCL_EXCL: nature of mask ('inclusive' or 'exclusive')
%   - MASK.BIN: 1 (create a binary mask) or 0 (do not create)
%     If MASK.BIN=1:
%     - MASK.BIN_THRESH: threshold value to create binary image:
%     i1+i2+...>MASK_BIN_THRESH, where iN are the images in MASK_IMGS_TEXT
% 
% - CLEANUP: 1 (remove some auxiliary files) or 0 (don't remove)
%   - The auxiliary files are the thresholded NIfTI files in each delay
%   folder and in the results folder, and the MAT file(s) in the results
%   folder.
%   - It will ignore the options SAVE_THR and SAVE_SLOVER.
% 
% 
%==========================================================================
%                                OUTPUT
%==========================================================================
% 
% 1) Preprocessing (PREPROC=1)
%    -------------------------
% - XXP-preproc.mat in the same folder as the raw files (XX is the subject
%   number)
% - Realigment:
%   - rp_SUBJECT_NAME_SCAN_X.txt for each functional session (run)
%   - SUBJECT_NAME_SCAN_X.mat for each session
%   - meanSUBJECT_NAME_SCAN_X.nii
%   - XXP-realign-coreg.ps
% - Slice timing correction (STC=1):
%   - aSUBJECT_NAME_SCAN_X.nii for each session
% - Segmentation (NORM.RUN=1):
%   - SUBJECT_NAME_T1W3D_X_seg_sn.mat
%   - SUBJECT_NAME_T1W3D_X_seg_inv_sn.mat
%   - wc1SUBJECT_NAME_T1W3D_X.nii
%   - mwc1SUBJECT_NAME_T1W3D_X.nii
%   - wc2SUBJECT_NAME_T1W3D_X.nii
%   - mwc2SUBJECT_NAME_T1W3D_X.nii
%   - wc3SUBJECT_NAME_T1W3D_X.nii
%   - mwc3SUBJECT_NAME_T1W3D_X.nii
%   - mSUBJECT_NAME_T1W3D_X.nii
% - Normalization (NORM.RUN=1):
%   - waSUBJECT_NAME_SCAN_X.nii for each session (STC=1)
%   -  wSUBJECT_NAME_SCAN_X.nii for each session (STC=0)
%   - wmeanSUBJECT_NAME_SCAN_X.nii
%   - wmSUBJECT_NAME_T1W3D_X.nii
% - Smoothing (SMOOTH.RUN=1):
%   - swaSUBJECT_NAME_SCAN_X.nii for each session (STC=1, NORM.RUN=1)
%   -  swSUBJECT_NAME_SCAN_X.nii for each session (STC=0, NORM.RUN=1)
%   - sarSUBJECT_NAME_SCAN_X.nii for each session (STC=1, NORM.RUN=0)
%   -  srSUBJECT_NAME_SCAN_X.nii for each session (STC=0, NORM.RUN=0)
% 
% 
% 2) Statistical analysis: Specify and estimate model (STAT_DES=1)
%    -------------------------------------------------------------
% - Inside the subject's folder: folder WORKSHEET-SPM_delay_DELAY with
%   WORKSHEET-design_matrix.mat and WORKSHEET-rp-conditions.png
% - If there is already a non-empty folder named WORKSHEET-SPM_delay_DELAY,
%   it will be renamed to WORKSHEET-SPM_delay_DELAY-backup (or
%   WORKSHEET-SPM_delay_DELAY-backup02, WORKSHEET-SPM_delay_DELAY-backup03,
%   ...)
% - In folder WORKSHEET-SPM_delay_DELAY:
%   - beta_XXXX.hdr/img, mask.hdr/img, ResMS.hdr/img, RPV.hdr/img, SPM.mat
% 
% 
% 3) Statistical analysis: Contrasts (STAT_CON=1)
%    --------------------------------------------
% - WORKSHEET-contrasts.mat: batch file in the same folder as SPM.mat
% - SPM.mat will be updated.
% - If STAT_DES=0, it is assumed that SPM.mat is inside folder
%   WORKSHEET-SPM_delay_DELAY
% 
% 
% 4) Statistical analysis: Results
%    -----------------------------
% - If STAT_DES=0, it is assumed that SPM.mat is inside folder
%   WORKSHEET-SPM_delay_DELAY
% 
% - If MASK=1 and MASK_BIN=1: mask_for_results.nii in the same folder as
% the raw files -- DISABLED
% 
% - If PS_SAVE=1:
%   - WORKSHEET_DELAYs_P0.XXXX_Vvox_unc.ps or SID_DELAYs_P0.XXXX_Vvox_FWE.ps
%     in the same folder as the SPM.mat file, where:
%     - 0.XXXX is the P value and V is the extend threshold
%   - Two pages for each contrast: the first one shows the glass brain view
%   and the p-values for the whole brain; the second one shows the
%   statistical parametric map over the UNDERLAY image
% 
% - If PS_SAVE=1 or SAVE_THR=1 or SAVE_SLOVER=1 or SLICES.SAVE=1 or
% GROUP.SAVE=1 or EVENT_SEQ.SAVE=1, in the same folder as
% SPM.mat:
%   - Files "THR_PREF"_SID_CONDNAME_TEST_P0.XXXX_Vvox_CORR_DELAY.hdr/img,
%     where:
%     - SID     : the subject ID (e.g., 'XXPy')
%     - CONDNAME: condition name
%     - TEST    : 'posBOLD', 'negBOLD' or 'Ftest'
%     - 0.XXXX  : P value
%     - V       : extent threshold
%     - CORR    : 'unc' or 'FWE'
%     - DELAY   : time in seconds from the stimulus onset to the beginning
%     of the hemodynamic response function
%   - File
%   PARENT_DIR\SBJ_FOLDER{XLS}{W}\FMRI_PATH\WORKSHEET{XLS}{W}-Results_DELAY
%   s\WORKSHEET{XLS}{W}-slover_obj-DDdelays.mat (DD is the number of
%   delays) "safe_res": structure containing the following fields:
%     If SAVE_SLOVER=1 or SLICES.SAVE=1:
%     - MN_MX: 2Cx4 matrix containing the minimum (columns 1 and 3) and
%     maximum (columns 2 and 4) values of the statistics. C is the number
%     of conditions. Odd rows: F test (columns 1 and 2 only); even rows: T
%     test (positive BOLD in columns 1 and 2; negative BOLD in columns 3
%     and 4). For every two rows, the conditions are changed
%     - OBJ: 2Cx1 cell array containing the slover object for every
%     condition. Odd rows: F test; even rows: T test
%     If SAVE_SLOVER=1 or SAVE_THR=1 or SLICES.SAVE=1 or GROUP.SAVE=1:
%     - THR: structure containing the following fields:
%       - FNAME   : 3Cx1 cell array, containing the full path of the
%       thresholded images. Rows 1,4,7,...: F test; rows 2,5,8,...:
%       positive T test; rows 3,6,9,...: negative T test
%       - ADD_STAT: 3Cx1 matrix, with the summation of the statistics value
%       for each image in THR.FNAME
%       - VOXELS  : 3Cx1 matrix, with the number of "active" voxels for
%       each image in THR.FNAME
% 
% - If GROUP.SAVE=1: folder "...\Results_DELAYS\Grouped_delays" with
%   the delays grouped.
% 
% - If EVENT_SEQ.SAVE=1: folder "...\Results_DELAYS\Event_seq" with:
%   - sequence of images for all conditions and tests
%   - file Event_seq-Summary.txt with the data for each condition and each
%   statistical test of the following: maximum value, sum of the statistics
%   for the whole image and number of voxels above the threshold
%   - folder "Plots" containing the graphics for Event_seq-Summary.txt
% 
% 
% DEFAULT VALUES
%==========================================================================
% basic.parent_dir   = --
% basic.fmri_path    = fullfile('fMRI','Analise','SPM8');
% basic.xlsfile      = --
% basic.worksheet    = --
% basic.sbj_folder   = --
% basic.eeg_fmri_exp = 1 (EEG-fMRI experiment) or 0 (only fMRI)
% opt_fmri.basic = basic;
% 
% 
% preproc.run = 1;
%   preproc.stc           = 1 (event-related design) or 0 (block design)
%   preproc.norm.run      = 1;
%   preproc.norm.struc_vx = [-1 -1 -1];
%   preproc.norm.func_vx  = [-1 -1 -1];
%   preproc.smooth.run    = 1;
%   preproc.smooth.fwhm   = [-2 -2 -2];
% opt_safe.preproc = preproc;
% 
% 
% stat_des.run = 1;
%   stat_des.pref_func      = 'swa' (event-related design) or 'sw' (block design)
%   stat_des.pref_struc     = 'wm';
%   stat_des.tag_struc      = '';
%   stat_des.use_rp         = 1;
%   stat_des.hpf            = 128;
%   stat_des.bf             = 'can_hrf';
%   stat_des.bf_length      = 32;
%   stat_des.bf_order       = 3;
%   stat_des.delay          = (-14:2:4) or (-12:2:6) if bf='gamma'
% opt_safe.stat_des = stat_des;
% 
% 
% opt_safe.stat_con = 1;
% 
% 
% stat_res.run = 1;
% 
% stat_res.thresh.adj = 'none';
% stat_res.thresh.P   = 0.001 (thresh.adj='none') or 0.05 (thresh.adj='FWE')
% stat_res.thresh.ET  = 0;
% 
% stat_res.ps_save = 1;
% 
% stat_res.save_slover = 0;
% stat_res.save_thr = 1;
% 
% stat_res.slices.save = 1;
%   stat_res.slices.img_ext = 'png';
%   stat_res.slices.step    = 3;
%   stat_res.slices.number  = 46;
%   stat_res.slices.slices  = --
%   stat_res.slices.orient  = {'axial'};
% 
% stat_res.group.save = 1;
%   stat_res.group.scale.adj   = 1;
%   stat_res.group.scale.ftest = 'max';
%   stat_res.group.scale.ttest = 'max';
%   stat_res.group.number      = 3;
% 
% stat_res.event_seq.save = 1;
%   stat_res.event_seq.scale.ftest   = 'max';
%   stat_res.event_seq.scale.posbold = 'max';
%   stat_res.event_seq.scale.negbold = 'max';
% 
% stat_res.mask.use = 0;
%   stat_res.mask.imgs_text  = 'mwc*';
%   stat_res.mask.incl_excl  = 'inclusive';
%   stat_res.mask.bin        = 1;
%   stat_res.mask.bin_thresh = 0.2;
% 
% stat_res.cleanup = 1;
% 
% opt_safe.stat_res = stat_res;
% 
% 
% EXAMPLE
%==========================================================================
% basic.parent_dir   = 'D:\Study\fMRI\SubjsData';
% basic.fmri_path    = fullfile('fMRI','Analise','SPM8');
% basic.xlsfile      = {'D:\Study\fMRI\SubjsData\fMRIcond1.xls';'D:\Study\fMRI\SubjsData\fMRIcond2.xls'};
% basic.worksheet    = {{'C01a';'C01b'},{'C100';'C101';'C102'}};
% basic.sbj_folder   = {{'C01-Abc';'C01-Abc'},{'C100-Def';'C101-Ghi';'C102-Jkl'}};
% basic.eeg_fmri_exp = 0;
% opt_safe.basic = basic;
% 
% 
% preproc.run = 1;
%   preproc.stc           = 1;
%   preproc.norm.run      = 1;
%   preproc.norm.struc_vx = [-1 -1 -1];
%   preproc.norm.func_vx  = [-1 -1 -1];
%   preproc.smooth.run    = 1;
%   preproc.smooth.fwhm   = [-2 -2 -2];
% opt_safe.preproc = preproc;
% 
% 
% stat_des.run = 1;
%   stat_des.pref_func  = 'swa'
%   stat_des.pref_struc = 'wm';
%   stat_des.tag_struc  = '';
%   stat_des.use_rp     = 1;
%   stat_des.hpf        = 128;
%   stat_des.bf         = 'gamma';
%   stat_des.bf_length  = 32;
%   stat_des.bf_order   = 3;
%   stat_des.delay      = (0:2:720);
% opt_safe.stat_des = stat_des;
% 
% 
% opt_safe.stat_con = 1;
% 
% 
% stat_res.run = 1;
% 
% stat_res.thresh.adj = 'none';
% stat_res.thresh.P   = 0.001 (thresh.adj='none') or 0.05 (thresh.adj='FWE')
% stat_res.thresh.ET  = 0;
% 
% stat_res.ps_save = 1;
% 
% stat_res.save_slover = 0;
% stat_res.save_thr = 1;
% 
% stat_res.slices.save = 1;
%   stat_res.slices.img_ext = 'png';
%   stat_res.slices.step    = 3;
%   stat_res.slices.number  = 46;
%   stat_res.slices.slices  = --
%   stat_res.slices.orient  = {'axial'};
% 
% stat_res.group.save = 0;
%   stat_res.group.scale.adj   = 1;
%   stat_res.group.scale.ftest = 'max';
%   stat_res.group.scale.ttest = 'max';
%   stat_res.group.number      = 3;
% 
% stat_res.event_seq.save = 1;
%   stat_res.event_seq.scale.ftest   = 'max';
%   stat_res.event_seq.scale.posbold = 'max';
%   stat_res.event_seq.scale.negbold = 'max';
% 
% stat_res.mask.use = 0;
%   stat_res.mask.imgs_text  = 'mwc*';
%   stat_res.mask.incl_excl  = 'inclusive';
%   stat_res.mask.bin        = 1;
%   stat_res.mask.bin_thresh = 0.2;
% 
% stat_res.cleanup = 1;
% 
% opt_safe.stat_res = stat_res;
% 
% 
% See also SAFE_PREPROC, SAFE_STAT_DESIGN, SAFE_STAT_CONTRASTS, SAFE_STAT_RESULTS
% 
%__________________________________________________________________________
% Copyright (C) 2012 Guilherme Coco Beltramini


% Guilherme Coco Beltramini - 2013-May-14


% Read options
%==========================================================================

    eeg_fmri_exp = opt_safe.basic.eeg_fmri_exp;
    
% 1) Parent directory for all subjects
%-------------------------------------
    parent_dir = opt_safe.basic.parent_dir;
    fmri_path  = opt_safe.basic.fmri_path;
    
% 2) XLS files & Worksheets
%--------------------------
    xlsfile   = opt_safe.basic.xlsfile(:);
    worksheet = opt_safe.basic.worksheet(:);
    
% 3) Subject folders
%-------------------
    sbj_folder = opt_safe.basic.sbj_folder(:);
    
% 4) Preprocessing
%-----------------
    preproc = opt_safe.preproc.run;
        stc          = opt_safe.preproc.stc;
        img_norm.run = opt_safe.preproc.norm.run;
        if img_norm.run
            img_norm.struc_vx = opt_safe.preproc.norm.struc_vx;
            img_norm.func_vx  = opt_safe.preproc.norm.func_vx;
        end
        img_smooth.run = opt_safe.preproc.smooth.run;
        if img_smooth.run
            img_smooth.fwhm = opt_safe.preproc.smooth.fwhm;
        end
        
% 5) Statistical analysis
%------------------------
    
    % 5.1) Design
    %------------
    stat_des = opt_safe.stat_des.run;
        pref_func  = opt_safe.stat_des.pref_func;
        pref_struc = opt_safe.stat_des.pref_struc;
        tag_struc  = opt_safe.stat_des.tag_struc;
        use_rp     = opt_safe.stat_des.use_rp;
        hpf        = opt_safe.stat_des.hpf;
        bf         = opt_safe.stat_des.bf;
        bf_length  = opt_safe.stat_des.bf_length;
        bf_order   = opt_safe.stat_des.bf_order;
        delay      = sort(opt_safe.stat_des.delay);
    
    % 5.2) Contrasts
    %---------------
    stat_con = opt_safe.stat_con;
    
    % 5.3) Results
    %-------------
    
    stat_res = opt_safe.stat_res.run;
    
    adj = opt_safe.stat_res.thresh.adj;
    P   = opt_safe.stat_res.thresh.P;
    ET  = opt_safe.stat_res.thresh.ET;
    
    ps_save = opt_safe.stat_res.ps_save;
    
    save_slover = opt_safe.stat_res.save_slover;
    save_thr    = opt_safe.stat_res.save_thr;
    thr_pref    = 'Thresh_';
    
    slices.save = opt_safe.stat_res.slices.save;
        slices.img_ext = opt_safe.stat_res.slices.img_ext;
        slices.step    = opt_safe.stat_res.slices.step;
        slices.number  = opt_safe.stat_res.slices.number;
        slices.slices  = opt_safe.stat_res.slices.slices;
        slices.orient  = opt_safe.stat_res.slices.orient;
    
    group.save = opt_safe.stat_res.group.save;
        if group.save, save_thr=1; end
        group.img_ext  = slices.img_ext;
        group.slc_step = slices.step;
        group.slc_nmbr = slices.number;
        group.slices   = slices.slices;
        group.orient   = slices.orient;
        group.scale    = opt_safe.stat_res.group.scale; % fields: adj, ftest, ttest
        group.number   = opt_safe.stat_res.group.number;
        group.thr_pref = thr_pref;
    
    event_seq.save = opt_safe.stat_res.event_seq.save;
    if event_seq.save
        event_seq.img_ext       = slices.img_ext;
        event_seq.thr_pref      = thr_pref;
        event_seq.scale.ftest   = opt_safe.stat_res.event_seq.scale.ftest;
        event_seq.scale.posbold = opt_safe.stat_res.event_seq.scale.posbold;
        event_seq.scale.negbold = opt_safe.stat_res.event_seq.scale.negbold;
    end    
    
    movie.save = 0; %opt_safe.stat_res.movie.save;
        if movie.save
            movie.sess          = opt_safe.stat_res.movie.sess(:)';
            movie.eventdur      = opt_safe.stat_res.movie.eventdur;
            movie.movie_ti      = opt_safe.stat_res.movie.movie_ti;
            movie.movie_tf      = opt_safe.stat_res.movie.movie_tf;
            movie.write_event   = opt_safe.stat_res.movie.write_event;
            movie.framedur      = opt_safe.stat_res.movie.framedur;
            movie.scale.adjust  = opt_safe.stat_res.movie.scale.adjust;
            movie.scale.Ftest   = opt_safe.stat_res.movie.scale.ftest;
            movie.scale.posBOLD = opt_safe.stat_res.movie.scale.posbold;
            movie.scale.negBOLD = opt_safe.stat_res.movie.scale.negbold;
            movie.orient        = slices.orient;
        else
            movie.sess     = 1;
            movie.eventdur = 6;
            movie.movie_ti = [];
            movie.movie_tf = [];
        end
    
    mask = 0; %opt_safe.mask.use;
    if mask
        mask_imgs_text = {'mwc1*','mwc2*'}; %opt_safe.mask.imgs_text;
        mask_incl_excl = 'inclusive';       %opt_safe.mask.incl_excl;
        mask_bin       = 0;                 %opt_safe.mask.bin;
        if mask_bin
            mask_bin_thresh = 0.5;          %opt_safe.mask.bin_thresh;
        end
    end
    
    cleanup = opt_safe.stat_res.cleanup;
    
    
    
%==========================================================================
%                  PREPROCESSING & STATISTICAL ANALYSIS
%==========================================================================

%curr_dir = pwd;
prev_dir = ' dummy previous directory';
moviedir = '';
spm_jobman('initcfg'); % Initialise jobs configuration and set MATLAB path accordingly

save_slover = movie.save || event_seq.save || save_slover;
if stat_res
    stat_res_part1 = ps_save || slices.save || save_thr || save_slover;
    stat_res_part2 = event_seq.save || group.save || movie.save;
else
    stat_res_part1 = 0;
    stat_res_part2 = 0;
end

if ~mask
    mask_suffix = [];
else
    mask_suffix = '_mask';
end


% Loop for all XLS files
%==========================================================================
for xls=1:length(xlsfile)

if isempty(xlsfile{xls})
    continue
end


% Loop for all worksheets
%==========================================================================
for ww=1:length(worksheet{xls})

if isempty(worksheet{xls}{ww})
    continue
end

subject_dir = sbj_folder{xls}{ww};
fpath = fullfile(parent_dir,subject_dir,fmri_path);


%==========================================================================
% Preprocessing
%==========================================================================
if preproc && ~strcmp(prev_dir,subject_dir) % the same data set will not be preprocessed more than once
    opt_pre.xlsfile   = xlsfile{xls};
    opt_pre.worksheet = worksheet{xls}{ww};
    opt_pre.path      = fpath;
    opt_pre.stc       = stc;
    opt_pre.norm      = img_norm;
    opt_pre.smooth    = img_smooth;
    safe_preproc(opt_pre,'run')
end


% Find the preprocessed structural image
%---------------------------------------
if ~stat_des
    
    if ps_save || save_slover || slices.save || event_seq.save || group.save || movie.save
        
        % Get all file names
        %underlay = dir([fpath filesep pref_struc '*' tag_struc]); this
        %does not guarantee that a file is returned
        underlay = dir(fpath);
        underlay = {underlay(~[underlay.isdir]).name}';
        
        % Find file
        tmp = regexp(underlay,['^' pref_struc '\w*' tag_struc]);
        tmp = ~cellfun(@isempty,tmp);
        
        % Get file name
        underlay = underlay(find(tmp,1,'first'));
        if ~isempty(underlay)
            fprintf('\nUsing %s as underlay image\n',underlay{1})
            underlay = fullfile(fpath,underlay{1});
        else
            error('Underlay file could not be found')
        end
        
    elseif save_thr
        underlay = '';
    end
    
end

% Statistical analysis
%=====================
if stat_des || stat_con || stat_res
	
    % Loop for the sessions
    %======================================================================
    for ss=1:length(movie.sess)
        
        % Output folders
        %---------------
        if ~movie.save
            path_stat = fpath;
            tmp       = sprintf('%s-', worksheet{xls}{ww});
        else
            path_stat = fullfile(fpath, ...
                sprintf('%s-Movie-Sess%.2d', worksheet{xls}{ww}, movie.sess(ss)));
            tmp       = '';
        end
        if length(delay)>1
            path_res = fullfile(path_stat, ...
                sprintf('%sResults_%.1fs_%.1fs%s', tmp, delay(1), delay(end), mask_suffix));
        else % only one delay
            path_res = fullfile(path_stat, ...
                sprintf('%sResults_%.1fs%s', tmp, delay, mask_suffix));
        end
        if movie.save
            moviedir = fullfile(path_res, 'Movie');
        end
        
        safe_res = cell(length(delay),1);
        safe_res = struct('thr',safe_res, 'mn_mx',safe_res, ...
            'obj',safe_res, 'movie_file',safe_res);
        delay_skipped = [];
        
    % Loop for the delays
    %======================================================================
    for dd=1:length(delay)
        
        SPM = fullfile(path_stat, ...
            sprintf('%s-SPM_delay_%.1f', worksheet{xls}{ww}, delay(dd)), ...
            'SPM.mat');
        
        %==================================================================
        % Design
        %==================================================================
        if stat_des
            opt_des.xlsfile       = xlsfile{xls};
            opt_des.worksheet     = worksheet{xls}{ww};
            opt_des.path          = fpath;
            opt_des.eeg           = eeg_fmri_exp;
            opt_des.delay         = delay(dd);
            opt_des.pref_func     = pref_func;
            opt_des.pref_struc    = pref_struc;
            opt_des.use_rp        = use_rp;
            opt_des.hpf           = hpf;
            opt_des.bf            = bf;
            opt_des.bf_length     = bf_length;
            opt_des.bf_order      = bf_order;
            opt_des.outdir        = path_stat;
            opt_des.moviesave     = movie.save;
            opt_des.moviesess     = movie.sess(ss);
            opt_des.movieeventdur = movie.eventdur;
            opt_des.movie_ti      = movie.movie_ti;
            opt_des.movie_tf      = movie.movie_tf;
            underlay = safe_stat_design(opt_des,'run');
            underlay = fullfile(fpath,underlay);
        end
        
        %==================================================================
        % Contrasts
        %==================================================================
        if stat_con
            opt_con.sID = worksheet{xls}{ww};
            safe_stat_contrasts(SPM, opt_con, 'run');
        end
        
        %==================================================================
        % Results
        %==================================================================
        if stat_res_part1
            opt_res.P        = P;
            opt_res.ET       = ET;
            opt_res.adj      = adj;
            opt_res.underlay = underlay;
            opt_res.mask     = mask;
            if mask
                mask_imgs_text = cellstr(mask_imgs_text);
                tmp = [];
                for ii=1:length(mask_imgs_text)
                    tmp = [tmp ; ...
                        strcat(select_files(fpath,mask_imgs_text{ii},'','','','files','path'),',1')];
                end
                if mask_bin
                    clear matlabbatch
                    matlabbatch{1}.spm.util.imcalc.input          = tmp;
                    matlabbatch{1}.spm.util.imcalc.output         = 'mask_for_results.nii';
                    matlabbatch{1}.spm.util.imcalc.outdir         = {fpath};
                    tmp = [];
                    for ii=1:size(matlabbatch{1}.spm.util.imcalc.input,1)
                        tmp = [tmp sprintf('i%d+',ii)];
                    end
                    tmp(end) = [];
                    tmp = [tmp '>' num2str(mask_bin_thresh)];
                    matlabbatch{1}.spm.util.imcalc.expression     = tmp;
                    matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0; % No - don't read images into data matrix
                    matlabbatch{1}.spm.util.imcalc.options.mask   = 0; % No implicit zero mask
                    matlabbatch{1}.spm.util.imcalc.options.interp = 1; % Trilinear
                    matlabbatch{1}.spm.util.imcalc.options.dtype  = 4; % INT16 - signed short
                    spm_jobman('run',matlabbatch)
                    opt_res.mask_imgs = {[fullfile(fpath,matlabbatch{1}.spm.util.imcalc.output) ',1']};
                else
                    opt_res.mask_imgs = tmp;
                end
                opt_res.mask_incl_excl = mask_incl_excl;
            end
            opt_res.save_slover = save_slover;
            opt_res.save_ps     = ps_save;
            opt_res.save_thr    = save_thr;
            opt_res.thr_pref    = thr_pref;
            opt_res.save_slc    = slices.save;
            opt_res.img_ext     = slices.img_ext;
            opt_res.save_movie  = movie.save;
            opt_res.moviedir    = moviedir;
            opt_res.slc_step    = slices.step;
            opt_res.slc_nmbr    = slices.number;
            opt_res.slc_slices  = slices.slices;
            opt_res.orient      = slices.orient;
            try
                safe_res(dd) = safe_stat_results(SPM, opt_res);
            catch
                safe_res(dd) = struct('thr',[], 'mn_mx',[], 'obj',[], 'movie_file',[]);
            end
            if isequal(safe_res(dd),struct('thr',[], 'mn_mx',[], 'obj',[], 'movie_file',[]))
                delay_skipped = [delay_skipped dd];
            end
        end
        
    end % loop for the delays
    
    
    % Keep only the useful delays
    %----------------------------
    delay_orig              = delay;
    delay(delay_skipped)    = [];
    safe_res(delay_skipped) = [];
    
    
    % Save results
    %=============
    if stat_res_part1 && ~isempty(delay)
        
        create_folder(path_res)
        tmp = fullfile(path_res,...
            sprintf('%s-slover_obj-%ddelays%s.mat',worksheet{xls}{ww},length(delay_orig),mask_suffix));
        save(tmp,'safe_res')
        fprintf('Written:\n')
        fprintf('    %s\n',tmp)
        
    end
    
    
    if stat_res_part2 && length(delay)>1
        
        % Load "safe_res" variable
        %-------------------------
        if ~stat_res_part1
            safe_res = fullfile(path_res,...
                sprintf('%s-slover_obj-%ddelays%s.mat',worksheet{xls}{ww},length(delay_orig),mask_suffix));
            load(safe_res)
            if wait_for_existence('safe_res','var',.2,5)
                warning('Could not load %s',safe_res)
                continue
            end
        end
        
        %==================================================================
        % Create event sequence
        %==================================================================
        if event_seq.save
            
            % Create folder
            %--------------
            %out_dir = fullfile(path_res,sprintf('%s-Event_sequence%s',worksheet{xls}{ww},mask_suffix));
            out_dir = fullfile(path_res,sprintf('Event_sequence%s',mask_suffix));
            create_folder(out_dir)
            
            % Sequence of images
            %-------------------
            opt_event_seq.path          = out_dir;
            opt_event_seq.img_ext       = event_seq.img_ext;
            opt_event_seq.thr_pref      = event_seq.thr_pref;
            opt_event_seq.scale.ftest   = event_seq.scale.ftest;
            opt_event_seq.scale.posbold = event_seq.scale.posbold;
            opt_event_seq.scale.negbold = event_seq.scale.negbold;
            safe_event_seq(safe_res,opt_event_seq)
            
            % Plot statistics for delays
            %---------------------------
            opt_plot_delays.delay    = delay;
            opt_plot_delays.safe_res = safe_res;
            opt_plot_delays.path_res = out_dir;
            opt_plot_delays.thr_pref = thr_pref;
            safe_plot_delays(opt_plot_delays)
            
        end
        
        %==================================================================
        % Combine thresholded images
        %==================================================================
        if group.save
            
            % Create folder
            %--------------
            %out_dir = fullfile(path_res,sprintf('%s-Grouped_delays_%.1fs_%.1fs',worksheet{xls}{ww},delay(1),delay(end)));
            %out_dir = fullfile(path_res,sprintf('%s-Grouped_delays%s',worksheet{xls}{ww},mask_suffix));
            out_dir = fullfile(path_res,sprintf('Grouped_delays%s',mask_suffix));
            create_folder(out_dir)
            
            % Group thresholded images
            %-------------------------
            if ~movie.save
                tmp = '';
            else
                tmp = sprintf('%s-Movie-Sess%.2d',worksheet{xls}{ww},movie.sess(ss));
            end
            tmp = fullfile(fpath,tmp);
            opt_group.path       = tmp;
            opt_group.sID        = worksheet{xls}{ww};
            opt_group.delay      = delay;
            opt_group.Ngroup     = group.number;
            opt_group.thr_pref   = group.thr_pref;
            opt_group.thr_dir    = out_dir;
            opt_group.save_slc   = group.save;
            opt_group.img_ext    = group.img_ext;
            opt_group.slc_step   = group.slc_step;
            opt_group.slc_nmbr   = group.slc_nmbr;
            opt_group.slc_slices = group.slices;
            opt_group.orient     = group.orient;
            opt_group.scale      = group.scale;
            opt_group.underlay   = underlay;
            safe_res_group = safe_groupthresh(opt_group);
            
            tmp = fullfile(out_dir,...
                sprintf('%s-stats-%dgroups%s.mat',worksheet{xls}{ww},group.number,mask_suffix));
            save(tmp,'safe_res_group')
            
            % Plot statistics for groups
            %---------------------------
            opt_plot_delays.delay     = [];
            opt_plot_delays.safe_res  = safe_res_group;
            opt_plot_delays.path_res  = out_dir;
            opt_plot_delays.thr_pref  = thr_pref;
            safe_plot_delays(opt_plot_delays)
            
        end
        
        %==================================================================
        % Create movie
        %==================================================================
        if movie.save
            create_folder(moviedir)
            opt_safe_create_movie.delay      = delay;
            opt_safe_create_movie.xlsfile    = xlsfile{xls};
            opt_safe_create_movie.event_name = movie.write_event;
            if movie.write_event
                opt_safe_create_movie.worksheet = worksheet{xls}{ww};
            end
            opt_safe_create_movie.session    = movie.sess(ss);
            opt_safe_create_movie.framedur   = movie.framedur;
            opt_safe_create_movie.scale      = movie.scale;
            opt_safe_create_movie.num_orient = length(movie.orient);
            opt_safe_create_movie.eeg        = eeg_fmri_exp;
            safe_create_movie(safe_res,opt_safe_create_movie)
        end
        
        
    elseif length(delay)<2 && (event_seq.save || movie.save || group.save)
        disp(' ')
        disp('There must be at least 2 delays to combine the images')
    end
    
    
    % Remove auxiliary files
    %=======================
    if cleanup && stat_res
        
        delete_thresh = 1;
        warning off MATLAB:DELETE:FileNotFound % disable warning
        
        % Load "safe_res" variable
        %-------------------------
        if ~ps_save && ~slices.save && ~group.save && ~event_seq.save && ~movie.save
            safe_res = fullfile(path_res,...
                sprintf('%s-slover_obj-%ddelays%s.mat',worksheet{xls}{ww},length(delay_orig),mask_suffix));
            if exist(safe_res,'file')==2
                load(safe_res)
                if wait_for_existence('safe_res','var',.2,5)
                    warning('Could not load %s',safe_res)
                    continue
                end
            else
                delete_thresh = 0;
            end
        end
        
        % Delete thresholded images in the DELAY folders
        %-----------------------------------------------
        if delete_thresh
            for dd=1:length(delay)
                if ~isfield(safe_res(dd).thr,'fname')
                    continue
                end
                for ff=1:length(safe_res(dd).thr.fname)
                    try
                        delete([safe_res(dd).thr.fname{ff} '*'])
                    catch
                        continue
                    end
                end
            end
            pause(1)
        end
        
        % Clean up results folder
        %------------------------
        tmp = fullfile(path_res,...
            sprintf('%s-slover_obj-%ddelays%s.mat',worksheet{xls}{ww},length(delay_orig),mask_suffix));
        if exist(tmp,'file')==2
            delete(tmp)
        end
        if length(dir(path_res))==2 % only '.' and '..'
            rmdir(path_res);
        else
            tmp = fullfile(path_res,'Grouped_delays',...
                sprintf('%s-stats-%dgroups%s.mat',worksheet{xls}{ww},group.number,mask_suffix));
            if exist(tmp,'file')==2 % get the file names and delete the files
                load(tmp)
                if wait_for_existence('safe_res_group','var',.2,5)
                    warning('Could not load %s',tmp)
                    continue
                end
                for gg=1:length(safe_res_group)
                    for ff=1:length(safe_res_group(gg).thr.fname)
                        delete(safe_res_group(gg).thr.fname{ff})
                    end
                end
                delete(tmp)
            else % find the thresholded files and delete them
                ff  = dir([path_res filesep 'Grouped_delays']);
                ff  = {ff(~[ff.isdir]).name}';
                tmp = regexp(ff,[thr_pref '*']);
                tmp = ~cellfun(@isempty,tmp);
                ff  = ff(tmp);
                ff  = strcat([path_res filesep 'Grouped_delays' filesep],ff);
                for tmp=1:length(ff)
                    delete(ff{tmp})
                end
            end
        end
        
        warning on MATLAB:DELETE:FileNotFound % restore warning state
        
        disp('Auxiliary files were removed')
        disp('-------------------------------------------------------------------')
        
    end
    
    end % loop for the sessions
    
end
prev_dir = subject_dir;
end % loop for the worksheets
end % loop for the XLS files

%cd(curr_dir)