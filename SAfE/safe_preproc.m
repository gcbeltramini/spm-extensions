function safe_preproc(opt_pre, run)
%
% SAFE_PREPROC preprocesses the fMRI raw files of one subject.
% 
% 
%USAGE
%-----
% safe_preproc(opt_pre, run)
% 
% 
%INPUT
%-----
% - OPT_PRE: structure with the following fields:
%   - XLSFILE  : XLS file name
%   - WORKSHEET: worksheet in the XLS file (must follow the correct pattern)
%   - PATH     : folder where the raw structural functional files are
%   - STC      : 1 (with slice timing correction) or 0 (without slice
%                timing correction)
%   - NORM: structure with the following fields:
%     - RUN: 1 (normalize) or 0 (don't normalize)
%     If RUN=1:
%     - STRUC_VX: voxel size of the normalized structural image
%       If negative numbers are used in some dimension, to get the voxel
%       size of the normalized structural image, their opposite sign, that
%       is, a positive number, is multiplied by the voxel size of the
%       structural image in that dimension. See the examples.
%       Default: [-1 -1 -1]
%     - FUNC_VX : voxel size of the normalized functional image
%       Analogous. Default: [-1 -1 -1]
%   - SMOOTH: structure with the following fields:
%     - RUN: 1 (smooth) or 0 (do not smooth)
%     If RUN=1:
%     - FWHM: FWHM of the Gaussian kernel
%     Analogous. Zero is allowed. Default: [-2 -2 -2]
% - RUN: 'run' or 'norun'
% 
%IMPORTANT
% - All the raw files must be in the same folder
% - Slice timing correction:
%   - Mandatory when DCM is performed; otherwise, it is optional
%   - If slice order is not acending, change the following parameters in
%     the script: matlabbatch{step}.spm.temporal.st.so
%                 matlabbatch{step}.spm.temporal.st.refslice
%   - Recommendation: for interleaved acquisitions or if subject movement
%     is moderate, slice time correct first; for sequential acquisitions or
%     if there is pronounced inter-slice movements, realign first. By
%     default, in SAfE: slice timing is done after realignment.
% - The SPM8 default values are always used except where "Customized" is
%   written in the code
% 
% 
%OUTPUT
%------
% - WORKSHEET-preproc.mat in the same folder as the raw files
% - If RUN='run', in the same directory:
%   - Realigment:
%     - rp_SUBJECT_NAME_SCAN_X.txt for each functional session (run)
%     - SUBJECT_NAME_SCAN_X.mat for each session
%     - meanSUBJECT_NAME_SCAN_X.nii
%     - WORKSHEET-realign-coreg.ps
%   - Slice timing correction (STC=1):
%     - aSUBJECT_NAME_SCAN_X.nii for each session
%   - Segmentation (NORM.RUN=1):
%     - SUBJECT_NAME_T1W3D_X_seg_sn.mat
%     - SUBJECT_NAME_T1W3D_X_seg_inv_sn.mat
%     - wc1SUBJECT_NAME_T1W3D_X.nii
%     - mwc1SUBJECT_NAME_T1W3D_X.nii
%     - wc2SUBJECT_NAME_T1W3D_X.nii
%     - mwc2SUBJECT_NAME_T1W3D_X.nii
%     - wc3SUBJECT_NAME_T1W3D_X.nii
%     - mwc3SUBJECT_NAME_T1W3D_X.nii
%     - mSUBJECT_NAME_T1W3D_X.nii
%   - Normalization (NORM.RUN=1):
%     - waSUBJECT_NAME_SCAN_X.nii for each session (STC=1)
%     -  wSUBJECT_NAME_SCAN_X.nii for each session (STC=0)
%     - wmeanSUBJECT_NAME_SCAN_X.nii
%     - wmSUBJECT_NAME_T1W3D_X.nii
%   - Smoothing (SMOOTH.RUN=1):
%     - swaSUBJECT_NAME_SCAN_X.nii for each session (STC=1, NORM.RUN=1)
%     -  swSUBJECT_NAME_SCAN_X.nii for each session (STC=0, NORM.RUN=1)
%     - sarSUBJECT_NAME_SCAN_X.nii for each session (STC=1, NORM.RUN=0)
%     -  srSUBJECT_NAME_SCAN_X.nii for each session (STC=0, NORM.RUN=0)
% 
% 
%EXAMPLES
%--------
% 1)
% opt_pre = struct(...
%     'xlsfile'      ,'D:\Study\fMRI\SubjsData\fMRI-conditions.xls',...
%     'worksheet'    ,'Motor01',...
%     'path'         ,'D:\Study\fMRI\SubjsData\Sbj01',...
%     'stc'          ,1,...
%     'norm.run'     ,1,...
%     'norm.struc_vx',[2 2 -3],...
%     'norm.func_vx' ,[3 3 5],...
%     'smooth.run'   ,1,...
%     'smooth.fwhm'  ,[8 -2 0]);
% safe_preproc(opt_pre,'run')
% 
% In this case, the voxel size of the normalized structural image will be:
% [2 2 3x(structural voxel size in the Z-dimension)]
% The FWHM will be: [8 2x(functional voxel size in the Y-dimension) 0]
% 
% 2)
% opt_pre.xlsfile    = 'D:\Study\EEG-fMRI\SubjsData\fMRI-conditions.xls';
% opt_pre.worksheet  = 'P09';
% opt_pre.path       = 'D:\Study\EEG-fMRI\SubjsData\Sbj03';
% opt_pre.stc        = 0;
% opt_pre.norm.run   = 0;
% opt_pre.smooth.run = 0;
% safe_preproc(opt_pre,'norun')
% 
% See also SAFE_STAT_DESIGN, SAFE_STAT_CONTRASTS, SAFE_STAT_RESULTS
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2013-Feb-20


% Dependency
%--------------------------------------------------------------------------
% 1 (use dependencies) or 0 (don't use dependencies)
% shouldn't matter which option you choose
dep = 1;


% Read options
%--------------------------------------------------------------------------
xlsfile    = opt_pre.xlsfile;
worksheet  = opt_pre.worksheet;
fpath      = opt_pre.path;
stc        = opt_pre.stc;
norm_run   = opt_pre.norm.run;
smooth_run = opt_pre.smooth.run;

if norm_run
    tmp               = opt_pre.norm.struc_vx>0;
    norm_struc_vx     = opt_pre.norm.struc_vx.*tmp;     % get the positive numbers
    norm_struc_factor = -opt_pre.norm.struc_vx.*(~tmp); % get the non-positive numbers
    % This way, the voxel size of the normalized image is:
    % norm_struc_vx + norm_struc_factor.*struc_vox
    
    tmp              = opt_pre.norm.func_vx>0;
    norm_func_vx     = opt_pre.norm.func_vx.*tmp;     % get the positive numbers
    norm_func_factor = -opt_pre.norm.func_vx.*(~tmp); % get the non-positive numbers
end

if smooth_run
    tmp           = opt_pre.smooth.fwhm>0;
    smooth_fwhm   = opt_pre.smooth.fwhm.*tmp;     % get the non-negative numbers
    smooth_factor = -opt_pre.smooth.fwhm.*(~tmp); % get the negative numbers
    % This way, the FWHM is:
    % smooth_fwhm + smooth_factor.*func_vox
end


% Close SPM graphics window
%--------------------------------------------------------------------------
spm_figure('Close','Graphics')


% "Run" option
%--------------------------------------------------------------------------
switch run
    case 'run'
        run = 1;
    case 'norun'
        run = 0;
    otherwise
        error('Invalid "run" option')
end


tmp = clock;
disp('-------------------------------------------------------------------')
fprintf('(%.2d:%.2d:%02.0f) Running %s.m...\n',...
    tmp(4),tmp(5),tmp(6),mfilename)


% Read the XLS file: TR, structural and functional scans
%==========================================================================
[xlsnum,xlstxt] = xlsread(xlsfile,worksheet);

% Remove empty columns and rows in the beginning of XLSNUM and XLSTXT
%--------------------------------------------------------------------
[xlsnum,xlstxt] = trim_array('beginning',xlsnum,xlstxt);

% Repetition time
%----------------
% The TR appears before two NaN's in the first column of "xlsnum" (for
% fMRI or EEG-fMRI, with the subject name being a number or not)
% for tmp=3:size(xlsnum,1)
%     if sum(isnan(xlsnum(tmp-1:tmp,1)),1)==2
%         break
%     end
% end
% TR = xlsnum(tmp-2,1);

if isnan(xlsnum(3,1)) % subject ID is not a number => the ID is xlstxt{1,2}
    TR = xlsnum(1,1);
else % subject ID is a number => xlstxt{1,2} is empty, and the ID is xlsnum(1,1)
    % in Linux, xlstxt{1,1} is also empty => create empty row
    if ~isempty(xlstxt{1,2})
        xlstxt = [repmat({''},1,size(xlstxt,2)) ; xlstxt];
    end
    TR = xlsnum(3,1);
end

% Structural scan
%----------------
struc_scan = fullfile(fpath,xlstxt{2,2});
if exist(struc_scan,'file')~=2
    error('Structural scan "%s" not found',struc_scan)
end
tmp        = spm_read_hdr(struc_scan);
struc_vox  = round(tmp.dime.pixdim(2:4)*10)/10;


% Functional scans - select files for Realign: Estimate & Reslice
%==========================================================================

new_sess  = find(strcmp(xlstxt(:,1),'fMRI:'));
sess      = size(new_sess,1); % number of sessions
func_scan = cell(1,sess);

% Read files and get number of sessions
%--------------------------------------------------------------------------
for s=1:sess
    
    func_scan_name = fullfile(fpath,xlstxt{new_sess(s),2});
    if exist(func_scan_name,'file')~=2
        error('File %s not found',func_scan_name)
    end

    % Get voxel size and check with the other sessions
    %-------------------------------------------------
    hdr = spm_read_hdr(func_scan_name);
    if s==1 % first session
        nslices  = hdr.dime.dim(4);
        func_vox = hdr.dime.pixdim(2:4);
    elseif ~isequal(hdr.dime.dim(4),nslices)
        warning(sprintf('SAfE:%s:DifferentSlices',mfilename),'The functional scans do not have the same number of slices')
        %error('The functional scans must have the same number of slices')
    elseif ~isequal(hdr.dime.pixdim(2:4),func_vox)
        %error('The functional scans must have the same voxel size')
        warning(sprintf('SAfE:%s:DifferentVoxels',mfilename),'The functional scans do not have the same voxel size')
    end

    % Get all the functional images from each session
    %------------------------------------------------
    dyn = hdr.dime.dim(5);
    func_scan{1,s} = cell(dyn,1);
    for tt=1:dyn
        func_scan{1,s}{tt} = [func_scan_name ',' num2str(tt)];
    end
    
end

% Round the voxel size
%---------------------
func_vox = round(func_vox*10)/10;


%==========================================================================
%                          START PREPROCESSING
%==========================================================================

% Change directory (realignment PS file will be saved in the current folder)
%--------------------------------------------------------------------------
old_dir = pwd;
cd(fpath)

clear matlabbatch
step = 0;


% Realign: Estimate & Reslice
%==========================================================================
step = step + 1;

% Data
%--------------------------------------------------------------------------
matlabbatch{step}.spm.spatial.realign.estwrite.data = func_scan;

% Estimation Options
%--------------------------------------------------------------------------
matlabbatch{step}.spm.spatial.realign.estwrite.eoptions.quality = .9;      % quality (0-1)
matlabbatch{step}.spm.spatial.realign.estwrite.eoptions.sep     = 4;       % separation [mm]
matlabbatch{step}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;       % smoothing (FWHM) [mm]
matlabbatch{step}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;       % num passes ("register to first" [0] or "register to mean" [1])
matlabbatch{step}.spm.spatial.realign.estwrite.eoptions.interp  = 2;       % interpolation ("2nd degree B-spline")
matlabbatch{step}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0]; % wrapping ("no wrap")
matlabbatch{step}.spm.spatial.realign.estwrite.eoptions.weight  = '';      % weighting

% Reslice Options
%--------------------------------------------------------------------------
if norm_run % customized
    matlabbatch{step}.spm.spatial.realign.estwrite.roptions.which  = [0 1]; % resliced images ("mean image only")
else
    matlabbatch{step}.spm.spatial.realign.estwrite.roptions.which  = [2 1]; % resliced images ("all images + mean image")
end
matlabbatch{step}.spm.spatial.realign.estwrite.roptions.interp = 4;       % interpolation ("4th degree B-spline")
matlabbatch{step}.spm.spatial.realign.estwrite.roptions.wrap   = [0 0 0]; % wrapping ("no wrap")
matlabbatch{step}.spm.spatial.realign.estwrite.roptions.mask   = 1;       % masking ("mask images" [1] or "dont mask images" [0])
matlabbatch{step}.spm.spatial.realign.estwrite.roptions.prefix = 'r';     % filename prefix

% Customized
%--------------------------------------------------------------------------
matlabbatch{step}.spm.spatial.realign.estwrite.eoptions.wrap = [0 1 0]; % wrap Y
matlabbatch{step}.spm.spatial.realign.estwrite.roptions.wrap = [0 1 0];


% Slice Timing
%==========================================================================
if stc
    
    step = step + 1;
    
    % Data - Sessions (realigned images)
    %----------------------------------------------------------------------
    if dep % dependency for sessions ("Realigned Images" or "Resliced Images")
        
        for s=1:sess
                matlabbatch{step}.spm.temporal.st.scans{s}(1)                      = cfg_dep;
                matlabbatch{step}.spm.temporal.st.scans{s}(1).tname                = 'Session';
                matlabbatch{step}.spm.temporal.st.scans{s}(1).tgt_spec{1}(1).name  = 'filter';
                matlabbatch{step}.spm.temporal.st.scans{s}(1).tgt_spec{1}(1).value = 'image';
                matlabbatch{step}.spm.temporal.st.scans{s}(1).tgt_spec{1}(2).name  = 'strtype';
                matlabbatch{step}.spm.temporal.st.scans{s}(1).tgt_spec{1}(2).value = 'e';
        end
        if norm_run
            for s=1:sess
                matlabbatch{step}.spm.temporal.st.scans{s}(1).sname        = sprintf('Realign: Estimate & Reslice: Realigned Images (Sess %d)',s);
                matlabbatch{step}.spm.temporal.st.scans{s}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{step}.spm.temporal.st.scans{s}(1).src_output   = substruct('.','sess', '()',{s}, '.','cfiles');
            end
        else
            for s=1:sess
                matlabbatch{step}.spm.temporal.st.scans{s}(1).sname        = sprintf('Realign: Estimate & Reslice: Resliced Images (Sess %d)',s);
                matlabbatch{step}.spm.temporal.st.scans{s}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{step}.spm.temporal.st.scans{s}(1).src_output   = substruct('.','sess', '()',{s}, '.','rfiles');
            end
        end
        
    else
        
        if norm_run
            matlabbatch{step}.spm.temporal.st.scans = func_scan;
        else
            matlabbatch{step}.spm.temporal.st.scans = addprefix(func_scan, 'r', 'sess');
        end
        
    end
    
    % Image parameters - customized
    %----------------------------------------------------------------------
    matlabbatch{step}.spm.temporal.st.nslices = nslices;       % number of slices
    matlabbatch{step}.spm.temporal.st.tr      = TR;            % TR
    matlabbatch{step}.spm.temporal.st.ta      = TR-TR/nslices; % TA
    matlabbatch{step}.spm.temporal.st.so      = 1:nslices;     % slice order
    % ascending: 1:nslices
    % descending: nslices:-1:1
    % interleaved bottom->up: [1:2:nslices 2:2:nslices]
    % interleaved top->down: [nslices:-2:1 nslices-1:-2:1]
    % interleaved middle->top - gives weird results:
    %    for k=1:nslices
    %        round( (nslices-k)/2 + rem(nslices-k,2)*(nslices-1)/2 ) + 1
    %    end
    matlabbatch{step}.spm.temporal.st.refslice = ceil(nslices/2); % reference slice
    % ascending/descending: ceil(nslices/2)
    % interleaved bottom->up: nslices-rem(nslices-1,2)
    % interleaved top->down: 1+rem(nslices-1,2)
    matlabbatch{step}.spm.temporal.st.prefix = 'a'; % filename prefix
    
end


% Coregister: Estimate
%==========================================================================
step = step + 1;

% Reference image (mean image)
%--------------------------------------------------------------------------
if ~dep
    [fpath,tmp,ext] = fileparts(func_scan{1,1}{1});
    matlabbatch{step}.spm.spatial.coreg.estimate.ref = {fullfile(fpath,['mean' tmp ext])};
else % dependency ("Reference Image")
    matlabbatch{step}.spm.spatial.coreg.estimate.ref(1)                      = cfg_dep;
    matlabbatch{step}.spm.spatial.coreg.estimate.ref(1).tname                = 'Reference Image';
    matlabbatch{step}.spm.spatial.coreg.estimate.ref(1).tgt_spec{1}(1).name  = 'filter';
    matlabbatch{step}.spm.spatial.coreg.estimate.ref(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{step}.spm.spatial.coreg.estimate.ref(1).tgt_spec{1}(2).name  = 'strtype';
    matlabbatch{step}.spm.spatial.coreg.estimate.ref(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{step}.spm.spatial.coreg.estimate.ref(1).sname                = 'Realign: Estimate & Reslice: Mean Image';
    matlabbatch{step}.spm.spatial.coreg.estimate.ref(1).src_exbranch         = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{step}.spm.spatial.coreg.estimate.ref(1).src_output           = substruct('.','rmean');
end

% Source image
%--------------------------------------------------------------------------
matlabbatch{step}.spm.spatial.coreg.estimate.source = {[struc_scan ',1']};

% Other images
%--------------------------------------------------------------------------
matlabbatch{step}.spm.spatial.coreg.estimate.other = {''};

% Estimation options
%--------------------------------------------------------------------------
matlabbatch{step}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi'; % objective function ("normalised mutual information")
matlabbatch{step}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2]; % separation [mm]
matlabbatch{step}.spm.spatial.coreg.estimate.eoptions.tol      = [.02 .02 .02 .001 .001 .001 .01 .01 .01 .001 .001 .001]; % tolerances
matlabbatch{step}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7]; % histogram smoothing


% Segment
%==========================================================================
if norm_run
    
    step = step + 1;

    % Data (coregistered image)
    %----------------------------------------------------------------------
    if ~dep
        matlabbatch{step}.spm.spatial.preproc.data = matlabbatch{step-1}.spm.spatial.coreg.estimate.source;
    else % dependency ("Coregistered Images")
        matlabbatch{step}.spm.spatial.preproc.data(1)                      = cfg_dep;
        matlabbatch{step}.spm.spatial.preproc.data(1).tname                = 'Data';
        matlabbatch{step}.spm.spatial.preproc.data(1).tgt_spec{1}(1).name  = 'filter';
        matlabbatch{step}.spm.spatial.preproc.data(1).tgt_spec{1}(1).value = 'image';
        matlabbatch{step}.spm.spatial.preproc.data(1).tgt_spec{1}(2).name  = 'strtype';
        matlabbatch{step}.spm.spatial.preproc.data(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{step}.spm.spatial.preproc.data(1).sname                = 'Coreg: Estimate: Coregistered Images';
        matlabbatch{step}.spm.spatial.preproc.data(1).src_exbranch         = substruct('.','val', '{}',{step-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{step}.spm.spatial.preproc.data(1).src_output           = substruct('.','cfiles');
    end

    % Output files
    %----------------------------------------------------------------------
    matlabbatch{step}.spm.spatial.preproc.output.GM      = [0 0 1]; % native space
    matlabbatch{step}.spm.spatial.preproc.output.WM      = [0 0 1];
    matlabbatch{step}.spm.spatial.preproc.output.CSF     = [0 0 0]; % none
    matlabbatch{step}.spm.spatial.preproc.output.biascor = 1;       % save bias corrected
    matlabbatch{step}.spm.spatial.preproc.output.cleanup = 0;       % don't do cleanup

    % Tissue probability maps
    %----------------------------------------------------------------------
    tmp = [fileparts(which('spm')) filesep 'tpm' filesep]; % SPM8 path
    matlabbatch{step}.spm.spatial.preproc.opts.tpm{1,1} = [tmp 'grey.nii'];
    matlabbatch{step}.spm.spatial.preproc.opts.tpm{2,1} = [tmp 'white.nii'];
    matlabbatch{step}.spm.spatial.preproc.opts.tpm{3,1} = [tmp 'csf.nii'];

    matlabbatch{step}.spm.spatial.preproc.opts.ngaus    = [2;2;2;4]; % Gaussians per class
    matlabbatch{step}.spm.spatial.preproc.opts.regtype  = 'mni';     % affine regularisation ("ICBM space template - European brains")
    matlabbatch{step}.spm.spatial.preproc.opts.warpreg  = 1;         % warping regularisation
    matlabbatch{step}.spm.spatial.preproc.opts.warpco   = 25;        % warp frequency cutoff
    matlabbatch{step}.spm.spatial.preproc.opts.biasreg  = .0001;     % bias regularisation ("very light regularisation")
    matlabbatch{step}.spm.spatial.preproc.opts.biasfwhm = 60;        % bias FWHM
    matlabbatch{step}.spm.spatial.preproc.opts.samp     = 3;         % sampling distance
    matlabbatch{step}.spm.spatial.preproc.opts.msk      = {''};      % masking image

    % Customized
    %----------------------------------------------------------------------
    matlabbatch{step}.spm.spatial.preproc.output.GM  = [1 1 1]; % modulated + unmodulated + native
    matlabbatch{step}.spm.spatial.preproc.output.WM  = [1 1 1];
    matlabbatch{step}.spm.spatial.preproc.output.CSF = [1 1 1];
    
    
    % Normalise: Write - structural image
    %======================================================================
    step = step + 1;

    % Data
    %----------------------------------------------------------------------

    %   Parameter file (subj -> MNI)
    %----------------------------------------------------------------------
    if ~dep
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname = {[struc_scan(1:end-4) '_seg_sn.mat']};
    else % dependency ("Norm Params Subj->MNI"):
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname(1)                      = cfg_dep;
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname(1).tname                = 'Parameter File';
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).name  = 'filter';
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).name  = 'strtype';
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname(1).sname                = 'Segment: Norm Params Subj->MNI';
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname(1).src_exbranch         = substruct('.','val','{}',{step-1},'.','val','{}',{1},'.','val','{}',{1});
        matlabbatch{step}.spm.spatial.normalise.write.subj.matname(1).src_output           = substruct('()',{1}, '.','snfile', '()',{':'});
    end

    %   Images to write (bias corrected image)
    %----------------------------------------------------------------------
    if ~dep
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample = matlabbatch{step-2}.spm.spatial.coreg.estimate.source;
    else % dependency ("Bias Corr Images")
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample(1)                      = cfg_dep;
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample(1).tname                = 'Images to Write';
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).name  = 'filter';
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).value = 'image';
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).name  = 'strtype';
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample(1).sname                = 'Segment: Bias Corr Images';
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample(1).src_exbranch         = substruct('.','val', '{}',{step-1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{step}.spm.spatial.normalise.write.subj.resample(1).src_output           = substruct('()',{1}, '.','biascorr', '()',{':'});
    end

    % Writing options
    %----------------------------------------------------------------------
    matlabbatch{step}.spm.spatial.normalise.write.roptions.preserve = 0;                       % preserve ("preserve concentrations")
    matlabbatch{step}.spm.spatial.normalise.write.roptions.bb       = [-78 -112 -50;78 76 85]; % bounding box
    matlabbatch{step}.spm.spatial.normalise.write.roptions.vox      = [2 2 2];                 % voxel sizes
    matlabbatch{step}.spm.spatial.normalise.write.roptions.interp   = 1;                       % interpolation
    matlabbatch{step}.spm.spatial.normalise.write.roptions.wrap     = [0 0 0];                 % wrapping
    matlabbatch{step}.spm.spatial.normalise.write.roptions.prefix   = 'w';                     % filename prefix

    % Customized
    %----------------------------------------------------------------------
    matlabbatch{step}.spm.spatial.normalise.write.roptions.vox  = norm_struc_vx + norm_struc_factor.*struc_vox;
    matlabbatch{step}.spm.spatial.normalise.write.roptions.wrap = [0 1 0]; % wrap Y


    % Normalise: Write - functional image
    %======================================================================
    step = step + 1;

    % Data
    %----------------------------------------------------------------------

    %   Parameter file (subj -> MNI)
    %----------------------------------------------------------------------
    matlabbatch{step}.spm.spatial.normalise.write.subj.matname = ...
        matlabbatch{step-1}.spm.spatial.normalise.write.subj.matname;

    %   Images to write (bias corrected image)
    %----------------------------------------------------------------------
    if ~dep
        if stc
            matlabbatch{step}.spm.spatial.normalise.write.subj.resample = addprefix(func_scan, 'a');
        else
            matlabbatch{step}.spm.spatial.normalise.write.subj.resample = addprefix(func_scan, '');
        end
    else % dependency for "sess" sessions ("Slice Timing: Slice Timing Corr. Images (Sess XXX)" [with slice timing] or
         % "Realign: Estimate & Reslice: Realigned Images (Sess XXX)" [no slice timing])
        for s=1:sess
            matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s)                      = cfg_dep;
            matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).tname                = 'Images to Write';
            matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).tgt_spec{1}(1).name  = 'filter';
            matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).tgt_spec{1}(1).value = 'image';
            matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).tgt_spec{1}(2).name  = 'strtype';
            matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).tgt_spec{1}(2).value = 'e';
            if stc
                matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).sname        = sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)',s);
                matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).src_exbranch = substruct('.','val','{}',{2},'.','val','{}',{1},'.','val', '{}',{1});
                matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).src_output   = substruct('()',{s},'.','files');
            else
                matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).sname        = sprintf('Realign: Estimate & Reslice: Realigned Images (Sess %d)',s);
                matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.', 'val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{step}.spm.spatial.normalise.write.subj.resample(s).src_output   = substruct('.','sess', '()',{s}, '.','cfiles');
            end
        end
    end

    % Writing options - customized (see above)
    %----------------------------------------------------------------------
    matlabbatch{step}.spm.spatial.normalise.write.roptions = ...
        matlabbatch{step-1}.spm.spatial.normalise.write.roptions;

    % Customized
    %----------------------------------------------------------------------
    matlabbatch{step}.spm.spatial.normalise.write.roptions.vox = norm_func_vx + norm_func_factor.*func_vox;

end


% Smoothing
%==========================================================================
if smooth_run
    
    step = step + 1;
    
    % Images to smooth (normalised functional images)
    %----------------------------------------------------------------------
    if ~dep
        
        if stc
            if norm_run
                matlabbatch{step}.spm.spatial.smooth.data = addprefix(func_scan, 'wa');
            else
                matlabbatch{step}.spm.spatial.smooth.data = addprefix(func_scan, 'ar');
            end
        else
            if norm_run
                matlabbatch{step}.spm.spatial.smooth.data = addprefix(func_scan, 'w');
            else
                matlabbatch{step}.spm.spatial.smooth.data = addprefix(func_scan, 'r');
            end
        end
        
    else
        
        if norm_run % dependency ("Normalised Images (Subj 1)"):
            matlabbatch{step}.spm.spatial.smooth.data(1)                      = cfg_dep;
            matlabbatch{step}.spm.spatial.smooth.data(1).tname                = 'Images to Smooth';
            matlabbatch{step}.spm.spatial.smooth.data(1).tgt_spec{1}(1).name  = 'filter';
            matlabbatch{step}.spm.spatial.smooth.data(1).tgt_spec{1}(1).value = 'image';
            matlabbatch{step}.spm.spatial.smooth.data(1).tgt_spec{1}(2).name  = 'strtype';
            matlabbatch{step}.spm.spatial.smooth.data(1).tgt_spec{1}(2).value = 'e';
            matlabbatch{step}.spm.spatial.smooth.data(1).sname                = 'Normalise: Write: Normalised Images (Subj 1)';
            matlabbatch{step}.spm.spatial.smooth.data(1).src_exbranch         = substruct('.','val', '{}',{step-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            matlabbatch{step}.spm.spatial.smooth.data(1).src_output           = substruct('()',{1}, '.','files');
        else
            for s=1:sess
                matlabbatch{step}.spm.spatial.smooth.data(s)                      = cfg_dep;
                matlabbatch{step}.spm.spatial.smooth.data(s).tname                = 'Images to Smooth';
                matlabbatch{step}.spm.spatial.smooth.data(s).tgt_spec{1}(1).name  = 'filter';
                matlabbatch{step}.spm.spatial.smooth.data(s).tgt_spec{1}(1).value = 'image';
                matlabbatch{step}.spm.spatial.smooth.data(s).tgt_spec{1}(2).name  = 'strtype';
                matlabbatch{step}.spm.spatial.smooth.data(s).tgt_spec{1}(2).value = 'e';
            end
            if stc
                for s=1:sess
                    matlabbatch{step}.spm.spatial.smooth.data(s).sname        = sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)',s);
                    matlabbatch{step}.spm.spatial.smooth.data(s).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
                    matlabbatch{step}.spm.spatial.smooth.data(s).src_output   = substruct('()',{s}, '.','files');
                end
            else
                for s=1:sess
                    matlabbatch{step}.spm.spatial.smooth.data(s).sname        = sprintf('Realign: Estimate & Reslice: Resliced Images (Sess %d)',s);
                    matlabbatch{step}.spm.spatial.smooth.data(s).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
                    matlabbatch{step}.spm.spatial.smooth.data(s).src_output   = substruct('.','sess', '()',{s}, '.','rfiles');
                end
            end
        end
        
    end
    
    matlabbatch{step}.spm.spatial.smooth.fwhm   = [8 8 8]; % FWHM
    matlabbatch{step}.spm.spatial.smooth.dtype  = 0;       % data type ("SAME")
    matlabbatch{step}.spm.spatial.smooth.im     = 0;       % implicit masking
    matlabbatch{step}.spm.spatial.smooth.prefix = 's';     % filename prefix
    
    % Customized
    %----------------------------------------------------------------------
    matlabbatch{step}.spm.spatial.smooth.fwhm = smooth_fwhm + smooth_factor.*func_vox;
    
end


% Save, run and finish
%==========================================================================
batch_file = fullfile(fpath,sprintf('%s-preproc.mat',worksheet));
save(batch_file,'matlabbatch')
fprintf('File %s was created\n',batch_file)

if run
    
    % Show the progress bar and the graphics
    %Fmenu  = spm('CreateMenuWin','off');
    Fmenu  = [];
    Finter = spm('CreateIntWin','off');
    Fgraph = spm_figure('Create','Graphics','Graphics','off');
    set([Fmenu,Finter,Fgraph],'Visible','on');
    
    % Run
    %spm_jobman('initcfg')
    spm_jobman('run',matlabbatch)
    
    % PS file
    s = movefile(sprintf('spm_%s.ps',datestr(now,'yyyymmmdd')),...
        sprintf('%s-realign-coreg.ps',worksheet));
    if ~s % didn't find the file - look for a file from the previous day
        % if the script was running near midnight
        movefile(sprintf('spm_%s.ps',datestr(now-1,'yyyymmmdd')),...
            sprintf('%s-realign-coreg.ps',worksheet));
    end
    
end

cd(old_dir)

tmp = clock;
fprintf('(%.2d:%.2d:%02.0f) %s.m done!\n',...
    tmp(4),tmp(5),tmp(6),mfilename)
disp('-------------------------------------------------------------------')

end


% =========================================================================
%                          AUXILIARY FUNCTIONS
% =========================================================================

function out = addprefix(func_img, prefix, varargin)
%
% ADDPREFIX(FUNC_IMG, PREFIX, OPT) adds a PREFIX to set of functional files
% FUNC_IMG, arranged in a 1xSESS cell array, where SESS is the number of
% sessions. Each cell is then a Nx1 cell array, where N is the number of
% scans in that session.
%
% If OPT is empty, the output is a Sx1 cell array, where S is the total
% number of scans.
%
% If OPT='sess', the cell organization of the input will be preserved.
%

% Guilherme Coco Beltramini - 2012-Jan-31, 03:07 pm

if nargin==2
    
    out = [];
    for s=1:length(func_img) % loop for sessions
        dyn = length(func_img{s});
        tmp = cell(dyn,1);
        for d=1:dyn % loop for scans (time)
            [fpath,file,ext] = fileparts(func_img{s}{d});
            tmp{d} = fullfile(fpath,[prefix file ext]);
        end
        out = [out ; tmp];
    end
    
else
    
    % Keep the cell structure of the input
    %----------------------------------------------------------------------
    out = func_img;
    for s=1:length(func_img) % loop for sessions
        for d=1:length(func_img{s}) % loop for scans (time)
            [fpath,file,ext] = fileparts(func_img{s}{d});
            out{s}{d}        = fullfile(fpath,[prefix file ext]);
        end
    end
    
end

end