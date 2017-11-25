function struc_scan = safe_stat_design(opt_des, run)
%
% SAFE_STAT_DESIGN creates the design matrix for an fMRI experiment.
% 
% 
%USAGE
%-----
% struc_scan = safe_stat_design(opt_des, run)
% 
% 
%INPUT
%-----
% - OPT_DES: structure with the following fields:
%   - XLSFILE   : XLS file name
%   - WORKSHEET : worksheet in the XLS file (must follow the correct pattern)
%   - PATH      : folder where all the preprocessed functional files are
%   - EEG       : 0 (fMRI experiment) or 1 (EEG-fMRI experiment)
%   - PREF_FUNC : prefix for the preprocessed functional files (default: 'swa')
%   - PREF_STRUC: prefix for the structural file (default: 'wm')
%   - DELAY     : time in seconds from the stimulus onset to the beginning
%                 of the hemodynamic response function. DELAY can be any
%                 number (positive, negative or zero)
%   - USE_RP    : 1 (use realignment parameters) or 0 (don't use) (default: 1)
%   - HPF       : high-pass filter [s] (default: 128)
%   - BF        : basis function: 'can_hrf', 'fourier', 'fourier_han',
%                 'gamma' or 'fir' (default: 'can_hrf')
%     - If BF = 'fourier', 'fourier_han', 'gamma', 'fir':
%       BF_LENGTH: post-stimulus window length [s] (default: 32)
%   - BF_ORDER  : number of basis functions (default: 3)
%   - OUTDIR    : output folder
% - RUN: 'run' or 'norun'
% 
% The use of delays was inspired by:
%  Bagshaw AP et al. (2004) "EEG-fMRI of Focal Epileptic Spikes: Analysis
%  With Multiple Haemodynamic Functions and Comparison With
%  Gadolinium-Enhanced MR Angiograms". Human Brain Mapping 22:179-192
% 
% 
%OUTPUT
%------
% - Inside OUTDIR folder WORKSHEET-SPM_delay_DELAY with files:
%   - WORKSHEET-design_matrix.mat
%   - WORKSHEET-rp-conditions.png
%   If RUN='run', in folder WORKSHEET-SPM_delay_DELAY:
%     - beta_XXXX.hdr/img, mask.hdr/img, ResMS.hdr/img, RPV.hdr/img, SPM.mat
% - STRUC_SCAN: name of the preprocessed structural scan file
% 
% N.B.: If there is already a non-empty folder named
%  WORKSHEET-SPM_delay_DELAY, it will be renamed to
%  WORKSHEET-SPM_delay_DELAY-backup (or WORKSHEET-SPM_delay_DELAY-backup02,
%  WORKSHEET-SPM_delay_DELAY-backup03, ...)
% 
% 
%EXAMPLE
%-------
% opt_des.xlsfile       = D:\Study\EEG-fMRI\SubjsData\Conditions.xls;
% opt_des.worksheet     = 'Subj01';
% opt_des.path          = 'D:\Study\EEG-fMRI\SubjsData';
% opt_des.eeg           = 1;
% opt_des.pref_func     = 'swa';
% opt_des.pref_struc    = 'wm';
% opt_des.delay         = -5;
% opt_des.use_rp        = 1;
% opt_des.hpf           = 128;
% opt_des.bf            = 'gamma';
% opt_des.bf_length     = 32;
% opt_des.bf_order      = 3;
% safe_stat_design(opt_des,'run')
% 
% See also SAFE_PREPROC, SAFE_STAT_CONTRASTS, SAFE_STAT_RESULTS
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini


% NOTES:
% - The SPM8 default values are always used except where "Customized" is
%   written in the code
% - Alternative: all preprocessed files must have the same number of slices
%   if the customized "microtime resolution" and "microtime onset" are used
%   (by using the default values, the number of slices can be different)


% Guilherme Coco Beltramini - 2013-Jul-03


% Read options
%==========================================================================
xlsfile    = opt_des.xlsfile;
worksheet  = opt_des.worksheet;
fpath      = opt_des.path;
eeg        = opt_des.eeg;
pref_func  = opt_des.pref_func;
pref_struc = opt_des.pref_struc;
delay      = opt_des.delay;
use_rp     = opt_des.use_rp;
hpf        = opt_des.hpf;
bf         = opt_des.bf;
if isfield(opt_des,'bf_length')
    bf_length = opt_des.bf_length;
end
bf_order   = opt_des.bf_order;
outdir     = opt_des.outdir;

% Movie options
%--------------
moviesave = 0; %opt_des.moviesave;
if moviesave
    moviesess     = opt_des.moviesess;
    movieeventdur = opt_des.movieeventdur;
    movie_ti      = opt_des.movie_ti;
    movie_tf      = opt_des.movie_tf;
end

% "Run" option
%-------------
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


% Read the XLS file
%==========================================================================
[xlsnum,xlstxt] = xlsread(xlsfile,worksheet);

% Remove empty columns and rows in the beginning of XLSNUM and XLSTXT
%--------------------------------------------------------------------
[xlsnum,xlstxt] = trim_array('beginning',xlsnum,xlstxt);

% Structural scan
%----------------
struc_scan = [pref_struc xlstxt{2,2}];

% Repetition time
%----------------
if isnan(xlsnum(3,1)) % subject ID is not a number => the ID is xlstxt{1,2}
    TR = xlsnum(1,1);
else % subject ID is a number => xlstxt{1,2} is empty, and the ID is xlsnum(1,1)
    % in Linux, xlstxt{1,1} is also empty => create empty row
    if ~isempty(xlstxt{1,2})
        xlstxt = [repmat({''},1,size(xlstxt,2)) ; xlstxt];
    end
    TR = xlsnum(3,1);
    xlsnum(1:2,:) = [];
end

% Get only onsets and durations and round to one decimal place
%-------------------------------------------------------------
if eeg
    tmp = 8;
    if size(xlsnum,2)<9
        eeg = 0;
        tmp = 1;
        warning('Index exceeds matrix dimensions in the XLS file. Assuming it is not an EEG-fMRI experiment.')
    end
else
    tmp = 1;
end
xlsnum = round(xlsnum(:,tmp:(tmp+1))*10)/10;

% Session information
%--------------------
new_sess = find(strcmp(xlstxt(:,1),'fMRI:'));

% Number of sessions
if ~moviesave
    sess = size(new_sess,1);
else
    sess = 1;
end

if use_rp
    rp = cell(1,sess);
end
func_scan = cell(1,sess);


% Not saving movie
%==========================================================================
if ~moviesave

s = 1; % counter for the useful sessions
while s<=sess % loop for sessions
    
    tmp = fullfile(fpath,[pref_func xlstxt{new_sess(s),2}]);
    if exist(tmp,'file')~=2
        new_sess(s) = [];
        sess = sess - 1;
        if sess==0
            error('No functional scan was found for spreadsheet %s',worksheet)
        end
        continue
    end
    
    % Number of scans
    %----------------
    dyn = spm_read_hdr(tmp);
    try
        dyn = dyn.dime.dim(5);
    catch ME
        error(ME.identifier, ME.message)
    end
    
    % Functional scans
    %-----------------
    func_scan{1,s} = strcat(tmp,',',...
        regexprep(cellstr(num2str((1:dyn)')),' ',''));
    % number to string -> string to cell -> remove spaces -> concatenate
    
    sess_row = new_sess(s);
    
    % Check if there are conditions
    %------------------------------
    if eeg
        tmp = strcmp(xlstxt{sess_row+3,1},'Type');
    else
        tmp = strcmp(xlstxt{sess_row+1,1},'Description');
    end
    
    % 1) There are conditions for this session
    %======================================================================
    if tmp
        
        % Counter for rows in the session
        %--------------------------------
        if eeg
            row0_txt = 5; % number of rows with the condition name after the identification of a new fMRI session
            row0_num = 3; % number of rows with the condition onset after the identification of a new fMRI session
            col_cond = 2; % column with the condition name
        else
            row0_txt = 2;
            row0_num = 0;
            col_cond = 1;
        end
        
        c       = 0; % counter for conditions
        row_cnt = 0;
        prev    = ''; % initialization for the previous condition
        c_name  = xlstxt{sess_row+row0_txt,col_cond}; % condition name
        
        while ~isempty(c_name) % loop for conditions
            onset_tmp = xlsnum(sess_row+row0_num+row_cnt,1) + delay;
            dur_tmp   = xlsnum(sess_row+row0_num+row_cnt,2);
            %if onset_tmp+dur_tmp>=-2*TR && onset_tmp<=(dyn-3)*TR
            % Because of the derivatives. The checking for this is now
            % incorporated in the contrast function.
            if onset_tmp+dur_tmp>=0 && onset_tmp<=(dyn-1)*TR
                % to guarantee that the onsets are inside the scan time
                
                % 1.1) Same condition as the previous
                %----------------------------------------------------------
                if isequal(c_name,prev)
                    cond{1,s}(c).onset    = [cond{1,s}(c).onset    ; onset_tmp];
                    cond{1,s}(c).duration = [cond{1,s}(c).duration ; dur_tmp];
                    
                % 1.2) Another condition
                %----------------------------------------------------------
                else
                    c = c + 1;
                    cond{1,s}(c).name     = c_name;    % name
                    cond{1,s}(c).onset    = onset_tmp; % onset
                    cond{1,s}(c).duration = dur_tmp;   % duration
                    cond{1,s}(c).tmod     = 0;         % time modulation
                    cond{1,s}(c).pmod     = struct('name',{},'param',{},'poly',{}); % parametric modulations
                    prev = cond{1,s}(c).name;
                    
                end
                
            end
            
            row_cnt = row_cnt + 1; % go to the next row
            if sess_row+row0_txt+row_cnt<=size(xlstxt,1)
                c_name = xlstxt{sess_row+row0_txt+row_cnt,col_cond};
            else % reached the end of "xlstxt"
                break
            end
        end
        
        % No useful condition for this session
        %-------------------------------------
        if c==0
            if use_rp % still have realignment parameters => keep session
                cond{1,s} = struct('name',{},...
                    'onset',{},...
                    'duration',{},...
                    'tmod',{},...
                    'pmod',{});
            else % discard session
                new_sess(s) = [];
                sess        = sess - 1;
                continue
            end
        end
        
        
    % 2) No conditions
    %======================================================================
    else
        
        % 2.1) Using realignment parameters
        %----------------------------------
        if use_rp
            cond{1,s} = struct('name',{},...
                'onset',{},...
                'duration',{},...
                'tmod',{},...
                'pmod',{});
        
        % 2.2) Not using realignment parameters
        %--------------------------------------
        else % discard session
            new_sess(s) = [];
            sess        = sess - 1;
            continue
        end
        
    end
    
    if use_rp
        rp{1,s} = fullfile(fpath,['rp_' xlstxt{new_sess(s),2}(1:end-3) 'txt']);
    end
    
    s = s + 1;
    
end

% Save movie
%==========================================================================
else
    
    % Initialize
    %-----------
    s    = 1;
    sess = 1;
    tmp  = fullfile(fpath,[pref_func xlstxt{new_sess(moviesess),2}]);
    
    if exist(tmp,'file')~=2
        error('File %s was not found',tmp)
    end
    
    
    % Number of scans
    %----------------
    dyn = spm_read_hdr(tmp);
    try
        dyn = dyn.dime.dim(5);
    catch ME
        error(ME.identifier,ME.message)
    end
    
    
    % Functional scans
    %-----------------
    func_scan{1,s} = strcat(tmp,',',...
        regexprep(cellstr(num2str((1:dyn)')),' ',''));
    % number to string -> string to cell -> remove spaces -> concatenate
    
    
    % Adjust the onset (delay)
    %-------------------------
    if delay<0
        movieeventdur = movieeventdur + delay;
        onset_tmp     = 0;
    else
        onset_tmp     = 0 + delay;
    end
    
    dur_tmp = movieeventdur;
    if onset_tmp+dur_tmp>=0 && onset_tmp<=(dyn-1)*TR
        % to guarantee that the onsets are inside the scan time
        
        % Get selected dynamic scans
        %---------------------------
        if isfinite(movie_ti) && isfinite(movie_tf) % both finite
            tmp = floor((onset_tmp-movie_ti)/TR+1) : ceil((onset_tmp+dur_tmp+movie_tf)/TR+1);
        elseif isfinite(movie_ti) % only "movie_ti" is finite
            tmp = floor((onset_tmp-movie_ti)/TR+1) : dyn;
        elseif isfinite(movie_tf) % only "movie_tf" is finite
            tmp = 1 : ceil((onset_tmp+dur_tmp+movie_tf)/TR+1);
        else % both infinite or NaN
            tmp = 1:dyn;
        end
        tmp(tmp<=0) = []; tmp(tmp>dyn) = [];
        
        % Conditions
        %-----------
        if ~isempty(tmp)
            
            % Adjust the onset
            %-----------------
            if onset_tmp>movie_ti
            % the same equation works for the other case because tmp(1)=1
                onset_tmp = onset_tmp - (tmp(1) - 1)*TR;
            end
            
            c = 1;
            cond{1,s}(c).name     = 'Event';   % name
            cond{1,s}(c).onset    = onset_tmp; % onset
            cond{1,s}(c).duration = dur_tmp;   % duration
            cond{1,s}(c).tmod     = 0;         % time modulation
            cond{1,s}(c).pmod     = struct('name',{},'param',{},'poly',{}); % parametric modulations
        else
            cond{1,s} = struct('name',{},...
                'onset',{},...
                'duration',{},...
                'tmod',{},...
                'pmod',{});
            tmp = 1:dyn; % select all the images (will not make a difference)
        end
        
    else
        cond{1,s} = struct('name',{},...
            'onset',{},...
            'duration',{},...
            'tmod',{},...
            'pmod',{});
        tmp = 1:dyn; % select all the images (will not make a difference)
    end
    
    % Select images
    %--------------
    selected_scans = tmp;
    func_scan{1,s} = func_scan{1,s}(selected_scans);
    
    % Realignment parameters
    %-----------------------
    if use_rp
        rp{1,s} = load(fullfile(fpath,['rp_' xlstxt{new_sess(moviesess),2}(1:end-3) 'txt']));
    end
    
end


% Create folder
%--------------------------------------------------------------------------
delay_folder = fullfile(outdir, sprintf('%s-SPM_delay_%.1f', worksheet, delay));
create_folder(delay_folder)


% Model specification
%==========================================================================

% Load defaults
%--------------------------------------------------------------------------
matlabbatch{1}.spm.stats.fmri_spec.dir            = {delay_folder};
matlabbatch{1}.spm.stats.fmri_spec.timing.units   = 'secs';     % units for design
matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = 0;          % interscan interval (TR) [s]
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;         % microtime resolution
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;          % microtime onset
matlabbatch{1}.spm.stats.fmri_spec.fact           = struct('name',{}, 'levels',{}); % factorial design
matlabbatch{1}.spm.stats.fmri_spec.volt           = 1;          % model interactions (Volterra)
matlabbatch{1}.spm.stats.fmri_spec.global         = 'None';     % global normalisation
matlabbatch{1}.spm.stats.fmri_spec.mask           = {''};       % explicit mask
matlabbatch{1}.spm.stats.fmri_spec.cvi            = 'AR(1)';    % serial correlations

% Customized
%--------------------------------------------------------------------------
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
if strcmpi(bf,'can_hrf')
    if bf_order==1
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % model derivatives (canonical HRF)
    elseif bf_order==2
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
    elseif bf_order>2
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
    end
elseif any(strcmpi(bf,{'fourier','fourier_han','gamma','fir'}))
    matlabbatch{1}.spm.stats.fmri_spec.bases.(lower(bf)) = struct(...
        'length',bf_length,...
        'order',bf_order);
else
    error('Invalid basis function')
end


for s=1:sess
    
    % Load defaults
    %----------------------------------------------------------------------
    matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans     = '<UNDEFINED>';              % scans
    matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond      = struct('name',{}, 'onset',{}, 'duration',{}, 'tmod',{}, 'pmod',{});
    matlabbatch{1}.spm.stats.fmri_spec.sess(s).multi     = {''};                       % multiple conditions
    matlabbatch{1}.spm.stats.fmri_spec.sess(s).regress   = struct('name',{},'val',{}); % regressors
    matlabbatch{1}.spm.stats.fmri_spec.sess(s).multi_reg = {''};                       % multiple regressors
    matlabbatch{1}.spm.stats.fmri_spec.sess(s).hpf       = 128;                        % high-pass filter [s]
    
    % Customized
    %----------------------------------------------------------------------
    matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans = func_scan{1,s};
    matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond  = cond{1,s};
    if use_rp
        if ~moviesave
            matlabbatch{1}.spm.stats.fmri_spec.sess(s).multi_reg = rp(1,s);
        else
            tmp = ['st';'nd';'rd';'th';'th';'th'];
            for mm=1:6
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(mm).name = sprintf('%d%s movement parameter',mm,tmp(mm,:));
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(mm).val = rp{1}(selected_scans,mm);
            end
        end
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(s).hpf = hpf;
    % Common approach specially if slice timing correction is applied
    % (doesn't seem to have much effect):
%     if isequal(s,1)
%         % Number of slices for the first scan:
%         nslices = getfield(getfield(spm_read_hdr(func_scan{1,s}{1}),'dime'),'dim',{4});
%         matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = nslices;
%         matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(nslices/2);
%     end

end


% Estimate GLM parameters
%==========================================================================
% Not using dependency:
%matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(delay_folder,'SPM.mat')};
% Using dependency:
matlabbatch{2}.spm.stats.fmri_est.spmmat(1)                      = cfg_dep;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname                = 'Select SPM.mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name  = 'filter';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name  = 'strtype';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname                = 'fMRI model specification: SPM.mat File';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch         = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output           = substruct('.','spmmat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical               = 1;


% Save and run
%==========================================================================
new_file = fullfile(delay_folder,[worksheet '-design_matrix.mat']);
save(new_file,'matlabbatch')
fprintf('The file %s was created\n',new_file)


no_rp = 0;
if ~all(cellfun(@isempty,cond)) % there is at least one condition
    if run
        try
            %spm_jobman('initcfg');
            spm_jobman('run',matlabbatch)
        catch
            warning('Model estimation was not performed. %s.m was skipped',mfilename)
            no_rp = 1;
        end
    end
else
    warning('Model design and estimation were not performed. %s.m was skipped',mfilename)
    no_rp = 1;
end

if ~no_rp
    
    % Conditions on the movement parameters
    %======================================================================
    rp_fmri_files = cell(size(func_scan,2),1);
    fpath         = fileparts(func_scan{1,1}{1});
    if size(func_scan,2)>1
        tmp = 2; % there will be a session separator
    else
        tmp = 1; % no session separator
    end
    rp_cond_time  = cell(tmp,1);
    last_sess     = 0;
    for ff=1:size(func_scan,2)
        
        % File name
        %----------
        tmp = strfind(func_scan{1,ff}{1},filesep);
        rp_fmri_files{ff} = func_scan{1,ff}{1}(tmp(end)+1+length(pref_func):end); % remove path and prefix
        
        % Remove the comma and the number that comes after
        tmp = strfind(rp_fmri_files{ff},',') - 1;
        if ~isempty(tmp) % found ","
            tmp = tmp(end); % get the last ","
            if tmp+2<=size(rp_fmri_files{ff},2) && ...
                    ( ~isreal(str2double(rp_fmri_files{ff}(tmp+2:end))) || ...
                    isnan(str2double(rp_fmri_files{ff}(tmp+2:end))) )
                % "," must be followed by a real number. If this doesn't happen
                % "," is in a weird place in the file name
                tmp = size(rp_fmri_files{ff},2);
            end
        else % couldn't find ","
            tmp = size(rp_fmri_files{ff},2);
        end
        rp_fmri_files{ff} = fullfile(fpath,rp_fmri_files{ff}(1:tmp));
        
        % Adjust beginning of the sessions
        %---------------------------------
        if ff>1
            last_sess = last_sess + size(func_scan{1,ff-1},1);
            rp_cond_time{2,1} = [rp_cond_time{2,1} ; last_sess last_sess+1];
        end
        
        % Get all the conditions
        %-----------------------
        if ~moviesave
            for cc=1:size(cond{1,ff},2)
                rp_cond_time{1,1} = [rp_cond_time{1,1};
                    [cond{1,ff}(1,cc).onset cond{1,ff}(1,cc).onset+cond{1,ff}(1,cc).duration]/TR + 1 + last_sess];
            end
        else
            rp_cond_time{1,1} = [delay delay+movieeventdur]/TR + 1;
        end

    end
    
    % Make event appear
    %------------------
    % If the event duration is last than 1 scan, it won't appear
    tmp = rp_cond_time{1,1}(:,2)-rp_cond_time{1,1}(:,1)<1;
    rp_cond_time{1,1}(tmp,2) = rp_cond_time{1,1}(tmp,1) + ones(sum(tmp,1),1);
    
    
    % Plot figure
    %------------
    [fighandle mov_pars] = fmri_plot_rp(rp_fmri_files,rp_cond_time,[[1 1 1]/3;[0 0 0]]);
    
    % Make black bars dividing sessions
    %----------------------------------
    tmp = findobj('Type','patch','FaceColor',[0 0 0]);
    for ff=1:length(tmp)
        set(tmp(ff),'FaceAlpha',1)
    end
    
    % Save figure
    %------------
    drawnow; pause(.1)
    saveas(fighandle,fullfile(delay_folder,[worksheet '-rp-conditions.png']))
    %fig_print(fighandle,...
    %    fullfile(delay_folder,[worksheet '-rp-conditions.png']),...
    %    [1200 900])
    % Same output, except for the grey background. White is better for
    % printing
    close(fighandle)
    
end


% Finish
%==========================================================================

% Close SPM windows
close(findobj('Tag','Graphics','Type','figure'))
close(findobj('Tag','Interactive','Type','figure'))

tmp = clock;
fprintf('(%.2d:%.2d:%02.0f) %s.m done!\n',...
    tmp(4),tmp(5),tmp(6),mfilename)
disp('-------------------------------------------------------------------')