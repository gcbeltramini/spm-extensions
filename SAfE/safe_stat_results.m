function safe_res = safe_stat_results(SPM, opt_res)
% 
% SAFE_STAT_RESULTS saves the statistical results.
% 
% 
%USAGE
%-----
% safe_stat_results(spm, opt_res)
% 
% 
%INPUT
%-----
% - SPM    : SPM.mat file (after SPM model estimation and contrasts creation)
% - OPT_RES: structure with the following fields:
%   - P       : P-value (any non-negative number smaller than 1)
%     (default: 0.001 [uncorrected] or 0.05 [FWE])
%   - ET      : extent threshold in voxels (ET >= 0)
%     (default: 0 [uncorrected or FWE])
%   - ADJ     : 'FWE' or 'none' (P-value adjustment) (default: 'none')
%   - UNDERLAY: underlay image (full path)
%   - MASK    : 1 (apply masking) or 0 (don't apply)
%     - If MASK=1:
%       - MASK_IMGS: Nx1 cell array containing the mask images (full path)
%       - MASK_INCL_EXCL: nature of mask: 'inclusive' or 'exclusive'
%   - SAVE_SLOVER: 1 (create slover object) or 0 (don't create)
%   - SAVE_PS    : 1 (create PS file) or 0 (don't create)
%   - SAVE_THR   : 1 (create thresholded image) or 0 (don't create)
%   - THR_PREF   : prefix for the thresholded images
%   - SAVE_SLC   : 1 (create image with slices) or 0 (don't create)
%     - If SAVE_SLC=1:
%       IMG_EXT: image extension (default: 'png')
%   If SAVE_SLC=1:
%     - SLC_STEP: spacing between slices (mm) (any positive number)
%     (default: 3)
%     - SLC_NMBR: number of slices (to use this option, SLC_STEP must be 0)
%     (default: 46)
%     - SLC_SLICES: coordinates of the slices (to use this option, SLC_STEP
%     and SLC_NMBR must be 0)
%     - ORIENT  : cell array of string, containig at least one of the
%     following options: 'axial', 'coronal', 'sagittal'
% 
% When SAVE_PS=1 or SAVE_SLC=1 or SAVE_THR=1, this function
% assumes that SPM.mat is inside a folder "...\..._delay_DELAY". If it is
% not, please rename it so that it contains the characters "_delay_"
% followed by the value of the delay in the end of the folder name.
% 
% 
%OUTPUT
%------
% - If SAVE_SLOVER=1 or SAVE_PS=1 or SAVE_THR=1 or SAVE_SLC=1:
%   - SAFE_RES: structure containing the following fields:
%     If SAVE_SLOVER=1 or SAVE_SLC=1:
%     - MN_MX: 2Cx4 matrix containing the minimum (columns 1 and 3) and
%     maximum (columns 2 and 4) values of the statistics. C is the number
%     of conditions. Odd rows: F test (columns 1 and 2 only); even rows: T
%     test (positive BOLD in columns 1 and 2; negative BOLD in columns 3
%     and 4). For every two rows, the conditions are changed
%     - OBJ: 2Cx1 cell array containing the slover object for every
%     condition. Odd rows: F test; even rows: T test
%     If SAVE_SLOVER=1 or SAVE_THR=1 or SAVE_SLC=1:
%     - THR: structure containing the following fields:
%       - FNAME   : 3Cx1 cell array, containing the full path of the
%       thresholded images. Rows 1,4,7,...: F test; rows 2,5,8,...:
%       positive T test; rows 3,6,9,...: negative T test
%       - ADD_STAT: 3Cx1 matrix, with the summation of the statistics value
%       for each image in THR.FNAME
%       - VOXELS  : 3Cx1 matrix, with the number of "active" voxels for
%       each image in THR.FNAME
% 
% - If SAVE_PS=1:
%   - SID_DELAYs_P0.XXXX_Vvox_CORR.ps in the same folder as the SPM.mat
%   file, where:
%     - SID   : subject ID (e.g., 'XXPy')
%     - DELAY : time in seconds from the stimulus onset to the beginning
%     of the hemodynamic response function
%     - 0.XXXX: P value
%     - V     : extent threshold
%     - CORR  : 'unc' or 'FWE'
%   - Two pages for each contrast: the first one shows the glass brain view
%   and the p-values for the whole brain; the second one shows the
%   statistical parametric map over the UNDERLAY image
% 
% - If SAVE_THR=1 or SAVE_SLOVER=1 or SAVE_SLC=1:
%   - In the same folder as SPM.mat:
%     files "THR_PREF"_SID_CONDNAME_TEST_P0.XXXX_Vvox_CORR_DELAY.hdr/img,
%     where:
%     - SID     : the subject ID (e.g., 'XXPy')
%     - CONDNAME: condition name
%     - TEST    : 'posBOLD', 'negBOLD' or 'Ftest'
%     - 0.XXXX  : P value
%     - V       : extent threshold
%     - CORR    : 'unc' or 'FWE'
%     - DELAY   : time in seconds from the stimulus onset to the beginning
%     of the hemodynamic response function
% 
% - If SAVE_SLC=1:
%   - SAVE_THR is an intermediate step.
%   - SID_CONDNAME_TEST_P0.XXXX_Vvox_CORR_DELAYs_ORI.IMG_EXT files in the
%   same folder as the SPM.mat file, where:
%     - SID     : the subject ID (e.g., 'XXPy')
%     - CONDNAME: condition name
%     - TEST    : 'Ftest' or 'pos-negBOLD'
%     - 0.XXXX  : P value
%     - V       : extent threshold
%     - CORR    : 'unc' or 'FWE'
%     - DELAY   : time in seconds from the stimulus onset to the beginning
%     of the hemodynamic response function
%     - ORI     : 'axi', 'cor' and/or 'sag', depending on what was chosen
%     in ORIENTATION
% 
% 
%EXAMPLES
%--------
% opt_res = struct(...
%     'P'          ,0.05,...
%     'ET'         ,0,...
%     'adj'        ,'FWE',...
%     'underlay'   ,'D:\EEG-fMRI\01P\SPM8\wm_01P_T13D.nii',...
%     'save_slover',1,...
%     'save_ps'    ,1,...
%     'save_thr'   ,1,...
%     'thr_pref'   ,'Thresh_',...
%     'save_slc'   ,0);
% safe_stat_results('D:\EEG-fMRI\01P\SPM8\SPM.mat',opt_res)
% 
% opt_res = struct(...
%     'P'          ,0.001,...
%     'ET'         ,20,...
%     'adj'        ,'none',...
%     'underlay'   ,'D:\fMRI\05P\SPM8\wm_05P_T13D.nii',...
%     'save_slover',1,...
%     'save_ps'    ,0,...
%     'save_thr'   ,0,...
%     'thr_pref'   ,'Thresh_',...
%     'save_slc'   ,1,...
%     'img_ext'    ,'png',...
%     'slc_step'   ,0,...
%     'slc_nmbr'   ,30,...
%     'slc_slices' ,0,...
%     'orient'     ,{'axial','sagittal'});
% safe_stat_results('D:\EEG-fMRI\05P\SPM8\SPM.mat',opt_res)
% 
% See also SAFE_PREPROC, SAFE_STAT_DESIGN, SAFE_STAT_CONTRASTS
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2014-Jul-14


max_time = 20; % maximum tolerance (in seconds) to calculate the height threshold


% Read options
%--------------------------------------------------------------------------
P        = opt_res.P;
ET       = opt_res.ET;
adj      = opt_res.adj;
underlay = opt_res.underlay;
if exist(underlay,'file')~=2
    error('Underlay image "%s" not found',underlay)
end
mask     = opt_res.mask;
if mask
    mask_imgs      = opt_res.mask_imgs(:);
    mask_incl_excl = opt_res.mask_incl_excl;
    switch mask_incl_excl
        case 'inclusive'
            mask_incl_excl = 0;
        case 'exclusive'
            mask_incl_excl = 1;
    end
else
    mask_imgs = [];
end
save_slover = opt_res.save_slover;
save_ps     = opt_res.save_ps;
save_slc    = opt_res.save_slc;
if save_slc
    img_ext = opt_res.img_ext;
end
save_thr   = opt_res.save_thr;
thr_pref   = opt_res.thr_pref;
save_movie = 0; %opt_res.save_movie;
if save_movie
    moviedir = opt_res.moviedir;
end
if save_slc || save_movie || save_slover
    slc_step   = opt_res.slc_step;
    slc_nmbr   = opt_res.slc_nmbr;
    slc_slices = opt_res.slc_slices;
    orient     = opt_res.orient;
    num_orient = length(orient);
end
visible = 'off'; % figure visibility: 'on' or 'off'


% Output
%--------------------------------------------------------------------------
safe_res = struct('thr',[], 'mn_mx',[], 'obj',[], 'movie_file',[]);


% P value and number of voxels
%--------------------------------------------------------------------------
if P<0 || P>1, error('P must be greater than 0 and less than 1'), end
if ET<0, error('Extent threshold (voxels) must be non-negative'), end


tmp = clock;
disp('-------------------------------------------------------------------')
fprintf('(%.2d:%.2d:%02.0f) Running %s.m...\n',...
    tmp(4),tmp(5),tmp(6),mfilename)


% To avoid warning with command: "spm_print(ps_file)"
%spm_jobman('initcfg'); % this is now run in safe_preproc_stat


% Load SPM.mat
%--------------------------------------------------------------------------
fpath = fileparts(SPM);
if isempty(fpath), fpath = pwd; end
SPM = fullfile(fpath, 'SPM.mat');
SPMmat = SPM;
try
    load(SPM)
catch
    warning('%s not found. The path is wrong, or you must specify and estimate the model first.',SPMmat)
    safe_stat_results_finish
    return
end


% Test for number of contrasts when creating a movie
%--------------------------------------------------------------------------
if save_movie && isfield(SPM,'xCon') && length(SPM.xCon)>4
    error('Too many contrasts to make movie. There should be one condition at most with F test, and positive and negative BOLD; and a movement parameters F test.')
end


% Delay between task onset and beginning of HRF
%--------------------------------------------------------------------------
tmp   = strfind(fpath,'_delay_');
tmp   = tmp(end);
delay = str2double(fpath(tmp+7:end));


% Get subject identification
%--------------------------------------------------------------------------
if ~isfield(SPM,'xCon') || isempty(SPM.xCon)
    warning('You must create the constrasts first. %s.m was skipped',mfilename)
    safe_stat_results_finish
    return
end
%sID = sscanf(SPM.xCon(1).name,'%s'); % get first contrast name and remove spaces
%sID = sID(1:regexp(sID,'-','once')-1);
tmp = strfind(SPM.xCon(1).name,' - ');
tmp = tmp(1) - 1;
sID = SPM.xCon(1).name(1:tmp);


% P-value adjustment
%--------------------------------------------------------------------------
switch lower(adj)
    case 'fwe'
        adj      = 'FWE';
        adj_name = 'FWE';
    case 'none'
        adj      = 'none';
        adj_name = 'unc';
    otherwise
        error('Invalid option for P value adjustment')
end


% Mask suffix
%--------------------------------------------------------------------------
% "_mask" must be in the end (other functions expect that)
if ~mask
    mask_suffix = [];
else
    mask_suffix = '_mask';
end


% PS file name
%--------------------------------------------------------------------------
if save_ps
    
    ps_file = sprintf('%s_%.1fs_P%g_%dvox_%s%s.ps',sID,delay,P,ET,adj_name,mask_suffix);
    
    % Preserve old files
    %-------------------
    if exist(fullfile(fileparts(SPMmat),ps_file),'file')==2
        ps_file = sprintf('%s_1.ps',ps_file(1:end-3));
    end
    tmp = 1;
    while exist(fullfile(fileparts(SPMmat),ps_file),'file')==2
        % PS file already exists
        tmp = tmp + 1;
        ps_file = sprintf('%s%d.ps',ps_file(1:end-floor(log10(tmp-1))-4),tmp);
    end

    % Rename old files  ---  Now the old files are preserved
    %-----------------
%     if exist(ps_file_full,'file')
%         if exist([ps_file_full(1:end-3) '-backup.ps'],'file')
%             tmp = 2;
%             while exist([ps_file_full(1:end-3) sprintf('-backup%.2d.ps',tmp)],'file')
%                 tmp = tmp + 1;
%             end
%             s = movefile(ps_file_full,[ps_file_full(1:end-3) sprintf('-backup%.2d.ps',tmp)]);
%         else
%             s = movefile(ps_file_full,[ps_file_full(1:end-3) '-backup.ps']);
%         end
%         if ~s
%             warning('%s already exists and could not be renamed',ps_file)
%             ps_file = [ps_file(1:end-3) '_new.ps'];
%         end
%     end

end


spm('defaults', 'FMRI') % because of the function spm_results_ui.m
if save_slc || save_movie || save_slover
    pos_neg   = 0;
    count_obj = 0;
    safe_res.mn_mx = zeros(length(SPM.xCon),4);
    safe_res.obj   = cell(num_orient*length(SPM.xCon),1);
end


% Options to save figure
%--------------------------------------------------------------------------
if save_slc
    opt_slvr_save.delay    = delay;
    opt_slvr_save.save_slc = save_slc;
end


if save_ps
    % turn off warning when creating a new spreadsheet
    warning off MATLAB:xlswrite:AddSheet
end

sheet_name = cell(length(SPM.xCon), 1);
for con = 1:length(SPM.xCon)
    
    % Identify the contrast
    %----------------------------------------------------------------------
    contr_name = SPM.xCon(con).name;
    contr_name = contr_name((length(sID)+4):end); % remove subject ID
    sgn        = strfind(contr_name,'baseline') - 2;
    
    % Type of contrast
    %----------------------------------------------------------------------
    if strcmp(contr_name(sgn-1:sgn), '<>')
        
        % F contrast
        %------------------------------------------------------------------
        contr_type = 'Ftest';
        name_end   = sgn - 3;
        % Colormap
        %spm_figure('ColorMap', 'gray-hot')
        tmp = hot(64 + 16);  tmp = tmp((1:64) + 16,:);
        colormap([gray(64); tmp]);
        
    elseif strcmp(contr_name(sgn-1:sgn), ' >')
        
        % Positive BOLD
        %------------------------------------------------------------------
        contr_type = 'posBOLD';
        name_end   = sgn - 2;
        tmp = hot(64 + 16);  tmp = tmp((1:64) + 16,:);
        colormap([gray(64); tmp]);
        
    elseif strcmp(contr_name(sgn-1:sgn), ' <')
        
        % Negative BOLD
        %------------------------------------------------------------------
        contr_type = 'negBOLD';
        name_end   = sgn - 2;
        %spm_figure('ColorMap', 'gray-cool')
        cool = [zeros(10,1) zeros(10,1) linspace(0.5,1,10)';
            zeros(31,1) linspace(0,1,31)' ones(31,1);
            linspace(0,1,23)' ones(23,1) ones(23,1) ];
        colormap([gray(64); cool]);
        
    else
        
        % Movement parameters
        %------------------------------------------------------------------
        contr_type = '';
        name_end   = 8; % 'Movement' has 8 characters
        tmp = hot(64 + 16);  tmp = tmp((1:64) + 16,:);
        colormap([gray(64); tmp]);
        
    end
    
    
    % Contrast name - remove invalid characters in a file name
    %----------------------------------------------------------------------
    contr_name = regexprep(contr_name(1:name_end),{'"','*',':','<','>','?','\'},'_');
    contr_name(strfind(contr_name,'|')) = '_';
    % apparently there is a bug with regexp, because it cannot find '|'
    
    
    % Print PS file
    %======================================================================
    if save_ps
        
        % Create variable xSPM
        %------------------------------------------------------------------
        xSPM.swd       = fpath;              % SPM.mat folder
        xSPM.Ic        = con;                % Contrast(s)
        xSPM.u         = P;                  % Threshold
        xSPM.Im        = mask_imgs;          % Masking
        if mask
            xSPM.Ex = mask_incl_excl; % exclusive or inclusive
        end
        xSPM.thresDesc = adj;                % Threshold type: 'none', 'FWE'
        xSPM.title     = SPM.xCon(con).name; % Results Title
        xSPM.k         = ET;                 % Extent (voxels)
        xSPM.units     = {'mm','mm','mm'};   % Data type: Volumetric (2D/3D)
        
        % Test time to calculate "height threshold" in spm_results_ui (from spm_uc_RF.m)
        %------------------------------------------------------------------
        stop = test_time(SPM, xSPM, max_time);
        if stop
            fprintf('%.1f s have elapsed. Aborting...\n', max_time)
            safe_stat_results_finish
            safe_res = struct('thr',[], 'mn_mx',[], 'obj',[], 'movie_file',[]);
            return
        end
        
        [hReg,xSPM,SPM] = spm_results_ui('Setup', xSPM); % show results
        %set(spm_figure('FindWin','Graphics'),'Visible',visible) % big window
        % it will show up again below
        set(spm_figure('FindWin','Interactive'), 'Visible', visible) % small window
        
        
        % Clusters for the whole brain
        %------------------------------------------------------------------
        TabDat = spm_list('List', xSPM, hReg); % show p-values for the whole brain
        %TabDat = get(findobj('Tag','TabDat'),'UserData'); % extract table data structure
        fig_graphic = spm_figure('FindWin','Graphics');
        set(fig_graphic, 'Visible', visible) % big window
        drawnow; pause(.1)
        spm_print(ps_file)                   % print PS file
        
        
        % Save statistics table
        %------------------------------------------------------------------
        if ~isempty(contr_type) % ignore movement contrast
            % 1) Print text table in Command Window:
            %    spm_list('TxtList',TabDat,3)
            % 2) Export table to Excel (depends on the version):
            %    export2excel(obj,evd,h) in spm_list
            %    spm_list('XLSList',TabDat)
            d          = [TabDat.hdr(1:2,:); TabDat.dat];
            xyz        = d(3:end,end);
            xyz        = num2cell([xyz{:}]');
            d(:,end+1) = d(:,end);
            d(:,end+1) = d(:,end);
            if size(TabDat.dat,1)>0
                d(3:end,end-2:end) = xyz;
            else  % no suprathreshold voxels
                d{3,1} = 'No suprathreshold clusters';
            end
            
            % Export to Excel
            %----------------
            % Remove invalid characters
            sheet_name{con} = regexprep([contr_name '-' contr_type],{'/','\','?','*','[',']',':'},'');
            if strcmp(sheet_name{con}(1),'''') % leading single quote is not allowed
                sheet_name{con}(1) = [];
            end
            if strcmp(sheet_name{con}(end),'''') % trailing single quote is not allowed
                sheet_name{con}(end) = [];
            end
            tmp = length(sheet_name{con});
            if tmp>31 % maximum allowed name length is 31 characters
                sheet_name{con} = sheet_name{con}(tmp-30:tmp);
            elseif tmp==0 % empty
                sheet_name{con} = ' ';
            end
            if any(strcmp(sheet_name{con},sheet_name(1:con-1))) % spreadsheet already exists
                tmp_name = [sheet_name{con} '1'];
                tmp = length(tmp_name);
                if tmp>31
                    tmp_name = tmp_name(tmp-30:tmp);
                end
                tmp = 2;
                while any(strcmp(tmp_name,sheet_name(1:con-1)))
                    tmp_name = sprintf('%s%d',sheet_name{con},tmp);
                    tmp = tmp + 1;
                end
                sheet_name{con} = tmp_name;
            end
            
            try
                xlswrite([ps_file(1:end-2) 'xls'],d,sheet_name{con});
            catch
                
                % Convert empty cells to NaN
                tmp = d(3:end,:);
                tmp(cellfun(@isempty,tmp)) = {NaN};
                d(3:end,:) = tmp;
                
                % Remove spaces from the header
                tmp = find(~cellfun(@ischar,d(:,1)),1,'first') - 1; % no. of lines in the header
                d(1:tmp,:) = strrep(d(1:tmp,:),' ','');
                
                fid = fopen([ps_file(1:end-3) '-' sheet_name{con} '.txt'],'w');
                
                % Write text
                tmp = repmat('%s\t',1,size(d,2));
                tmp = [tmp(1:end-1) 'r\n'];
                for rr=1:2
                    fprintf(fid,tmp,d{rr,:});
                end
                
                % Write numbers
                tmp = repmat('%g\t',1,size(d,2));
                tmp = [tmp(1:end-1) 'r\n'];
                for rr=3:size(d,1)
                    fprintf(fid,tmp,d{rr,:});
                end
                
                fclose(fid);
                
            end
        end
        
        
        % Go to global maximum
        %------------------------------------------------------------------
        if ~isempty(xSPM.Z) % there are suprathreshold voxels
            % otherwise:
            % spm('alert!','No suprathreshold voxels to jump to!',mfilename,0);
            spm_mip_ui('Jump', spm_mip_ui('FIndMIPax', fig_graphic), 'glmax'); % go to global maxima
        end
        
        
        % Show functional map overlaid on the chosen image
        %------------------------------------------------------------------
        spm_sections(xSPM, hReg, underlay)
        set(fig_graphic, 'Visible', visible) % big window
        drawnow; pause(.1)
        spm_print(ps_file) % print PS file
        
    end
    
    
    % Save thresholded SPM as image
    %======================================================================
    if (save_thr || save_slc || save_movie || save_slover) && ~isempty(contr_type)
        % Ignore movement contrast
        % Alternative: isempty(strfind(SPM.xCon(i).name,' - Movement'))
        
        
        thr_name = sprintf('%s%s_%s_%s_P%g_%dvox_%s_%.1fs%s', ...
            thr_pref, sID, contr_name, contr_type, ...
            P, ET, adj_name, ...
            delay, mask_suffix);
        safe_res.thr.fname{con,1} = fullfile(fpath, thr_name);
        
        
        % Write thresholded image
        %------------------------------------------------------------------
        %if ~exist([safe_res.thr.fname{con,1} '.img'],'file')
        % thresholded image does not exist
        % It is better to overwrite previous analyses
            
            if ~save_ps
            % variable "xSPM" does not exist
                xSPM.swd       = fileparts(SPMmat);  % SPM.mat folder
                xSPM.Ic        = con;                % Contrast(s)
                xSPM.u         = P;                  % Threshold
                xSPM.Im        = mask_imgs;          % Masking
                if mask
                    xSPM.Ex = mask_incl_excl; % exclusive or inclusive
                end
                xSPM.thresDesc = adj;                % Threshold type: 'none', 'FWE'
                xSPM.title     = SPM.xCon(con).name; % Results Title
                xSPM.k         = ET;                 % Extent (voxels)
                xSPM.units     = {'mm','mm','mm'};   % Data type: Volumetric (2D/3D)
                
                % Test time to calculate "height threshold" in spm_results_ui (from spm_uc_RF.m)
                %----------------------------------------------------------
                stop = test_time(SPM,xSPM,max_time);
                if stop
                    fprintf('%.1f s have elapsed. Aborting...\n',max_time)
                    safe_stat_results_finish
                    safe_res = struct('thr',[],'mn_mx',[],'obj',[],'movie_file',[]);
                    return
                end
                
                [hReg, xSPM, SPM] = spm_results_ui('Setup', xSPM); % show results
                set(spm_figure('FindWin', 'Graphics'), 'Visible', visible) % big window
                set(spm_figure('FindWin', 'Interactive'), 'Visible', visible) % small window
            end
            
            if (save_thr || save_slc || save_movie || save_slover)
                spm_write_filtered_mod(xSPM.Z, xSPM.XYZ, xSPM.DIM, xSPM.M,...
                    sprintf('SPM{%c}-filtered: u = %5.3f, k = %d', xSPM.STAT, xSPM.u, xSPM.k),...
                    thr_name);
            end
            
        %end
        
        
        % Add statistics and count the number of voxels above the threshold
        %------------------------------------------------------------------
        %safe_res.thr.contr_type{con,1} = contr_type;
        
        tmp = nifti([safe_res.thr.fname{con,1} '.img']);
        tmp = tmp.dat(:,:,:);
        tmp(isnan(tmp)) = 0; % NaN -> 0
        safe_res.thr.add_stat(con,1) = sum(sum(sum(tmp,1),2),3);
        
        tmp = tmp>0;
        safe_res.thr.voxels(con,1) = sum(sum(sum(tmp,1),2),3);
        
    end
    
    
    % Save: Slices  and/or  Movies
    %======================================================================
    if (save_slc || save_movie || save_slover) && ~isempty(contr_type)
        
        % For all conditions, the order is: F test (if available), positive
        % T-test and negative T-test.
        
        % F test
        %------------------------------------------------------------------
        if strcmp(contr_type, 'Ftest')
            
            
            % Create slover object
            %--------------------------------------------------------------
            opt_slvr_create.imgs       = {underlay, fullfile(fpath, [thr_name '.img'])};
            opt_slvr_create.type       = {'Structural', 'Blobs'};
            opt_slvr_create.cmap       = {'', 'hot'};
            opt_slvr_create.slc_step   = slc_step;
            opt_slvr_create.slc_nmbr   = slc_nmbr;
            opt_slvr_create.slc_slices = slc_slices;
            
            
            % Loop for all orientations
            %--------------------------------------------------------------
            for oo=1:num_orient
                
                count_obj = count_obj + 1;
                
                
                % Create slover object
                %----------------------------------------------------------
                opt_slvr_create.orient = orient{oo};
                [mn_mx, slvr_obj] = safe_slover_create(opt_slvr_create);
                
                
                % Save figure
                %----------------------------------------------------------
                if save_slc
                    opt_slvr_save.save_movie = 0;
                    opt_slvr_save.slc_output = fullfile(fpath,...
                        sprintf('%s_%s_Ftest_P%g_%dvox_%s_%.1fs_%s%s.%s',...
                        sID,contr_name,P,ET,adj_name,delay,orient{oo}(1:3),mask_suffix,img_ext));
					if save_movie
						opt_slvr_save.write_txt = 0;
					else
						opt_slvr_save.write_txt = 1;
						opt_slvr_save.txt       = [sID ' ' contr_name];
					end
                    safe_slover_save(slvr_obj, opt_slvr_save)
                end
                
                
                % Get data to create the movie
                %----------------------------------------------------------
                if save_movie
%                     safe_res.movie_file{count_obj,1} = fullfile(moviedir,...
%                         sprintf('%s_%s_Ftest_P%g_%dvox_%s_%s.gif',sID,contr_name,P,ET,adj_name,orient{oo}(1:3)));
                     safe_res.movie_file{count_obj,1} = fullfile(moviedir,...
                         sprintf('%s_%s_Ftest_P%g_%dvox_%s_%s.avi',sID,contr_name,P,ET,adj_name,orient{oo}(1:3)) );
                end
                
                safe_res.obj{count_obj} = slvr_obj;
                
            end
            
            % Color scale
            safe_res.mn_mx(count_obj/num_orient,:) = [mn_mx(2,:) 0 0];
            
            
        % Positive or negative BOLD
        %------------------------------------------------------------------
        else
            
            pos_neg = pos_neg + 1;
            if strcmp(contr_type, 'posBOLD')
                pos_neg_img{1} = fullfile(fpath, [thr_name '.img']);
            else % 'negBOLD'
                pos_neg_img{2} = fullfile(fpath, [thr_name '.img']);
            end
            
            if pos_neg==2 % created thresholded images for positive and negative BOLD
                
                pos_neg = 0;
                
                
                % Create slover object
                %----------------------------------------------------------
                opt_slvr_create.imgs       = {underlay, pos_neg_img{1}, pos_neg_img{2}};
                opt_slvr_create.type       = {'Structural','Blobs','Negative blobs'};
                opt_slvr_create.cmap       = {'', 'hot', 'winter'};
                opt_slvr_create.slc_step   = slc_step;
                opt_slvr_create.slc_nmbr   = slc_nmbr;
                opt_slvr_create.slc_slices = slc_slices;
                
                
                % Loop for all orientations
                %----------------------------------------------------------
                for oo=1:num_orient
                    
                    count_obj = count_obj + 1;
                    
                    
                    % Create slover object
                    %------------------------------------------------------
                    opt_slvr_create.orient = orient{oo};
                    [mn_mx, slvr_obj] = safe_slover_create(opt_slvr_create);
                    
                    
                    % Save figure
                    %------------------------------------------------------
                    if save_slc
                        opt_slvr_save.save_movie = 0;
                        opt_slvr_save.slc_output = fullfile(fpath,...
                            sprintf('%s_%s_pos-negBOLD_P%g_%dvox_%s_%.1fs_%s%s.%s',...
                            sID,contr_name,P,ET,adj_name,delay,orient{oo}(1:3),mask_suffix,img_ext));
                        %opt_slvr_save.write_txt = 0;
                        opt_slvr_save.write_txt = 1;
                        opt_slvr_save.txt       = [sID ' ' contr_name];
                        safe_slover_save(slvr_obj, opt_slvr_save)
                    end
                    
                    
                    % Get data to create the movie
                    %------------------------------------------------------
                    if save_movie
%                         safe_res.movie_file{count_obj,1} = fullfile(moviedir,...
%                             sprintf('%s_%s_pos-negBOLD_P%g_%dvox_%s_%s.gif',sID,contr_name,P,ET,adj_name,orient{oo}(1:3)));
                        safe_res.movie_file{count_obj,1} = fullfile(moviedir,...
                            sprintf('%s_%s_pos-negBOLD_P%g_%dvox_%s_%s.avi',sID,contr_name,P,ET,adj_name,orient{oo}(1:3)));
                    end
                    
                    
                    safe_res.obj{count_obj} = slvr_obj;
                    
                end
                
                % Color scale
                safe_res.mn_mx(count_obj/num_orient,:) = [mn_mx(2,:) mn_mx(3,:)];
                
            end
        end
    end
end

if save_ps
    warning on MATLAB:xlswrite:AddSheet % restore warning state
    xls_delete_sheets([ps_file(1:end-2) 'xls'],[1 2 3]) % delete the first 3 sheets
    xls_protect_sheets([ps_file(1:end-2) 'xls'],'all','protect','') % protect all sheets
end

% Delete extra elements
%--------------------------------------------------------------------------
if save_slc || save_movie || save_slover
    safe_res.obj((count_obj+1):end) = [];
    safe_res.mn_mx((count_obj/num_orient+1):end,:) = [];
end


if exist('safe_res','var')~=1 % only saved the PS file
    safe_res = 'PS file was saved'; % return a dummy value
end

safe_stat_results_finish

end


% =========================================================================
%                           AUXILIARY FUNCTIONS
%==========================================================================

function stop = test_time(SPM,xSPM,max_time)
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

end

%==========================================================================

function safe_stat_results_finish

% Close SPM windows
tmp = findobj('Tag','Graphics');
close(tmp)
tmp = findobj('Tag','Interactive');
close(tmp)

tmp = clock;
fprintf('(%.2d:%.2d:%02.0f) %s.m done!\n',...
    tmp(4),tmp(5),tmp(6),mfilename)
disp('-------------------------------------------------------------------')

end

%==========================================================================

function Vo = spm_write_filtered_mod(Z, XYZ, DIM, M, descrip, F)
% 
% Modified version of spm_write_filtered.m. Using zeros instead of NaN's
% because the slice overlay (slover) image has fewer and smaller blobs with
% NaN's.

% 
% Writes the filtered SPM as an image
% FORMAT V0 = spm_write_filtered(Z,XYZ,DIM,M,descrip,F)
%
% Z       - {1 x ?} vector point list of SPM values for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (voxel coordinates)
% DIM     - image dimensions {voxels}
% M       - voxels -> mm matrix [default: spm_matrix(-(DIM+1)/2)]
% descrip - description string [default: 'SPM-filtered']
% F       - output file's basename [default: user query]
%
% FORMAT V0 = spm_write_filtered(xSPM)
%
% xSPM    - SPM results structure from spm_getSPM
%
% Vo      - output image volume information
%__________________________________________________________________________
%
% spm_write_filtered takes a pointlist image (parallel matrices of
% co-ordinates and voxel intensities), and writes it out into an image
% file.
%
% It is intended for writing out filtered SPM's from the results
% section of SPM, but can be used freestanding.
%__________________________________________________________________________
% Copyright (C) 1996-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_write_filtered.m 4351 2011-06-13 17:18:25Z ged $

%-Get filename
%--------------------------------------------------------------------------
F = spm_str_manip(F, 'sd');
if isempty(F), F = 'output'; end
F = [F '.img'];
spm('Pointer', 'Watch')

%-Set up header information
%--------------------------------------------------------------------------
Vo = struct(...
    'fname'  , F,...
    'dim'    , DIM',...
    'dt'     , [spm_type('float32') spm_platform('bigend')],...
    'mat'    , M,...
    'descrip', descrip);
if all(Z==1) % binary map
    Vo.dt(1) = spm_type('uint8');
elseif all(ismember(Z,0:max(Z))) % n-ary map
    Vo.dt(1) = spm_type('uint16');
end

%-Reconstruct (filtered) image from XYZ & Z pointlist
%--------------------------------------------------------------------------
Y      = zeros(DIM(1:3)'); %nan(DIM(1:3)');
OFF    = XYZ(1,:) + DIM(1)*(XYZ(2,:)-1 + DIM(2)*(XYZ(3,:)-1));
Y(OFF) = Z.*(Z > 0);
    
%-Write the reconstructed volume
%--------------------------------------------------------------------------
Vo = spm_write_vol(Vo, Y);
spm('alert"', {'Written:',['    ',spm_select('CPath',F)]},mfilename,1);

%-End
%--------------------------------------------------------------------------
spm('Pointer', 'Arrow');

end
