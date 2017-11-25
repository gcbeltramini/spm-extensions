function safe_res_group = safe_groupthresh(opt_group)
% 
% SAFE_GROUPTHRESH creates combined images from thresholded functional maps
% 
% 
%USAGE
%-----
% safe_res_group = safe_groupthresh(opt_group)
% 
% 
%INPUT
%-----
% OPT_GROUP: structure with the following fields:
% - PATH    : subject folder where the analysis for all delays are
%   Folder where the thresholded images must be:
%   "PATH"\"SID"-SPM_delay_"DELAY"
% - SID     : subject ID
% - DELAY   : vector of any real number (positive, negative or zero)
% - NGROUP  : number of delays per group (any natural number) (default: 3)
% - THR_PREF: prefix of the thresholded images (default: 'Thresh_')
% - THR_DIR : directory where the combined images will be saved
% - SAVE_SLC: 1 (save slices) or 0 (don't save)
% If SAVE_SLC=1:
%   - IMG_EXT : image extension, e.g., 'bmp', 'fig', 'jpg', 'png', 'tif'
%     (default: 'png')
%   - UNDERLAY: underlay image
%   - SLC_STEP: spacing between slices (mm) (any positive number)
%     (default: 3)
%   - SLC_NMBR: number of slices (to use this option, SLC_STEP must be 0)
%     (default: 46)
%   - SLC_SLICES: coordinates of the slices (to use this option, SLC_STEP
%     and SLC_NMBR must be 0)
%   - SCALE.ADJ: 1 (adjust scale) or 0 (don't adjust)
%   If SCALE.ADJ=1:
%     - SCALE.FTEST: 'max' (the first value will be 10^(-10) and the last
%       value will be the maximum value from all delays); or 1x2 matrix
%       containing the first and last value in the colormap (the first
%       value should be 10^(-10)) (default: 'max')
%     - SCALE.TTEST: 'max' (the first value will be 10^(-10) and the last
%       value will be the maximum value from all delays); or  1x4 matrix
%       containing the first and last value in the colormap for positive
%       (SCALE.TTEST([1 2])) and negative (SCALE.TTEST([3 4])) BOLD
%       (SCALE.TTEST(1) and SCALE.TTEST(3) should be 10^(-10)) (default: 'max')
% 
% 
%OUTPUT
%------
% - Inside THR_DIR:
%   - Thresholded images: THR_PREF_SID_CONDNAME_TEST_..._groupXX, where
%   CONDNAME is the condition name, TEST is "Ftest" or "posBOLD" or
%   "negBOLD", and XX is the group number (01, 02, ...). These temporary
%   files are deleted.
%   - If SAVE_SLC=1:
%     - Images: SID_CONDNAME_TEST_...groupXX_ORI.IMG_EXT, where TEST is
%     "Ftest" or "posnegBOLD", and ORI is the orientation ("axi", "sag" or
%     "cor")
% - SAFE_RES_GROUP: structure array to be used with SAFE_PLOT_DELAYS
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2013-Jul-01


zero_thresh = 10^(-10);


% Read options
%==========================================================================
fpath    = opt_group.path;
sID      = opt_group.sID;
delay    = opt_group.delay(:);
Ngroup   = opt_group.Ngroup;
thr_pref = opt_group.thr_pref;
thr_dir  = opt_group.thr_dir;
save_slc = opt_group.save_slc;
if save_slc
    underlay   = opt_group.underlay;
    img_ext    = opt_group.img_ext;
    slc_step   = opt_group.slc_step;
    slc_nmbr   = opt_group.slc_nmbr;
    slc_slices = opt_group.slc_slices;
    orient     = opt_group.orient;
    num_orient = length(orient);
    scale      = opt_group.scale;
end


if exist(thr_dir,'dir')~=7
    fprintf('%s does not exist. Nothing will be done.\n',thr_dir)
    return
end


tmp = clock;
disp('-------------------------------------------------------------------')
fprintf('(%.2d:%.2d:%02.0f) Running %s.m...\n',...
    tmp(4),tmp(5),tmp(6),mfilename)


% Number of groups and group size
%==========================================================================
Ndelays = length(delay);
if Ngroup<Ndelays
    if Ngroup>Ndelays/2
        group_n = 2;
        group_sz = [Ngroup Ndelays-Ngroup];
    else
        group_n  = floor( (Ndelays+1)/Ngroup ); % number of groups
        group_sz = [repmat(Ngroup,1,group_n-1) Ndelays-(group_n-1)*Ngroup]; % size of each group
        if group_sz(end)==0
            group_n       = group_n - 1;
            group_sz(end) = [];
        end
    end
else
    group_n  = 1;
    group_sz = Ndelays;
end
% The code below doesn't work in the case: delay=(-9:2:9), Ngroup=3,
% because it results in: group_n=3 & group_sz=[3 3 3 1]
% tmp     = Ndelays/Ngroup;
% group_n = floor(tmp); % number of groups
% if group_n==tmp
%     group_sz = repmat(Ngroup,1,floor(tmp)); % size of each group
% else
%     group_sz = [repmat(Ngroup,1,floor(tmp)) Ndelays-floor(tmp)*Ngroup];
% end


% Get names of folders for each group
%==========================================================================
gr_folders = cell(group_n,1);
counter    = 0;
for gg=1:group_n % loop for the groups
    gr_folders{gg} = cell(group_sz(gg),1);
    for gsz=1:group_sz(gg) % loop for the delays inside each group
        counter = counter + 1;
        gr_folders{gg}{gsz} = fullfile(fpath,sprintf('%s-SPM_delay_%.1f',sID,delay(counter)));
    end
end


% Get all thresholded images
%==========================================================================
% imgs{d}{T(d),[1 2]}:
% - d: delay
% - T(d): number of thresholded images as a function of the delay
% - imgs{d}{T(d),1}: thresholded image name
% - imgs{d}{T(d),2}: condition name
curr_dir = pwd;
imgs     = cell(length(delay),1);
count    = 0;
for gg=1:group_n
    for gsz=1:group_sz(gg)
        
        cd(gr_folders{gg}{gsz})
        count = count + 1;
        
        % Image name
        %-----------
        tmp = dir(sprintf('%s*_mask.img',thr_pref));
        if isempty(tmp) % masking was not applied
            mask = 0;
            mask_suffix = [];
            tmp = dir(sprintf('%s*s.img',thr_pref));
        else
            mask = 1;
            mask_suffix = '_mask';
        end
        if size(tmp,1)==0 % no thresholded image
            continue
        end
        imgs{count} = cell(size(tmp,1),4);
        % col. 1: filename; col. 2: filename without delay; col. 3: type of test; col. 4: condition
        [imgs{count}{:,1}] = deal(tmp.name);
        % names will come in alphabetic order: first condition name, then
        % test ('Ftest', 'negBOLD', 'posBOLD')
        
        % Condition name & Part of filename
        %----------------------------------
        for contr=1:size(imgs{count},1)
            
            tmp = strfind(imgs{count}{contr,1},'_');
            if ~mask
                tmp = tmp(end) - 1;
            else
                tmp = tmp(end-1) - 1;
            end
            imgs{count}{contr,2} = [imgs{count}{contr,1}(1:tmp) mask_suffix];
            
            tmp = strfind(imgs{count}{contr,1},'_Ftest_');
            if ~isempty(tmp) % F test
                tmp = tmp(end);
                imgs{count}{contr,3} = 'Ftest';
            else
                tmp = strfind(imgs{count}{contr,1},'_negBOLD_');
                if ~isempty(tmp) % T test (negative BOLD)
                    tmp = tmp(end);
                    imgs{count}{contr,3} = 'Tneg';
                else
                    tmp = strfind(imgs{count}{contr,1},'_posBOLD_');
                    if ~isempty(tmp) % T test (positive BOLD)
                        tmp = tmp(end);
                        imgs{count}{contr,3} = 'Tpos';
                    else
                        % something wrong happened
                    end
                end
            end
            imgs{count}{contr,4} = imgs{count}{contr,1}((length(thr_pref)+length(sID)+2):(tmp-1));
        end
    end
end


% Create thresholded group images
%==========================================================================
curr_delay = 1;
if save_slc
    count_imgs  = 0;
    thresh_imgs = cell(length(imgs{1})*group_n,4);
    % col. 1: file name; col. 2: type of test; col. 3: condition; col. 4: text for the image
    % initialize "thresh_imgs" with an approximate number of cells
end

safe_res_group = cell(group_n,1);
safe_res_group = struct('thr',safe_res_group, 'mn_mx',safe_res_group);
for gg=1:group_n
    
    if gg>1
        curr_delay = group_sz(gg-1) + curr_delay; 
    end
    
    % Find the first delay with conditions for this group
    %----------------------------------------------------
    for dd = curr_delay : (group_sz(gg)+curr_delay-1) % loop for all thresholded images across delays
        if size(imgs{dd},1)~=0 % there are conditions for this delay
            break
        else
            curr_delay   = curr_delay + 1;
            group_sz(gg) = group_sz(gg) - 1;
        end
    end
    if group_sz(gg)==0 % last delay is empty => go to the next group
        continue
    end
    
    
    % Text for the description
    %-------------------------
    descrip = 'Combined image for delays (s): ';
    for dd = curr_delay : (group_sz(gg)+curr_delay-1)
        descrip = sprintf('%s%d, ',descrip,delay(dd));
    end
    descrip(end-1:end) = [];
    
    
    % Get maximum data and initialize the combined image header
    %----------------------------------------------------------
    Nmax = cell(size(imgs{curr_delay},1),1);
    Ncmb = Nmax;
    safe_res_group(gg).thr.fname    = cell(size(imgs{curr_delay},1),1);
    safe_res_group(gg).thr.add_stat = zeros(size(imgs{curr_delay},1),1);
    safe_res_group(gg).thr.voxels   = zeros(size(imgs{curr_delay},1),1);
    safe_res_group(gg).mn_mx        = zeros(size(imgs{curr_delay},1),4);
    mn_mx_row = 1;
    for contr=1:size(imgs{curr_delay},1) % loop for all contrasts
        
        % Initialize with the first delay
        %--------------------------------
        
        % Find the first non-empty delay for this group
        tmp = curr_delay - sum(group_sz(1:gg-1));
        tmp = fullfile(gr_folders{gg}{tmp},imgs{curr_delay}{contr,1});
        if exist(tmp,'file')~=2 % file does not exist
            continue
        end
        N = nifti(tmp);
        Ncmb{contr} = N; % combined image
        Ncmb{contr}.dat.fname = fullfile(thr_dir,sprintf('%s_group%.2d.nii',imgs{curr_delay}{contr,2},gg));
        Nmax{contr} = N.dat(:,:,:); % maximum image
        Nmax{contr}(isnan(Nmax{contr})) = 0; % NaN -> 0
                
        % Update the maximum values
        %--------------------------
        for dd = curr_delay+1 : (group_sz(gg)+curr_delay-1)
            
            % Find the same contrast for each delay
            %--------------------------------------
            tmp = strfind(imgs{dd}(:,2),imgs{curr_delay}{contr,2});
            if isempty(tmp)
                continue
            end
            tmp = imgs{dd}{~cellfun(@isempty,tmp),1};
            % in case of more than 1 match, get only the first one
            
            N = nifti(fullfile(gr_folders{gg}{dd-curr_delay+1},tmp));
            N = N.dat(:,:,:);
            N(isnan(N)) = 0; % NaN -> 0
            tmp = N>Nmax{contr}; % tmp = abs(N)>abs(Nmax{contr});
            Nmax{contr} = N.*tmp + Nmax{contr}.*(~tmp);
        end
        
        % Create image with the maximum values
        %-------------------------------------
        Ncmb{contr}.descrip = descrip;
        create(Ncmb{contr})
        Ncmb{contr}.dat(:,:,:) = Nmax{contr};
        
        % Output
        %-------
        safe_res_group(gg).thr.fname{contr}    = Ncmb{contr}.dat.fname;
        safe_res_group(gg).thr.add_stat(contr) = sum(sum(sum(Nmax{contr},1),2),3);
        safe_res_group(gg).thr.voxels(contr)   = sum(sum(sum(Nmax{contr}>0,1),2),3);
        % Because of the command "dir", the names are in alphabetic order:
        % first condition name, then test ('Ftest', 'negBOLD', 'posBOLD').
        % So when it is a negative T-test, the following contrast is the
        % positive T-test for the same condition name. That's why the line
        % in "mn_mx" is only skipped after an F or positive T-test.
        switch imgs{curr_delay}{contr,3}
            case 'Ftest'
                safe_res_group(gg).mn_mx(mn_mx_row,2) = max(max(max(Nmax{contr},[],1),[],2),[],3);
                mn_mx_row = mn_mx_row + 1;
            case 'Tneg'
                safe_res_group(gg).mn_mx(mn_mx_row,4) = max(max(max(Nmax{contr},[],1),[],2),[],3);
            case 'Tpos'
                safe_res_group(gg).mn_mx(mn_mx_row,2) = max(max(max(Nmax{contr},[],1),[],2),[],3);
                mn_mx_row = mn_mx_row + 1;
        end
        
        if save_slc
            count_imgs = count_imgs + 1;
            thresh_imgs{count_imgs,1} = Ncmb{contr}.dat.fname;
            thresh_imgs{count_imgs,2} = imgs{curr_delay}{contr,3};
            thresh_imgs{count_imgs,3} = imgs{curr_delay}{contr,4};
            
            % Text in the figure
            %-------------------
            tmp = sum(group_sz(1:(gg-1)));
            tmp = delay((tmp+1):(tmp+group_sz(gg))); % delay must be a column vector
            tmp = num2str(tmp);
            txt = [];
            for ii=1:size(tmp,1) % remove excess of spaces (problem with non-integer numbers)
                txt = [txt strtrim(tmp(ii,:)) ' ']; % or txt = [txt strtrim(tmp(ii,:)) ','];
            end
            txt = txt(1:(end-1)); % remove last space
            txt = [thresh_imgs{count_imgs,3} ' (' txt ' s)'];
            thresh_imgs{count_imgs,4} = txt;
        end
        
    end
    safe_res_group(gg).mn_mx(mn_mx_row:end,:) = [];
    
end


% Create group images
%==========================================================================
if save_slc
    
    % Clear unused space in "thresh_imgs"
    %------------------------------------
    thresh_imgs(cellfun(@isempty,thresh_imgs(:,1)),:) = [];
    
    
    % Initialize
    %-----------
    count_obj  = 0;
    sz_thr_dir = size(thr_dir,2);
    tmp        = cell(size(thresh_imgs,1)*num_orient,1);
    output     = struct(...
        'cond_name',tmp,...
        'test'     ,tmp,...
        'mn_mx'    ,tmp,...
        'slvr_obj' ,tmp,...
        'slc_name' ,tmp);
     % Could multiply size by a factor 2/3 (F, posBOLD, negBOLD --> F,
     % posnegBOLD) or 1/2 (posBOLD, negBOLD --> posnegBOLD). Anyway, the
     % extra elements will be removed in the end
    
    
    % Save figure - options
    %----------------------
    opt_slvr_save.save_slc   = 1;
    opt_slvr_save.save_movie = 0;
    opt_slvr_save.write_txt  = 1;
    
    
    % Loop for the thresholded images
    %================================
    for img=1:size(thresh_imgs,1)
        
        thr_name = thresh_imgs{img,1};
        
        % Same issue regarding the alphabetic order, as pointed out above.
        switch thresh_imgs{img,2}
            
            case 'Ftest'
                
                % Create slover object - options
                %-------------------------------
                opt_slvr_create.imgs       = {underlay,thr_name};
                opt_slvr_create.type       = {'Structural','Blobs'};
                opt_slvr_create.cmap       = {'','hot'};
                opt_slvr_create.slc_step   = slc_step;
                opt_slvr_create.slc_nmbr   = slc_nmbr;
                opt_slvr_create.slc_slices = slc_slices;
                
                % Remove prefix
                %--------------
                thr_name = [thr_name(1:(sz_thr_dir+1)) ...
                    thr_name((sz_thr_dir+2+size(thr_pref,2)):(end-4))];
                
                % Loop for all orientations
                %--------------------------
                for oo=1:num_orient
                    
                    count_obj = count_obj + 1;
                    
                    % Create slover object
                    %---------------------
                    opt_slvr_create.orient = orient{oo};
                    [mn_mx,slvr_obj] = safe_slover_create(opt_slvr_create);
                    
                    output(count_obj).cond_name = thresh_imgs{img,3};
                    output(count_obj).test      = 'F';
                    output(count_obj).mn_mx     = mn_mx;
                    output(count_obj).slvr_obj  = slvr_obj;
                    output(count_obj).slc_name  = [thr_name '_' orient{oo}(1:3) '.' img_ext];
                    output(count_obj).txt       = thresh_imgs{img,4};
                    
                end
                
            case 'Tneg'
                
                thr_name_neg = thresh_imgs{img};
                
            case 'Tpos'
                
                % Create slover object - options
                %-------------------------------
                opt_slvr_create.imgs       = {underlay,thr_name,thr_name_neg};
                opt_slvr_create.type       = {'Structural','Blobs','Negative blobs'};
                opt_slvr_create.cmap       = {'','hot','winter'};
                opt_slvr_create.slc_step   = slc_step;
                opt_slvr_create.slc_nmbr   = slc_nmbr;
                opt_slvr_create.slc_slices = slc_slices;
                
                % Replace "posBOLD" with "posnegBOLD" and remove the prefix
                %----------------------------------------------------------
                tmp = strfind(thr_name,'_posBOLD_');
                tmp = tmp(end);
                thr_name = [thr_name(1:tmp) 'posneg' thr_name((tmp+4):(end-4))];
                thr_name = [thr_name(1:(sz_thr_dir+1)) ...
                    thr_name((sz_thr_dir+2+size(thr_pref,2)):end)];
                
                % Loop for all orientations
                %--------------------------
                for oo=1:num_orient
                    
                    count_obj = count_obj + 1;
                    
                    % Create slover object
                    %---------------------
                    opt_slvr_create.orient = orient{oo};
                    [mn_mx,slvr_obj] = safe_slover_create(opt_slvr_create);
                    
                    output(count_obj).cond_name = thresh_imgs{img,3};
                    output(count_obj).test      = 'T';
                    output(count_obj).mn_mx     = mn_mx;
                    output(count_obj).slvr_obj  = slvr_obj;
                    output(count_obj).slc_name  = [thr_name '_' orient{oo}(1:3) '.' img_ext];
                    output(count_obj).txt       = thresh_imgs{img,4};
                    
                end
            
        end
        
    end
    output(count_obj+1:end) = []; % remove unused elements
    
    
    % Adjust scale
    %======================================================================
    if scale.adj
        
        % Sort "output" by condition name
        %--------------------------------
        tmp       = cell(size(output,1),1);
        [tmp{:}]  = deal(output(:).cond_name);
        [tmp,idx] = sort(tmp);
        output    = output(idx);
        
        % Initialize variables for the loop
        %----------------------------------
        condF = '';
        condT = condF;
        maxF  = [zero_thresh 0]; % add small number to show the underlay image
        maxTp = [zero_thresh 0];
        maxTn = [zero_thresh 0];
        imgF  = [];
        imgT  = [];
        
        for img=1:count_obj
            
            switch output(img).test
                
                case 'F' % F test
                    
                    if strcmp(scale.ftest,'max') % get the maximum value
                        if strcmp(condF,output(img).cond_name) % same condition
                            maxF = max([output(img).mn_mx(2,:) ; maxF],[],1);
                            imgF = [imgF img];
                        else % another condition
                            for ii=imgF
                                output(ii).slvr_obj.img(2).range = maxF;
                            end
                            imgF  = img;
                            condF = output(img).cond_name;
                            maxF  = output(img).mn_mx(2,:) + [zero_thresh 0];
                        end
                    else
                        output(img).slvr_obj.img(2).range = scale.ftest;
                    end
                    
                case 'T' % positive and negative BOLD
                    
                    if strcmp(scale.ttest,'max') % get the maximum value
                        if strcmp(condT,output(img).cond_name) % same condition
                            maxTp = max([output(img).mn_mx(2,:) ; maxTp],[],1);
                            maxTn = max([output(img).mn_mx(3,:) ; maxTn],[],1);
                            imgT = [imgT img];
                        else % another condition
                            for ii=imgT
                                output(ii).slvr_obj.img(2).range = maxTp;
                                output(ii).slvr_obj.img(3).range = maxTn;
                            end
                            condT = output(img).cond_name;
                            imgT  = img;
                            maxTp = output(img).mn_mx(2,:) + [zero_thresh 0];
                            maxTn = output(img).mn_mx(3,:) + [zero_thresh 0];
                        end
                    else
                        output(img).slvr_obj.img(2).range = scale.ttest([1 2]); % positive BOLD
                        output(img).slvr_obj.img(3).range = scale.ttest([3 4]); % negative BOLD
                    end
                    
            end
        end
        
        % Assign the maximum value to the last condition
        if strcmp(scale.ftest,'max')
            for ii=imgF
                output(ii).slvr_obj.img(2).range = maxF;
            end
        end
        if strcmp(scale.ttest,'max')
            for ii=imgT
                output(ii).slvr_obj.img(2).range = maxTp;
                output(ii).slvr_obj.img(3).range = maxTn;
            end
        end
        
    end
    
    
    % Save figure
    %----------------------------------------------------------------------
    wb = waitbar(0,'Saving group images...',...
        'Name','Union of delays',...
        'CloseRequestFcn','');
    % will beep if 'WindowStyle' = 'modal' ("beep off" won't solve)
    percent = 0;
    fprintf('Saving group images...\n')
    for ii=1:count_obj
        opt_slvr_save.slc_output = output(ii).slc_name;
        opt_slvr_save.txt        = output(ii).txt;
        safe_slover_save(output(ii).slvr_obj,opt_slvr_save)
        percent = percent + 1;
        waitbar(percent/count_obj,wb)
    end
    %close(wb) % won't work because 'CloseRequestFcn' = ''
    delete(wb)
    
    % Delete thresholded images
    %--------------------------
%     for ii=1:size(thresh_imgs,1)
%         delete(thresh_imgs{ii,1})
%     end
    
end


% Finish
%==========================================================================
cd(curr_dir)

% Close SPM windows
tmp = findobj('Tag','Graphics');
close(tmp)
tmp = findobj('Tag','Interactive');
close(tmp)

tmp = clock;
fprintf('(%.2d:%.2d:%02.0f) %s.m done!\n',...
    tmp(4),tmp(5),tmp(6),mfilename)
disp('-------------------------------------------------------------------')