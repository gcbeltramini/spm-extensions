function safe_event_seq(safe_res, opt_event_seq)
% 
% SAFE_EVENT_SEQ creates figures of the conditions for various delays.
% 
% 
%USAGE
%-----
% safe_event_seq(safe_res, opt_event_seq)
% 
% 
%INPUT
%-----
% - SAFE_RES: full path for the file that is output of safe_stat_results
% for various delays
% - OPT_EVENT_SEQ: structure containing the following fields
%   - PATH         : where the figures will be saved (full path)
%   - IMG_EXT      : image extension ('bmp', 'fig', 'jpg', 'png', 'tif')
%   - THR_PREF     : prefix for the thresholded images (used to edit the
%   file names)
%   - SCALE.FTEST  : 'max' (the first value will be 10^(-10) and the last
%   value will be the maximum value from all delays); or 1x2 matrix
%   containing the first and last value in the colormap (the first value
%   should be 10^(-10))
%   - SCALE.POSBOLD: analogous
%   - SCALE.NEGBOLD: analogous
% 
% 
%OUTPUT
%------
% - Inside OPT_EVENT_SEQ.PATH: images for every slover object in SAFE_RES,
% that is: SID_COND_TEST_P0.XXXX_Vvox_CORR_DD_DELAY.IMG_EXT, where
%   - SID   : subject ID
%   - COND  : condition name
%   - TEST  : "Ftest" or "posnegBOLD"
%   - 0.XXXX: P value
%   - V     : extent threshold
%   - CORR  : 'unc' or 'FWE'
%   - DD    : number of the delay (01, 02, 03, ...)
%   - DELAY : time in seconds from the stimulus onset to the beginning of
%     the hemodynamic response function
%   - IMG_EXT: chosen image extension
% 
% 
% Default values:
% --------------
% safe_res                    = --
% opt_event_seq.path          = --
% opt_event_seq.img_ext       = 'png';
% opt_event_seq.thr_pref      = 'Thresh_';
% opt_event_seq.scale.ftest   = 'max';
% opt_event_seq.scale.posbold = 'max';
% opt_event_seq.scale.negbold = 'max';
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2014-Jul-14


zero_thresh = 10^(-10);


% Read options
%--------------------------------------------------------------------------
fpath         = opt_event_seq.path;
img_ext       = opt_event_seq.img_ext;
thr_pref      = opt_event_seq.thr_pref;
scale.ftest   = opt_event_seq.scale.ftest;
scale.posbold = opt_event_seq.scale.posbold;
scale.negbold = opt_event_seq.scale.negbold;
visible       = 'off'; % figure visibility: 'on' or 'off'


fig = spm_figure('GetWin','Graphics'); % SPM graphics window
bg_color = [0 0 0];
set(fig, 'Color'  , bg_color)          % set background to black
set(fig, 'Visible', visible)           % figure visibility


DD      = length(safe_res);        % no. of delays
OO      = length(safe_res(1).obj); % no. of slover objects
num_ori = size(safe_res(1).obj,1)/size(safe_res(1).mn_mx,1); % no. of orientations


% Adjust color scale in the case of 'max'
%--------------------------------------------------------------------------
if any([strcmp(scale.ftest,'max'),strcmp(scale.posbold,'max'),strcmp(scale.negbold,'max')])
    max_stat = zeros(size(safe_res(1).mn_mx)) + zero_thresh;
    for dd=1:DD % loop for the delays
        try
            max_stat = max(safe_res(dd).mn_mx,max_stat);
        catch
            % in case there is no field mn_mx or this field is empty or
            % with the wrong dimensions
            continue
        end
    end
end


wb        = waitbar(0,'Saving images...','Name','Sequence of delays');
percent   = 0;
num_digit = numel(num2str(DD));
for dd = 1:DD % loop for the delays
    for oo = 1:OO % loop for the slover objects
        
        % Adjust file name
        %------------------------------------------------------------------
        [tmp,fname] = fileparts(safe_res(dd).obj{oo}.img(2).vol.fname);
        % "fname" will correspond only to F or positive T test
        tmp         = strfind(fname,'_');
        if ~strcmp(fname((end-4):end),'_mask') % masking was not applied
            tmp = tmp(end);
        else
            tmp = tmp(end-1);
        end
        % insert delay number to show files in the correct order
        fname = [fname(1:tmp) sprintf('ori%d_',rem(oo-1,num_ori)+1) ...
            sprintf(sprintf('%%.%dd_',num_digit),dd) fname((tmp+1):end)];
        
        
        % Adjust color scale
        %------------------------------------------------------------------
        
        tmp = strfind(fname,'_Ftest_');
        idx = floor((oo-1)/num_ori)+1;
            
        if ~isempty(tmp) % F test
            
            % File name
            %----------
            %tmp = tmp(end);
            fname = [fname((size(thr_pref,2)+1):end) '.' img_ext];
            
            if strcmp(opt_event_seq.scale.ftest,'max')
                safe_res(dd).obj{oo}.img(2).range = max_stat(idx,1:2);
            else
                safe_res(dd).obj{oo}.img(2).range = opt_event_seq.scale.ftest;
            end
            
        else % T test
            
            % File name
            %----------
            tmp   = strfind(fname,'_posBOLD_');
            tmp   = tmp(end);
            fname = [fname((size(thr_pref,2)+1):tmp) 'posneg' fname((tmp+4):end) '.' img_ext];
            
            if strcmp(opt_event_seq.scale.posbold,'max')
                safe_res(dd).obj{oo}.img(2).range = max_stat(idx,1:2);
            else
                safe_res(dd).obj{oo}.img(2).range = opt_event_seq.scale.posbold;
            end
            
            if strcmp(opt_event_seq.scale.negbold,'max')
                safe_res(dd).obj{oo}.img(3).range = max_stat(idx,3:4);
            else
               safe_res(dd).obj{oo}.img(3).range = opt_event_seq.scale.negbold;
            end
            
        end
        
        
        % Create figure
        %------------------------------------------------------------------
        safe_res(dd).obj{oo}.figure = fig;
        paint(safe_res(dd).obj{oo});
        set(fig, 'Visible', visible) % figure visibility
        %set(findobj('Parent',fig,'Tag','cbar'),...
        %    'FontSize',.1,'FontWeight','bold') % colorbar font
        
        
        % Save figure
        %------------------------------------------------------------------
        %set(fig,'PaperPositionMode','auto') % this line corrects the fact that sometimes the output figure is cut
        drawnow; pause(.1)
        %saveas(fig, fullfile(fpath, fname))
        fig_print(fig, fullfile(fpath,fname), [1080 1557], 1, 1, bg_color)
        
        
        % Update the progress bar
        %------------------------------------------------------------------
        percent = percent + 1;
        waitbar(percent/OO/DD, wb)
        
    end
end

close(wb)
close(fig)