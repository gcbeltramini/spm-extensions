function safe_plot_delays(opt_plot_delays)
% 
% SAFE_PLOT_DELAYS creates the figures of some image statistics as a
% function of the delays or groups.
% 
% 
%USAGE
%-----
% safe_plot_delays(opt_plot_delays)
% 
% 
%INPUT
%-----
% OPT_PLOT_DELAYS: structure with the following fields:
% - DELAY: time in seconds from the stimulus onset to the beginning of the
%   hemodynamic response function. DELAY is a vector of any real number
%   (positive, negative or zero). If DELAY=[], SAFE_PLOT_DELAYS will create
%   the graphics for groups and not delays.
% - SAFE_RES: output of safe_stat_results.m. The necessary fields are, for
%   delay DD and condition CC:
%   - safe_res(dd).thr.fname{cc}: thresholded file
%   - safe_res(dd).thr.add_stat : (no. conditions)x1 vector with the sum of
%   the statistics for the whole image
%   - safe_res(dd).thr.voxels   : (no. conditions)x1 vector with the number
%   of suprathreshold voxels
%   - safe_res(dd).mn_mx        : for a certain row, if the 4th column has
%   a zero value, the row contains in the 1st column the minimum and the
%   2nd column the maximum value for an F test; otherwise, the 1st and 2nd
%   column are for a positive T test, and the 3rd and 4th column for a
%   negative T test
% - PATH_RES: path where the results will be saved
% - THR_PREF: prefix for the thresholded images (used to edit the file names)
% 
% 
%OUTPUT
%------
% Inside PATH_RES:
% - Event_seq-Summary.txt or Group-Summary.txt: Condition name, Type of
%   data, Statistics, X and Y data for all the graphics (see below)
% - Folder "Plots" with:
%  XX_max.fig , XX_max_scaled.fig , XX_stat_sum.fig , XX_stat_sum_scaled.fig , XX_vox.fig , XX_vox_scaled.fig
%  XX_max.png , XX_max_scaled.png , XX_stat_sum.png , XX_stat_sum_scaled.png , XX_vox.png , XX_vox_scaled.png
%   If the bar plot option is enabled (via script change), "Plot" will
%   contain the folders "Line_plots" and "Bar_plots", both wilh the
%   corresponding figures.
% 
% 
% - XX     : "subject ID"_"condition name"
% - XX_max*: Delay/Group vs. Maximum value for the statistical test (F & T)
% - XX_vox*: Delay/Group vs. No. of suprathreshold voxels for the statistical test
% - XX_stat_sum*: Delay/Group vs. Sum of the statistics for all suprathreshold voxels
% 
% 
%EXAMPLE
%-------
% opt_plot_delays.delay    = (-12:2:6);
% opt_plot_delays.safe_res = safe_res;
% opt_plot_delays.path_res = 'D:\fMRI\Data\Sbj23\Analysis\23-Results_-12.0s_6.0s\23-Event_sequence';
% opt_plot_delays.thr_pref = 'Thresh_';
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme C. Beltramini - 2013-Jul-02


% Read options
%==========================================================================
delay    = opt_plot_delays.delay(:);
safe_res = opt_plot_delays.safe_res;
path_res = opt_plot_delays.path_res;
thr_pref = opt_plot_delays.thr_pref;
visible  = 'off'; % figure visibility: 'on' or 'off'
if length(delay)>1
    bar_plot = 0;
else
    bar_plot = 1;
end


tmp = clock;
disp('-------------------------------------------------------------------')
fprintf('(%.2d:%.2d:%02.0f) Running %s.m...\n',...
    tmp(4),tmp(5),tmp(6),mfilename)


curr_dir = pwd;
cd(path_res)
create_folder('Plots')
cd('Plots')
if bar_plot
    create_folder('Line_plots')
    create_folder('Bar_plots')
end

movie = isfield(safe_res,'movie_file') && ~isempty(safe_res(1).movie_file);


% Initialize variables
%==========================================================================

% Number of delays/groups
%------------------------
num_delays = length(safe_res);
if num_delays==0 % something wrong happened
    return
end

% Groups or delays
%-----------------
if isempty(delay)
    group = 1;
    delay = 1:num_delays;
else
    group = 0;
end

% Number of conditions (F/T tests and different conditions)
%----------------------------------------------------------
try
    num_cond = size(safe_res(1).thr.fname,1);
catch ME % something wrong happened in the results
    error('Variable "safe_res": %s',ME.message)
end

plot_add_stat = zeros(num_cond,num_delays);
plot_voxels   = plot_add_stat;
plot_mn_mx    = plot_add_stat;


% Read the conditions
%==========================================================================
plot_fname = cell(num_cond,1);
conds      = cell(num_cond,3);
% 1st col.: condition name; 2nd col.: test (F, positive T or negative T);
% 3rd col.: color for the graphic

for cc=1:num_cond % loop for the conditions
    
    % Figure name
    %------------
    tmp = strfind(safe_res(1).thr.fname{cc},filesep);
    % not using "fileparts" because of the '.' in the delay
    tmp = tmp(end);
    plot_fname{cc} = safe_res(1).thr.fname{cc}(tmp+1:end); % get only file name
    plot_fname{cc} = plot_fname{cc}(length(thr_pref)+1:end); % remove prefix
    tmp = strfind(plot_fname{cc},'_');
    if ~strcmp(plot_fname{cc}((end-4):end),'_mask') % masking was not applied
        tmp = tmp(end);
    else
        tmp = tmp(end-1);
    end
    plot_fname{cc} = plot_fname{cc}(1:tmp-1); % remove delay/group and mask suffix
    
    % Conditions
    %-----------
    tmp = strfind(plot_fname{cc},'_Ftest_');
    if ~isempty(tmp)
        
        % F test
        %-------
        conds{cc,1} = plot_fname{cc}(1:tmp(end)-1);
        conds{cc,2} = 'F test';
        conds{cc,3} = 'k';
        
    else
        tmp = strfind(plot_fname{cc},'_posBOLD_');
        if ~isempty(tmp)
            
            % T test (positive)
            %------------------
            conds{cc,1} = plot_fname{cc}(1:tmp(end)-1);
            conds{cc,2} = 'T (positive BOLD)';
            conds{cc,3} = 'r';
            
        else
            tmp = strfind(plot_fname{cc},'_negBOLD_');
            if ~isempty(tmp)
                
                % T test (negative)
                %------------------
                conds{cc,1} = plot_fname{cc}(1:tmp(end)-1);
                conds{cc,2} = 'T (negative BOLD)';
                conds{cc,3} = 'b';
                
            else
                % something wrong happened
                return
            end
        end
    end
    
end


% Read the data
%==========================================================================
for dd=1:num_delays % loop for the delays/groups
    
    plot_add_stat(:,dd) = safe_res(dd).thr.add_stat;
    plot_voxels(:,dd)   = safe_res(dd).thr.voxels;
    
    mn_mx_row = 1; % counter for the rows in mn_mx
    pos_neg   = 0; % counter for positive and negative T-tests
    for cc=1:num_cond
        switch conds{cc,2}
            case 'F test'
                tmp = safe_res(dd).mn_mx(mn_mx_row,2);
                mn_mx_row = mn_mx_row + 1;
            case 'T (negative BOLD)'
                tmp = safe_res(dd).mn_mx(mn_mx_row,4);
                pos_neg = pos_neg + 1;
            case 'T (positive BOLD)'
                tmp = safe_res(dd).mn_mx(mn_mx_row,2);
                pos_neg = pos_neg + 1;
        end
        if tmp==1 % no suprathreshold voxels 
            %tmp = NaN;
            tmp = 0;
        end
        
        plot_mn_mx(cc,dd) = tmp;
        
        if pos_neg==2 % positive and negative BOLD for the same condition
            mn_mx_row = mn_mx_row + 1;
            pos_neg   = 0;
        end
        
    end
end


% Plot graphics & Save the results
%==========================================================================

% Initialize variables
%---------------------
curr_cond  = {'' 0}; % condition name and first row for the condition
movie_data = delay;  % initialize data to write in the file

% The first 6 figures are line plots; the other 3 are bar plots:
fig_names = {'vox','vox_scaled',...
    'max','max_scaled',...
    'stat_sum','stat_sum_scaled',...
    'vox','max','stat_sum'};
if ~bar_plot
    fig_names(7:9) = [];
end
num_figs  = length(fig_names); % number of figures
line_plot = [ones(6,1) ; zeros(num_figs-6,1)];
figs      = zeros(num_figs,1);

% Prepare summary file
%---------------------

% Create file
if ~group
    fid = fopen(['..' filesep 'Event_seq-Summary.txt'],'wt');
else
    fid = fopen(['..' filesep 'Group-Summary.txt'],'wt');
end

% Header
if ~movie
    fprintf(fid,'Condition\tType\tStatistics');
    for dd=1:length(delay)
        if ~group
            fprintf(fid, '\t%.1fs', delay(dd));
        else
            fprintf(fid, '\tGroup%.2d', delay(dd));
        end
    end
else
    if ~group
        fprintf(fid, 'Delay');
    else
        fprintf(fid,'Group');
    end
end

for cc=1:num_cond % loop for the conditions
    
% New condition
%==========================================================================
if ~strcmp(conds{cc,1},curr_cond{1})
    
    % Bar plot
    %---------
    if bar_plot
        plot_voxels_bar.val     = zeros(num_delays,3) - 1; % 1 column for each test
        plot_mn_mx_bar.val      = plot_voxels_bar.val;
        plot_add_stat_bar.val   = plot_voxels_bar.val;
        plot_voxels_bar.color   = cell(1,3);
        plot_mn_mx_bar.color    = plot_voxels_bar.color;
        plot_add_stat_bar.color = plot_voxels_bar.color;
        bar_col = 1;
    end
    
    
    % 1) Number of voxels above threshold
    %----------------------------------------------------------------------
    
    % Line plot
    %----------
    figs(1) = figure('Name','Voxels','Visible',visible);
    plot(delay,plot_voxels(cc,:),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    title(conds{cc,1},'FontSize',14,'Interpreter','none')
    ylabel('Number of suprathreshold voxels', 'FontSize', 12)
    hold on
    
    figs(2) = figure('Name','Voxels_scaled','Visible',visible);
    plot(delay,plot_voxels(cc,:)/max(plot_voxels(cc,:)),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    title(conds{cc,1},'FontSize',14,'Interpreter','none')
    ylabel('Suprathreshold voxels (scaled)', 'FontSize', 12)
    hold on
    
    % Bar plot
    %---------
    if bar_plot
        figs(7) = figure('Name','Voxels','Visible',visible);
        set(figs(7),'UserData','Suprathreshold voxels')
        plot_voxels_bar.val(:,bar_col) = plot_voxels(cc,:);
        plot_voxels_bar.color{bar_col} = conds{cc,3};
    end
    
    
    % 2) Maximum value for the statistical test
    %----------------------------------------------------------------------
    figs(3) = figure('Name','Maximum','Visible',visible);
    plot(delay,plot_mn_mx(cc,:),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    title(conds{cc,1},'FontSize',14,'Interpreter','none')
    ylabel('Maximum value','FontSize',12)
    hold on
    
    figs(4) = figure('Name','Maximum_scaled','Visible',visible);
    plot(delay,plot_mn_mx(cc,:)/max(plot_mn_mx(cc,:)),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    title(conds{cc,1},'FontSize',14,'Interpreter','none')
    ylabel('Maximum value (scaled)','FontSize',12)
    hold on
    
    % Bar plot
    %---------
    if bar_plot
        figs(8) = figure('Name','Maximum','Visible',visible);
        set(figs(8),'UserData','Maximum value')
        plot_mn_mx_bar.val(:,bar_col) = plot_mn_mx(cc,:);
        plot_mn_mx_bar.color{bar_col} = conds{cc,3};
    end
    
    
    % 3) Sum of the statistics for the whole image
    %----------------------------------------------------------------------
    figs(5) = figure('Name','Sum','Visible',visible);
    plot(delay,plot_add_stat(cc,:),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    title(conds{cc,1},'FontSize',14,'Interpreter','none')
    ylabel('Sum of the statistics for the whole image','FontSize',12)
    hold on
    
    figs(6) = figure('Name','Sum_scaled','Visible',visible);
    plot(delay,plot_add_stat(cc,:)/max(plot_add_stat(cc,:)),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    title(conds{cc,1},'FontSize',14,'Interpreter','none')
    ylabel('Sum of the statistics for the whole image (scaled)','FontSize',12)
    hold on
    
    % Bar plot
    %---------
    if bar_plot
        figs(9) = figure('Name','Sum','Visible',visible);
        set(figs(9),'UserData','Sum of the statistics for the whole image')
        plot_add_stat_bar.val(:,bar_col) = plot_add_stat(cc,:);
        plot_add_stat_bar.color{bar_col} = conds{cc,3};
    end
    
    
    curr_cond{1} = conds{cc,1};
    curr_cond{2} = cc;
    if bar_plot
        bar_col = bar_col + 1;
    end
    
    
% Same condition
%==========================================================================
else
    
    % 1) Number of voxels above threshold
    %----------------------------------------------------------------------
    figure(figs(1)), set(figs(1),'Visible',visible)
    plot(delay,plot_voxels(cc,:),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    hold on
    
    figure(figs(2)), set(figs(2),'Visible',visible)
    plot(delay,plot_voxels(cc,:)/max(plot_voxels(cc,:)),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    hold on
    
    % Bar plot
    %---------
    if bar_plot
        %figure(figs(7)), set(figs(7),'Visible',visible)
        plot_voxels_bar.val(:,bar_col) = plot_voxels(cc,:);
        plot_voxels_bar.color{bar_col} = conds{cc,3};
    end
    
    
    % 2) Maximum value for the statistical test
    %----------------------------------------------------------------------
    figure(figs(3)), set(figs(3),'Visible',visible)
    plot(delay,plot_mn_mx(cc,:),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    hold on
    
    figure(figs(4)), set(figs(4),'Visible',visible)
    plot(delay,plot_mn_mx(cc,:)/max(plot_mn_mx(cc,:)),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    hold on
    
    % Bar plot
    %---------
    if bar_plot
        %figure(figs(8)), set(figs(8),'Visible',visible)
        plot_mn_mx_bar.val(:,bar_col) = plot_mn_mx(cc,:);
        plot_mn_mx_bar.color{bar_col} = conds{cc,3};
    end
    
    
    % 3) Sum of the statistics for the whole image
    %----------------------------------------------------------------------
    figure(figs(5)), set(figs(5),'Visible',visible)
    plot(delay,plot_add_stat(cc,:),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    hold on
    
    figure(figs(6)), set(figs(6),'Visible',visible)
    plot(delay,plot_add_stat(cc,:)/max(plot_add_stat(cc,:)),conds{cc,3},...
        'LineWidth',2,'Marker','o','MarkerFaceColor',conds{cc,3})
    hold on
    
    % Bar plot
    %---------
    if bar_plot
        %figure(figs(9)), set(figs(9),'Visible',visible)
        plot_add_stat_bar.val(:,bar_col) = plot_add_stat(cc,:);
        plot_add_stat_bar.color{bar_col} = conds{cc,3};
    end
    
    
    if bar_plot
        bar_col = bar_col + 1;
    end
    
end


% Condition will change
%==========================================================================
if (cc<num_cond && ~strcmp(conds{cc+1,1},curr_cond{1})) || cc==num_cond
    
    if bar_plot
        
        % Bar plot: Remove unused columns
        %------------------------------------------------------------------
        if bar_col<=3
            plot_voxels_bar.val(:,bar_col:end)   = [];
            plot_mn_mx_bar.val(:,bar_col:end)    = [];
            plot_add_stat_bar.val(:,bar_col:end) = [];
            plot_voxels_bar.color(bar_col:end)   = [];
            plot_mn_mx_bar.color(bar_col:end)    = [];
            plot_add_stat_bar.color(bar_col:end) = [];
        end
        bar_col = bar_col - 1;
        
        
        % Plot bars
        %------------------------------------------------------------------
        
        
        % 1) Number of voxels above threshold
        %------------------------------------
        figure(figs(7)), set(figs(7),'Visible',visible)
        if length(delay)>1
            bar_handle = bar(delay,plot_voxels_bar.val,'hist');
        else
            bar_handle = bar([delay;delay+1],...
                [plot_voxels_bar.val;plot_voxels_bar.val],'hist');
            xlim([delay-0.5 delay+0.5])
        end
        for bb=1:length(bar_handle) % change bar color
            set(bar_handle(bb),'FaceColor',plot_voxels_bar.color{bb})
        end
        
        % 2) Maximum value for the statistical test
        %------------------------------------------
        figure(figs(8)), set(figs(8),'Visible',visible)
        if length(delay)>1
            bar_handle = bar(delay,plot_mn_mx_bar.val,'hist');
        else
            bar_handle = bar([delay;delay+1],...
                [plot_mn_mx_bar.val;plot_mn_mx_bar.val],'hist');
            xlim([delay-0.5 delay+0.5])
        end
        for bb=1:length(bar_handle) % change bar color
            set(bar_handle(bb),'FaceColor',plot_mn_mx_bar.color{bb})
        end
        
        % 3) Sum of the statistics for the whole image
        %---------------------------------------------
        figure(figs(9)), set(figs(9),'Visible',visible)
        if length(delay)>1
            bar_handle = bar(delay,plot_add_stat_bar.val,'hist');
        else
            bar_handle = bar([delay;delay+1],...
                [plot_add_stat_bar.val;plot_add_stat_bar.val],'hist');
            xlim([delay-0.5 delay+0.5])
        end
        for bb=1:length(bar_handle) % change bar color
            set(bar_handle(bb),'FaceColor',plot_add_stat_bar.color{bb})
        end
        
    end
    
    
    % Edit plot, write data to file, save and close
    %----------------------------------------------------------------------
    for ff=1:num_figs
        
        % Edit plot
        %----------
        figure(figs(ff)), set(figs(ff),'Visible',visible)
        if ~group
            xlabel('Delay (s)','FontSize',12)
        else
            xlabel('Group','FontSize',12)
            set(gca,'XTick',delay)
        end
        grid on
        %axis tight
        if ~line_plot(ff) && bar_plot
            title(conds{cc,1},'FontSize',14,'Interpreter','none')
            ylabel(get(figs(ff),'UserData'),'FontSize',12)
        end
        legend(conds(curr_cond{2}:cc,2))
        set(gca,'FontSize',12) % increase size of numbers and legend font
        
        % Write data to file
        %-------------------
        if line_plot(ff) && isempty(strfind(get(figs(ff),'Name'),'_scaled'))
        % save only when the figure is a line plot and the values are not scaled
            for ii=1:3
                switch ii
                    case 1
                        stat = 'F';
                        tmp = findobj(figs(ff),'DisplayName','F test');
                    case 2
                        stat = 'Tneg';
                        tmp = findobj(figs(ff),'DisplayName','T (negative BOLD)');
                    case 3
                        stat = 'Tpos';
                        tmp = findobj(figs(ff),'DisplayName','T (positive BOLD)');
                end
                if isempty(tmp), continue, end
                if ~movie
                    fprintf(fid,'\n%s\t%s\t%s',conds{cc,1},fig_names{ff},stat);
                    fprintf(fid,'\t%.4f',get(tmp,'YData'));
                else
                    fprintf(fid,'\t%s-%s',fig_names{ff},stat); % header
                    movie_data = [movie_data get(tmp,'YData')'];
                end
            end
        end
        
        % Save and close figure
        %----------------------
        if line_plot(ff)
            if bar_plot
                tmp = ['Line_plots' filesep];
            else
                tmp = '';
            end
            fname = [tmp conds{cc,1} '_' fig_names{ff}];
        elseif bar_plot
            fname = ['Bar_plots' filesep conds{cc,1} '_' fig_names{ff}];
        end
        drawnow; pause(.1)
        saveas(figs(ff), [fname '.png'])
        %fig_print(figs(ff), [fname '.png'], [1200 900])
        % Same output, except for the grey background. White is better for
        % printing
        saveas(figs(ff), [fname '.fig'])
        close(figs(ff))
    end
    
end

end % loop for the conditions


% Finish
%==========================================================================
fprintf('Graphics were saved in %s\n',fullfile(path_res,'Plots'))
if ~group
    fprintf('Data summary was saved in %s\n',fullfile(path_res,'Event_seq-Summary.txt'))
else
    fprintf('Data summary was saved in %s\n',fullfile(path_res,'Group-Summary.txt'))
end
if movie
    tmp = repmat('%.4f\t',1,size(movie_data,2)-1);
    tmp(end-1:end) = [];
    fprintf(fid,['\n%.1f\t' tmp],movie_data');
end
fclose(fid);
pause(0.5) % give some time for the figures to be saved
cd(curr_dir)

tmp = clock;
fprintf('(%.2d:%.2d:%02.0f) %s.m done!\n',...
    tmp(4),tmp(5),tmp(6),mfilename)
disp('-------------------------------------------------------------------')