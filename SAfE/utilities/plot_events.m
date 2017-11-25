function plot_events(fig, cond_time, color)
% 
% 
%USAGE
%-----
% plot_events(fig, cond_time)
% plot_events(fig, cond_time, color)
% 
% 
%INPUT
%-----
% - FIG      : figure handle
% - COND_TIME: Cx1 cell array, where C is the number of conditions.
%   COND_TIME{cc} is an EVx2 matrix with the onsets in column 1 and the end
%   of the stimulus in column 2 (EV is the number of events)
%   N.B.: The units must be the same as the x-axis in FIG
% - COLOR    : Cx3 double array with the color in RGB for condition C
%   - [1 1 1]/3 is the default gray
% 
% 
%OUTPUT
%------
% - FIG will have rectangles indicating COND_TIME (useful for fMRI to show
%   the conditions) with the corresponding colors
% 
%__________________________________________________________________________
% Copyright (C) 2013-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2012-Oct-18, 05:21 pm
% 
% Changelog
%--------------------------------------------------------------------------
% 2012-Nov-21, 04:20 pm: corrected getting y-axis limits when there is no
%   legend
% 2013-Feb-23, 01:07 pm: different conditions are allowed and their color
%   can be different


if isempty(cond_time)
    return
end
if ~iscell(cond_time)
    cond_time = {cond_time};
end
if nargin<3
    color = repmat([1 1 1]/3, length(cond_time), 1);
end


% Get the Y-axis limits
%----------------------
figure(fig)
tmp   = gca(fig);
ylims = get(tmp,'YLim');


% Draw the rectangles
%--------------------
for cc=1:length(cond_time)
    tmp = [repmat(ylims(1),2,size(cond_time{cc},1)) ; repmat(ylims(2),2,size(cond_time{cc},1))];
    patch([cond_time{cc}(:,1)' ; cond_time{cc}(:,2)' ; cond_time{cc}(:,2)' ; cond_time{cc}(:,1)'],...
        tmp,color(cc,:),'FaceAlpha',0.5,'EdgeAlpha',0)
end
drawnow; pause(.1)