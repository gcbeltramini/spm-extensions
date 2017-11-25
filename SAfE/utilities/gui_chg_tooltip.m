function gui_chg_tooltip(obj1_handle, obj2_handle, tooltip_string)
% 
% Change the tooltip string
% 
% 
%USAGE
%-----
% chg_tooltip(obj1_handle, obj2_handle, tooltip_string)
% 
% 
%INPUT
%-----
% - OBJ1_HANDLE   : handle of the edit box object
% - OBJ2_HANDLE   : handle of the text object describing OBJ1_HANDLE (can
%   be empty)
% - TOOLTIP_STRING: tooltip string (can be empty)
% 
% 
%OUTPUT
%------
% Tooltip string of:
% - OBJ1_HANDLE: TOOLTIP_STRING if the text in OBJ1_HANDLE fits the box;
%   the text, otherwise
% - OBJ2_HANDLE: empty if the text in OBJ1_HANDLE fits the box;
%   TOOLTIP_STRING, otherwise
% 
%__________________________________________________________________________
% Copyright (C) 2013-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2013-May-15, 01:47 am


oldunits = get(obj1_handle, 'Units');    % get the original units
set(obj1_handle, 'Units', 'characters')
width    = get(obj1_handle, 'Position'); % get the width in character units
width    = width(3);
set(obj1_handle, 'Units', oldunits)      % reset units

tmp = get(obj1_handle, 'String');
tmp = regexprep(tmp, '\s+', ' '); % remove repeated spaces

if length(tmp)<=width % content is smaller than the box width
    set(obj1_handle, 'TooltipString', tooltip_string)
    set(obj2_handle, 'TooltipString', '')
else
    if length(tmp)>55 % too big for the tooltip
        tmp = [tmp(1:50) '...'];
    end
    set(obj1_handle, 'TooltipString', tmp)
    set(obj2_handle, 'TooltipString', tooltip_string)
end
