function varargout = trim_array(action,varargin)
% 
% Remove columns and rows that are empty or that contain only NaN's in the
% input arrays, in the beginning, in the end, in both extremities or in the
% entire array
% 
% 
%USAGE
%-----
% [out1,out2,...] = trim_array(action,array1,array2,...)
% output          = trim_array(action,array1,array2,...)
% 
% 
%INPUT
%-----
% - ACTION: 'beginning', 'end', 'edges', or 'all'
%   - 'beginning': remove only the first rows and columns
%   - 'end'      : remove only the last rows and columns
%   - 'edges'    : remove the both the first and last rows and columns
%   - 'all'      : remove all rows and columns with NaN or empty elements
% - ARRAY : cell or numeric array
% 
% 
%OUTPUT
%------
% - OUTPUT: If N input arrays are given, with N outputs, each output is the
%   corresponding trimmed array. If one output is given, it is a cell array
%   containing the trimmed arrays in each element.
% 
%__________________________________________________________________________
% Copyright (C) 2013-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2013-Feb-20, 11:48 am

num_out = nargout;
num_in  = nargin - 1;

if num_out~=num_in && num_out~=1 && num_out~=0
    error('Invalid number of outputs')
end

output = varargin;
if isempty(output)
    error('Invalid number of inputs')
end

for vv=1:num_in % loop for all input variables
    
    
    % Find NaN's and empty elements
    %==============================
    tmp = output{vv}; % to preserve the input
    if iscell(tmp)
        test = cellfun(@isempty,tmp);
        tmp(test) = {0}; % put 0 in the empty cells
        tmp(cellfun(@ischar,tmp)) = {0}; % put 0 where the cell is a character array
        test = test + cell2mat(cellfun(@isnan,tmp,'UniformOutput',0));
    else
        test = isnan(tmp);
    end
    
    
    % Remove rows
    %============
    remove = sum(test,2)==size(output{vv},2);
    tmp = [];
    switch action
        case 'beginning'
            if remove(1) % there are empty rows in the beginning
                tmp = 1:find(diff(remove),1,'first');
            end
        case 'end'
            if remove(end) % there are empty rows in the end
                tmp = (find(diff(remove),1,'last')+1):size(output{vv},1);
            end
        case 'edges'
            if remove(1) % there are empty rows in the beginning
                tmp = 1:find(diff(remove),1,'first');
            end
            if remove(end) % there are empty rows in the end
                tmp = [tmp ...
                    (find(diff(remove),1,'last')+1):size(output{vv},1)];
            end
        case 'all'
            if any(remove)
                tmp = remove;
            end
        otherwise
            error('Invalid "action"')
    end
    output{vv}(tmp,:) = [];
    
    
    % Remove columns
    %===============
    remove = sum(test,1) == (size(output{vv},1)+nnz(tmp));
    tmp = [];
    switch action
        case 'beginning'
            if remove(1) % there are empty columns in the beginning
                tmp = 1:find(diff(remove),1,'first');
            end
        case 'end'
            if remove(end) % there are empty columns in the end
                tmp = (find(diff(remove),1,'last')+1):size(output{vv},2);
            end
        case 'edges'
            if remove(1) % there are empty columns in the beginning
                tmp = 1:find(diff(remove),1,'first');
            end
            if remove(end) % there are empty columns in the end
                tmp = [tmp ...
                    (find(diff(remove),1,'last')+1):size(output{vv},2)];
            end
        case 'all'
            if any(remove)
                tmp = remove;
            end
    end
    output{vv}(:,tmp) = [];
    
    
    % Output
    %=======
    if num_out>1 || num_in==1
        varargout{vv} = output{vv};
    end
    
end

% Output
%=======
if num_out<=1 && num_in~=1
    varargout{1} = output;
end