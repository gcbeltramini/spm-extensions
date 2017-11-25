function data_mean = block_average(data, start_end, shift)
% 
% Block-wise data average
% 
% 
%USAGE
%-----
% data_mean = block_average(data, start_end)
% data_mean = block_average(data, start_end, shift)
% 
% 
%INPUT
%-----
% - DATA     : NxM matrix with M signals along the rows
% - START_END: Nblocks x 2 matrix with the first and last elements to take
%   the average of the signal
% - SHIFT    : maximum number of data points to shift the elements in
%   START_END to find the best match (greatest correlation coefficient).
%   Only the first column in DATA is optimized. The same shift will be
%   applied to all other columns. Default: 0
% 
% 
%OUTPUT
%------
% - DATA_MEAN: average of DATA for the elements in START_END
% 
% 
%EXAMPLE
%-------
% >> data_mean = block_average(data,[5 35;35 65;65 95],3);
% 
%__________________________________________________________________________
% Copyright (C) 2013-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2013-Feb-25, 12:10 pm


% Input
%==========================================================================

if isvector(data)
    data = data(:);
end

if nargin<3
    shift = 0;
elseif ceil(shift)~=floor(shift)
    error('"SHIFT" must be integer')
else % valid input for "shift"
    shift = abs(shift);
end

if any(start_end(:,1)>start_end(:,2))
    error('The numbers in the 2nd column in "start_end" must be greater than the ones in the 1st column')
end

if ~isequal(floor(start_end), ceil(start_end))
    error('The indices must be positive integers')
end

if any(start_end(:,1)<1)
    error('Invalid element. The smallest number must be 1')
    %start_end(start_end(:,1)<1,1) = 1;
end

% Number of data points must be the same
if any(diff(start_end(:,2)-start_end(:,1)))
    error('The number of data points is not the same for all intervals')
end


% Initialize
%==========================================================================
N    = size(data, 1);      % number of points for each channel
Nint = size(start_end, 1); % number of intervals
if any(start_end(:,2)>N) % the greatest number must be N
    %error('Invalid element. The greatest number must be %d', N)
    [tmp, idx] = sort(start_end, 1);
    tmp(:, 2) = tmp(:, 1) + ( N - tmp(Nint, 1) ); % adjust the maximum element
    
    % Undo sort
    start_end = zeros(size(tmp));
    start_end(idx(:, 1), 1) = tmp(:, 1);
    start_end(idx(:, 1), 2) = tmp(:, 2);
    %for ii = 1:2 % loop for the columns in "start_end"
    %    start_end(idx(:, ii), ii) = tmp(:, ii);
    %end
    
end


% Find the best intervals
%==========================================================================
Nval      = start_end(1,2) - start_end(1,1) + 1; % number of values for each channel
data_mean = zeros(Nval, size(data,2), Nint);     % initialize block averaged data
% dimension 1: values for each channel
% dimension 2: channels
% dimension 3: intervals/blocks
data_mean(:,:,1) = data(start_end(1,1):start_end(1,2), :); % first interval for all channels

ch = 1; % use only the first channel to optimize the shift
corr_tmp = data_mean(:,ch,1) - sum(data_mean(:,ch,1), 1)/Nval; % zero-mean
corr_tmp = corr_tmp/sqrt(sum(abs(corr_tmp).^2, 1)); % normalize (abs guarantees a real result)

for ii=2:Nint % loop for the intervals
    
    % Look for a better interval by shifting the given one back and forth
    best = -Inf;
    for ss=-shift:shift % loop for the shifts
        
        block = (start_end(ii,1):start_end(ii,2)) + ss;
        if any(block<1) || any(block>N) % indices out of bounds
            continue
        end
        
        if shift>0
            % Correlation coefficient
            tmp = data(block, ch) - sum(data(block, ch), 1)/Nval; % zero-mean
            tmp = tmp/sqrt(sum(abs(tmp).^2, 1));                  % normalize
            tmp = sum(corr_tmp.*tmp, 1);                          % correlation
            % Equivalent method using Matlab's implementation:
            %tmp = corrcoef(data_mean(:,ch,1), data(block, ch)); tmp = tmp(2);
            if tmp>best % the correlation is greater => use this block
                best = tmp;
                data_mean(:,:,ii) = data(block, :);
            end
        else
            data_mean(:,:,ii) = data(block, :);
        end
        
    end
end


% Calculate the average across the intervals
%==========================================================================
data_mean = sum(data_mean, 3)/Nint;