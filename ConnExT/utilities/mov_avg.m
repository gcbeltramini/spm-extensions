function yavg = mov_avg(y, window)
% 
% Calculate the simple moving average for each column. Note that this only
% works well when the points in the independent variable are uniformly
% sampled.
% 
%USAGE
%-----
% mov_avg(y)
% mov_avg(y, window)
% 
% 
%INPUT
%-----
% - Y     : NxS data matrix of S series of N points
% - WINDOW: window size (default: 3)
%   Note that WINDOW < N
% 
% 
%OUTPUT
%------
% - YAVG: NxS data matrix of S series averaged with the moving average
%         filter
% 
% 
% Because of the endpoints, this function is better than:
% >> conv(y, ones(1,window)/window, 'same');
%    (this is better than "filter", but has problem with the endpoints)
% >> filter(ones(1,window)/window, 1, y);
%    (the averaged curve is shifted to the right)
% 
% 
% Based on Matlab's "smooth" function from the Curve Fitting Toolbox

% Guilherme Coco Beltramini - 2014-Jun-01, 06:36 pm


% Input data matrix
%------------------
if isvector(y)
    y = y(:);
end
N = size(y,1);


% Window size
%------------
if nargin<2 || window<1
    window = 3;
elseif window==1
    yavg = y;
    return
else
    window = min(floor(window), N); % the window size cannot be greater than N
    window = window - 1 + mod(window, 2); % force it to be odd
end


% Moving average
%---------------
yavg = filter(ones(window,1)/window, 1, y);


% Adjust the endpoints
%---------------------

% When "y" is a vector:
%avgi = cumsum(y(1:window-2));
%avgi = avgi(1:2:end) ./ (1:2:(window-2))';
%avgf = cumsum(y(N:-1:N-window+3));
%avgf = avgf(end:-2:1) ./ (window-2:-2:1)';
%yavg = [avgi; yavg(window:end); avgf];

% When "y" is a matrix:
avgi = cumsum(y(1:window-2, :), 1);
avgi = bsxfun(@rdivide, avgi(1:2:end, :), (1:2:(window-2))');
avgf = cumsum(y(N:-1:N-window+3, :), 1);
avgf = bsxfun(@rdivide, avgf(end:-2:1, :), (window-2:-2:1)');
yavg = [avgi; yavg(window:end, :); avgf];