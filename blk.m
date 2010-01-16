function y = blk(varargin)
% BLK_ROW extract block
% y = blk(x, d, a = 1, b = 1, step = 1)
%

y = blk_row(varargin{:});
d = varargin{2};
y = y(:, 1:d);
