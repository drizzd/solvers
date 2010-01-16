function y = blk_row(x, d, a, m, step)
% BLK_ROW extract block rows
% y = blk_row(x, d, a = 1, m = 1, step = 1)
%

if nargin < 5
	step = 1;
end
if nargin < 4
	m = 1;
end

n = size(x, 1);
if mod(n, d) ~= 0
	error('size of x is not a multiple of d')
end

m = m * abs(step);

if step < 0
	a = a - m + 1;
	u = m;
	v = 1;
else
	u = 1;
	v = m;
end

idx = reshape((a-1)*d+1:(a+m-1)*d, d, m);
idx = idx(:, u:step:v);

y = x(idx, :);
