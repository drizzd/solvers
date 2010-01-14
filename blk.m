function y = blk(x, d, a, b, step)
% BLK extract block
% y = blk(x, d, a = 1, b = 1, step = 1)
%

if nargin < 5
	step = 1;
end
if nargin < 4
	b = 1;
end

n = size(x, 1);
if mod(n, d) ~= 0
	error('size of x is not a multiple of d')
end

if step ~= -1 && step ~= 1
	error('invalid step size')
end

if step == -1
	a = a - b + 1;
end

y = x((a-1)*d+1:(a+b-1)*d, 1:min(size(x, 2), d));

if step == -1
	yr = zeros(size(y));
	S = size(y, 1);
	if mod(S, d) ~= 0
		error('this should not happen');
	end
	S = S/d;
	for k = 1:S
		yr((k-1)*d+1:k*d, :) = y((S-k)*d+1:(S-k+1)*d, :);
	end
	y = yr;
end
