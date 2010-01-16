function x = solve_lev_block(A, y, d)
% SOLVE_LEV_BLOCK block Levinson solver
% x = solve_lev_block(A, b, d)
%
% Solves `A x = b' using Block Levinson algorithm if A is a d x d block
% Toeplitz matrix.
%
% `d'	block matrix size (default 1)
%

if nargin < 3
	d = 1;
end

n = size(A, 1);
if size(A, 2) ~= n
	error('A must be square')
end
if size(y, 1) ~= n
	error('dimensions of A and y do not match')
end
if mod(n, d) ~= 0
	error('dimensions of A must be a multiple of the block size')
end

if n == d
	x = A\y;
	return;
end

b = blk(A, d, 1);
c = b;
a = -b\blk(A.', d, 2).';
v = -c\blk(A, d, 2);
x = b\blk_row(y, d, 1);
Y = a;
W = v;
N = n/d;
for k = 1:N-1
	b = b - b * a * v;
	c = c - c * v * a;
	u = blk_row(y, d, k+1);
	for l = 2:k+1
		u = u - blk(A, d, l) * blk_row(x, d, k - l + 2);
	end
	mu = c\u;
	x = [x + blk(Y, d, k, k, -1) * mu; mu];
	if k < N - 1
		Y_ = Y;
		u = blk(A.', d, k+2).';
		for l = 2:k+1
			u = u + blk(A.', d, l).' * blk(Y, d, k - l + 2);
		end
		a = -b\u;
		Y = [Y + blk(W, d, k, k, -1) * a; a];
		u = blk(A, d, k+2);
		for l = 2:k+1
			u = u + blk(A, d, l) * blk(W, d, k - l + 2);
		end
		v = -c\u;
		W = [W + blk(Y_, d, k, k, -1) * v; v];
	end
end
