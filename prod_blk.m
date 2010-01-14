function y = prod_blk(f, A, x, n)
% PROD_BLK block product
% y = prod_blk(f, A, x, n)
%
% Compute A x by applying f to n x n subblocks of A and x. E.g., f = @multiply
% will always work. However, if subblocks of A have structure, there may be
% more efficient methods of multiplication, such as FFT based convolution for
% block toeplitz matrices.
%
% See also multiply.
%

L = length(x);
if mod(L, n) ~= 0
	error('length of x is not a multiple of n');
end
m = L/n;
if size(A) ~= [L, L]
	error('matrix-vector dimensions do not match');
end

y = zeros(size(x));
for k = 1:m:m*n
	for l = 1:m:m*n
		y(k:k+m-1) = y(k:k+m-1) + f(A(k:k+m-1, l:l+m-1), x(l:l+m-1));
	end
end
