function A = randtoep(n, d)
% RANDTOEP random toeplitz matrix
% A = randtoep(n, d)

if nargin < 2
	d = 1;
end

H = [];
for j = 1:d
	H = [H; convmtx(complex(rand(1, n), rand(1, n)), n)];
	%H = [H; toeplitz([1; rand(n-1, 1)], [1; rand(n-1, 1)])];
end
A = H * H';
