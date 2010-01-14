function [x, R] = solve_schur(A, b, d)
% SOLVE_SCHUR block Schur solver
% x = solve_schur(A, b, d)
%
% Solves `A x = b' using Schur algorithm. A must be Hermitian d x d block
% Toeplitz matrix. Only first block column of A is used.
%
% `d'	block matrix size (default 1)
%

if nargin < 3
	d = 1;
end

n = size(A, 1);
if size(A, 2) < d
	error('A must contain at least one block column');
end
if size(b, 1) ~= n
	error('dimensions of A and b do not match')
end
if mod(n, d) ~= 0
	error('dimension of A must be a multiple of the block size')
end

N = n/d;

% QR factorization
T = A(:, 1:d)';
L = chol(T(1:d, 1:d), 'lower');
T(:, 1:d) = L';
T(:, 1+d:N*d) = L^-1 * T(:, 1+d:N*d);

R = zeros(n, n);
R(1:d, :) = T;

T = repmat(T, 2, 1);

T(1+d:end, 1:d) = zeros(d);
W = diag([ones(1, d) -ones(1, d)]);
Z = zeros(n);
Z(1:end-d, 1+d:end) = eye(n - d);
assert(norm(A - Z'*A*Z - T'*W*T) < 1e-9);

for j = 1:N-1
	T = [T(1:d, 1:end-d); T(d+1:2*d, 1+d:end)];
	for k = 1:d
		x = T(k:end, k);
		W = diag([ones(d-k+1, 1); -ones(d,1)]);
		v = W * x;
		alpha = sqrt(x' * W * x);
		v(1) = v(1) + x(1)/abs(x(1)) * alpha;
		U = eye(2*d);
		beta = v'*W*v;
		U(k:end, k:end) = W - 2/beta * v * v';
		T = U * T;
	end
	R(j*d+1:(j+1)*d, j*d+1:end) = T(1:d, :);
end

% back substitution
y = R'\b;
x = R\y;
