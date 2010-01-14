function [x, r] = solve_pcg(C, A, b, x, nr_iter, prod_matvec)
% SOLVE_PCG preconditioned conjugate gradient solver
% [x, r] = solve_pcg(C, A, b, x, nr_iter, prod_matvec)
%
% Solves `C^-1 A x = C^-1 b' using the preconditioned conjugate gradient
% method.
%
% `nr_iter'	number of iterations (default: length(x))
% `prod_matvec'	optionally supply custom matrix vector multiplication (default:
%		@multiply builtin)
%

if nargin < 4
	x = zeros(size(b));
end
if nargin < 5
	nr_iter = size(A, 1);
end
if nargin < 6
	prod_matvec = @multiply;
end

n = size(A, 1);
if size(A, 2) ~= n
	error('A must be square')
end
if size(b, 1) ~= n || size(x, 1) ~= n
	error('dimensions of A, b or x do not match')
end
M = size(b, 2);

if nr_iter > n
	nr_iter = n;
end

rz = zeros(M, 1);
pAp = zeros(M, 1);

r = b - prod_matvec(A, x);
for k = 1:nr_iter
	z = C\r;
	rz = real(mdot(r, z));
	if all(rz == 0)
		break;
	end
	if k == 1
		p = z;
	else
		beta = rz./rz_prev;
		beta(rz == 0) = 0;
		p = z + p * diag(beta);
	end
	Ap = prod_matvec(A, p);
	pAp = real(mdot(p, Ap));
	alpha = rz./pAp;
	alpha(rz == 0) = 0;
	x = x + p * diag(alpha);
	r = r - Ap * diag(alpha);
	rz_prev = rz;
end
