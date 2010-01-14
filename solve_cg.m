function [x, r] = solve_cg(A, b, x, nr_iter, prod_matvec)
% SOLVE_CG conjugate gradient solver
% [x, r] = solve_cg(A, b, x, nr_iter, prod_matvec)
%
% Solves `A x = b' using the conjugate gradient method.
%
% `nr_iter'	number of iterations (default: length(x))
% `prod_matvec'	optionally supply custom matrix vector multiplication (default:
%		@multiply builtin)
%

if nargin < 3
	x = zeros(size(b));
end
if nargin < 4
	nr_iter = size(A, 1);
end
if nargin < 5
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

rr = zeros(M, 1);
pAp = zeros(M, 1);

r = b - prod_matvec(A, x);
for k = 1:nr_iter
	rr = mnorm(r).^2;
	if all(rr == 0)
		break;
	end
	if k == 1
		p = r;
	else
		beta = rr./rr_prev;
		beta(rr == 0) = 0;
		p = r + p * diag(beta);
	end
	Ap = prod_matvec(A, p);
	pAp = real(mdot(p, Ap));
	alpha = rr./pAp;
	alpha(rr == 0) = 0;
	x = x + p * diag(alpha);
	r = r - Ap * diag(alpha);
	rr_prev = rr;
end
