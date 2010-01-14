function x = solve_sd(A, b, x, nr_iter, prod_matvec)
% SOLVE_SD steepest descent solver
% x = solve_sd(A, b, x, nr_iter, prod_matvec)
%
% Solves `A x = b' using the steepest descent method.
%
% `nr_iter'	number of iterations (default: length(x))
% `prod_matvec'	optionally supply custom matrix vector multiplication (default:
%		builtin)
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
if length(b) ~= n || length(x) ~= n
	error('dimensions of A, b or x do not match')
end

if nr_iter > n
	nr_iter = n;
end

r = b - prod_matvec(A, x);
for k = 1:nr_iter
	rr = norm(r)^2;
	if rr == 0
		break;
	end
	Ar = prod_matvec(A, r);
	alpha = rr/real(r'*Ar);
	x = x + alpha * r;
	r = r - alpha * Ar;
end
