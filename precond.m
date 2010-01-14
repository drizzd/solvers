function c = precond(A, type)
% PRECOND first column of preconditioner
% c = precond(A, type)
%
% 'type' is on of { 'strang', 'chan' }
%
% See also solve_pcg.
%

N = size(A, 1);
if size(A, 2) ~= N
	error('A must be square')
end

switch type
case 'strang'
	n = ceil(N/2);
	m = 1-mod(N, 2);
	x = (A(n+1:n+m, 1) + A(1, n+1:n+m).')/2;
	c = [A(1:n, 1); x; A(1, n:-1:2).'];
case 'chan'
	a = [0:1/N:1].';
	c = a(end:-1:2) .* A(:, 1) + a(1:end-1) .* [0; A(1, end:-1:2).'];
otherwise
	error(['unknown preconditioner: ' repr(type)]);
end
