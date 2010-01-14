function B = blkfun(f, A, n)
% BLKFUN block function
% B = blkfun(f, A, n)
%

if size(A, 1) ~= size(A, 2)
	error('A is not square');
end
if mod(size(A, 1), n) ~= 0
	error('size of A is not a multiple of n');
end
N = size(A, 1);
M = N/n;

F = size(f(A(1:M, 1:M)));
B = zeros(n * F);
for k = 1:n
	for l = 1:n
		B(1+(k-1)*F(1):k*F(1), 1+(l-1)*F(2):l*F(2)) = ...
				f(A(1+(k-1)*M:k*M, 1+(l-1)*M:l*M));
	end
end
