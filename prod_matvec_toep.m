function y = prod_matvec_toep(A, x)
% PROD_MATVEC_TOEP toeplitz matrix vector multiplication
% y = prod_matvec_toep(A, x)
%

if size(A) ~= [1 1] * length(x)
	error('matrix-vector dimensions do not match');
end

y = circconv([A(:, 1); A(1, end:-1:2).'], [x; zeros(length(x) - 1, 1)]);
y = y(1:length(x));
