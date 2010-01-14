function s = mdot(x, y)
% MDOT multi-column vector dot product
% s = mdot(x, y)
%
% s(k) = x(:, k)' * y(:, k)
%

s = shiftdim(sum(conj(x) .* y, 1), 1);
