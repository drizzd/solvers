function xx = mnorm(x, t)
% MNORM multi-column vector norm
% a = mnorm(x, t)
%
% a(k) = norm(x(:, k), t)
%

if nargin < 2
	t = 2;
end

S = size(x);
xx = zeros([S(2:end) 1]);
n = numel(xx);
for k = 1:n
	xx(k) = norm(x(:, k), t);
end
