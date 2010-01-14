function A = randcorr(n)
% RANDCORR random correlation matrix
% A = randcorr(n)
%

if nargin < 2
	reg = 0;
end

nb = ceil(n/2);
ev = [zeros(n-nb, 1); ones(nb, 1)*n/nb];
H = complex(gallery('randcorr', ev), gallery('randcorr', ev));
A = H * H';
