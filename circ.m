function C = circ(c)
% CIRC Circular matrix.
% C = circ(c)
%
% c is the first column of the circular matrix C.
%

C = toeplitz(c, [c(1); c(end:-1:2)]);
