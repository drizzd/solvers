function s = repr(x)
% REPR string representation
% s = repr(x)
%

if ischar(x(1))
	s = ['''' x(:)' ''''];
else
	s = mat2str(x);
end
