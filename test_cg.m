function A = test_cg(n, iter, reg)
% TEST_CG test conjugate gradient algorithms
% A = test_cg(n, iter, reg)
%

if nargin < 3
	reg = 0;
end
if nargin < 2
	iter = n;
end

%A = randtoep(n);
A = randcorr(n);
% regularization
A = A + reg * eye(n);
b = complex(rand(n, 1), rand(n, 1));
clf;
hold on;
plot(abs(fft(A(:,1))), 'g');
cg = norm(b - A * solve_cg(A, b, zeros(n, 1), iter));
fprintf('                k       2       1       inf     F\n');
fprintf('CG:    %-8.2g %-7.2g\n', [cg, cond(A)]);
cg = norm(b - A * pcg(A, b, 0, iter));
fprintf('pcg:   %-8.2g\n', cg);
cg = norm(b - A * bicg(A, b, 0, iter));
fprintf('bicg:   %-8.2g\n', cg);
C = circ(precond(A, 'strang'));
plot(abs(fft(C(:,1))), 'r');
E = C-A;
strang = norm(b - A * solve_pcg(C, A, b, zeros(n, 1), iter));
fprintf('1 PCG: %-8.2g %-7.2g %-7.2g %-7.2g %-7.2g %-7.2g\n', ...
	[strang, cond(C), norm(E)/n, norm(E, 1)/n, ...
	norm(E, inf)/n, norm(E, 'fro')/n]);
C = circ(precond(A, 'chan'));
plot(abs(fft(C(:,1))), 'b');
E = C-A;
chan = norm(b - A * solve_pcg(C, A, b, zeros(n, 1), iter));
fprintf('F PCG: %-8.2g %-7.2g %-7.2g %-7.2g %-7.2g %-7.2g\n', ...
	[chan, cond(C), norm(E)/n, norm(E, 1)/n, ...
	norm(E, inf)/n, norm(E, 'fro')/n]);
chan = norm(b - A * pcg(A, b, 0, iter, C));
fprintf('F pcg: %-8.2g\n', chan);
