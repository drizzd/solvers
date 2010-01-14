function solve = get_solver(type, Q, Nd)
% GET_SOLVER
% solve = get_solver(type, Q, Nd)
%

swap_solve = @(solve, x, y)swap_vec(solve(swap_blk(x, Q), swap_vec(y, Q), Nd), Nd);
switch type
case 'pseudo'
	solve = @(x, y)(pinv(x)*y);
case 'matlab'
	warning('off', 'MATLAB:singularMatrix');
	solve = @(x, y)(x\y);
case 'lev'
	% Levinson solver
	solve = @(x, y)swap_solve(@solve_lev_block, x, y);
case 'schur'
	% Schur algorithm
	solve = @(x, y)swap_solve(@solve_schur, x, y);
case {'cg', 'sd', 'pcg-strang', 'pcg-chan'}
	switch type
	case 'cg'
		solve = @solve_cg;
	case 'sd'
		solve = @solve_sd;
	case {'pcg-strang', 'pcg-chan'}
		prec = @(A)blkfun(@(B)circ(precond(B, type(5:end))), A, Nd);
		solve = @(A, b, x, nr_iter, prod_matvec)(...
			solve_pcg(prec(A), A, b, x, nr_iter, prod_matvec));
	end
otherwise
	error(['unknown solver ' repr(type)]);
end
