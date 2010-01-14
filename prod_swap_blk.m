function y = prod_swap_blk(f, A, x, n)
% PROD_SWAP_BLK inside out block product
% y = prod_swap_blk(f, A, x, n)
%
% Apply block product to A, x inside out.
%
% See also prod_blk, swap_blk, swap_vec.
%

N = length(x);
y = swap_vec(prod_blk(f, swap_blk(A, n), swap_vec(x, n), n), N/n);
