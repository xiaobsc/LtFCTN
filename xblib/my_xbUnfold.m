function W_out = my_xbUnfold(W, dim, i)
n = length(dim);
perm_order = [i, 1:i-1, i+1:n-1, n];
W = permute(W, perm_order);
mid_dim = prod(dim([1:i-1, i+1:n-1]));
W_out = reshape(W, dim(i), mid_dim, dim(n));
end