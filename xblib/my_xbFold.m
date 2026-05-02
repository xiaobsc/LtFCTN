function W = my_xbFold(W, dim, i)
n = length(dim);
m= [i, 1:i-1, i+1:n-1, n];
W = reshape(W, dim(m));
W = ipermute(W,m);
end