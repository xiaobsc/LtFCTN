function I = xbteye(rho,n,n3)
I = zeros(n,n,n3);
I(:,:,1) = rho*eye(n);
end
