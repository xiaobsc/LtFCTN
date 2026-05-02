function Out = tensor_xbcontraction(X, Y, Sx, Sy, n, m)
Nx = ndims(X);
Ny = ndims(Y);
Lx = size(X);
Ly = size(Y);
tp=Lx(:,end);
if  Nx < Sx
    tempLx = ones(1,Sx-Nx);
    Lx = [Lx,tempLx];
end
if  Ny < Sy
    tempLy = ones(1,Sy-Ny);
    Ly = [Ly,tempLy];
end

indexx = 1:Sx;
indexy = 1:Sy;
indexx(n) = [];
indexy(m) = [];
tempX = permute(X,[indexx,n,Nx]); 
tempXX=reshape(tempX,[prod(Lx(indexx)),prod(Lx(n)),tp]);
tempY = permute(Y,[m,indexy,Ny]);
tempYY=reshape(tempY,[prod(Ly(m)),prod(Ly(indexy)),tp]);
tempOut = tprod(tempXX,tempYY);
Out     = reshape(tempOut,[Lx(indexx),Ly(indexy),tp]);
end
