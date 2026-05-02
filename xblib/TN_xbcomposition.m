function Out = TN_xbcomposition(G)
N = length(G);
m = 2; n = 1;
Out = G{1};
M = N;
M1=M-1;
N1=N-1;
for i = 1:N-2
    Out = tensor_xbcontraction(Out,G{i+1},M1,N1,m,n);
    M1 = M1+N1-2*i;
    n = [n,1+i];
    tempm = 2+i*(N1-i);
    if i>1
        m(2:end)=m(2:end)-[1:i-1];
    end
    m   = [m, tempm];
end
end