function Out = sub_xbTN(G,k)

N = length(G);
a = [k+1:N-1, 1:k,N];
for i = [1:k-1,k+1:N]
    G{i}      = permute(G{i},a);
end
m = 2; n = 1;
Out = G{a(1)};
M = N;
M1=M-1;
N1=N-1;
for i = 1:N1-2
    Out = tensor_xbcontraction(Out,G{a(i+1)},M1,N1,m,n);
    M1 = M1+N1-2*i;
    n = [n,1+i];
    tempm = 2+i*(N1-i);
    if i>1
        m(2:end)=m(2:end)-[1:i-1];
    end
    m   = [m, tempm];
end
szeOut=size(Out);
In=length(szeOut);
p = zeros(1,2*(N1-k));

for i =1:(N1-k)
    p(2*i-1)= 2*i;
    p(2*i)  = 2*i-1;
end
Out = permute(Out,[p,2*(N1-k)+1:2*(N1-1),In]);
Out = permute(Out,[2*(N1-k)+1:2*(N1-1),1:2*(N1-k),In]);
end