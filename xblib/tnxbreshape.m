function Out = tnxbreshape(Grest,N,i)
szeGrest=size(Grest);
tp=szeGrest(:,end);
In=length(szeGrest);
Nway = size(Grest);
if  length(Nway) < 2*(N-1)
    Nway = [Nway,ones(1,2*(N-1)-length(Nway))];
end

m = zeros(N-2,1);   n = zeros(N-2,1);
for k=1:N-2
    if k<i
        m(k)=2*k;n(k)=2*k-1;
    else
        m(k)=2*k-1;n(k)=2*k;
    end
end
tempG = permute(Grest,[m;n;In]);
Out = reshape(tempG,[prod(Nway(m)),prod(Nway(n)),tp]);
end