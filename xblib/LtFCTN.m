function [X,G,Out] = LtFCTN(F,Omega,optsLtFCTN)
if isfield(optsLtFCTN, 'tol');         tol   = optsLtFCTN.tol;              end
if isfield(optsLtFCTN, 'maxit');       maxit = optsLtFCTN.maxit;            end
if isfield(optsLtFCTN, 'Xtrue');       XT    = optsLtFCTN.Xtrue;            end
if isfield(optsLtFCTN, 'rho');         rho   = optsLtFCTN.rho;              end
if isfield(optsLtFCTN, 'mu');          mu    = optsLtFCTN.mu;               end
if isfield(optsLtFCTN, 'P');           P     = optsLtFCTN.P;                end
if isfield(optsLtFCTN, 'R');           R     = optsLtFCTN.R;                end


F2=permute(F,[1,2,4,3]);
XT=permute(XT,[1,2,4,3]);
N = ndims(F2);
X = F2;
Tp=P(:,end);
[N1,N2,N3,N4]= size(F2);
R1=P(1);R2=P(2);R3=P(3);R4=P(4);
E=rand(R1,R2,R3,R4);
NwayE      = size(E);

%%%%%%%%%%%%%%%%% 
DCT_matrix1 = dctmtx(N1);
P1 = DCT_matrix1(1:R1, :);
DCT_matrix2 = dctmtx(N2);
P2 = DCT_matrix2(1:R2, :);
DCT_matrix3 = dctmtx(N3);
P3 = DCT_matrix3(1:R3, :);
DCT_matrix4 = dctmtx(N4);
P4 = DCT_matrix4(1:R4, :);


%%%%%%%%%%%%%%%%%
U = {P1, P2, P3, P4};
R11= R(1:N-1,1:N-1);
TPdim = diag(NwayE(1:N-1))+R11+R11';
tempdim=zeros(N,N);
tempdim(1:N-1,1:N-1)=TPdim;
tempdim(N, :) = R(N, :);
tempdim(:, N) = R(:, N);

G = cell(1,N);
for i = 1:N
    G{i} = rand(tempdim(i,:));
end

AllMode = 1:N;         
Out.RSE = [];

for k = 1:maxit
    Xold = X;
    Eold = E;
    Uold = U;
    UT = cellfun(@(x) x', U, 'UniformOutput', false);
    % Update C
    for i = 1:N-1
        Ei = my_xbUnfold(E,NwayE,i);
        Gi = my_xbUnfold(G{i},tempdim(i,:),i);
        Girest = tnxbreshape(sub_xbTN(G,i),N,i);
        tempC = tprod(Ei,tran(Girest))+tprod(Gi,xbteye(rho,size(Gi,2),Tp));
        tempA = tprod(Girest,tran(Girest))+xbteye(rho,size(Gi,2),Tp);
        TGG=tprod(tempC ,tinv(tempA));
        G{i}  = my_xbFold(TGG,tempdim(i,:),i);
    end

    % Update Y
    XALLD=tmodexb(Xold, Uold(AllMode), AllMode);
    E = (TN_xbcomposition(G)+mu*XALLD+rho*Eold)/(1+mu+rho);

    % Update X
    EALLD=tmodexb(E, UT(AllMode), AllMode);
    X = ( mu*EALLD+rho*Xold)/(mu+rho);
    X(Omega) = XT(Omega);
    % Update E
    for j = 1:N
        OthMode = setdiff(AllMode, j); 
        Q = tmodexb(E, UT(OthMode), OthMode); 
        X_mat = unfold(X, j);  
        Q_mat = unfold(Q, j);  
        PT = Q_mat * X_mat' + Uold{j};
        [U_svd, ~, V_svd] = svd(PT, 'econ');
        U{j} = U_svd * V_svd';
    end



    %% check the convergence
    rse=norm(X(:)-Xold(:))/norm(Xold(:));
    Out.RSE = [Out.RSE,rse];
     if mod(k, 100) == 0  ||   k == 1 
       fprintf('LtFCTN: iter = %d   RSE=%f   \n', k, rse);
     end
    if rse < tol
        break;
    end

end
X=ipermute(X,[1,2,4,3]);
end




function X_unf = unfold(X, mode)
sz = size(X);
order = [mode, 1:mode-1, mode+1:numel(sz)];  
X_perm = permute(X, order);                  
X_unf = reshape(X_perm, sz(mode), []);     
end


function Y = tmodexb(X, V, dims)
Y = X;
for i = 1:numel(dims)
    Y = tmodexb_single(Y, V{i}, dims(i));  
end
end


function Y = tmodexb_single(X, V, n)
sz = size(X);
perm_order = [n, 1:n-1, n+1:ndims(X)];   
X_perm = permute(X, perm_order);         
X_reshaped = reshape(X_perm, sz(n), []); 
Y_data = V * X_reshaped;                 
new_size = [size(V, 1), sz([1:n-1, n+1:end])];  
Y = ipermute(reshape(Y_data, new_size), perm_order); 
end

