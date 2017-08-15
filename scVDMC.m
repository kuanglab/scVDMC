function [U, V, B, sortB, obj] = scVDMC(X, d, k, w, lambda, alpha, U_ini, V_ini, max_iter)
% Objective function as follows
% 1/2 sum_d ||diag(sqrt(B)) (X{d} - U{d} V{d}') ||_F^2 - w sum_d B'*Var(U_d)
% +alpha sum_{i,j} B_i Var(Y{i,j})
% s.t.
% sum(V(d)(i,:))=1 and V(d) is binary matrix
% sum(B) = lambda and B is binary vector
%
%
% m: number of features;
% n(d): number of samples in domain d;
% k: number of clusters;
% X(d): m by n_d matrix for domain d;
% U(d): m by k matrix
% V(d): n(d) by k matrix
% B: m by 1 binary vector
% Var(U(d)): vector with variance of U_d(i,:)
% U_ini, V_ini: initialization of U and V
% max_iter: max number of iteration
% obj: objective function value
% sortBeta:

%% initialization
n = zeros(d, 1);
for dd = 1:d
    [m, n(dd)] = size(X{dd});
end;

% initial of U
if ~exist('U_ini', 'var') || isempty(U_ini)
    % set U to be centers of k-means clusters
    IX = cell(d,1);
    U = cell(d,1);
    for dd = 1:d
        [~, SCORE, ~] = pca(X{dd}');
        IX{dd} = kmeans(SCORE, k, 'Distance', 'correlation', 'Replicates', 20);
        for kk = 1:k
            U{dd}(:, kk)  = mean(X{dd}(:, IX{dd} == kk), 2);
        end;
    end;
else
    U = U_ini;
end;

% initial V
if ~exist('V_ini', 'var') || isempty(V_ini)
    % set V to be all zeros
    V = cell(d,1);
    for dd = 1:d
        V{dd} = zeros(n(dd),k);
    end;
else
    V = V_ini;
end;

% inital B
B = ones(m,1) * lambda/m;

% avoid empty class
for dd = 1:d
    [U{dd}, V{dd}] = clearmty(U{dd}, V{dd}, X{dd});
end;


%% algorithm
obj = zeros(max_iter, 4);
V_old = cell(d,1);
U_old = cell(d,1);
for iter = 1: max_iter
    % construct Y
    Y=zeros(m,k);
    for i=1:m
        for j=1:k
            vec=zeros(d,1);
            for L=1:d
                vec(L)=U{L}(i,j);
            end
            Y(i,j)=var(vec,1);
        end
    end
    Y_vec=sum(Y,2); 
    obj(iter,4)=alpha*B'*Y_vec;
    % solve B
    B_old = B;
    var_vector = zeros(m, 1);
    for dd = 1:d
        var_vector = var_vector + var(U{dd}, 1, 2);
    end;
    A_vector = zeros(m,1);
    for dd = 1:d
        A = X{dd} - U{dd}*V{dd}';
        A_vector = A_vector + (diag(A*A'));
    end;
    
    
    thef = 0.5*A_vector-w*var_vector+alpha*Y_vec;
    
    [~, ix] = sort(thef);
    tops = floor(lambda);
    % handle case when B is not integer
    left = lambda - tops;
    B = zeros(m, 1);
    B(ix(1:tops)) = 1;
    B(ix(tops+1)) = left;
    
    if left > 0
        sortB = ix(1:(tops+1));
    else
        sortB = ix(1:tops);
    end;
    
    Uold=U;
    for dd = 1:d
        V_old{dd} = V{dd};
        U_old{dd} = U{dd};
        
        % solve V
        V{dd} = SolveV(X{dd}, U{dd}, V{dd}, B);
        
        % avoid empty V class
        [U{dd}, V{dd}] = clearmty(U{dd}, V{dd}, X{dd});
        
        % solve U
        U{dd} = SolveU(X{dd}, V{dd}, Uold, w, alpha, B ,dd);
    end;
    
    % Calculate objective function
    disp(['iter:' num2str(iter) ' finished']);
    [obj(iter,1), obj(iter,2), obj(iter,3)] = objfunction(X, U, V, w, B);
    % check convergence
    U_cov = zeros(d,1);
    V_cov = zeros(d,1);
    for dd = 1:d
        U_cov(dd) = norm(U{dd} - U_old{dd}, 'fro')/norm(U_old{dd}, 'fro');
        V_cov(dd) = norm(V{dd} - V_old{dd}, 'fro')/norm(V_old{dd}, 'fro');
    end;
    Beta_cov = norm(B - B_old)/norm(B_old);
    con = max([max(U_cov), max(V_cov), Beta_cov]);
    disp(['residue: ' num2str(con)]);
    if con < 1e-3
        break;
    end;
end;

if iter == max_iter
    disp(['algo didn''t converge in ' int2str(max_iter) ' iterations.']);
else
    disp(['algo converge in ' int2str(iter) ' iterations.']);
end;

end

function Ud = SolveU(Xd, Vd, U, w, alpha, B ,dd)
Ud=U{dd};
[~, K] = size(Ud);
d=length(U);
phi=eye(d)-ones(d)/d;
ker1 = inv(Vd'*Vd-2*w/K*(eye(K)-ones(K)/K)+2*alpha*K*(1-1/d)/d*eye(K));
IX = find(B > 0);
% only solve selected features
for ix = 1:length(IX)
    r = IX(ix);
    ker2=zeros(K,1);
    for L=1:d
        ker2=ker2+phi(dd,L)*U{L}(r,:)';
    end
    ker2=ker2-phi(dd,dd)*Ud(r,:)';
    ker2=2*alpha*K/d*ker2;
    Ud(r, :) = ker1*(Vd'*Xd(r,:)'-ker2);
end
end

function V = SolveV(X, U, V, Beta)
[n, k] = size(V);
V = zeros(n, k);

for nn = 1 : n
    MSE = zeros(k,1);
    for kk = 1:k
        MSE(kk) = norm(sqrt(Beta) .* X(:, nn) - sqrt(Beta) .* U(:, kk));
    end;
    [~, IX] = min(MSE);
    V(nn, IX) = 1;
end
end

function [obj_all, recon_err, var_val] = objfunction(X, U, V, w, Beta)
% reconstruction error
recon_err = 0;
% variance
var_val = 0;
d = length(X);
for dd = 1:d
    recon_err = recon_err + 1/2*norm(diag(sqrt(Beta))*(X{dd} - U{dd}*V{dd}'), 'fro')^2;
    var_val = var_val + w * var(U{dd}, 1, 2)'* Beta;
end;
obj_all = recon_err - var_val;
end

function [U, V] = clearmty(U, V, X)
clu_size = sum(V);
IX = find(clu_size == 0);
while ~isempty(IX)
    disp('empty cluster correction!');
    [~, big_clu] = max(clu_size);
    SID = find(V(:,big_clu));
    newIX = kmeans(X(:, SID)', 2, 'Distance', 'correlation', 'Replicates', 20);
    % split SID into 2 cluster
    ept_clu = IX(1);
    V(SID(newIX == 2), ept_clu) = 1;
    V(SID(newIX == 2), big_clu) = 0;
    
    U(:, ept_clu) = mean(X(:,SID(newIX == 2)),2);
    U(:, big_clu) = mean(X(:,SID(newIX == 1)),2);
    
    clu_size = sum(V);
    IX = find(clu_size == 0);
end;
end
