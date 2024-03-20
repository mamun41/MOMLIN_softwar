function Penalty = get_connectivity(Data,alpha)
% data : expression (n*p)
% alpha : 2

% X1 = corr(Data,'Type','Spearman').^alpha;  % correlation mat
X1 = corr(Data).^alpha;  % correlation mat
diag_X1=diag(diag(X1));  % corr. diag. elements

X1=X1-diag_X1;           % set diagonal elements 0
L1=diag(sum(X1,2));      % degree mat > rowSums

Penalty=L1-X1;           % Laplacian computation
end