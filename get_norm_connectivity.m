function Penalty = get_norm_connectivity(Data,alph)
% data : expression (n*p)

X1 = corr(Data).^alph;  % correlation mat
diag_X1=diag(diag(X1));  % corr. diag. elements

% Adjency mat
X1=X1-diag_X1;           % set diagonal elements 0

% degree diagonal mat
D=diag(sum(X1,2));      % degree mat > rowSums

% Laplacian computation
L =D-X1;           

% Compute the inverse square root of D
D_inv_sqrt = D^(-1/2);

% Handle possible divide-by-zero in D^(-1/2) for isolated nodes
D_inv_sqrt(isinf(D_inv_sqrt)) = 0;

% Compute the normalized Laplacian matrix L_norm
Penalty = D_inv_sqrt * L * D_inv_sqrt;

end