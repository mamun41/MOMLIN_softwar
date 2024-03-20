function [U, V, W] = momlin_main(data, opts)

global scale ;


% Multi-task Sparse Canonical Correlation Analysis
for i = 1:size(data.Y,1) 
    [~, p(i)] = size(data.X{i,1});  % mRNA data
    [~, q(i)] = size(data.Y{i,1});  % miRNA, Methy & CNV
end

% disease/diagnosis class
n_class = data.n_class;

% mRNA weight matrix
U = ones(p(1), n_class);


% miRNA, Methy, & CNV weights matrix
for j = 1:size(data.Y,1)
    V{j,1} = ones(q(j), n_class);
end


W = ones(1, n_class);


% Calculate connectivity penalty
for j = 1:size(data.Y,1)
    Dgn_u{j,1} = get_connectivity(data.X{j,1},2);   % GN penalty for X
end

% class balancing by oversampling
[X, Y, Z] = do_oversample_2(data);       


% initial covariance & scale matrix
for c = 1 : n_class        % n_class, disease condition

    for i = 1:size(Y,1)    % i, different omics type
        % normalization
        % ops: norm, rnorm, center, std, rstd, minmax, median, robust
        X{i,1}{c} = normalize(X{i,1}{c}, scale); 
        Y{i,1}{c} = normalize(Y{i,1}{c}, scale);
        Z{i,1}{c} = Z{i,1}{c};  
        
        % Pre-calculate covariance
        XX{i,1}{c} = X{i,1}{c}' * X{i,1}{c};
        XY{i,1}{c} = X{i,1}{c}' * Y{i,1}{c};
        YY{i,1}{c} = Y{i,1}{c}' * Y{i,1}{c};
        YX{i,1}{c} = XY{i,1}{c}';
        
        ZZ{i,1}{c} = Z{i,1}{c}' * Z{i,1}{c};
        XZ{i,1}{c} = X{i,1}{c}' * Z{i,1}{c};
        YZ{i,1}{c} = Y{i,1}{c}' * Z{i,1}{c};
        ZX{i,1}{c} = XZ{i,1}{c}';
        ZY{i,1}{c} = YZ{i,1}{c}';
        
        % Scale U and V
        U(:, c) = U(:, c) / norm(X{i,1}{c} * U(:, c));
        V{i,1}(:, c) = V{i,1}(:, c) / norm(Y{i,1}{c} * V{i,1}(:, c));
        W(:, c) = W(:, c) / norm(Z{i,1}{c} * W(:, c));
        
        % initialize adaptive(/loss) weights
        Xu = X{i,1}{c} * U(:, c); Yv = Y{i,1}{c} * V{i,1}(:, c); Zw = Z{i,1}{c} * W(:, c);
        omegaij{i,1}{c} = 1 / norm(Xu - Yv);
        omegaiz{i,1}{c} = 1 / norm(Xu - Zw);
        omegajz{i,1}{c} = 1 / norm(Yv - Zw);
        
    end  
end


% Set tuned parameters
% for mRNA
lambda_u = opts.lambda_u;   % L1 norm
% Other modality
lambda_v{1} = opts.lambda_v1; % L1 norm
lambda_v{2} = opts.lambda_v2; % L1 norm
lambda_v{3} = opts.lambda_v3; % L1 norm
lambda_v{4} = opts.lambda_v4; % L1 norm

beta = opts.beta;  % laplacian constrain


% Set iterative termination criterion
max_Iter = 100;
t = 0;
tol = 1e-5;
tu = inf;
tv = inf;

% Iteration
while (t < max_Iter && (tu > tol || tv > tol))
    
    t = t + 1;
    
    U_old = U;
    V_old = V;   

    % Update U
    
    for c = 1 : n_class
        Du_c = updateD2(U(:, c)); % L1-norm
        
        for i = 1:size(Y,1) 
            % solve each u_c
            Fu = lambda_u * beta* Du_c + lambda_u * (1-beta) * Dgn_u{i,1} + XX{i,1}{c};
            bu = omegaij{i,1}{c} * XY{i,1}{c} * V{i,1}(:, c) + omegaiz{i,1}{c} * XZ{i,1}{c};

            U(:, c) = Fu \ bu;

            % scale each u_c
            U(:, c) = U(:, c) / norm(X{i,1}{c} * U(:, c));
        end
    end
    
    
    
    % Update V
    
    for c = 1 : n_class      % c: task, T1,T2,T3,T4 & T5 
        
        for i = 1:size(Y,1)  % i : diff modality (CNV, Methy, miRNA)
            Dv_ic = updateD2(V{i,1}(:,c));   % L1-norm
           
            Fv = lambda_v{i}*Dv_ic + YY{i,1}{c}; 
            bv = omegaij{i,1}{c} * YX{i,1}{c} * U(:, c) + omegajz{i,1}{c} * YZ{i,1}{c};

            V{i,1}(:, c) = Fv \ bv;

            % scale each v_c
            V{i,1}(:, c) = V{i,1}(:, c) / norm(Y{i,1}{c} * V{i,1}(:, c));
        end
    end
    
    
    
    % Update W
    for c = 1 : n_class
        for i = 1:size(Y,1) 
            % solve each w_c
            Fw = ZZ{i,1}{c};
            bw = ZX{i,1}{c} * U(:, c) + ZY{i,1}{c}* V{i,1}(:, c);
            W(:, c) = Fw \ bw;
            % scale each u_c
            W(:, c) = W(:, c) / norm(Z{i,1}{c} * W(:, c));
        end
    end
    
      
    % Iteration termination

    tu = max(max(abs(U - U_old)));
    for i = 1:size(Y,1)
        tv(i) = max(max(abs(V{i,1} - V_old{i,1})));
    end
    tv = max(tv);
    
end






