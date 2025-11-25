function D = updateD2(W)

if nargin == 1
    % for L2,1-norm & L1,1-norm
    d = 1 ./ sqrt(sum(W .^ 2, 2) + eps);  % sum(w.^2,2) --> indicated row sums
else
    error('Error type.');
end

D = diag(0.5 * d);