function CCCs = calcCCC_2_test(Data, U, V)

% calculate Canonical Correlation Coefficient
global scale;

X = Data.X;
Y = Data.Y;

for c = 1 : Data.n_class
    for j = 1:size(Y,1)
        
%         CCCs{j,1}(:,c) = abs(corr(zscore(X{j,1}) * U(:, c), zscore(Y{j,1}) * V{j,1}(:, c)));
        CCCs{j,1}(:,c) = abs(corr(normalize(X{j,1}, scale) * U(:, c), normalize(Y{j,1},scale) * V{j,1}(:, c)));
    end
end

