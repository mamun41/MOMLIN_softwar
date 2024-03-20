function [X, Y, Z] = do_oversample_2(rawData)

% rawData = trainData_m2scca
% Oversample for class balancing
for ii = 1:size(rawData.Y,1)

    X1 = rawData.X{ii,1}; Y1 = rawData.Y{ii,1}; Z1 = rawData.Z;
    % X = trainData.X; Y = trainData.Y; Z = trainData.Z;

    for c = 1 : size(Z1, 2)
        idx_p = find(Z1(:, c) == 1);
        num_p = length(idx_p);     % Num. of positive samples
        idx_n = find(Z1(:, c) == 0);
        num_n = length(idx_n);     % Num. of negative samples
    
        if num_p > num_n
            idx_oversample = idx_n(randi([1, num_n], num_p - num_n, 1));
        elseif num_n > num_p
            idx_oversample = idx_p(randi([1, num_p], num_n - num_p, 1));
        else
            idx_oversample = [];
        end
    
        tempX{c} = [X1; X1(idx_oversample, :)];
        tempY{c} = [Y1; Y1(idx_oversample, :)];
        tempZ{c} = [Z1(:, c); Z1(idx_oversample, c)];
    end

X{ii,:} = tempX; Y{ii,:} = tempY; Z{ii,1} = tempZ;

end
