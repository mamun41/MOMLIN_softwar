% ----------------------------------------
% Author: Md Mamunur Rashid <mamun.stat92@gmail.com>
% Date: SEP 25, 2023
% ----------------------------------------

% close all
clear all
clc

% add data path
addpath('data/');
addpath('SupervisedPCA');
% Turns Warning Off
warning('off','all');

%% attach necessary folders: data, models, others
% ====================
ResultsFile = 'momlin_out';
if ~isdir(ResultsFile)
    mkdir(ResultsFile);
end


% =======
global scale ;

scale = "std";

%% Load BRCA triple omics data 

% ====  mRNA: X  (log2 TPM)

mRNA = importdata('01_rna_tpm.csv');    % 
X1 = mRNA.data;   
mRNA_name = mRNA.textdata(1,2:end);    


% ===  Mutation Y1

mutation = importdata('02_dna_count.csv');           % 
Y1 = mutation.data;    
mutation_name = mutation.textdata(1, 2:end); 


% ===  Clinical Y2
clinical = importdata('03_clinical_pheno.csv');
Y2 = clinical.data;    
clinical_name = clinical.textdata(1, 2:end); 

% hormon
SubtypeL = readtable('00_class_metadata.csv');
hormon_status = table2array(SubtypeL(:,"ERHER2_status"));
ER_HER2_status = string(hormon_status)=='ER- HER2-';
ER_HER2_status = double(ER_HER2_status);
ER_HER2_status(ER_HER2_status == 0) = -1;

Y2 = [Y2, ER_HER2_status];
clinical_name = [clinical_name, {'"ER-HER2-"'}];




% ===  TiME feature: Y3

TiME = importdata('04_TiME_score.csv');
Y3 = TiME.data;    
TiME_name = TiME.textdata(1, 2:end); 


% ===  Pathgway activity: Y4

path = importdata('05_p.acivity.gsva.csv');    % 
Y4= path.data;   
path_name = path.textdata(1, 2:end);   
 

% ===  sample label  : metadata preparation
SubtypeL = readtable('00_class_metadata.csv');
Slabel = table2array(SubtypeL(:,"RCB_category"));
n_class = length(unique(Slabel)); 



% Categorical phenotype (One-hat encoding: T1,..T4)                      
Class_b = [string(Slabel)=='pCR' string(Slabel)=='RCB-I'...
    string(Slabel)=='RCB-II' string(Slabel)=='RCB-III'];

% 1: pCR; 2: RCB-I; 3: RCB-II; 4: RCB-III;
clasID = double(Class_b(:,1) + 2*Class_b(:,2) + 3*Class_b(:,3) + 4*Class_b(:,4));
Z = double(Class_b);


% assigned datasets
XX1 = X1; % rna
YY1 = Y1; % mutation 
YY2 = Y2; % clinical
YY3 = Y3; % TiME
YY4 = Y4; % Pathways

% Variable Feature selection in high dim data
X1raw = XX1;     %
Y3raw = YY3;     % 
Y4raw = YY4;     % 

% calculate median absolute deviation (mad) of raw data
X1var = flipud(sortrows([mad(X1raw,1,1)',(1:size(X1raw,2))'])); 
Y3var = flipud(sortrows([mad(Y3raw,1,1)',(1:size(Y3raw,2))']));
Y4var = flipud(sortrows([mad(Y4raw,1,1)',(1:size(Y4raw,2))']));

% quantile(Y4var(:,1),.25)

X1indices2keep = floor(0.1*size(X1raw,2));  % Find top 40% of highly variant mRNAs 
Y3indices2keep = floor(1.0*size(Y3raw,2));  % 
Y4indices2keep = floor(1.0*size(Y4raw,2));  % 

XX1 = X1raw(:,sortrows(X1var(1:X1indices2keep,2)));   
YY3 = Y3raw(:, sortrows(Y3var(1:Y3indices2keep,2))); 
YY4 = Y4raw(:, sortrows(Y4var(1:Y4indices2keep,2)));   

mRNA_name = mRNA_name(sortrows(X1var(1:X1indices2keep,2)));
TiME_name = TiME_name(sortrows(Y3var(1:Y3indices2keep,2)));
path_name = path_name(sortrows(Y4var(1:Y4indices2keep,2)));


% Initiation 
[n, p] = size(XX1);        % n, p
[~, q1] = size(YY1);       % q 
[~, q2] = size(YY2);       % z 
[~, q3] = size(YY3); 
[~, q4] = size(YY4); 



%% Set tuned parameters

% weight :U for mRNA
opts.lambda_u = 0.80;            % L1-norm for X

% weight :v1_3 for DNA, Clinic, TiME, Path
opts.lambda_v1 = 0.60;          % L1-norm for Y1 
opts.lambda_v2 = 0.40;          % L1-norm for Y2
opts.lambda_v3 = 0.60;          % L1-norm for Y3
opts.lambda_v4 = 0.70;          % L1-norm for Y4


% Adjust co-expressed module
opts.beta = 0.5;            % GN-norm for X dependency (3)


% #diagnosis class
trainData.n_class = n_class;
testData.n_class = n_class;



%% ::::::::::::::: Run Main Algorithm M2SCCA  ::::::::::::::::::::::::::

% initiations
nModels = 5;      % Number of MOMLIN modules to run
k_fold = 3;       % Cross-validation number for MOMLIN
dOut_AUCs = {};
    

for ii = 1:1:nModels
% Kfold cross validation based on 70% train Data

indices = cvpartition(clasID,'KFold',k_fold,'Stratify',true);

for k = 1 : indices.NumTestSets
    fprintf('[conduct fold %d ', k);

    % Split training data and test data
    idx_test = indices.test(k);
    idx_train = ~idx_test;    

    % save reapative random indices for further use
    rCVF_train(:,k,ii) = idx_train;
    rCVF_test(:,k,ii) = idx_test;

    trainData.n_class = n_class;
    testData.n_class = n_class;

    % training sets
    trainData.X{1,1} = XX1(idx_train, :);    % mRNA
    trainData.X{2,1} = XX1(idx_train, :);    % 
    trainData.X{3,1} = XX1(idx_train, :);    % 
    trainData.X{4,1} = XX1(idx_train, :);    % 

    trainData.Y{1,1} = YY1(idx_train, :);    % mutation
    trainData.Y{2,1} = YY2(idx_train, :);    % clinic
    trainData.Y{3,1} = YY3(idx_train, :);    % TiME  
    trainData.Y{4,1} = YY4(idx_train, :);    % pathways 
    trainData.Z = Z(idx_train, :);

    % testing sets
    testData.X{1,1} = XX1(idx_test, :);   % mRNA
    testData.X{2,1} = XX1(idx_test, :);
    testData.X{3,1} = XX1(idx_test, :);
    testData.X{4,1} = XX1(idx_test, :);

    testData.Y{1,1} = YY1(idx_test, :);  % mutation
    testData.Y{2,1} = YY2(idx_test, :);  % clinic
    testData.Y{3,1} = YY3(idx_test, :);  % TiME 
    testData.Y{4,1} = YY4(idx_test, :);  % pathways 
    testData.Z = Z(idx_test, :);


    %% Train model
    tic;

    [U(:,:,k), V] = momlin_main(trainData, opts);

    V1(:,:,k) = V{1};
    V2(:,:,k) = V{2};
    V3(:,:,k) = V{3};
    V4(:,:,k) = V{4};

    time(k, 1) = toc;

    %% Calculate canonical correlation coefficients (CCCs)
    CCCs_train = calcCCC_2(trainData, U(:, :, k), V);
    CCCs_train1(:,:,k) = CCCs_train{1}; % RNA-Mutatiion
    CCCs_train2(:,:,k) = CCCs_train{2}; % RNA-clinical
    CCCs_train3(:,:,k) = CCCs_train{3}; % RNA-TiME
    CCCs_train4(:,:,k) = CCCs_train{4}; % RNA-path

    CCCs_test = calcCCC_2(testData, U(:, :, k), V);
    CCCs_test1(:,:,k) = CCCs_test{1}; % mRNA-mutation
    CCCs_test2(:,:,k) = CCCs_test{2}; % mRNA-clinical
    CCCs_test3(:,:,k) = CCCs_test{3}; % mRNA-molecule
    CCCs_test4(:,:,k) = CCCs_test{4}; % mRNA-molecule

    fprintf('(%.2fs)]\n', time(k));

end 

U_mean_i(:,:,ii) = mean(U,3);
V1_mean_i(:,:,ii) = mean(V1,3);
V2_mean_i(:,:,ii) = mean(V2,3);
V3_mean_i(:,:,ii) = mean(V3,3);
V4_mean_i(:,:,ii) = mean(V4,3);


mean_CCCs_train_i(1,:,ii) = mean(CCCs_train1,3);
mean_CCCs_train_i(2,:,ii) = mean(CCCs_train2,3);
mean_CCCs_train_i(3,:,ii) = mean(CCCs_train3,3);
mean_CCCs_train_i(4,:,ii) = mean(CCCs_train4,3);

mean_CCCs_test_i(1,:,ii) = mean(CCCs_test1,3);
mean_CCCs_test_i(2,:,ii) = mean(CCCs_test2,3);
mean_CCCs_test_i(3,:,ii) = mean(CCCs_test3,3);
mean_CCCs_test_i(4,:,ii) = mean(CCCs_test4,3);

end



% ------ Calculate Canonical Weights -----------
U_mean = sum(U_mean_i, 3)./sqrt(sum(sum(U_mean_i, 3).^2,1)); % mRNA
% ---------
V1_mean = sum(V1_mean_i, 3)./sqrt(sum(sum(V1_mean_i, 3).^2,1)); % Dana
% ---------
V2_mean = sum(V2_mean_i, 3)./sqrt(sum(sum(V2_mean_i, 3).^2,1)); % Clinical
% ---------
V3_mean = sum(V3_mean_i, 3)./sqrt(sum(sum(V3_mean_i, 3).^2,1)); % TiME 
% ---------
V4_mean = sum(V4_mean_i, 3)./sqrt(sum(sum(V4_mean_i, 3).^2,1)); % Pathways


% Correlation
% training
CCCs_mean_train = mean(mean_CCCs_train_i,3);  % mRNA-DNA
% testing
CCCs_mean_test = mean(mean_CCCs_test_i,3);  % mRNA-DNA

% save MOMLIN loadings
% save(strcat(ResultsFile,"/result_momlin.mat"))



    %%  canonical weight heatmap vis. 
    figure(); colormap(jet);
    caxis_range = 0.50;
    subplot(5, 1, 1); imagesc(U_mean', [-caxis_range caxis_range]);
    colorbar('Ticks', [-caxis_range 0 caxis_range]);
    % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.YTick = 1:n_class;
    % Set TickLabels;
    ax.YTickLabel = {'pCR','RCB-I','RCB-II','RCB-III'};
    xlabel('mRNA', 'Color', 'k', 'FontSize',14);
    ylabel('Response class', 'Color', 'k', 'FontSize',14);
    yL.FontSize = 20; xL.FontSize = 20;
    %title('Proteome canonical weights','Color','k', 'FontSize',20)

    %caxis_range = 0.9;
    colormap(jet);
    subplot(5, 1, 2); imagesc(V1_mean', [-caxis_range caxis_range]);
    colorbar('Ticks', [-caxis_range 0 caxis_range]);
    % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.YTick = 1:n_class;
    % Set TickLabels;
    ax.YTickLabel = {'pCR','RCB-I','RCB-II','RCB-III'};
    xlabel('DNA', 'Color', 'k', 'FontSize',14);
    ylabel('Response class', 'Color', 'k', 'FontSize',14);
    yL.FontSize = 20; xL.FontSize = 20;
    %title('Metabolome canonical weights','Color','k', 'FontSize',20)

    %caxis_range = 0.9;
    colormap(jet);
    subplot(5, 1, 3); imagesc(V2_mean', [-caxis_range caxis_range]);
    colorbar('Ticks', [-caxis_range 0 caxis_range]);
    % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.YTick = 1:n_class;
    % Set TickLabels;
    ax.YTickLabel = {'pCR','RCB-I','RCB-II','RCB-III'};
    xlabel('Clinical', 'Color', 'k', 'FontSize',14);
    ylabel('Response class', 'Color', 'k', 'FontSize',14);
    yL.FontSize = 20; xL.FontSize = 20;

    %caxis_range = 0.9;
    colormap(jet);
    subplot(5, 1, 4); imagesc(V3_mean', [-caxis_range caxis_range]);
    colorbar('Ticks', [-caxis_range 0 caxis_range]);
    % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.YTick = 1:n_class;
    % Set TickLabels;
    ax.YTickLabel = {'pCR','RCB-I','RCB-II','RCB-III'};
    xlabel('TiME', 'Color', 'k', 'FontSize',14);
    ylabel('Response class', 'Color', 'k', 'FontSize',14);
    yL.FontSize = 20; xL.FontSize = 20;


    %caxis_range = 0.9;
    colormap(jet);
    subplot(5, 1, 5); imagesc(V4_mean', [-caxis_range caxis_range]);
    colorbar('Ticks', [-caxis_range 0 caxis_range]);
    % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.YTick = 1:n_class;
    % Set TickLabels;
    ax.YTickLabel = {'pCR','RCB-I','RCB-II','RCB-III'};
    xlabel('Pathways', 'Color', 'k', 'FontSize',14);
    ylabel('Response class', 'Color', 'k', 'FontSize',14);
    yL.FontSize = 20; xL.FontSize = 20;



