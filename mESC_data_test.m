clc;
clear;
close all;

load mESC_data.mat;

%% prepare data, convert label to matrix
% working on sub labels in serum
serum_data = log10(data(:, 455:end)+1);
serum_label = SID_class_label(455:end,:);
[m, n] = size(serum_data);
k = 3;
a = perms(1:k);
d = 3;
max_iter = 40;
alpha=0.1;
V_True = zeros(n, k);
labels = unique(serum_label(:,7));
for i = 1:k
    V_True(strcmp(labels{i}, serum_label(:,7)), i) = 1;
end;

%% separate data and normalize
sp_data = cell(d,1);
sp_label = cell(d,1);
sp_data{1} = serum_data(:, cell2mat(serum_label(:,6)) == 1);
sp_data{2} = serum_data(:, cell2mat(serum_label(:,6)) == 2);
sp_data{3} = serum_data(:, cell2mat(serum_label(:,6)) == 3);
sp_data_orig = sp_data;
sp_label{1} = serum_label(cell2mat(serum_label(:,6)) == 1,:);
sp_label{2} = serum_label(cell2mat(serum_label(:,6)) == 2,:);
sp_label{3} = serum_label(cell2mat(serum_label(:,6)) == 3,:);
nn = cellfun(@length, sp_label);
V_True_sp = cell(d,1);
for dd = 1:d
    V_True_sp{dd} = zeros(nn(dd), k);    
    for i = 1:k        
        V_True_sp{dd}(strcmp(labels{i}, sp_label{dd}(:,7)), i) = 1;
    end;
end;

% remove non variance feature
nanix = [];
for dd = 1:d
    mD = mean(sp_data{dd}, 2);
    vD = std(sp_data{dd}, [], 2);
%     sp_data{dd} = sp_data{dd} ./ repmat(vD, 1, nn(dd));
    sp_data{dd} = (sp_data{dd} - repmat(mD, 1, nn(dd))) ./ repmat(vD, 1, nn(dd));
%     sp_data{dd} = sp_data{dd} - repmat(mD, 1, nn(dd));
%     sp_data{dd} = sp_data{dd} ./ repmat(sqrt(diag(sp_data{dd}*sp_data{dd}')), 1, nn(dd));
    nanix = union(nanix, find(vD == 0));
end;

for dd = 1:d
    sp_data{dd}(nanix,:) = [];
    sp_data_orig{dd}(nanix,:) = [];
end;
GID_clean = GID;
GID_clean(nanix) = [];
m_clean = length(GID_clean);


%% pool_mean method
pool_data = serum_data;
pool_data(nanix,:) = [];
IX = kmeans(pool_data', k, 'Distance', 'correlation', 'Replicates', 20);   
V_PKMS = zeros(n, k);
for i = 1:k
        V_PKMS(IX == i, i) = 1;
end;  
V_PKMS = V_PKMS(:, a(ix(1),:));

%% using kmeans as our start point
wpool = [0 0.1 0.5 1];
lambdapool = [10:10:100];
err_TF = zeros(d, length(wpool), length(lambdapool));
U_TF = cell(d, length(wpool), length(lambdapool));
V_TF = cell(d, length(wpool), length(lambdapool));
TF_sorted_gene = cell(length(wpool), length(lambdapool));
TF_sorted_IX = cell(length(wpool), length(lambdapool));

V_ini = cell(d,1);
U_ini = cell(d,1);
V_ini{1} = V_PKMS(cell2mat(serum_label(:,6)) == 1,:);
V_ini{2} = V_PKMS(cell2mat(serum_label(:,6)) == 2,:);
V_ini{3} = V_PKMS(cell2mat(serum_label(:,6)) == 3,:);
for dd = 1:d
    for kk = 1:k
        U_ini{dd}(:, kk)  = mean(sp_data{dd}(:, V_ini{dd}(:,kk) == 1), 2);
    end;
end;

%% scVDMC algorithm
for wi = 1:length(wpool)
    w = wpool(wi);
    for li = 1:length(lambdapool)
        lambda = lambdapool(li);
        [U_TF_once, V_TF_once, Beta_TF_once, sortBeta_TF_once, Obj_TF_once] = ...
            scVDMC(sp_data, d, k, w, lambda,alpha,U_ini, V_ini, max_iter);
        
        for dd = 1:d            
            err = zeros(size(a,1), 1);
            for i = 1:size(a,1)
                err(i) = length(find(sum(V_True_sp{dd} ~= V_TF_once{dd}(:, a(i,:)), 2)));
            end;
            [err_TF(dd,wi,li), ix] = min(err);
            V_TF{dd,wi,li} = V_TF_once{dd}(:, a(ix(1),:));
            U_TF{dd,wi,li} = U_TF_once{dd}(:, a(ix(1),:));            
        end;
        TF_sorted_gene{wi, li} = GID_clean(sortBeta_TF_once);
        TF_sorted_IX{wi, li} = sortBeta_TF_once;
    end;
end;
save('mESC_scVDMC.mat',...
    'wpool', 'lambdapool',...
    'TF_sorted_gene', 'TF_sorted_IX', ...
    'err_TF', 'U_TF', 'V_TF');
FigHandle = figure('Position', [100, 100, 1000, 600]);
plot(squeeze(sum(err_TF, 1))', 'linewidth', 2, 'LineStyle', '-.');
set(gca, 'xtick', 1:length(lambdapool), 'XTickLabel', lambdapool);
title('mESC data');
xlabel('\lambda');
ylabel('err');  



