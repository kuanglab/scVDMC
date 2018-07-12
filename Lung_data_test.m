clc;
clear;
close all;

load Lung_data;

%% prepare data
d = 3;
n = zeros(d,1);
sp_data = cell(d,1);
sp_label = cell(d,1);
for dd = 1:d
    sp_data{dd} = cleandata(CellRep == dd, :)';
    [m, n(dd)] = size(sp_data{dd});
    sp_label{dd} = Celltype(CellRep == dd);
end

%% sparate data
% remove ciliated since only shown in data3
% all data have 4 labels
ix = find(strcmp(sp_label{3}, 'ciliated'));
sp_data{3}(:,ix) = [];
sp_label{3}(ix) = [];
n(3) = n(3) - length(ix);

k = 4;
a = perms(1:k);
% sort sample by cell type
for dd = 1:d
    [sp_label{dd}, IX] = sort(sp_label{dd});
    sp_data{dd} = sp_data{dd}(:,IX);
end
%% get unuseable features
% unnormalized separate data
sp_data_orig = sp_data;

nanix = [];
for dd = 1:d
    mD = mean(sp_data{dd}, 2);
    vD = std(sp_data{dd}, [], 2);
    sp_data{dd} = (sp_data{dd} - repmat(mD, 1, n(dd))) ./ repmat(vD, 1, n(dd));
    nanix = union(nanix, find(vD == 0));
end

for dd = 1:d
    sp_data{dd}(nanix,:) = [];
    %sp_data_orig{dd}(nanix,:) = [];
end

m_clean = m - length(nanix);
genelist_clean = genelist;
genelist_clean(nanix, :) = [];

%% get pool, unnormalized
% pool data
pool_data = cell2mat(sp_data_orig');
pool_label = [];
for dd = 1:d
    pool_label = [pool_label; sp_label{dd}];
end

%% original data, check the center and means
V_True_sp = cell(d, 1);
labels = unique(pool_label);
for dd = 1:d
    V_True_sp{dd} = zeros(n(dd), k);
    for i = 1:k
        V_True_sp{dd}(strcmp(labels{i}, sp_label{dd}), i) = 1;
    end
end

V_True_pool = cell2mat(V_True_sp);
%% pool_mean
IX = kmeans(pool_data', k, 'Distance', 'correlation', 'Replicates', 20);
V_PKMS = zeros(sum(n), k);
 for i = 1:k
        V_PKMS(IX == i, i) = 1;
 end
err = zeros(size(a,1), 1);
    for i = 1:size(a,1)
        err(i) = length(find(sum(V_True_pool ~= V_PKMS(:, a(i,:)), 2)));
    end
[err_PKMS, ix] = min(err);
V_PKMS = V_PKMS(:, a(ix(1),:));       
%% using kmeans as our start point
err_bag = zeros(d, 1);
U_bag = cell(d, 1);
V_bag = cell(d, 1);
V_ini = cell(d,1);
U_ini = cell(d,1);
iix = [0; cumsum(n)];
V_ini{1} = V_PKMS(iix(1)+1:iix(2),:);
V_ini{2} = V_PKMS(iix(2)+1:iix(3),:);
V_ini{3} = V_PKMS(iix(3)+1:iix(4),:);
for dd = 1:d
    for kk = 1:k
        U_ini{dd}(:, kk)  = mean(sp_data{dd}(:, V_ini{dd}(:,kk) == 1), 2);
    end
end


%% run scVDMC algorithm
w = 1; % hyperparameters
lambda = 40;
alpha=1;
max_iter = 50;
[U, V, Beta, sortBeta, Obj]...
            = scVDMC(sp_data, d, k, w, lambda, alpha, U_ini, V_ini, max_iter);
for dd = 1:d            
    err = zeros(size(a,1), 1);
    for i = 1:size(a,1)
        err(i) = length(find(sum(V_True_sp{dd} ~= V{dd}(:, a(i,:)), 2)));
    end
    [err_bag(dd), ix] = min(err);
    V_bag{dd} = V{dd}(:, a(ix(1),:));
    U_bag{dd} = U{dd}(:, a(ix(1),:));            
end
sorted_markergenes = genelist_clean(sortBeta);


