clc;
clear;
close all;

%%
% This file contains single cell RNA-seq expression data (log2(FPKM) values)
% for all 80 lung epithelial cells at E18.5 together with the putative cell type 
% of each cell in a .txt file.
T = readtable('nature13173-s4.txt', 'delimiter', '\t');
genelist_full = T.Properties.VariableNames(5:end)';
Cellname = table2cell(T(1:80, 1));
Celltype = table2cell(T(1:80, 4));

% replicate information
CellRep = zeros(80,1);
for i = 1:80    
    tmp = regexp(Cellname{i}, '_', 'split');
    CellRep(i) = str2double(tmp{2});
end;

Celltype_list = unique(Celltype);
data = table2array(T(1:80, 5:end));

%% keep high mean features, keep rare but high value features
genelist = genelist_full((sum(data > 1) > 2));
cleandata = data(:, (sum(data > 1) > 2));

save Lung_data data cleandata genelist Cellname Celltype Celltype_list CellRep;