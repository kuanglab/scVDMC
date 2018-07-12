clc;
clear;
close all;

%%
% This file contains single cell RNA-seq expression data (log2(FPKM) values)
% for all 80 lung epithelial cells at E18.5 together with the putative cell type 
% of each cell in a .txt file.
fileID = fopen('nature13173-s4.txt');
tline = fgets(fileID);
head = strsplit(tline);
genelist_full = strrep(head(5:end-1)','"','');
fmt = [repmat('%s', 1, 4), repmat('%f', 1, length(genelist_full))];
C = textscan(fileID,fmt);
fclose(fileID);
Cellname =  strrep(C{1}(1:80),'"','');
Celltype = strrep(C{4}(1:80),'"','');
data = cell2mat(C(5:end));
data = data(1:80,:);

% replicate information
CellRep = zeros(80,1);
for i = 1:80    
    tmp = regexp(Cellname{i}, '_', 'split');
    CellRep(i) = str2double(tmp{2});
end

Celltype_list = unique(Celltype);

%% keep high mean features, keep rare but high value features
genelist = genelist_full((sum(data > 1) > 2));
cleandata = data(:, (sum(data > 1) > 2));

save Lung_data data cleandata genelist Cellname Celltype Celltype_list CellRep;