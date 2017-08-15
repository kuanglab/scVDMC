clc;
clear;
close all;

%%
T1 = readtable('cell_states_conditions.txt', 'ReadVariableNames', 0, 'delimiter', ' ');
T1 = table2cell(T1(:,2:3));
cellName = T1(:,1);
tmp = regexp(T1(:,1), '_', 'split');
tmp1 = [];
for i = 1:length(tmp)
    tmp1 = [tmp1; {tmp{i}{3} tmp{i}{4} tmp{i}{5}(1:end-7)}];
end;

addlabel = [tmp1 T1(:,2)];

%%
T = readtable('nCountGenesBatchAdjusted.csv', 'ReadRowNames', 1);

SID = T.Properties.VariableNames';
GID = T.Properties.RowNames;
data = table2array(T);

%%
tmp = regexp(SID, '_', 'split');
tmp1 = [];
for i = 1:length(tmp)
    tmp1 = [tmp1; tmp{i}];
end;
tmp1 = [tmp1(:,1) tmp1(:,1) tmp1(:,2)];
for i = 1:length(tmp1)
    tmp1{i,1} = tmp1{i,1}(1:end-1);
    tmp1{i,2} = tmp1{i,2}(end);
end;

SID_class_label = zeros(length(tmp1), 3);

% 
SID_class_label(1:82, 1) = 1;
SID_class_label(83:141, 1) = 2;
SID_class_label(142:213, 1) = 3;
SID_class_label(214:295, 1) = 4;
SID_class_label(296:388, 1) = 5;
SID_class_label(389:454, 1) = 6;
SID_class_label(455:533, 1) = 7;
SID_class_label(534:623, 1) = 8;
SID_class_label(624:704, 1) = 9;

% cell type
SID_class_label(1:82, 2) = 1;
SID_class_label(83:141, 2) = 1;
SID_class_label(142:213, 2) = 1;
SID_class_label(214:295, 2) = 1;
SID_class_label(296:388, 2) = 2;
SID_class_label(389:454, 2) = 2;
SID_class_label(455:533, 2) = 3;
SID_class_label(534:623, 2) = 3;
SID_class_label(624:704, 2) = 3;

% batch
SID_class_label(1:82, 3) = 2;
SID_class_label(83:141, 3) = 3;
SID_class_label(142:213, 3) = 4;
SID_class_label(214:295, 3) = 4;
SID_class_label(296:388, 3) = 2;
SID_class_label(389:454, 3) = 3;
SID_class_label(455:533, 3) = 3;
SID_class_label(534:623, 3) = 2;
SID_class_label(624:704, 3) = 1;

SID_class_label = [tmp1 num2cell(SID_class_label)];
SID_class_label = [SID_class_label, addlabel([1:454 626:704 536:625 455:535],4)];
%%
save mESC_data data GID SID_class_label cellName;
