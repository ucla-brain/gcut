clear all;
clc;

if ismac || isunix
    path_1 = './demo_data/'; %the path of Sample data
    path_2 = './demo_result/'; %the path to store the result
elseif ispc
    path_1 = '.\demo_data\'; %the path of Sample data
    path_2 = '.\demo_result\'; %the path to store the result
else
    disp('Platform not supported')
end


soma_index = 0; 
%soma location format. 
%0 - sequence number of soma in SWC file.
%1 - x, y, z location of soma
Neuron_split_sparse_beta(path_1, path_2, soma_index, 's10', '', [], [], []); 