clear all;
clc;
path_1 = '.\Sample\data\'; %the path of Sample data
path_2 = '.\Sample\result\'; %the path to store the result
soma_index = 0; 
%soma location format. 
%0 - sequence number of soma in SWC file.
%1 - x, y, z location of soma
Neuron_split_sparse_beta(path_1, path_2, soma_index); 