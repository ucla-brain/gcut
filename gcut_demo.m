clear all;
clc;

[gcut_dir, script_name, script_ext] = fileparts(mfilename('fullpath'));

if ismac || isunix
    path_1 = [gcut_dir '/demo_data/'];       %the path of Sample data
    path_2 = [gcut_dir '/demo_result/'];    %the path to store the result
elseif ispc
    path_1 = [gcut_dir '\demo_data\'];       %the path of Sample data
    path_2 = [gcut_dir '\demo_result\'];    %the path to store the result
else
    disp('Platform not supported')
end

soma_index = 0; 
%soma location format. 
%0 - sequence number of soma in SWC file.
%1 - x, y, z location of soma
if exist(path_1, 'dir') == 7 && exist(path_2, 'dir') == 7
    Neuron_split_sparse_beta(path_1, path_2, soma_index, 's10', '', [], [], []); 
else
    disp('sample data or sample output directories not found')
end    