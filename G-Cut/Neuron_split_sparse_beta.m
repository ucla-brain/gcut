function [] = Neuron_split_sparse_beta(path_1, path_2, soma_option)

% This script will run the neuron segmentation algorithm in parallel mode. 
% A minimun of 8 CPU cores are assumed available.
% [] = Neuron_split_spare_sever_version(path_1, path_2, soma_option) 
% path_1: input data folder. must contain neuron cluster in swc format, and 
%         corresponding text file with soma information. The swc file and text
%         file should be named [cluster_name].swc and [cluster_name]_soma_ind.txt
%         respectively. Multiple neuron clusters can exist in the input data
%         folder, as long as their corresponding soma information files are
%         present .
% path_2: output folder. SWC file for segemented individual neurons are stored
%         in this folder as [cluster_name]_[number].swc, where number ranges 
%         between 1 and total number of neurons in the cluster.
% soma_option - input soma information file type:
% 0 - soma_index: The node number of somas in SWC file
%                 if node 5 and 23 are somas, the file should be written as:
%                 5
%                 23
% 1 - soma_location: The x,y,z locations of somas
%                    if node 5 has coordinates (1, 2, 3) and node 23 has 
%                    coordinates (4, 5, 6), the file should be written as:
%                    1 2 3
%                    4 5 6

if soma_option ~= 0 && soma_option ~=1
        disp('Please input a valid soma index!');
        return;
end
file_list = dir(fullfile(path_1,'*.swc'));
file_n = length(file_list);

for file_i = 1:1:file_n
    file_name = file_list(file_i).name;
    [path, or_name, ext] = fileparts(file_name);
    A_1 = swc_read(strcat(path_1,file_name));
    index_import = soma_option;
    if index_import == 0
        index_x_1 = importdata(strcat(path_1 ,or_name,'_soma_ind.txt'));
        index_x_1 = index_x_1';
    else
         soma_location = importdata(strcat(path_1,or_name,'_soma_ind.txt'));
        [lo_m lo_n] = size(soma_location);
        soma_index = zeros(lo_m,1);
        ma_location = A_1(:,3:5);
        for lo_ind = 1:1:lo_m
            location_1 = soma_location(lo_ind,:);
            dist_ma = pdist2(ma_location,location_1);
            [~, soma_ind_1] = min(dist_ma);
            soma_index(lo_ind) = soma_ind_1;
        end
        index_x_1 = soma_index';
    end
    %---------------no soma merge-----------------------------------
    A = A_1;
    index_x = index_x_1;
    %----------------------------------------------------------------
    index_x = sort(index_x);
    tic;
    [m,n] = size(A);
%-----------------------------------------------------------
    Parent_list = ones(m,2);
    Parent_list(:,1) = A(:,1);
    Parent_list(:,2) = A(:,7);
    Child_list = ones(m,2);
    Child_list(:,1) = A(:,7);
    Child_list(:,2) = A(:,1);
     Child_list((Child_list(:,1) == -1),:) = [];
    elem_A = A(:,1);
    [count_elem, node_ind] = hist(Child_list(:,1),elem_A);
    leaf_node = node_ind(count_elem == 0);
    branch_node = node_ind(count_elem > 1);
    [neurite_ma] = split_neurite_delta(branch_node,leaf_node,Parent_list, Child_list);
    Set_index = zeros(m,1);
    Set_index(index_x) = 1;
    [era_i era_j] = size(neurite_ma);
    for pre_i =1:1:era_i
        era = find(neurite_ma(pre_i,:)==0);
        neurite_ma(pre_i,era) = NaN;
    end
    disp('split complete');
    toc;
    [in_neuron, rec_matrix] = connect_neurite_beta_2(A,neurite_ma,Set_index,m);
    for i = 1:1:length(rec_matrix)
        save_file_std_swc(strcat(path_2, or_name, '_split_',num2str(i),'.swc'),rec_matrix{i});
    end
    disp(strcat(num2str(file_i),'_complete!'));
end

%matlabpool close;
