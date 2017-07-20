function [] = Neuron_cluster_segmentation(path_1, path_2, soma_lo_ind, options)
%This is a matlab function used to segment neuron cluster contain
%multi-interweaving neurons. The accept file format is SWC with
%reconstructed single neuron cluster. And the soma index file should be
%named as "(neuron_cluster_swc_file)_soma_index.txt) and under the same
%folder of SWC file.
%The function have four inputs:
%path_1: the load folder which have the neuron cluster need to segment. One
%folder can have several SWC files.
%
%path_2: the save folder to save the segmented individual neurons
%
%soma_lo_ind: '0' ---------soma location is represented by index in SWC file 
%               (default setting)
%             '1'----------soma location is represented by real location in
%             raw image stack 
%options: 'merge' ---------merge the neurite which very close to soma. This
%                 could reduce the deviation of stem number caused by 
%                 tracing method 


if nargin <2
    help Neuron_cluster_segmenation
    return;
end

if nargin <3
        soma_lo_ind = 0;
elseif nargin <4

        merge = 0;
else
        merge = 1;
end

if soma_lo_ind ~= 0 && soma_lo_ind ~=1
    disp('invalid soma index option!');
    return;
end

file_list = dir(fullfile(path_1,'*.swc'));
file_n = length(file_list);
for file_i = 1:1:file_n
    file_name = file_list(file_i).name;
    [path, or_name, ext] = fileparts(file_name);
    A_1 = swc_read(strcat(path_1,file_name));
    if soma_lo_ind == 0
        index_x_1 = importdata(strcat(path_1 ,or_name,'_soma_ind.txt'));
        soma_index = index_x_1';
    elseif soma_lo_ind ==1
        index_x_1 = importdata(strcat(path_1 ,or_name,'_soma_ind.txt'));
        soma_index = location_find(index_x_1,A_1);
    end
    if merge == 0
        neu = A_1;
    else
        [soma_index, neu] = soma_merge_beta(A_1,soma_index);
    end
    tic;
    
    [m, n] = size(neu);
    Parent_list = ones(m,2);
    Parent_list(:,1) = neu(:,1);
    Parent_list(:,2) = neu(:,7);
    Child_list = ones(m,2);
    Child_list(:,1) = neu(:,7);
    Child_list(:,2) = neu(:,1);
    Child_list((Child_list(:,1) == -1),:) = [];
    elem_neu = neu(:,1);
    [count_elem, node_ind] = hist(Child_list(:,1),elem_neu);
    leaf_node = node_ind(count_elem == 0);
    branch_node = node_ind(count_elem > 1);
    [neurite_ma] = split_neurite_delta(branch_node,leaf_node,Parent_list, Child_list);
    Set_index = zeros(m,1);

    Set_index(soma_index) = 1;
    [era_i era_j] = size(neurite_ma);
    for pre_i =1:1:era_i
        era = neurite_ma(pre_i,:)==0;
        neurite_ma(pre_i,era) = NaN;
    end
    disp('split complete');
    toc;
    
    [in_neuron, rec_matrix] = connect_neurite_beta_2(neu,neurite_ma,Set_index,m);
end
        