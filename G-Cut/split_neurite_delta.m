function [neurite_matrix] = split_neurite_delta(branch_matrix, leaf_matrix, parent_list, child_list)
n_id = 1;
if isempty(find(branch_matrix ==1))
    branch_matrix = [1;branch_matrix];
end
while isempty(leaf_matrix)~=1
    node_id = 1;
    index_1 = leaf_matrix(1);
    index_2 = index_1;
    while(isempty(find(branch_matrix == index_2)==1))
        index_1 = index_2;
        Neu_matrix(n_id,node_id) = index_1;
        index_2 = parent_list(parent_list(:,1)==index_1,2);
        child_list((child_list(:,2) == index_1),:) = [];
        node_id = node_id + 1;
    end
    Neu_matrix(n_id,node_id) = index_2;
    child_list((child_list(:,2) == index_1 ),:) = [];
    [l_m l_n] = size(leaf_matrix);
    if isempty(find(child_list(:,1) == index_2))&&(index_2~=branch_matrix(1))
        leaf_matrix = [leaf_matrix; index_2];
        branch_matrix((branch_matrix == index_2)) = [];
    end
    leaf_matrix(1) = [];
    n_id = n_id + 1;
end
neurite_matrix = Neu_matrix;