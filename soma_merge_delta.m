function [new_soma_set, merge_matrix ] = soma_merge_delta(raw_matrix, soma_set)
%This function is used to merge the nodes which are close to the soma node
%in spacial location.
[m n] = size(soma_set);
[m1 n1] = size(raw_matrix);
dist_matrix = raw_matrix(:,3:5);
del_ma = [];
parent_list = [raw_matrix(:, 1),raw_matrix(:, 7)];
child_list = [raw_matrix(:, 7), raw_matrix(:, 1)];

for i = 1:1:n
    or_ind = soma_set(i);
    radius = 1.5 * raw_matrix(or_ind,6);
    node_1 = or_ind;
    candi_par = [];
    candi_child = [];
    
    while(parent_list(node_1, 2) ~= -1)
        par_node = parent_list(node_1, 2);
        dist = sqrt(sum((raw_matrix(par_node, 3:5) - raw_matrix(or_ind, 3:5)).^2));
        ra_1 = max(raw_matrix(par_node, 6) + raw_matrix(or_ind, 6), radius);
        if dist<ra_1
            candi_par = [candi_par; par_node];
            node_1 = par_node;    
        else
            break;
        end
    end
    node_1 = raw_matrix(node_1, 7);
    
    node_2 = or_ind;
    while(~isempty(child_list(ismember(child_list(:, 1), node_2), 2)))
        ind_1 = ismember(child_list(:, 1), node_2);
        child_node = child_list(ind_1, 2);
        dist = pdist2(raw_matrix(child_node, 3:5), raw_matrix(or_ind, 3:5));
        child_dist = max([repmat(radius, length(child_node), 1), raw_matrix(child_node, 6) + repmat(raw_matrix(or_ind, 6), length(child_node), 1)], [], 2);
        dist_1 = child_dist - dist;
        candi_1 = find(dist_1 > 0);
        if ~isempty(candi_1)
            candi_child = [candi_child; child_node(candi_1)];
            node_2 = child_node(candi_1);
            raw_matrix(ismember(raw_matrix(:, 7), child_node(candi_1)), 7) = or_ind;
        else
            break;
        end
    end
    node_2 = child_list(ismember(child_list(:, 1), node_2), 2);
    child_2 = child_list(ismember(child_list(:, 1), candi_par), 2);
    raw_matrix(child_2, 7) = node_1;
    raw_matrix(node_2, 7) = or_ind;
    if ~isempty(find((raw_matrix(:,1) - raw_matrix(:, 7)) == 0))
        aaaa = 1;
    end
    del_ma = [del_ma; candi_par; candi_child];
    disp(num2str(i));
end

raw_matrix(del_ma,:) = [];
for i_1 = 1:1:n
    so_ind = find(raw_matrix(:,1)==soma_set(i_1));
    soma_set(i_1) = so_ind;
end
matrix_1 = tree_resort(raw_matrix);
merge_matrix = matrix_1;
new_soma_set = soma_set;