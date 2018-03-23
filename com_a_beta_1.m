function [child_ma count_i sum_node sum_penaty ] = com_a_beta_1 (child_ma, count_i,Children_matrix, index_1, node_weight) %迭代计算该结点的所有child的penaty以及数量(不包含结点本身）
%This function is used for computing the total penalty of all the child 
%node for a starting node
sum_penaty = 0;
sum_node = 0;
if isempty(find(Children_matrix(index_1,:)>0)) 
    parent = find(Children_matrix(:,index_1)>0);
    sum_penaty = sum_penaty + (Children_matrix(parent,index_1)*node_weight(index_1,1));
    sum_node = sum_node + 0;
    child_ma(count_i) = index_1;
    count_i = count_i + 1;
else
    child = find(Children_matrix(index_1,:)>0); 
    parent = Children_matrix(:,index_1)>0;
    sum_node = sum_node + sum(Children_matrix(index_1,:));
    [cm, cn] = size(child);    
    child_ma(count_i) = index_1;
    count_i = count_i + 1;

        penaty = (Children_matrix(parent,index_1)*node_weight(index_1,1)); 
    sum_penaty = sum_penaty + penaty;
    for i = 1:1:cn
        [child_ma count_i sum_node_1 sum_penaty_1] = com_a_beta_1(child_ma, count_i,Children_matrix,child(i), node_weight); 
        sum_node = sum_node + sum_node_1;
        sum_penaty = sum_penaty + sum_penaty_1;
    end
end