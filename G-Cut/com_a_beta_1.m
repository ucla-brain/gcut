function [child_ma count_i sum_node sum_penaty ] = com_a_beta_1 (child_ma, count_i,Children_matrix, index_1, node_weight) %迭代计算该结点的所有child的penaty以及数量(不包含结点本身）
%This function is used for computing the total penalty of all the child 
%node for a starting node
sum_penaty = 0;
sum_node = 0;
if isempty(find(Children_matrix(index_1,:)>0)) %如果该结点为叶结点，则迭代结束，返回
    parent = find(Children_matrix(:,index_1)>0);
    sum_penaty = sum_penaty + (Children_matrix(parent,index_1)*node_weight(index_1,1));
    sum_node = sum_node + 0;
    child_ma(count_i) = index_1;
    count_i = count_i + 1;
else
    child = find(Children_matrix(index_1,:)>0); %如果为分支结点，则找到其相应的所有子结点
    parent = Children_matrix(:,index_1)>0;
%    sum_fitness = sum_fitness + angle_ma(index_1); %该分支结点的矢量角
    sum_node = sum_node + sum(Children_matrix(index_1,:));
    [cm, cn] = size(child);    
    child_ma(count_i) = index_1;
    count_i = count_i + 1;
    %if isempty(parent)
        %penaty = 0;
   % else
        penaty = (Children_matrix(parent,index_1)*node_weight(index_1,1)); %该树枝的惩罚度为适应度(父节点到该结点）*该树枝的平均直径
    %end
    sum_penaty = sum_penaty + penaty;
    for i = 1:1:cn
        [child_ma count_i sum_node_1 sum_penaty_1] = com_a_beta_1(child_ma, count_i,Children_matrix,child(i), node_weight); %迭代地计算该分支结点的所有子节点的角度和数量总和
        sum_node = sum_node + sum_node_1;
        sum_penaty = sum_penaty + sum_penaty_1;
    end
end