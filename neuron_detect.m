function [ Parent_list, Child_list, branch_nodes, leaf_nodes ] = neuron_detect( raw_matrix )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    A = raw_matrix;
    [m, ~] = size(A);
    Parent_list = ones(m,2);
    Parent_list(:,1) = A(:,1);
    Parent_list(:,2) = A(:,7);
    Child_list = ones(m,2);
    Child_list(:,1) = A(:,7);
    Child_list(:,2) = A(:,1);
     Child_list((Child_list(:,1) == -1),:) = [];
    elem_A = A(:,1);
    [count_elem, node_ind] = hist(Child_list(:,1),elem_A);
    leaf_nodes = node_ind(count_elem == 0);
    branch_nodes = node_ind(count_elem > 1);
end

