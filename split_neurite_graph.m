function [ neurite_matrix ] = split_neurite_graph( cluster_graph )
%SPLIT_NEURITE_GRAPH Summary of this function goes here
%   Detailed explanation goes here

    %data_1 = cluster_graph;
    
    %edge_ma = cluster_graph(:, 1:3);
    
    %[r_1, ~] = find(edge_ma == -1);
    
    %edge_ma(r_1,:) = [];
    
    %graph_1 = edgeL2adj(edge_ma);
    
    graph_1 = cluster_graph;
    
    graph_1 = symmetrize(graph_1);
    
    neigh_n = sum(graph_1, 2);
    
    branch_node = find(neigh_n > 2);
    
    con_node = find(neigh_n == 2);
    
    leaf_node = find(neigh_n == 1);
    
    ind_ma = [ (1:length(neigh_n))', neigh_n];
    
%    neurite_matrix = [];

    i_c = 1;
    
    branch_ma = branch_node;
    
    while(~isempty(branch_ma))
        
        s_node = branch_ma(1);
        
        neurite_1 = s_node;
        
        neigh_1 = find(graph_1(s_node, :) == 1);
        
        %neigh_1(ind_ma(neigh_1, 2) == 0) = [];
        
        neigh_1(neigh_1 == s_node) = [];
        
        neigh_1 = neigh_1(1);
        
        graph_1(s_node, neigh_1) = 0;
        
        graph_1(neigh_1, s_node) = 0;
        
        neurite_1 = [neurite_1, neigh_1];
        
        while (isempty(find(branch_node == neigh_1, 1)) && isempty(find(leaf_node == neigh_1, 1)));
            
            neigh_2 = find(graph_1(neigh_1, :) == 1);
            
            %neigh_2(ind_ma(neigh_2, 2) == 0) = [];
            
            neigh_2(neigh_2 == s_node) = [];
            
            s_node = neigh_1;
            
%             if (find(neigh_2 == 4115))
%                 
%                 aaaa = 1;
%                 
%             end
            
%             disp(num2str(neigh_2));
            
            neigh_1 = neigh_2(1);
            
            graph_1(s_node, neigh_1) = 0;
        
            graph_1(neigh_1, s_node) = 0;
            
            neurite_1 = [neurite_1, neigh_1];
            
        end
        
%         if length(neurite_1) > 2
%             
%             ind_ma(neurite_1(2:end-1), 2) = ind_ma(neurite_1(2:end-1), 2) - 2;
%             
%         end
%         
%         ind_ma(neurite_1([1, end]), 2) = ind_ma(neurite_1([1, end]), 2) - 1;

        %del_branch = ind_ma(branch_ma, 2) == 0;
        
        del_branch = find(sum(graph_1(branch_ma, :), 2) == 0);
        
        branch_ma(del_branch) = [];
        
        for j_c = 1:1:length(neurite_1)
            
            neurite_matrix(i_c, j_c) = neurite_1(j_c);
            
        end
        
        i_c = i_c + 1;
        
    end
        
        


end

