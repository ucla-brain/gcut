function [ recon_matrix ] = topo_recon( neurite_matrix, leaf_node, raw_matrix, recon_distance )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%recon_matrix_p = neurite_matrix_p;
recon_matrix = neurite_matrix;
% dist_1 = raw_matrix(neurite_matrix_p(:,1), 3:5);
% dist_2 = raw_matrix(neurite_matrix_p(:,2), 3:5);
% dist = pdist2(dist_1, dist_2);
% dist(eye(length(neurite_matrix_p(:, 1))) == 1) = inf;
% [r_ind, c_ind] = find(dist < 10);
% for dd = 1:1:length(r_ind)
%     
%     recon_matrix_p(r_ind(dd), 1) = neurite_matrix_p(c_ind(dd), 2);
%     recon_matrix(r_ind(dd), 1) = neurite_matrix_p(c_ind(dd), 2);
%     
% end
dist_1 = raw_matrix(leaf_node, 3:5);
dist = pdist2(dist_1, dist_1);
dist(eye(length(leaf_node)) == 1) = inf;
dist = triu(dist);
dist(dist < recon_distance & dist > 0) = 1;
dist(dist~=1) = 0;
comp = find_conn_comp(dist);
%leaf_pair = [r_ind, c_ind];
for dd = 1:1:length(comp)
    
    %recon_matrix_p(recon_matrix_p == leaf_node(r_ind(dd))) = leaf_node(c_ind(dd));
    for dd_1 = 2:1:length(comp{dd})
        
        recon_matrix(recon_matrix == leaf_node(comp{dd}(dd_1))) = leaf_node(comp{dd}(1));
    
    end

end

