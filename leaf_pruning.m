function [prun_matrix, prun_neurite ] = leaf_pruning( raw_matrix, neurite_ma, leaf_node, prun_c, prun_mode )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

     m = length(leaf_node);
%     [~, n] = size(neurite_ma);
%     leaf_neurite = size(m, n);
    leaf_ind = ismember(neurite_ma(:, 1), leaf_node);
    leaf_ind = find(leaf_ind == 1);
    leaf_neurite = neurite_ma(leaf_ind, :);
    leaf_length = zeros(m, 1);
    if prun_mode == 1
        %----------------------
        for i = 1:1:length(leaf_neurite(:, 1))
            neu_1 = leaf_neurite(i, :);
            neu_1(neu_1 == 0) = [];
            eu = 0;
            for j = 1:1:length(neu_1) - 1
                eu = eu + sqrt(sum((raw_matrix(neu_1(j + 1), 3:5) - raw_matrix(neu_1(j), 3:5)).^2));
            end
            leaf_length(i) = eu;
        end
        length_1 = sort(leaf_length);
        thre_1 = length_1( ceil(length(length_1)*prun_c)); %根据树枝长度来进行剪枝
        %------------------------------
    elseif prun_mode == 2
        for i = 1:1:length(leaf_neurite(:, 1))
            neu_1 = leaf_neurite(i, :);
            neu_1(neu_1 == 0) = [];
            total_an = 0;
            for j = 1:1:length(neu_1) - 1
                neu_an_1 = raw_matrix(neu_1(j), 3:5) - raw_matrix(neu_1(j+1), 3:5);
                std_an_1 = raw_matrix(neu_1(j+1), 3:5) - raw_matrix(1, 3:5);
                an_1 = acos(sum(neu_an_1.* std_an_1)/(sqrt(sum(neu_an_1.^2))*sqrt(sum(std_an_1.^2))));
                if isnan(an_1)
                    an_1 = pi;
                end
                total_an = total_an + an_1;
            end
            av_an = total_an/(length(neu_1) - 1);
            leaf_length(i) = av_an;
        end
        length_1 = sort(leaf_length);
       % thre_1 = length_1( ceil(length(length_1)*prun_c)); %根据angle来进行剪枝（考虑转化为GOF）
       thre_1 = prun_c;
    end
    prun_ind = leaf_length > thre_1;
%    prun_ma = neurite_ma(leaf_ind(prun_ind), :);
    prun_neurite = neurite_ma;
    prun_neurite(leaf_ind(prun_ind), :) = [];
    rest_node = prun_neurite(prun_neurite ~=0);
    rest_node = union(rest_node, rest_node);
    prun_matrix = raw_matrix(rest_node, :);
    %prun_matrix(prun_node, :) = [];
end

