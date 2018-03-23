function [ neu_matrix, soma_location ] = find_soma_location( neu_ma, location_path, location_index )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    A_1 = neu_ma;
    if location_index == 0
        index_x_1 = importdata(location_path);
        index_x_1 = index_x_1';
    elseif location_index == 1
        soma_location = importdata(location_path);
        [lo_m lo_n] = size(soma_location);
        soma_index = zeros(lo_m,1);
        ma_location = A_1(:,3:5);
        for lo_ind = 1:1:lo_m
            location_1 = soma_location(lo_ind,:);
            dist_ma = pdist2(ma_location,location_1);
            candi_soma = A_1(dist_ma < (min(dist_ma) + 10), 1);
            [~, soma_ind] = max(A_1(candi_soma, 6));
             soma_ind_1 = A_1(candi_soma(soma_ind), 1);
             soma_index(lo_ind) = soma_ind_1;
             disp(num2str(lo_ind));
        end
        index_x_1 = soma_index';
        A_1(index_x_1, 6) = 10;
    end
    neu_matrix = A_1;
    soma_location = index_x_1;


end

