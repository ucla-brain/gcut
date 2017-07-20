function [new_soma_set, merge_matrix ] = soma_merge_beta(raw_matrix, soma_set)
[m n] = size(soma_set);
[m1 n1] = size(raw_matrix);
dist_matrix = raw_matrix(:,3:5);
del_ma = [];
% con_ma = sparse(m1,m1);
% for i = m1:-1:1
%     if raw_matrix(i,7) == -1
%         continue;
%     else
%         index_1 = raw_matrix(i,7);
%         con_ma(i,index_1) = 1;
%         con_ma(index_1,i) = 1;
%     end
% end

for i = 1:1:n
    or_ind = soma_set(i);
    %dist_soma = 0;
    radius = 1.5 * raw_matrix(or_ind,6);
    or_ma = dist_matrix(or_ind,:);
    dist_ma = (dist_matrix - (or_ma(ones(m1,1),:)));
    dist_ma = sqrt(dist_ma(:,1).^2 + dist_ma(:,2).^2 + dist_ma(:,3).^2);
    candi_node = find(dist_ma<radius);
    candi_node(find(candi_node==or_ind)) = [];
    fault_node = [];
    for i_1 = 1:1:length(candi_node)
        j = candi_node(i_1);
        fault_flag = 1;
        if j> or_ind
            pa_ind = raw_matrix(j,7);
            while(norm(dist_matrix(pa_ind,:)-dist_matrix(or_ind,:))<radius)
                if isequal(pa_ind,or_ind)
                    fault_flag = 0;
                    break;
                else
                    pa_ind = raw_matrix(pa_ind,7);
                end
            end
        else
            pa_ind = raw_matrix(or_ind,7);
            while(norm(dist_matrix(pa_ind,:)-dist_matrix(or_ind,:))<radius)
                if isequal(pa_ind,j)
                    fault_flag = 0;
                    break;
                else
                    pa_ind = raw_matrix(pa_ind,7);
                end
            end
        end
        if ~isequal(fault_flag,0)
            fault_node = [fault_node;i_1];
        end
    end
    candi_node(fault_node) = [];
%         path = dijkstra(con_ma,j,or_ind);
%         path = path';
%         dist_ma_1 = (dist_matrix(path,:) - (or_ma(ones(length(path),1),:)));
%         dist_ma_1 = sqrt(dist_ma_1(:,1).^2 + dist_ma_1(:,2).^2 + dist_ma_1(:,3).^2);
%         fault_node = find(dist_ma_1>radius);
            for i_1 = 1:1:length(candi_node)
                j = candi_node(i_1);
                in_1 = raw_matrix(:,7) == j;
%                 raw_matrix(in_1,7) = or_ind;
%                 raw_matrix(or_ind,6) = max(raw_matrix(or_ind,6),raw_matrix(j,6));
%                 raw_matrix(or_ind,7) = min(raw_matrix(or_ind,7),raw_matrix(j,7));
                if j>or_ind
                	raw_matrix(in_1,7) = or_ind;
                    raw_matrix(or_ind,6) = max(raw_matrix(or_ind,6),raw_matrix(j,6));
                    raw_matrix(or_ind,7) = min(raw_matrix(or_ind,7),raw_matrix(j,7));
                    del_ma = [del_ma; raw_matrix(j,1)];
                else 
                    raw_matrix(or_ind,6) = max(raw_matrix(or_ind,6),raw_matrix(j,6));
                    raw_matrix(or_ind,7) = min(raw_matrix(or_ind,7),raw_matrix(j,7));
                    raw_matrix(j,2:7) = raw_matrix(or_ind,2:7);
                    del_ma = [del_ma; or_ind];
                    or_ind = j;
                end
                
            end
    soma_set(i) = or_ind;
end
raw_matrix(del_ma,:) = [];
[m2 n2] = size(raw_matrix);
for i_1 = 1:1:n
    so_ind = find(raw_matrix(:,1)==soma_set(i_1));
    soma_set(i_1) = so_ind;
end
for i =1:1:m2
    in_2 = raw_matrix(:,7)==raw_matrix(i,1);
    raw_matrix(in_2,7) = i;
    raw_matrix(i,1) = i;
end
merge_matrix = raw_matrix;
new_soma_set = soma_set;