function [recon_Neuron ] = neurite_recon_ver_2(Neuron_matrix, raw_matrix, node_number, soma_set)
soma_1 = [];
[m n] = size(Neuron_matrix);
recon_Neuron = cell(m,n);
for i = 1:1:n
    recon_ma = Neuron_matrix{i};
    con_ma = zeros(node_number,7);
    [m_1 n_1] = size(recon_ma);
    con_ma(1,1) = soma_set(i);
    con_ma(1,7) = -1;
    con_ma(1,2:6) = raw_matrix(soma_set(i),2:6);
    branch_set = soma_set(i);
    index_i = 2;
    loop_error_flag = 0;
    while(isempty(recon_ma)~=1 && loop_error_flag <10)
        branch_set_1 = [];
        delet_ma = [];
        branch_delt =[];
        for re_i = 1:1:length(recon_ma(:,1))
            member_index = ismember(branch_set,recon_ma(re_i,:));
            member_set = find(member_index ==1);
            if isempty(member_set)==1
                continue;
            else
                po_index = find(recon_ma(re_i,:)==branch_set(member_set));
                ma_1 = recon_ma(re_i,:);
                ma_1(isnan(ma_1))=[];
                if po_index ==1
                    for re_j = 2:1:length(ma_1)
                        con_ma(index_i,1) = recon_ma(re_i,re_j);
                        con_ma(index_i,2:6) = raw_matrix(recon_ma(re_i,re_j),2:6);
                        ind = find(con_ma(:,1) == recon_ma(re_i,re_j-1));
                        if isempty(ind)
                            disp(strcat('error!_',num2str(recon_ma(re_i,re_j-1))));
                        end
                        con_ma(index_i,7) = ind;
                        index_i = index_i + 1;
%                         con_ma(recon_ma(re_i,re_j),1) = recon_ma(re_i,re_j);
%                         con_ma(recon_ma(re_i,re_j),2) = recon_ma(re_i,re_j-1);
                    end
                   % branch_delt = [branch_delt;find(branch_set==recon_ma(re_i,1))];
                    branch_set_1 = union(branch_set_1,ma_1);
                elseif po_index == length(ma_1)
                    for re_j = length(ma_1)-1:-1:1
                        con_ma(index_i,1) = recon_ma(re_i,re_j);
                        con_ma(index_i,2:6) = raw_matrix(recon_ma(re_i,re_j),2:6);
                        ind = find(con_ma(:,1) == recon_ma(re_i,re_j+1));
                        if isempty(ind)
                            disp(strcat('error!_',num2str(recon_ma(re_i,re_j+1))));
                        end
                        con_ma(index_i,7) = ind;
                        index_i = index_i + 1;
%                         con_ma(recon_ma(re_i,re_j),1) = recon_ma(re_i,re_j);
%                         con_ma(recon_ma(re_i,re_j),2) = recon_ma(re_i,re_j+1);
                    end
                  %  branch_set(branch_set==recon_ma(re_i,length(ma_1))) = [];
                    branch_set_1 = union(branch_set_1,ma_1);
                else
                    for re_j = po_index+1:1:length(ma_1)
                        con_ma(index_i,1) = recon_ma(re_i,re_j);
                        con_ma(index_i,2:6) = raw_matrix(recon_ma(re_i,re_j),2:6);
                        ind = find(con_ma(:,1) == recon_ma(re_i,re_j-1));
                        if isempty(ind)
                            disp(strcat('error!_',num2str(recon_ma(re_i,re_j-1))));
                        end
                        con_ma(index_i,7) = ind;
                        index_i = index_i + 1;
%                         con_ma(recon_ma(re_i,re_j),1) = recon_ma(re_i,re_j);
%                         con_ma(recon_ma(re_i,re_j),2) = recon_ma(re_i,re_j-1);
                    end
                    for re_j = po_index-1:-1:1
                        con_ma(index_i,1) = recon_ma(re_i,re_j);
                        con_ma(index_i,2:6) = raw_matrix(recon_ma(re_i,re_j),2:6);
                        ind = find(con_ma(:,1) == recon_ma(re_i,re_j+1));
                        if isempty(ind)
                            disp(strcat('error!_',num2str(recon_ma(re_i,re_j+1))));
                        end
                        con_ma(index_i,7) = ind;
                        index_i = index_i + 1;
%                         con_ma(recon_ma(re_i,re_j),1) = recon_ma(re_i,re_j);
%                         con_ma(recon_ma(re_i,re_j),2) = recon_ma(re_i,re_j+1);
                    end
                 %   branch_set(branch_set==recon_ma(re_i,po_index)) = [];
                    branch_set_1 = union(branch_set_1,ma_1);
                end
                delet_ma = [delet_ma;re_i];
            end
        end
        if isempty(delet_ma)
            loop_error_flag = loop_error_flag + 1;
        end
        branch_set = branch_set_1;
        branch_set = union(branch_set,branch_set);
        recon_ma(delet_ma,:)=[];
        
    end
    ma_1 = [1:node_number]';
    conec_ma = con_ma;
    conec_ma(:,1) = ma_1;
    conec_ma(find(conec_ma(:,7)==0),:)=[];
    conec_ma(find(conec_ma(:,7)==-1),2) = 1;
    recon_Neuron{i} = conec_ma;
end