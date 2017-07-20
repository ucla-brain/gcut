function [fitness_matrix com_order_matrix connection_matrix soma_connect] = fit_compute_beta_2(fitness_matrix, neurite_matrix, location_matrix, sort_neurite_matrix, branch_order_index, in_neuron_set, soma_set, soma_index, connection_matrix, soma_connect, neurite_weight)
%This function is used to compute the fitness of each neurite for a give
%soma and construct a weighted DAG.

[m1 n1] = size(fitness_matrix); %the size of fitness matrix
com_order_matrix = zeros(m1,n1); %the path orientation matrix
no_assin_ma = (find(fitness_matrix == -inf));
mm = length(no_assin_ma);
for sub_i = 1:1:mm 
        i = no_assin_ma(sub_i);
        con_c = ismember(neurite_matrix(i,:),in_neuron_set); 
        con_count = con_c==1; 
        con_set = neurite_matrix(i,con_count); 
        
        if isempty(con_set)==1 
            fitness_matrix(i) = -inf;
 
        elseif ~isempty(con_set)&&sort_neurite_matrix(i)~=0; 
            fitness_matrix(i) = inf;
            soma_connect(soma_index,sort_neurite_matrix(i))=1;

        else
            if length(con_set)>1
                disp(con_set);
            end
            [r c] = find(in_neuron_set==con_set); 
            index_neu = find_vector_beta(neurite_matrix, in_neuron_set(r,:)); 
            if  con_set == min(neurite_matrix(i,:))  %the neurite orientation is from right to left
                pro_ma = neurite_matrix(i,:);
                pro_ma(isnan(pro_ma))=[];  
                [m_2 n_2] = size(pro_ma);
                total_an = 0;
                for j = 1:1:n_2-1
                    neu_ve_1 = (location_matrix(neurite_matrix(i,j),2)-location_matrix(neurite_matrix(i,j+1),2));
                    neu_ve_2 = (location_matrix(neurite_matrix(i,j),3)-location_matrix(neurite_matrix(i,j+1),3));
                    neu_ve_3 = (location_matrix(neurite_matrix(i,j),4)-location_matrix(neurite_matrix(i,j+1),4));
                    std_ve_1 = (location_matrix(neurite_matrix(i,j+1),2)-location_matrix(soma_set(soma_index),2)); 
                    std_ve_2 = (location_matrix(neurite_matrix(i,j+1),3)-location_matrix(soma_set(soma_index),3));
                    std_ve_3 = (location_matrix(neurite_matrix(i,j+1),4)-location_matrix(soma_set(soma_index),4));
                    com_an = acos(((neu_ve_1*std_ve_1 + neu_ve_2*std_ve_2 + neu_ve_3*std_ve_3)/...
                            (sqrt(neu_ve_1^2 + neu_ve_2^2 + neu_ve_3^2)*sqrt(std_ve_1^2 + std_ve_2^2 + std_ve_3^2))));
                        if isnan(com_an)
                            com_an = 0;
                        end
                    total_an = total_an + com_an;
                end
                av_an = total_an/(n_2-1);
                neurite_order = branch_order_index(con_set); 
                com_order_matrix(i) = 1;    
                fitness_matrix(i) = Gassi_com_fit_beta(av_an,neurite_order); %compute the fitness of the neurite for a given soma
                if fitness_matrix(i)==inf
                    disp('inf here!');
                    disp(strcat(num2str(av_an),',',num2str(neurite_order)));
                end
                for i_1 = 1:1: length(index_neu)
                    connection_matrix(index_neu(i_1),i) = fitness_matrix(i)*2*(1+(neurite_weight(index_neu(i_1),2)-neurite_weight(i,2))/(neurite_weight(index_neu(i_1),2)+neurite_weight(i,2))); %从胞体出发，与胞体集合相连的该树枝的适应度定义矩阵中由胞体边到该树枝的值
                end
            elseif con_set == max(neurite_matrix(i,:)) %from left to right
                pro_ma = neurite_matrix(i,:);
                pro_ma(isnan(pro_ma))=[];  
                [m_2 n_2] = size(pro_ma);
                total_an = 0;
                for j = n_2-1:-1:1
                    neu_ve_1 = (location_matrix(neurite_matrix(i,j+1),2)-location_matrix(neurite_matrix(i,j),2));
                    neu_ve_2 = (location_matrix(neurite_matrix(i,j+1),3)-location_matrix(neurite_matrix(i,j),3));
                    neu_ve_3 = (location_matrix(neurite_matrix(i,j+1),4)-location_matrix(neurite_matrix(i,j),4));
                    std_ve_1 = (location_matrix(neurite_matrix(i,j+1),2)-location_matrix(soma_set(soma_index),2)); 
                    std_ve_2 = (location_matrix(neurite_matrix(i,j+1),3)-location_matrix(soma_set(soma_index),3));
                    std_ve_3 = (location_matrix(neurite_matrix(i,j+1),4)-location_matrix(soma_set(soma_index),4));
                    com_an = acos(((neu_ve_1*std_ve_1 + neu_ve_2*std_ve_2 + neu_ve_3*std_ve_3)/...
                            (sqrt(neu_ve_1^2 + neu_ve_2^2 + neu_ve_3^2)*sqrt(std_ve_1^2 + std_ve_2^2 + std_ve_3^2))));
                        if isnan(com_an)
                            com_an = 0;
                        end
                    total_an = total_an + com_an;
                end
                av_an = total_an/(n_2-1);
                neurite_order = branch_order_index(con_set);
                com_order_matrix(i) = 0;               
                fitness_matrix(i) = Gassi_com_fit_beta(av_an,neurite_order);
                if fitness_matrix(i)==inf
                    disp('inf here!!');
                    disp(strcat(num2str(av_an),',',num2str(neurite_order)));
                end
                for i_1 = 1:1: length(index_neu)
                    connection_matrix(index_neu(i_1),i) = fitness_matrix(i)*2*(1+(neurite_weight(index_neu(i_1),2)-neurite_weight(i,2))/(neurite_weight(index_neu(i_1),2)+neurite_weight(i,2))); %考虑加入树枝直径的变化
                end
            end
        end
end