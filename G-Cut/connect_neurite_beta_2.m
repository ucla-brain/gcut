function [Neuron, re_Neuron] = connect_neurite_beta_2(raw_matrix, neurite_matrix, Set_index, node_num) 
%This function compute the fitness of each neurite for each soma. Then it
%use topological constraints to determine the assignment of each neurite.
%After that, this function assign the neurite to the certain soma 
%according to the neurite assignment. Then construct individual neuron
%by iteratively joining the neurite to the given soma.

[raw_i raw_j] = size(raw_matrix);
location_matrix = zeros(raw_i,4);
location_matrix(:,1) = raw_matrix(:,1);
location_matrix(:,2:4) = raw_matrix(:,3:5);
radia_matrix = raw_matrix(:,6);
[m n] = size(neurite_matrix); 
soma_set = find(Set_index==1); 
[m1 n1] = size(soma_set);
Neuron = cell(n1,m1);
orin_soma_node = cell(n1,m1);
Neu_con_ma = cell(n1,m1); 
fitness_ma = zeros(m,m1); 
sort_neurite_ma = zeros(m,1); 
branch_order_index = zeros(node_num,m1); 
[node_weight ] = neu_length(neurite_matrix, location_matrix, radia_matrix); %
branch_order = 1; 
for soma_i = 1:1:m1    %initialize the graph
    for neu_i = 1:1:m
        ne_ind = (find(neurite_matrix(neu_i,:)==soma_set(soma_i)));
        ne_num = neurite_matrix(neu_i,ne_ind);
            if ~isempty(ne_ind)
                if ne_num == max(neurite_matrix(neu_i,:))
                    branch_order_index(neurite_matrix(neu_i,ne_ind),soma_i) = branch_order;
                    branch_order_index(min(neurite_matrix(neu_i,:)),soma_i) = branch_order_index(neurite_matrix(neu_i,ne_ind),soma_i) + 1;
                    sort_neurite_ma(neu_i) = soma_i; 
                elseif ne_num == min(neurite_matrix(neu_i,:))
                    branch_order_index(neurite_matrix(neu_i,ne_ind),soma_i) = branch_order;
                    branch_order_index(max(neurite_matrix(neu_i,:)),soma_i) = branch_order_index(neurite_matrix(neu_i,ne_ind),soma_i) + 1;
                    sort_neurite_ma(neu_i) = soma_i;
                else
                    branch_order_index(neurite_matrix(neu_i,ne_ind),soma_i) = branch_order;
                    branch_order_index(min(neurite_matrix(neu_i,:)),soma_i) = branch_order_index(neurite_matrix(neu_i,ne_ind),soma_i) + 1;
                    branch_order_index(max(neurite_matrix(neu_i,:)),soma_i) = branch_order_index(neurite_matrix(neu_i,ne_ind),soma_i) + 1;
                    sort_neurite_ma(neu_i) = soma_i;
                end
                sort_neurite_ma(neu_i) = soma_i;
                set_neu_ma = neurite_matrix(neu_i,:);
                set_neu_ma(isnan(set_neu_ma))=[];
                Set_index(set_neu_ma) = 1;
                Neuron{soma_i} = [Neuron{soma_i};neurite_matrix(neu_i,:)]; 
            end
    end
end

soma_connect = zeros(m1,m1);
for or_i = 1:1:m1
    for or_j = 1:1:length(Neuron{or_i}(:,1))
    orin_soma_node{or_i}(or_j) = find_vector_beta(neurite_matrix,Neuron{or_i}(or_j,:));
    end
end



for so_i = 1:1:m1 %We contruct a weighted DAG for each soma
    sub_neuron = Neuron{so_i};
    sub_fitness_ma = ones(m,1).*-inf;
    sub_con_ma = zeros(m,m);
    sub_sort_neurite_ma = sort_neurite_ma;
    orin_core = 0;
    add_core = -1;
    while(add_core ~= orin_core) 
        add_core = orin_core;
        [sub_fitness_ma com_order_ma sub_con_ma soma_connect] = fit_compute_beta_2(sub_fitness_ma, neurite_matrix,location_matrix, sub_sort_neurite_ma, branch_order_index(:,so_i), sub_neuron, soma_set, so_i, sub_con_ma, soma_connect, node_weight); %Compute the fitness of each neurite for a certain soma.
        fit_ass_neu = find(sub_fitness_ma~=-inf & sub_sort_neurite_ma == 0);
        if (isempty(fit_ass_neu))
            orin_core = add_core;
        else
            sub_sort_neurite_ma(fit_ass_neu) = so_i;
            sub_neuron = [sub_neuron; neurite_matrix(fit_ass_neu,:)];
            order_1 = com_order_ma(fit_ass_neu)==0;
            order_2 = com_order_ma(fit_ass_neu)==1;
            branch_order_index(min(neurite_matrix(fit_ass_neu(order_1),:),[],2),so_i) = branch_order_index(max(neurite_matrix(fit_ass_neu(order_1),:),[],2),so_i) + 1;
            branch_order_index(max(neurite_matrix(fit_ass_neu(order_2),:),[],2),so_i) = branch_order_index(min(neurite_matrix(fit_ass_neu(order_2),:),[],2),so_i) + 1;
            orin_core = orin_core + length(fit_ass_neu);
        end
    end
    Neu_con_ma{so_i} = sub_con_ma;
    Neuron{so_i} = sub_neuron;
end

disp('connection graph construction complete!');
toc;

[Neuron] = path_extract(Neu_con_ma, Neuron, neurite_matrix, soma_set, soma_connect, orin_soma_node, node_weight); 

disp('path assign complete!');
toc;
[re_Neuron] = neurite_recon_ver_2(Neuron, raw_matrix, node_num, soma_set);