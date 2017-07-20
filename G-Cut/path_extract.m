function [Neuron ] = path_extract(Neu_con_ma, Neuron, neurite_matrix, soma_set, soma_connect, orin_soma_node, node_weight)
%This function is used to find the common path between somas and use the
%topological constraints to find the neurite assignment for each soma.

[m n] = size(soma_set);
soma_connect = soma_connect - eye(m);
logic_Neu_con_ma = Neu_con_ma;
[m2 n2] = size(Neu_con_ma);
for i = 1:1:n2
    logic_Neu_con_ma{i}(logic_Neu_con_ma{i}>0) = 1; 
end
whole_logic_con_ma = logic_Neu_con_ma{1};
for i = 1:1:n2
    whole_logic_con_ma = whole_logic_con_ma | logic_Neu_con_ma{i}; 
end

circle_logic_con = whole_logic_con_ma & whole_logic_con_ma';
whole_com_path = mean(circle_logic_con,2);
whole_path_node = find(whole_com_path >0); %find the common path set
cell_ind = 1;
path_soma_set_cell = cell(1,m*m);
in_com_path_node_cell = cell(1,m*m);
while(~isempty(whole_path_node)) %find each common path 
    com_node = whole_path_node(1);
    in_com_path_node = [];
    path_soma_set = zeros(1,m); 
    for j = 1:1:m
        [d, P] = dijkstra(whole_logic_con_ma, orin_soma_node{j}, com_node );
        if d == inf
            path_soma_set(j) = 0;
        else
            path_soma_set(j) = 1;
            P(ismember(P,orin_soma_node{j}))=[];
            in_com_path_node = union(in_com_path_node, P);
        end
    end
    path_soma_set_cell{cell_ind} = path_soma_set;
    in_com_path_node_cell{cell_ind} = in_com_path_node;
    whole_path_node(ismember(whole_path_node,in_com_path_node))=[];
    cell_ind = cell_ind + 1;
end

id = cellfun('length',in_com_path_node_cell);
in_com_path_node_cell(id==0)=[];
[p_m p_n] = size(in_com_path_node_cell);
path_assin_struct_cell = cell(p_n,1);

for i = 1:1:p_n
    in_path_soma_con = path_soma_set_cell{i};
    in_path_node = in_com_path_node_cell{i};
    [h v] = size(in_path_node);
    if h<v
        in_path_node = in_path_node';
    end
    soma_con_ind = find(in_path_soma_con==1);
    soma_con_num = length(soma_con_ind);
    neuron_con_ma_cell = cell(1,soma_con_num);
    logic_con_ma_cell = cell(1,soma_con_num);
    in_soma_orin_cell = cell(1,soma_con_num);
    Neuron_in = cell(1,soma_con_num);
    for ce_i = 1:1:soma_con_num %store the information of common path
        neuron_con_ma_cell{ce_i} = Neu_con_ma{soma_con_ind(ce_i)};
        in_soma_orin_cell{ce_i} =  orin_soma_node{soma_con_ind(ce_i)};
        Neuron_in{ce_i} =  Neuron{soma_con_ind(ce_i)};
        logic_con_ma_cell{ce_i} = logic_Neu_con_ma{soma_con_ind(ce_i)};
    end
    [Neuron_in in_path_child_cell in_fitness_ma in_connect_ind new_path_node] = neurite_route_sort_delta_3(neuron_con_ma_cell, logic_con_ma_cell, Neuron_in, neurite_matrix, soma_con_num, in_soma_orin_cell, in_path_node, node_weight );
    for ce_i = 1:1:soma_con_num 
        Neuron{soma_con_ind(ce_i)} = Neuron_in{ce_i};
    end
    %defint the structure to store the processed common path information
    path_struct.fitness = in_fitness_ma;
    path_struct.connect_ind = in_connect_ind;
    path_struct.path_node = new_path_node;
    path_struct.connect_soma = soma_con_ind;
    path_struct.connect_soma_num = length(soma_con_ind);
    path_struct.path_child_cell = in_path_child_cell;

    path_assin_struct_cell{i} = path_struct;
end
    
for i = 1:1:p_n %Use linear programming with topological constraints to determine the neurite assignment
    [assign_ma] = linear_path_assign_pre(path_assin_struct_cell{i});
    m = path_assin_struct_cell{i}.connect_soma_num;
    neuron_id = path_assin_struct_cell{i}.connect_soma;
    path_node = path_assin_struct_cell{i}.path_node;
    path_child_cell = path_assin_struct_cell{i}.path_child_cell;
    [max_Y max_I] = max(assign_ma,[],2);
    for ii = 1:1:m 
        assign_node = (max_I==ii);
        recon_node = path_node(assign_node);
        [new_Neuron ] = path_reconnect(path_child_cell{ii}, recon_node, Neuron{neuron_id(ii)}, neurite_matrix);
        Neuron{neuron_id(ii)} = new_Neuron; 
    end
end