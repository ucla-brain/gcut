function [assign_ma] = linear_path_assign_pre(path_struct)

fitness_matrix = -path_struct.fitness;
connect_ind = path_struct.connect_ind;
path_node_matrix = path_struct.path_node;
connect_soma_num = path_struct.connect_soma_num;

[ assign_ma ] = linear_path_assign( fitness_matrix, connect_ind, path_node_matrix, connect_soma_num );