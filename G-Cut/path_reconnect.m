function [new_Neuron ] = path_reconnect(path_child_matrix, path_node_soma, Neuron, neurite_matrix)

soma_path_neu = path_node_soma;
path_child_matrix_1 = path_child_matrix;
[m_5] = length(path_node_soma);
for i_10 = 1:1:m_5
node_index_5 = path_node_soma(i_10);
node_index_6 = find(path_child_matrix(:,1)==node_index_5);
if ~isempty(node_index_6)
    node_child = path_child_matrix_1(node_index_6,:);
    node_child(node_child==0)=[];
    soma_path_neu = union(soma_path_neu,node_child);
end
end


new_Neuron = Neuron;

for i_10 = 1:1:length(soma_path_neu)
neurite_index = soma_path_neu(i_10);
add_neu = neurite_matrix(neurite_index,:);
new_Neuron = [new_Neuron;add_neu];
end