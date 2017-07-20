function [path ] = common_path_extract(neurite_matrix, soma_set,orin_soma_node)

%to solve loop, incompleted


[m n] = size(neurite_matrix);
[m_1 n_1] = size(orin_soma_node);
logic_con = zeros(m,m);
for i = 1:1:m
    neurite = neurite_matrix(i,:);
    neurite(isnan(neurite)) = [];
    neurite(ismember(neurite,soma_set)) = [];
    if neurite(1)~=
    [r c] = find(neurite_matrix == neurite(1));
    [r_1 c_1] = find(neurite_matrix == neurite(end));
    con_set = union(r,r_1);
    con_set(con_set ==i ) = [];
    logic_con(i,con_set) = 1;
end

for so_i = 1:1:m_1
    logic_con(orin_soma_node{so_i},orin_soma_node{so_i}) = 0;
end

path = cell(1,m_1);
soma_connect = zeros(m_1,m_1);

for soma_i = 1:1:m_1
    soma_j = 1:m_1;
    soma_j(soma_j == soma_i) = [];
    ind_1 = 1;
    while(~isempty(soma_j))
        [d P ] = dijkstra(logic_con, orin_soma_node{soma_i}, orin_soma_node{soma_j(ind_1)});
        if d == inf
            soma_connect(soma_i, soma_j(ind_1)) = 0;
            soma_j(ind_1) = [];
        else
            path_in = P;
            [d P] = dijkstra(logic_con, path_in, orin_soma_node{soma_j});
        end
    end
end