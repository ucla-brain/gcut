function [index ] = find_vector_beta(Matrix, vector)
index = [];
id_1 = find(Matrix(:,1)==vector(1));
id_2 = find(Matrix(:,2)==vector(2));
id = intersect(id_1,id_2);
if ~isempty(id)
    index = [index; id];
end