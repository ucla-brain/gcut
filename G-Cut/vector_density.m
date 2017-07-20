function [density_matrix] = vector_density(vector_matrix,n)
density_matrix = zeros(n,3);
x_min = min(vector_matrix(:,2));
x_max = max(vector_matrix(:,2));
length_x = (x_max - x_min)/n;
i = 1;
for len = x_min:length_x:(x_max-length_x)
    num_a = intersect(find(vector_matrix(:,2)>len),find(vector_matrix(:,2)<(len+length_x)));
    to_num = length(num_a);
    x_density = sum(vector_matrix(num_a,5))/to_num;
    y_density = sum(vector_matrix(num_a,6))/to_num;
    z_density = sum(vector_matrix(num_a,7))/to_num;
    density_matrix(i,:) = [x_density,y_density,z_density];
    i = i+1;
end
