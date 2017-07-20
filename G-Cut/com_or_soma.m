function [branch_order_s ] =  com_or_soma(index, Parent_matrix, Children_matrix, soma_node)
branch_order_s  = 1;
while(~isequal(index,soma_node))

        child = find(Children_matrix(index,:)==1); %如果为分支结点，则找到其相应的所有子结点
        [cm, cn] = size(child);      
        if cn>1
            branch_order_s = branch_order_s+1;
        end  
        index = (find(Parent_matrix(index,:)==1));
end
   