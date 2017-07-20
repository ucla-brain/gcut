function [Neuron path_child_cell fitness_matrix connect_ind new_path_node] = neurite_route_sort_delta_3(Neu_con_ma, logic_con_ma, Neuron, neurite_matrix, soma_con_num, orin_soma_node, path_node, node_weight )
%This function is used to search all the child neurite in the common path
%and simplify the segmentation problem.

m = soma_con_num;
[m2 n2] = size(Neu_con_ma);
glob_ma = cell(m2,n2);
path_child_cell = cell(m2,n2);
connect_ind = cell(m2,n2);
delect_union_cell = cell(1,m);



for i = 1:1:m
    Children_matrix = Neu_con_ma{i};
    Parent_matrix = Children_matrix';
    [l_m l_n] = size(logic_con_ma{i});
    index_2 =[];    
    [m3 n3] = size(path_node);
    add_path_node = zeros(l_m,1);
    add_node_child = zeros(l_m,l_m);
    add_node_fit_ma = zeros(l_m,m+n3);
    add_pre_ma = zeros(l_m,2);
    add_j = 1;
    
    %--------------------------------
    sparse_ma_previous = zeros(m3,1);
    sparse_ma_next = zeros(m3,1);
    
    %----------------------------------
    path_matrix = zeros(m3,m+n3);  
    path_matrix(:,1) = path_node;
    com_node = path_node;
    new_path_node = path_node;
    [soma_node  path_child_node] = find_connection(com_node,orin_soma_node{i},logic_con_ma{i});
    
        
    add_node_id = find(Children_matrix(soma_node,:)>0);
    
    add_node_id(ismember(add_node_id,path_child_node))=[];
    
    check_node = [];
    for node_check_i = 1:1:m
        [re_path_node re_path_add] = find_connection(add_node_id, path_child_node,logic_con_ma{node_check_i});
        check_node = union(check_node, re_path_add);
    end
    
    add_node_id = intersect(add_node_id, check_node);
    
    union_set_1 =[];
    if isempty(add_node_id)
        union_set_1 =[];
    else
        for add_i = 1:1:length(add_node_id)
            child_ma = zeros(l_m,1);
            count_i = 1;
            [child_ma count_i sum_node sum_penaty] = com_a_beta_1(child_ma, count_i,Children_matrix, add_node_id(add_i), node_weight); %compute all the penalty for the child node
            child_ma(child_ma==0) = [];
            add_path_node(add_j) = add_node_id(add_i);
            if length(child_ma)>1
                for iii = 1:1:length(child_ma)
                    add_node_child(add_j,iii) = child_ma(iii);
                end
            end
            add_node_fit_ma(add_j,n3) = add_node_id(add_i);
            add_node_fit_ma(add_j,i+n3) = sum_penaty;
            add_pre_ma(add_j,1) = add_node_id(add_i);
            add_pre_ma(add_j,2) = soma_node;
            add_j = add_j + 1;
            union_set_1 = union(union_set_1,child_ma);
        end
    end
    for child_i = 1:1:length(path_child_node) 
        
       
       path_index_2 = (find(path_node == path_child_node(child_i)));
        
       sparse_ma_next(path_index_2) = path_child_node(child_i);
       
       sparse_ma_previous(path_index_2) = soma_node;
       
        path_matrix(path_index_2,n3+i) = Children_matrix(soma_node,path_child_node(child_i))*node_weight(path_child_node(child_i),1); %Use the length of neurite as weight
    
    end

    index_2 = [index_2;path_child_node];
    
    if isempty(index_2)
        disp('error!!');
    end
    
    index_1 = index_2;
    
    com_node(ismember(com_node,index_1)) = [];
    
    while(~isempty(com_node))
        
        index_2 = []; 
        index_3 = [];
        
        for child_i = 1:1:length(index_1)
        	[orin_node path_child_node] = find_connection(com_node,index_1(child_i),logic_con_ma{i});
            if isempty(orin_node)
                orin_node = index_1(child_i);
            end
             add_node_id = find(Children_matrix(orin_node,:)>0);
             
             add_node_id(ismember(add_node_id,path_child_node))=[];
 
             for add_i = 1:1:length(add_node_id)
             
                 count_i = 1;
                 child_ma = zeros(l_m,1);
                 [child_ma count_i sum_node sum_penaty ] = com_a_beta_1 (child_ma, count_i,Children_matrix, add_node_id(add_i), node_weight); 
                 child_ma(child_ma==0)=[];       
                 add_path_node(add_j) = add_node_id(add_i);
                if length(child_ma)>1
                    for iii = 1:1:length(child_ma)
                        add_node_child(add_j,iii) = child_ma(iii);
                    end
                end
                 add_node_fit_ma(add_j,n3) = add_node_id(add_i);
                 add_node_fit_ma(add_j,i+n3) = sum_penaty;
                 add_pre_ma(add_j,1) = add_node_id(add_i);
                 add_pre_ma(add_j,2) = orin_node;
                 add_j = add_j + 1;
                 union_set_1 = union(union_set_1,child_ma);
             end
            
            
            
            for child_j = 1:1:length(path_child_node) 
                   
                   path_index_2 = (find(path_node == path_child_node(child_j)));
                   
                   sparse_ma_next(path_index_2) = path_child_node(child_j);
                   
                   sparse_ma_previous(path_index_2) = orin_node;
                   
                   path_matrix(path_index_2,n3+i) = Children_matrix(orin_node,path_child_node(child_j))*node_weight(path_child_node(child_j),1); %
                    
                   [orin_node_1 child_node] = find_connection(com_node,path_child_node(child_j),logic_con_ma{i});
                    
                   if isempty(orin_node_1)||isempty(child_node) 
                        
                        add_node_id = find(Children_matrix(path_child_node(child_j),:)>0);
                        for add_i = 1:1:length(add_node_id)
                            count_i = 1;
                            child_ma = zeros(l_m,1);
                            [child_ma count_i sum_node sum_penaty ] = com_a_beta_1 (child_ma, count_i,Children_matrix, add_node_id(add_i), node_weight);
                            child_ma(child_ma==0) = [];
                            add_path_node(add_j) = add_node_id(add_i);
                            if length(child_ma)>1
                                for iii = 1:1:length(child_ma)
                                    add_node_child(add_j,iii) = child_ma(iii);
                                end
                            end
                            add_node_fit_ma(add_j,n3) = add_node_id(add_i);
                            add_node_fit_ma(add_j,i+n3) = sum_penaty;
                            add_pre_ma(add_j,1) = add_node_id(add_i);
                            add_pre_ma(add_j,2) = path_child_node(child_j);
                            add_j = add_j + 1;
                            union_set_1 = union(union_set_1,child_ma);
                        end
                            
                            index_3 = [index_3;path_child_node(child_j)];
                    end
                    
            end
            index_2 = [index_2;path_child_node];
            

        end
        
        index_1 = index_2; 
        
        com_node(ismember(com_node,index_1)) = []; 
        
        index_1(ismember(index_1,index_3))=[]; 
    end
    add_path_node(add_path_node==0)=[];
    add_len = length(add_path_node);
    add_node_child(add_node_child(:,1)==0,:)=[];
    add_node_fit_ma(add_node_fit_ma(:,1)==0,:)=[];
    add_pre_ma(add_pre_ma(:,1)==0,:)=[];
    
    path_child_matrix = add_node_child;
    
    new_path_node = [new_path_node; add_path_node];
    path_matrix = [path_matrix;add_node_fit_ma];
    
    connect_ind_ma = [sparse_ma_next,sparse_ma_previous];
    connect_ind_ma = [connect_ind_ma; add_pre_ma];
    
    union_set_1 =new_path_node;
    [s_m s_n] = size(add_node_child);
    for un_i = 1:1:s_m
        union_add = add_node_child(un_i,:);
        union_add(union_add==0)=[];
        union_set_1 = union(union_set_1,union_add);
    end

    de_connect_ma = ismember(connect_ind_ma(:,2), soma_node);
    connect_ind_ma(de_connect_ma,:) = [];
    
    for de_i = 1:1:m
        de_connect_ma = ismember(connect_ind_ma(:,1),orin_soma_node{de_i});
        connect_ind_ma(de_connect_ma,:)=[];
    end
    
    connect_ind{i} = connect_ind_ma;
    
    glob_ma{i} = path_matrix;
    
    path_child_cell{i} = path_child_matrix;
    
    m1 = length(union_set_1);
    
    union_index_1 = zeros(m1,1);
    
    for set_i = 1:1:m1
        union_id = union_set_1(set_i);
        
        union_neurite = neurite_matrix(union_id,:);
        
        back_id_1 = find(Neuron{i}(:,1)==union_neurite(1));
        back_id_2 = find(Neuron{i}(:,2)== union_neurite(2));
        
        union_index_1(set_i) = intersect(back_id_1,back_id_2);
        
        if length(union_index_1(set_i))~=1
            disp('error');
        end
        
    end
   delect_union_cell{i} = union_index_1; 
    
        
    
end



[g_m g_n] = size(glob_ma{1});
glob_matrix = zeros(g_m, m+n3);
glob_matrix(:,1) = sort(glob_ma{1}(:,1));


for i =1:1:m
    Neuron{i}(delect_union_cell{i},:) = []; 

    for mg_i = 1:1:g_m
        mg_ind_1 = find(glob_ma{i}(:,1)==glob_matrix(mg_i,1));
        if isempty(mg_ind_1)
            disp('warning!');
            glob_matrix(mg_i,i+1) = 0;
        else
            glob_matrix(mg_i,i+1) = glob_ma{i}(mg_ind_1,i+1);
        end

    end

end


new_path_node = glob_matrix(:,1);

fitness_matrix = glob_matrix(:,2:m+1);
