clear all;
clc;
path_1 = 'F:\MCP_new\Neuron_connet\random_connect\test\';
path_2 = 'F:\MCP_new\Neuron_connet\random_connect\test\';

%parpool('local',2);

file_list = dir(fullfile(path_1,'*.swc'));
file_n = length(file_list);
for file_i = 1:1:file_n
    file_name = file_list(file_i).name;
    [path, or_name, ext] = fileparts(file_name);
    A_1 = swc_read(strcat(path_1,file_name));
    index_import = 0;
    if index_import == 0
        index_x_1 = importdata(strcat(path_1 ,or_name,'_soma_ind.txt'));
        index_x_1 = index_x_1';
    else
        %index_x_1 = [1,130,1301,1634,2632,3611,4580,5256,6031,7385,8666,10229,11638,12023,13326,14114];
        index_x_1 = [1,56,1801];
    end
    %index_x_1 = [1,200,1650,2249,4445,6746,7488,8381];
    %index_x_1 = [1,141,96,737,743,1628,2262,2735];

    %index_x_1 = [1,514,1149,1690,1757,3066,3737,4459,5573,6428,9588,9910,10860,12419,12760,12902,14672,15600];
    %index_x_1 = [1,299,1760];
    %index_x_1 = [1,124];
    %index_x_1 = [1,306,1898];
    %index_x_1 = [1,56,1801];
    %-----------------------soma merge------------------------------
    %[index_x A] = soma_merge_beta(A_1,index_x_1);
    %---------------no soma merge-----------------------------------
    A = A_1;
    index_x = index_x_1;
    %----------------------------------------------------------------
    index_x = sort(index_x);
    tic;
    [m,n] = size(A);
    branch_length = 0;
    surface = 0;
    Volume = 0;
    Diameter = 0;
    B = A(:,3:5);
    body_size = 0;
    Matrix_path = sparse(zeros(m,m));
    %Matrix_path = sparse(Matrix_path);
    vector_matrix = zeros(m,9); %矢量矩阵，记录结点上每个点的书序（1），空间坐标（2:4），该点与父节点的方向矢量（5:7），该点的父节点（8），以及该点矢量与标准矢量的夹角（9）。
    vector_matrix(:,1) = A(:,1);
    vector_matrix(:,2:4) = B;
    vector_matrix(:,8) = A(:,7);
    standard_vector = zeros(m,3); %标准矢量，记录从soma到该点的矢量
%     [YY II] = sort(A(:,1));
%     for i = 1:1:length(index_x)
%         index_x(i) = find(A(:,1)==index_x(i));
%     end
%     A = A(II,:);
    for i =m:-1:1
        if (A(i,7) == -1)

            continue;

        else
            index_1 = A(i,7);
%             r_index = find(A(:,1)==index_1);
%             if isempty(r_index)
%                 disp(num2str(index_1));
%             end
%             Matrix_path(i,r_index) = 1;
%             Matrix_path(r_index,i) = 1;
%             A(r_index,1) = r_index;
%             A(i,1) = i;
%             A(i,7) = r_index;
            Matrix_path(i,index_1) = 1;  %parents路径
            Matrix_path(index_1,i) = 1;  %Children 路径
        end
    end

    Ma_Mask_low = sparse(triu(ones(m,m),0));
    %Ma_Mask_low = sparse(Ma_Mask_low);
    %Matrix_path_childe = zeros(m,m);
    Matrix_path_childe = sparse(m,m);
    % for i = 1:1:m
    %     for j = 1:1:m
    %         Matrix_path_childe(i,j) = Matrix_path(i,j) * Ma_Mask_low(i,j);
    %     end
    % end
    Matrix_path_childe = Matrix_path .* Ma_Mask_low;
    %得到Children表
    Matrix_pa = mean(Matrix_path_childe,2).*m;
    leaf_node = find(Matrix_pa==0);
    branch_node = find(Matrix_pa>1);

    density_matrix = vector_density(vector_matrix,6);
    stem = Matrix_pa(1);
    Ma_Mask_up = sparse(tril(ones(m,m),0));
    %Ma_Mask_up = sparse(Ma_Mask_up);
    %Matrix_path_parent = zeros(m,m);
    Matrix_path_parent = sparse(m,m);
    % for i = 1:1:m
    %     for j = 1:1:m
    %         Matrix_path_parent(i,j) = Matrix_path(i,j) * Ma_Mask_up(i,j);
    %     end
    % end
    Matrix_path_parent = Matrix_path .* Ma_Mask_up;
    %得到parents表
    [neurite_ma] = split_neurite_beta(branch_node,leaf_node,Matrix_path_parent, Matrix_path_childe);
    %[prun_m,branch_an_vector,branch_child] = branch_pruning(vector_matrix,vector_matrix(:,9),Matrix_path_childe,Matrix_path_parent,branch_node,leaf_node,1); %进行剪枝[m1,n1] = size(leaf_node);
    Set_index = zeros(m,1);
    %index_x = [1,299,1760];

    Set_index(index_x) = 1;
    [era_i era_j] = size(neurite_ma);
    for pre_i =1:1:era_i
        era = find(neurite_ma(pre_i,:)==0);
        neurite_ma(pre_i,era) = NaN;
    end
    disp('split complete');
    toc;
    [in_neuron, rec_matrix] = connect_neurite_beta_1(A,neurite_ma,Set_index,m);
    for i = 1:1:length(rec_matrix)
        save_file_std_swc(strcat(path_2, or_name, '_split_',num2str(i),'.swc'),rec_matrix{i});
    end
    disp(strcat(num2str(file_i),'_complete!'));
end

%matlabpool close;