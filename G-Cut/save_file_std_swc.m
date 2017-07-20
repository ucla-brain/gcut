function [] = save_file(Save_name, matrix)
  fid=fopen(Save_name,'wt');%写入文件路径  
            [m,n]=size(matrix);  
            for save_i=1:1:m     
                for save_j=1:1:n        
                    if save_j==n           
                        fprintf(fid,'%d\n',matrix(save_i,save_j));       
                    else
                        if save_j==1 || save_j==2 || save_j==7
                            fprintf(fid,'%d ',matrix(save_i,save_j));
                        else
                            fprintf(fid,'%g ',matrix(save_i,save_j));
                        end
                    end
                end
            end
            fclose(fid);
