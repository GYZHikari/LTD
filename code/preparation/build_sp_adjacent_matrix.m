function sp_adj_mat = build_sp_adjacent_matrix(label,sp_num)
% $Description:
%    -compute the adjacent matrix
% $Agruments
% Input;
%    -M: superpixel label matrix
%    -N: superpixel number 
% Output:
%    -spAdjcMat: adjacent matrix


sp_adj_mat = zeros(sp_num,sp_num);
[m n] = size(label);

for i = 1:m-1
    for j = 1:n-1
        if(label(i,j)~=label(i,j+1))
            sp_adj_mat(label(i,j),label(i,j+1)) = 1;
            sp_adj_mat(label(i,j+1),label(i,j)) = 1;
        end;
        if(label(i,j)~=label(i+1,j))
            sp_adj_mat(label(i,j),label(i+1,j)) = 1;
            sp_adj_mat(label(i+1,j),label(i,j)) = 1;
        end;
        if(label(i,j)~=label(i+1,j+1))
            sp_adj_mat(label(i,j),label(i+1,j+1)) = 1;
            sp_adj_mat(label(i+1,j+1),label(i,j)) = 1;
        end;
        if(label(i+1,j)~=label(i,j+1))
            sp_adj_mat(label(i+1,j),label(i,j+1)) = 1;
            sp_adj_mat(label(i,j+1),label(i+1,j)) = 1;
        end;
    end;
end;    

    