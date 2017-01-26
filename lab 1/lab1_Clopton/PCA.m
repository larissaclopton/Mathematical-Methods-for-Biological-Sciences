function [PCs,Eigs] = PCA(A,dim,k)

    % calculate the covariance matrix along the specified dimension
    % dim 1 -> perform PCA on rows
    if dim == 1
        cov_matrix = cov(A');
    % dim 2 -> perform PCA on columns
    elseif dim == 2
        cov_matrix = cov(A);
    else
        error('PCA: Invalid dimension.');
    end
   
    % diagonalize the covariance matrix
    [eigvecs,eigvals] = eig(cov_matrix);

    % return k max eigenvalues in descending sorted order
    [diagonals, indexes] = sort(diag(eigvals),'descend');
    Eigs = diagonals(1:k)/sum(diagonals)';
    indexes = indexes(1:k);
    
    % return eigenvectors corresponding to eigenvalues
    PCs = eigvecs(:,indexes);
    
end



