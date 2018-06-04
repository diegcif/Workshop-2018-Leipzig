eigenVecLargestEigenval = X -> (
    if rank X == 0 then return {}; -- X degenerate matrix
    (eigvals, eigvecs) = eigenvectors X;
    if any(eigvals, x -> x < 0) then return {}; -- X not PSD
    return flatten entries eigvecs^{maxPosition eigvals};
)

