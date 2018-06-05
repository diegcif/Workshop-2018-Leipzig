eigenVecLargestEigenval = X -> (
    if rank X == 0 then return {}; -- X degenerate matrix
    (eigvals, eigvecs) = eigenvectors X;
    if any(eigvals, x -> x < 0) then return {}; -- X not PSD
    return flatten entries eigvecs^{maxPosition eigvals};
)

minimizePoly = p -> (
    -- take the ambient ring and adjoin variable t
    -- now assume it is there
    (Q, mon, X, tval) = solveSOS(f-t,{t},-t);
    return (tval, eigenVecLargestEigenval Q);
)
