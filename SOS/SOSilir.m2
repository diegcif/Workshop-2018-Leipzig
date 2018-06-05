eigenVecLargestEigenval = X -> (
    if rank X == 0 then return {}; -- X degenerate matrix
    (eigvals, eigvecs) = eigenvectors X;
    if any(eigvals, x -> x < 0) then return {}; -- X not PSD
    return flatten entries eigvecs^{maxPosition eigvals};
)

minimizePoly = (p, tvar) -> (
    ringp = ring p;
    if member(tvar, gens ringp) then error "Please provide a new variable name";
    coeffp = coefficientRing ringp;
    newR := coeffp[gens ringp | {tvar}];
    F = map(newR, ringp);
    newp = F(p); 
    use newR;
    tpoly = t;
    (Q, mon, X, tval) = solveSOS(newp-tpoly,{tpoly},-tpoly, RndTol=>12);
    return (tval, eigenVecLargestEigenval Q);
)
