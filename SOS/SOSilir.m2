loadPackage"SOS"
eigenVecLargestEigenval = (X, rndErr)  -> (
    if X == null then return {}; -- X null
    if rank X == 0 then return {}; -- X degenerate matrix
    (eigvals, eigvecs) = eigenvectors X;
    if any(eigvals, x -> x < 0) then return {}; -- X not PSD
    retval = flatten entries eigvecs^{maxPosition eigvals};
    return apply(retval, v -> round(rndErr, realPart(v)));
)

minimizePoly = (p, tvar, rndErr) -> (
    ringp = ring p;
    if member(tvar, gens ringp) then error "Please provide a new variable name";
    coeffp = coefficientRing ringp;
    newR := coeffp[gens ringp | {tvar}];
    F = map(newR, ringp);
    newp = F(p); 
    use newR;
    tpoly = t;
    (Q, mon, X, tval) = solveSOS(newp-tpoly,{tpoly},-tpoly, RndTol=>12);
    return (tval, eigenVecLargestEigenval(X, rndErr));
)

--Test
R=QQ[x,y];
f = x^2+y^2+1;
r = minimizePoly(f, t, 3);
print(r);
