needsPackage( "SOS", Configuration=>{"CSDPexec"=>"CSDP/csdp"} )
eigenVecLargestEigenval = (X, rndErr)  -> (
    if X === null then return; -- X null
    n := numRows X;
    (e, V) := eigenvectors(X,Hermitian=>true);
    if e#0 < 0 then return; -- X not PSD
    if e#(n-2) > 1e-3 then return; -- not rank one
    v := sqrt abs(e#(n-1)) * V_{n-1};
    retval := flatten entries v;
    return apply(retval, u -> round(rndErr, u));
)

minimizePoly = method(
     Options => {RndTol => -3, Solver=>"M2", Verbose => false} )
minimizePoly(RingElement,ZZ) := o -> (p, rndErr) -> (
    ringp := ring p;
    tvar := symbol t;
    coeffp := coefficientRing ringp;
    newR := coeffp[gens ringp | {tvar}];
    F := map(newR, ringp);
    newp := F(p); 
    tpoly := last gens newR;
    (Q, mon, X, tval) := solveSOS(newp-tpoly,{tpoly},-tpoly, o);
    opt := first tval;
    x := eigenVecLargestEigenval(X, rndErr);
    dic := if x===null then {}
        else for i to numRows mon-1 list (
            y := mon_(i,0);
            if first degree y!=1 then continue;
            y => x_i );
    return (opt, dic);
)

--Test
R=QQ[x];
--f = x^2+y^2+1;
f = (x-1)^2 + (x+1)^2;
r = minimizePoly(f, 3, RndTol=>12, Solver => "CSDP");
print(r);
