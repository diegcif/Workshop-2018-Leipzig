needsPackage( "SOS", Configuration=>{"CSDPexec"=>"CSDP/csdp"} )

rank1factor = method(
     Options => {EigTol => 1e-4} )
rank1factor(Matrix) := o -> (X) -> (
    n := numRows X;
    (e, V) := eigenvectors(X,Hermitian=>true);
    if e#0 < 0 then return; -- X not PSD
    if e#(n-2) > o.EigTol then return; -- not rank one
    v := sqrt abs(e#(n-1)) * V_{n-1};
    retval := flatten entries v;
    return retval;
    )

lowerBound = method(
     Options => {RndTol => -3, Solver=>"M2", Verbose => false, EigTol => 1e-4} )

lowerBound(RingElement,List) := o -> (f,pars) -> (
    -- sos lower bound for the polynomial f
    o' := new OptionTable from 
        {RndTol=>o.RndTol, Solver=>o.Solver, Verbose=>o.Verbose};
    R := ring f;
    tvar := symbol t;
    coeffp := coefficientRing R;
    newR := coeffp[gens R | {tvar}];
    phi := map(newR, R);
    tpoly := last gens newR;
    (Q, mon, X, tval) := solveSOS(phi(f)-tpoly,{tpoly}|phi\pars,-tpoly, o');
    if Q===null then return (,);
    bound := first tval;
    x := if X=!=null then rank1factor(X,EigTol=>o.EigTol) else null;
    sol := if x===null then {}
        else for i to numRows mon-1 list (
            y := mon_(i,0);
            if first degree y!=1 then continue;
            y => x_i );
    return (bound, sol);
    )

lasserreHierarchy = method(
     Options => {RndTol => -3, Solver=>"M2", Verbose => false, EigTol => 1e-4} )
lasserreHierarchy(RingElement,List,ZZ) := o -> (f,h,D) -> (
    -- Lasserre hierarchy for the problem
    -- min f(x) s.t. h(x)=0
    if odd D then error "D must be even";
    if all(h|{f}, isHomogeneous) then
        error "problem is homogeneous";
    (H,p) := genericCombination(h, D, false);
    S := ring H;
    f = sub(f,S);
    (bound,sol) := lowerBound(f+H, p, o);
    return (bound,sol);
    )


R=QQ[x];
--f = x^2+y^2+1;
f = (x-1)^2 + (x+3)^2;
r = lowerBound(f, {}, RndTol=>12, Solver => "CSDP");
print(r);

-- minimize f(x) subject to h(x) = 0
R=QQ[x,y,z];
f = 10 - x^2 - y
h = {x^2+y^2+z^2-1}
r = lasserreHierarchy(f,h,2, RndTol=>12, Solver => "CSDP")
print(r);
