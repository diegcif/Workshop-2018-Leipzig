needsPackage ("SOS" , Configuration => { "CSDPexec" => "CSDP/csdp"})
needsPackage ("NumericalAlgebraicGeometry")

sortBy = (L, fun) -> 
    last \ sort for l in L list (fun l, l);

zeroSolve = f -> (
    maxEntry := v -> max \\ abs \ flatten entries v;
    G := squareUp polySystem f;
    sol := solveSystem G;
    pts := for p in sol list
        if not isRealPoint p then continue
        else new Point from {Coordinates => realPart\coordinates p};
    F := matrix {f};
    for p in pts do p.Residual = maxEntry evaluate(F,p);
    pts = sortBy(pts, p -> p.Residual);
    return pts;
    )

roundPoints = (d,pts) -> (
    -- d is the precision in decimal digits
    for p in pts do 
        p.Coordinates = round_d \ p.Coordinates;
    n := #pts;
    recur := {};
    for i to n-1 do
        for j from i+1 to n-1 do
            if coordinates pts#i == coordinates pts#j then
                recur = append(recur,j);
    recur = set recur;
    pts = for i to n-1 list 
        if member(i,recur) then continue else pts#i;
    return pts;
    )

realZeros = method(
     Options => {RndTol => -3, Solver=>"CSDP", Verbose => false, CleanTol => 1e-3, ResTol => 1e-2} )

realZeros(List,ZZ) := o -> (h,d) -> (
    if #h==0 then error "list of polynomials is empty";
    R := ring h#0;
    if coefficientRing R === QQ then (
        R = changeRingField(RR,R);
        h = toRing_R \ h; 
        );
    (s,mult) := sosInIdeal(h,d,RndTol=>o.RndTol,Solver=>o.Solver,Verbose=>o.Verbose);
    if s===null then
        error "SOS polynomial not found. Try increasing the degree.";
    s' := cleanSOS(s,o.CleanTol);
    h' := h | gens s';
    pts := zeroSolve(h');
    pts' := select(pts, p -> p.Residual<o.ResTol);
    return pts';
    )


-- One polynomial
R = QQ[x,y]
h1 = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1 --Motzkin
--h1 = (x^6 + y^6 + 1 ) + 3*x^2*y^2 - (x^4*y^2 + x^4 + x^2*y^4 + x^2 + y^4 + y^2) -- robinson
h = {h1}
pts = realZeros(h,8, CleanTol=>1e-3, ResTol=>1e-2)
print roundPoints(4,pts)


-- Simple example (one zero)
--R = QQ[x,y,z]
--h1 = (x-1)^2 + (y-2)^2
--h2 = (y-2)^2 + z^2
--h3 = z*(z-1)
--h = {h1,h2,h3}
--pts = realZeros(h,4, CleanTol=>1e-3, ResTol=>1e-2)
--print roundPoints(4,pts)


-- More complicated (eight zeros)
--R = QQ[x,y,z]
--h1 = 5*x^9 - 6*x^5*y + x*y^4 + 2*x*z
--h2 = -2*x^6*y + 2*x^2*y^3 + 2*y*z
--h3 = x^2 + y^2 - 17/64
--h = {h1,h2,h3}
--pts = realZeros(h,10, CleanTol=>1e-3, ResTol=>1e-2)
--print roundPoints(4,pts)
