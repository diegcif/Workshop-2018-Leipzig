needsPackage ("SOS" , Configuration => { "CSDPexec" => "CSDP/csdp"})
needsPackage ("NumericalAlgebraicGeometry")

sortBy := (L, fun) -> 
    last \ sort for l in L list (fun l, l);

zeroSolve = f -> (
    G := squareUp polySystem f;
    sol := solveSystem G;
    pts := for p in sol list
        if not isRealPoint p then continue
        else new Point from {Coordinates => realPart\coordinates p};
    F := matrix {f};
    for p in pts do p.Residual = norm evaluate(F,p);
    pts = sortBy(pts, p -> p.Residual);
    return pts;
    )

roundPoints = (pts,d) -> (
    -- d is the precision in digits
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

R = QQ[x,y,z]
h1 = (x-1)^2 + (y-2)^2
h2 = (y-2)^2 + z^2
h3 = z*(z-1)
h = {h1,h2,h3}

R = QQ[x,y,z]
h1 = 5*x^9 - 6*x^5*y + x*y^4 + 2*x*z
h2 = -2*x^6*y + 2*x^2*y^3 + 2*y*z
h3 = x^2 + y^2 - 265625/1000000
h = {h1,h2,h3}

(h',g,d) = sosIdeal(h,10)
(g',d') = cleanSOS(g,d,1e-3)
pts = select(zeroSolve(g'|h'), p -> p.Residual<1e-3)
print roundPoints(pts,4)

end


R = QQ[x,y]
h1 = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1 --Motzkin
--h1 = (x+y+1)*(x^4+y^4+2)
h = {h1}

R = QQ[x,y,z]
a = 1; b = 3 ; c = -1
h1 = a*(x^6 + y^6 + z^6 ) + b*x^2*y^2*z^2 + c*(x^4*y^2 + x^4*z^2 + x^2*y^4 + x^2*z^4 + y^4*z^2 + y^2*z^4 ) -- robinson
h = {h1}

(f,p) = genericCombination(h, 8)
--out = solveSOS (f, p, Solver=>"CSDP",Verbose=>false)
(Q,mon,X,tval) = solveSOS (f, p, Solver=>"CSDP",Verbose=>false)
(g,d) = sosdec(Q,mon)

S = (ring Q)[x,y,z]
g = for i to #g-1 list 
    if d_i<1e-5 then continue else sub(g_i,S)

--N = 100
--L = for i to N-1 list
--    #solveSystem g;
--print tally L
--L = for i to N-1 list
--    #zeroSolve g;
--print tally L
