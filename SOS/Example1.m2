restart
needsPackage ("SOS" , Configuration => { "CSDPexec" => "CSDP/csdp"})

-- R = QQ[x1,x2]
-- h1 = x2^4*x1+3*x1^3-x2^^4 -3*x1^2
-- h2 = x1^2*x2 - 2*x1^2
-- h3 = 2*x2^4*x1 - x1^3 -2*x2^4 +x1^2
-- 
-- t = 4
-- 
-- -- Even simpler

R = QQ[x,y,z, g1, g2]
h1 = x^2 + y^2
h2 = y^2 + z^2

f = g1*h1 + g2*h2

-- (Q,mon, X) = solveSOS (f, {g1,g2}, Solver=>"CSDP")


makePPoly = (d, R, p, i) -> (
    -- d is the desired degree of the generic polynomial in R.  p is a
    -- symbol used to make variables for the coefficients and i is the
    -- first index of these variables.
    mon := flatten entries basis (d, R);
    newvars := p_(i,0) .. p_(i,#mon-1);
    return (newvars, mon)
    )


formPoly = (coeff, mons) -> (
    sum apply (coeff, mons, (i,j)->i*j))

makePPoly (2, R, symbol p, 0)


constructf = (h, t) -> (
    -- h is a list of polynomials
    -- 2t is a maximumd degree
    R := ring h#0;

    -- compute complementary degrees:
    dbar := for hi in h list 2*t - first degree hi;

    p := symbol p;
    i := -1;
    L := for d in dbar list (
	i = i + 1;
	makePPoly (d, R, p, i)
	);
    pvars := flatten for l in L list l#0;
    Q := newRing (R, Variables=> gens R|pvars);
    g := for i to #h-1 list (
	formPoly ( apply (L#i#0, m->Q_m) , apply (L#i#1, m -> sub (m, Q)) )
	);
    sum apply (h,g, (I,J)->sub(I,Q)*J)
    )

Q = constructf ({h1,h2},2)
    
    
    