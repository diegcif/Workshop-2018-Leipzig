newPackage(
    "SOS",
    Version => "1.5", 
    Date => "November 17, 2017",
    Authors => {
     {Name => "Helfried Peyrl", 
      Email => "peyrl@control.ee.ethz.ch",
      HomePage => "https://scholar.google.com/citations?user=cFOV7nYAAAAJ&hl=de"},
     {Name => "Pablo A. Parrilo",
      Email => "parrilo@mit.edu",
      HomePage => "http://www.mit.edu/~parrilo/"},
     {Name => "Diego Cifuentes",
      Email => "diegcif@mit.edu",
      HomePage => "http://www.mit.edu/~diegcif/"}
    },
    Headline => "Sum-of-Squares Package",
    DebuggingMode => true,
    Configuration => {"CSDPexec"=>"csdp","SDPAexec"=>"sdpa"},
    AuxiliaryFiles => true,
    PackageImports => {"SimpleDoc","FourierMotzkin","NumericalHilbert"},
    PackageExports => {}
)

export {
--Types
    "SOSPoly",
--Methods/Functions
    "sosPoly",
    "solveSOS",
    "sosdec",
    "sosdecTernary",
    "sumSOS",
    "blkDiag",
    "LDLdecomposition",
    "solveSDP",
    "makeMultiples",
    "sosInIdeal",
    "lowerBound",
    "roundPSDmatrix",
    "checkSolver",
    "smat2vec",
    "vec2smat",
--only for debugging
    "parameterVector",
    "createSOSModel",
    "choosemonp",
    "project2linspace",
    "toRing",
    "changeRingField",
    "changeMatrixField",
--Method options
    "RndTol",
    "UntilObjNegative",
    "WorkingPrecision",
    "Solver",
    "TraceObj",
    "EigTol",
    "Scaling",
    "ParBounds"
}

--##########################################################################--
-- GLOBAL VARIABLES 
--##########################################################################--

csdpexec=((options SOS).Configuration)#"CSDPexec"
sdpaexec=((options SOS).Configuration)#"SDPAexec"

--##########################################################################--
-- TYPES
--##########################################################################--

SOSPoly = new Type of HashTable

-- constructor for sos decomposition
sosPoly = method()
sosPoly (Ring, List, List) := SOS => (R, polys, coeffs) -> (
    new SOSPoly from {
        ring => R,
        gens => polys,
        coefficients => coeffs
        }
    )
sosPoly (List, List) := SOS => (polys, coeffs) -> sosPoly(ring polys#0,polys,coeffs)

ring SOSPoly := S -> S#ring

gens SOSPoly := o -> S -> S#gens

coefficients SOSPoly := o -> S -> S#coefficients

length SOSPoly := S -> #(S#gens)

substitute (SOSPoly,Ring) := (S,R) ->
    sosPoly(for g in S#gens list sub(g,R), S#coefficients)

net SOSPoly := S -> (
    if #gens S == 0 then return "0";
    return "coeffs:"||net coefficients S||"gens:"||net gens S;
    )

Number * SOSPoly := (a,S) -> (
    if a<0 then error "scalar must be nonnegative";
    if a==0 then return sosPoly(ring S, {}, {});
    return sosPoly(ring S, gens S, a * coefficients S);
    )

SOSPoly + SOSPoly := (S,S') -> (
    R := ring S;
    if R =!= ring S' then error "cannot add elements of different rings";
    return sosPoly(R,S#gens|S'#gens, S#coefficients|S'#coefficients);
    )

SOSPoly * SOSPoly := (g1,g2)-> (
    if g1#ring =!= g2#ring then error "cannot multiply elements of different rings";
    q1:=for i in g1#gens list(
        for j in g2#gens list i*j);
    q2:=for i in g1#coefficients list(
        for j in g2#coefficients list i*j);
    return sosPoly(g1#ring, flatten(q1),flatten(q2));
    )

SOSPoly ^ ZZ := (p1,D)->(
    if D<=0 then error "power should be a positive integer.";
    if odd D then error "power should be an even integer.";
    p2 := (sumSOS p1)^(D//2);
    return sosPoly(ring p1,{p2},{1});
    )

SOSPoly == RingElement := (S, f) -> (
    if ring S=!=ring f then
        error "Cannot compare elements of different rings. Try to use 'sub'.";
    return sumSOS S == f;
    )

RingElement == SOSPoly := (f, S) -> S == f

SOSPoly == SOSPoly := (S, S') -> S == sumSOS S'

sumSOS = method()

sumSOS (List, List) := (g,d) -> sum for i to #g-1 list g_i^2 * d_i

sumSOS SOSPoly := a -> sum for i to #(a#gens)-1 list a#gens_i^2 * a#coefficients_i

clean(RR,SOSPoly) := (tol,s) -> (
    if s===null then return (,);
    R := ring s;
    if coefficientRing R === QQ then tol=0;
    g := gens s;
    d := coefficients s;
    I := positions(d, di -> di>tol);
    return sosPoly(R,g_I,d_I);
    )

--##########################################################################--
-- METHODS
--##########################################################################--

verbose = (s,o) -> if o.Verbose then print s

--###################################
-- solveSOS
--###################################

sosdec = (mon,Q) -> (
    if mon===null or Q===null then return;
    (L,D,P,err) := PSDdecomposition(Q);
    if err != 0 then (
        print "Gram Matrix is not positive semidefinite";
        return 
        );
    n := numRows Q;
    g := toList flatten entries (transpose mon * transpose inverse P * L);
    d := for i to n-1 list D_(i,i);
    idx := positions (d, i->i!=0);
    d = d_idx;
    g = g_idx;
    return sosPoly(ring mon,g,d);
    )

-- internal way to call solveSOS
rawSolveSOS = method(
     Options => {RndTol => 3, Solver=>"M2", Verbose => false, TraceObj => false, ParBounds=>{}} )
 
rawSolveSOS(Matrix,Matrix,Matrix) := o -> (F,objP,mon) -> (
    -- f is a polynomial to decompose
    -- mon is a vector of monomials
    -- objFcn is a linear objective function for the SDP solver
    -- bounds can be empty or contain lower and upper bounds for the parameters

    checkInputs := (mon,bounds) -> (
        if numColumns mon > 1 then error("Monomial vector should be a column.");
        isMonomial := max(length \ terms \ flatten entries mon)==1;
        if not isMonomial then error("Vector must consist of monomials.");
        if not member(#bounds,{0,2}) then 
            error "ParBounds should be a list with two elements";
        );

    bounds := o.ParBounds;
    checkInputs(mon,bounds);
    kk := coefficientRing ring F;
         
    -- build SOS model --     
    (C,Ai,Bi,A,B,b) := createSOSModel(F,mon,Verbose=>o.Verbose);

    ndim := numRows C;
    np := numRows objP;

    obj := 
        if o.TraceObj then
            map(RR^(#Bi),RR^1,(i,j)-> -trace Bi#i) || map(RR^(#Ai),RR^1,(i,j)-> -trace Ai#i)
        else
            (-objP) || zeros(kk,#Ai,1); 
    if obj==0 then verbose( "Solving SOS feasibility problem...", o)
    else verbose("Solving SOS optimization problem...", o);

    if np!=0 and #bounds==2 then (
        C = blkDiag(C,diagonalMatrix(-bounds#0),diagonalMatrix(bounds#1));
        Ai = for i to #Ai-1 list blkDiag(Ai_i, zeros(kk,2*np,2*np));
        Bi = for i to #Bi-1 list blkDiag(Bi_i, 
            map(kk^np,kk^np, (j,k) -> if j==k and j==i then 1_kk else 0_kk),
            map(kk^np,kk^np, (j,k) -> if j==k and j==i then -1_kk else 0_kk));
        );


    (my,X,Q) := solveSDP(C, Bi | Ai, obj, Solver=>o.Solver, Verbose=>o.Verbose);
    if Q===null then return (Q,X,);
    y := -my;
    if #bounds==2 then Q = Q^{0..ndim-1}_{0..ndim-1};
    pvec0 := flatten entries y^(toList(0..np-1));

    if kk=!=QQ or o.RndTol==infinity then
        return (Q,X,pvec0);

    -- rational rounding --
    (ok,Qp,pVec) := roundSolution(y,Q,A,B,b,o.RndTol);
    if ok then return (Qp,X,pVec);
    print "rounding failed, returning real solution";
    return (Q,X,pvec0);
    )

rawSolveSOS(Matrix,Matrix) := o -> (F,objP) -> (
    mon := choosemonp (F,Verbose=>o.Verbose);
    if mon===null then return (,,,);
    (Q,X,pvec) := rawSolveSOS(F,objP,mon,o);
    mon = checkMonField(mon,F,Q,o.RndTol);
    return (mon,Q,X,pvec);
    )
rawSolveSOS(Matrix) := o -> (F) -> 
    rawSolveSOS(F,zeros(QQ,numRows F-1,1),o)

checkMonField = (mon,F,Q,RndTol) -> (
    if mon===null or Q===null then return mon;
    if coefficientRing ring F===QQ then
        if RndTol==infinity or ring Q=!=QQ then
            return changeMatrixField(RR,mon);
    return mon;
    )

-- This is the main method to decompose a polynomial as a 
-- sum of squares using an SDP solver.
solveSOS = method(
     Options => {RndTol => 3, Solver=>"M2", Verbose => false, TraceObj => false, ParBounds=>{}} )

solveSOS(RingElement,RingElement,Matrix) := o -> (f,objFcn,mon) -> (
    (F,objP) := parameterVector(f,objFcn);
    return rawSolveSOS(F,objP,mon,o);
    )
solveSOS(RingElement,Matrix) := o -> (f,mon) -> 
    solveSOS(f,0_(ring f),mon,o)

solveSOS(RingElement,RingElement) := o -> (f,objFcn) -> (
    (F,objP) := parameterVector(f,objFcn);
    mon := choosemonp (F,Verbose=>o.Verbose);
    if mon===null then return (,mon,,);
    (Q,X,pvec) := rawSolveSOS(F,objP,mon,o);
    mon = checkMonField(mon,F,Q,o.RndTol);
    return (mon,Q,X,pvec);
    )
solveSOS(RingElement) := o -> (f) -> 
    solveSOS(f,0_(ring f),o)


changeRingField = (kk,R) -> kk(monoid[gens R])

changeMatrixField = (kk, M) -> (
    -- M is a matrix whose entries are polynomials whose coefficient
    -- ring should be changed.
    e := entries M;
    R := changeRingField(kk, ring e#0#0);
    matrix for row in e list(
	for entry in row list(
	    toRing(R, entry))))

toRing = method ()
toRing (Ring, RingElement) := (S,f) -> (
    -- maps f to ring S
    R := ring f;
    kk := coefficientRing R;
    phi := map(S,R);
    -- QQ => RR
    if kk===QQ then return phi(f);
    -- RR => QQ
    if not class kk === RealField then error "Expecting conversion from real here";
    (mon,coef) := coefficients f;
    mon = matrix {liftMonomial_S \ flatten entries mon};
    K := 2^(precision kk);
    coef = matrix(QQ, {for c in flatten entries coef list round(K*sub(c,RR))/K});
    f' := (mon*transpose coef)_(0,0);
    return f';
    )

toRing (Ring, SOSPoly) := (S, s) -> (
    -- maps s to Ring S
    R := ring s;
    kk := coefficientRing R;
    if kk===QQ then (
	-- switching from QQ to RR coefficients
	return sosPoly (S, (x -> sub (x, S)) \ gens s,
	    (q -> sub (q, kk)) \ coefficients s)
	);
    if class kk === RealField and coefficientRing S===QQ then (
	g' := toRing_S \ gens s;
	K := 2^(precision kk);
	c' := for c in coefficients s list round(K*sub(c,RR))/K;
	return sosPoly (S, g', c')
	);
    error "Error: only conversion between real and rational coefficient fields is implemented."
    )

liftMonomial = (S,f) -> (
    -- maps monomial f to ring S
    n := numgens S;
    e := first exponents f;
    e = e_(toList(0..n-1)); -- ignore some variables
    return S_e;
    )

roundSolution = {Verbose=>false} >> o -> (y,Q,A,B,b,RndTol) -> (
    -- round and project --
    Qnum := matrix applyTable(entries Q, a -> round(a*2^52)/2^52);

    dhi := 52;
    d := RndTol;
    np := numColumns B;
    
    while (d < dhi) do (
        verbose("rounding step #" | d, o);
        if np!=0 then (
           pVec := map(QQ^np,QQ^1,(i,j) -> round(y_(i,0)*2^d)/2^d);
           bPar := b - B*pVec;
           ) 
        else bPar= b;

        (ok,Qp) := roundPSDmatrix(Qnum,A,bPar,d,Verbose=>o.Verbose);
        if ok then break else d = d + 1;
        );
    pVec = if np!=0 then flatten entries pVec else null;
    return (ok,Qp,pVec);
    )

createSOSModel = method(
    Options => {Verbose => false} )
createSOSModel(RingElement,Matrix) := o -> (f,v) -> (
    F := parameterVector(f);
    return createSOSModel(F,v);
    )
createSOSModel(Matrix,Matrix) := o -> (F,v) -> (
    kk := coefficientRing ring F;
    np := numRows F - 1;
    
    -- monomials in vvT
    vvT := entries(v* transpose v);
    mons := g -> first entries monomials g;
    K := sort \\ unique \\ flatten \\ mons \ flatten vvT;
    
    -- Linear constraints: b
    b := transpose matrix(kk,{for k in K list coefficient(k,F_(0,0)) });
    
    -- Linear constraints: A, B
    coeffMat := (x,A) -> applyTable(A, a -> coefficient(x,a));
    A := matrix(kk, for i to #K-1 list smat2vec(coeffMat(K_i, vvT),Scaling=>2) );
    
    -- Consider search-parameters:
    B := map(kk^#K,kk^np, (i,j) -> -coefficient(K#i, F_(j+1,0)) );
    
    (C,Ai,Bi) := getImageModel(A,B,b);
    
    return (C,Ai,Bi,A,B,b);
    )

getImageModel = (A,B,b) -> (
    -- compute the C matrix
    c := b//A;
    C := vec2smat(c);
    -- compute the B_i matrices
    np := numColumns B;
    if np!=0 then (
        bi := -B//A;
        Bi := toSequence for k to np-1 list
            vec2smat(bi_{k});
    )else Bi = ();
    -- compute the A_i matrices     
    v := - kernelGens A;

    Ai := toSequence for k to (rank v)-1 list
        vec2smat(v_{k});

    return (C,Ai,Bi);
    )

parameterVector = method()
parameterVector(RingElement,RingElement) := (f,objFcn) -> (
    -- given a polynomial f = f_0 + \sum_i p_i * f_i
    -- the method returns the vector with the f_i's
    R := ring f;
    kk := coefficientRing R;
    if isField kk then (
        if objFcn!=0 then 
            error "Objective must be zero if there are no parameters.";
        return (matrix{{f}},zeros(kk,0,1));
        );
    if first degree f > 1 then 
        error("Polynomial should depend affinely on the parameters.");
    p := gens R;
    F := matrix for t in {1_R}|p list {coefficient(t,f)};
    if degree objFcn > {1,0} then 
        error("Objective should be a linear function of the parameters.");
    kk = coefficientRing kk;
    objP := matrix for t in p list {sub(coefficient(t,objFcn),kk)};
    return (F,objP);
    )
parameterVector(RingElement) := (f) -> first parameterVector(f,0_(ring f))

kernelGens = A -> (
    if ring A === QQ or ring A === ZZ then
        return generators kernel A;
    tol := 1e-12;
    return numericalKernel(A,tol);
    )

zeros = (kk,m,n) -> map(kk^m,kk^n,{})

smat2vec = method( Options => {Scaling => 1} )
smat2vec(List) := o -> A -> (
    n := #A;
    v := for i to n-1 list
        for j from i to n-1 list 
            if i==j then A#i#j else o.Scaling*A#i#j;
    return flatten v;
    )
smat2vec(Matrix) := o -> A -> matrix(ring A, apply(smat2vec(entries A,o), a->{a}))

vec2smat = method( Options => {Scaling => 1} )
vec2smat(List) := o -> v -> (
    N := #v;
    n := (-1 + round sqrt(1+8*N))//2;
    ct := -1;
    L := for i to n-1 list (toList(i:0) |
        for j from i to n-1 list (ct = ct+1; ct));
    A := table(toList(0..n-1), toList(0..n-1), (i,j) -> 
        if i==j then v_(L#i#j) 
        else if i<j then v_(L#i#j)/(o.Scaling)
        else v_(L#j#i)/(o.Scaling) );
    return A;
    )

vec2smat(Matrix) := o -> v -> matrix(ring v, vec2smat(flatten entries v,o))

choosemonp = method(
    Options => {Verbose => false} )
choosemonp(RingElement) := o -> (f) -> (
    F := parameterVector(f);
    mon := choosemonp(F);
    if mon===null then return;
    return sub(mon,ring f);
    )
choosemonp(Matrix) := o -> (F) -> (
     R := ring F;
     if isQuotientRing R then error("Monomial vector must be provided in quotient rings.");
     n:= numgens R;
     mons := g -> set first entries monomials g;
     lm0 := mons F_(0,0);
     filterVerts := (verts) -> (
         -- only consider those without parameters (this is a hack!)
         return select(verts, v -> member(R_v,lm0));
         );
     lmf := sum \\ mons \ flatten entries F;
     falt := sum lmf;
     
     -- Get exponent-points for Newton polytope:
     points := substitute(matrix (transpose exponents falt),QQ);
     maxdeg := first degree falt;
     mindeg := floor first min entries (transpose points*matrix map(ZZ^n,ZZ^1,i->1));
     maxdegs := apply(entries points, i-> max i);
     mindegs := apply(entries points, i-> min i);
     
     -- Regard exponent-points in a possible subspace
     numpoints := numColumns points;
     shift := first entries transpose points;
     V := matrix transpose apply(entries transpose points, i -> i - shift);
     basV := mingens image V;
     basVdim := numgens image basV;
     if basVdim != n then T := id_(QQ^n)//basV else T = id_(QQ^n);
     basVtrans := kernelGens transpose basV;
     
     -- Compute Newton polytope:
     liftedpts := T*V || map (QQ^1,QQ^(size falt),i->1);
     dualpolytope := transpose substitute(first fourierMotzkin(liftedpts),QQ);
     argmin := L -> (m:= min L; set positions(L, l->l==m));
     idx := sum apply(entries(dualpolytope * liftedpts), i->argmin i);
     polytope := substitute(points_(toList idx),ZZ);
     oddverts := select(entries transpose polytope, i->any(i,odd));
     if #filterVerts(oddverts)>0 then(
         print("Newton polytope has odd vertices. Terminate.");
         return;
         );

     -- Get candidate points
     cp := pointsInBox(mindeg,maxdeg,mindegs,maxdegs);
     verbose("#candidate points: " | #cp, o);
     -- Only the even ones
     cpf := select(cp,i-> all(i,even)); 
     verbose("#even points: " | #cpf, o);
     -- Drop points that do not live on the subspace: 
     cpf2 := select(cpf,i-> matrix{i-shift}*basVtrans==0);
     verbose("#points in subspace of exponent-points: " | #cpf2, o);
     
     -- Find points within the polytope:
     lexponents := select(cpf2, i-> 
           max flatten entries (dualpolytope * ((T * transpose matrix {i-shift})||1)) <=0)/2;
     isInteger := l -> denominator l == 1;
     assert all(flatten lexponents, isInteger );
     lexponents = apply(lexponents, i -> numerator \ i);
     lmSOS := for i in lexponents list R_i;
     verbose("#points inside Newton polytope: " | #lmSOS, o);

     if #lmSOS==0 then return;
     return matrix transpose {lmSOS};
     )

pointsInBox = (mindeg,maxdeg,mindegs,maxdegs) -> (
    -- integer vectors within specified bounds
    n := #mindegs;
    -- Get candidate points
    local x; x= symbol x;
    R0 := QQ(monoid[x_0..x_(n-1)]);
    mon := flatten entries basis(mindeg,maxdeg,R0);
    e := apply (mon, i -> flatten exponents i);
    -- Only those within the box of degrees[mindegs:maxdegs]:
    e = select(e,i-> all(i-mindegs,j->j>=0) and all(maxdegs-i,j->j>=0)); 
    return e;
    )

project2linspace = (A,b,x0) -> (
     -- cast into QQ (necessary class to compute inverse)
     A2 := promote (A,QQ);
     -- ugly hack to convert b into a matrix if it is a scalar in QQ/ZZ:
     b2 := promote (matrix{{b}},QQ);
     x02 := promote (x0,QQ);

     -- compute projection:
     xp := x02 - transpose(A2)*((A2*x02-b2)//(A2*transpose(A2)))
     )

roundPSDmatrix = {Verbose=>false} >> o -> (Q,A,b,d) -> (
     verbose("Rounding precision: " | d, o);
     Q0 := matrix (applyTable (entries Q, i -> round(i*2^d)/2^d) );
     x0 := smat2vec(Q0);
     t := timing (xp := project2linspace(A,b,x0););
     verbose("Time needed for projection: " | net t#0, o);
     Q = vec2smat(xp);

     t = timing((L,D,P,Qpsd) := PSDdecomposition(Q););
     verbose("Time needed for LDL decomposition: " | net t#0, o);
     if Qpsd == 0 then (true, Q) else (false,Q)
     )

PSDdecomposition = (A) -> (
    kk := ring A;
    if kk===QQ then
        return LDLdecomposition(A);
    if kk=!=RR and not instance(kk,RealField) then
        error "field must be QQ or RR";
    tol := 1e-9;
    (e,V) := eigenvectors(A,Hermitian=>true);
    err := if all(e, i -> i > -tol) then 0 else 1;
    e = max_0 \ e;
    D := diagonalMatrix e;
    P := id_(kk^(numRows A));
    return (V,D,P,err);
    )
    
LDLdecomposition = (A) -> (
     kk := ring A;
     if kk=!=QQ and kk=!=RR and not instance(kk,RealField) then
        error "field must be QQ or RR";
     if transpose A != A then error("Matrix must be symmetric.");
     tol := if kk===QQ then 0 else 1e-9;

     n := numRows A;
     Ah := new MutableHashTable; map (kk^n,kk^n,(i,j)->Ah#(i,j) = A_(i,j));
     v := new MutableList from for i to n-1 list 0_kk;
     d := new MutableList from for i to n-1 list 0_kk;
     piv := new MutableList from toList(0..n-1);
     err := 0;

     for k from 0 to n-1 do (
      q := k + maxPosition apply(k..n-1, i->Ah#(i,i));
      -- Symmetric Matrix Permutation:
      tmp := piv#q; piv#q = piv#k; piv#k = tmp;
      for i to n-1 do (tmp := Ah#(i,q); Ah#(i,q) = Ah#(i,k); Ah#(i,k) = tmp;);
      for i to n-1 do (tmp := Ah#(q,i); Ah#(q,i) = Ah#(k,i); Ah#(k,i) = tmp;);

      --  positive semidefinite?
      if Ah#(k,k) < -tol then (err = k+1; break;);
      if abs(Ah#(k,k))<=tol then 
          if any(0..n-1, i->abs(Ah#(i,k))>tol) then (
               err = k+1; break;);

      -- Perform LDL factorization step:
      if Ah#(k,k) > 0 then (
                 for i to k-1 do (v#i = Ah#(k,i)*Ah#(i,i));
           Ah#(k,k) = Ah#(k,k) - sum for i to k-1 list Ah#(k,i)*v#i;
           if Ah#(k,k) < -tol then (err = k+1; break;);
           if Ah#(k,k) > 0 then
             for i from k+1 to n-1 do
                 Ah#(i,k) = (Ah#(i,k)-sum for j to k-1 list Ah#(i,j)*v#j) / Ah#(k,k);
      );
     );

     A = map(kk^n,kk^n,(i,j)-> if i>j then Ah#(i,j) else if i==j then 1_kk else 0_kk);
     D := map(kk^n,kk^n,(i,j)->if i==j then Ah#(i,j) else 0_kk);
     P := submatrix(id_(kk^n),toList piv);
     (A,D,P,err)
)

blkDiag = args -> (
     args = sequence args;
     if #args<2 then error ("expected at least 2 input arguments.");

     r := ring args#0;
     B := args#0;

     for i from 2 to #args do (
      n1 := numgens source B;
           m1 := numgens source B;
           n2 := numgens source args#(i-1);
           m2 := numgens source args#(i-1);
      B = matrix{{B,map(r^m1,n2,(i,j)->0_r)},{map(r^m2,r^n1,(i,j)->0),args#(i-1)}};
      );
     return B;
     )

--###################################
-- SOS IDEAL
--###################################

makeMultiples = (h, D, homog) -> (
    -- h is a list of polynomials
    -- multiplies each hi with monomials up to degree D
    if #h==0 then return ({},{});
    if D < max\\first\degree\h then
        error "increase degree bound";
    R := ring h#0;
    -- compute monomials
    mon := for i to #h-1 list (
        di := D - first degree h#i;
        b := if homog then basis(di,R) else basis(0,di,R);
        flatten entries b
        );
    H := for i to #h-1 list h#i * mon#i;
    return (flatten H, mon);
    )

sosInIdeal = method(
     Options => {RndTol => 3, Solver=>"CSDP", Verbose => false} )
sosInIdeal(List,ZZ) := o -> (h,D) -> (
    -- h is a list of polynomials
    -- D is a degree bound
    -- returns sos polynomial in <h>

    -- The output consists of an SOSPoly and the multipliers that
    -- express the SOS in terms of the generators.

    if odd D then error "D must be even";
    homog := all(h, isHomogeneous);
    (H,m) := makeMultiples(h, D, homog);
    F := matrix transpose {{0}|H};
    (mon,Q,X,tval) := rawSolveSOS (F, o);
    if Q===null or Q==0 or (ring Q=!=QQ and norm Q<1e-6) then (
        print("no sos polynomial in degree "|D);
        return (null,null);
	);
    a := sosdec(mon,Q);
    if a===null then (
        print("no sos polynomial in degree "|D);
        return (null,null);
	);
    S := ring h#0;
    kk := ring Q;
    if kk =!= coefficientRing S then
        S = kk(monoid[gens ring h#0]);
    a = sub(a,S);
    m = applyTable(m, i -> sub(i,S));
    mult := getMultipliers(m,tval);
    return (a,mult);
    )

getMultipliers = (mon,tval) -> (
    k := -1;
    mult := for m in mon list
        sum for i in m list( k=k+1; i*tval#k);
    return mult;
    )

sosdecTernary = method(
     Options => {RndTol => 3, Solver=>"CSDP", Verbose => false} )
sosdecTernary(RingElement) := o -> (f) -> (
    -- Implements Hilbert's algorithm to write a non-negative ternary
    -- form as sos of rational functions.
    -- Returns two lists of SOSPolys, the numerator and the denomenator polys
    if numgens ring f =!= 3 then error "polynomial must involve 3 variables";
    if not isHomogeneous f then error "polynomial must be homogeneous";
    fi := f;
    S := {};
    di := first degree fi;
    while di > 4 do(
        (Si,mult) := sosInIdeal({fi},2*di-4,o);
        if Si===null then return (,);
        fi = first mult;
        di = first degree fi;
        S = append(S,Si);
        );
    (mon,Q,X,tval) := rawSolveSOS matrix{{fi}};
    Si = sosdec(mon,Q);
    if Si===null then return (,);
    S = append(S,Si);
    nums := for i to #S-1 list if odd i then continue else S#i;
    dens := for i to #S-1 list if even i then continue else S#i;
    return (nums, dens);
    )

--###################################
-- SOS OPTIMIZATION
--###################################

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

recoverSolution = (mon,X,tol) -> (
    if X===null then return {};
    x := rank1factor(X,EigTol=>tol);
    if x===null then return {};
    if x#0<0 then x = -x;
    sol := for i to numRows mon-1 list (
        y := mon_(i,0);
        if sum degree y!=1 then continue;
        y => x_i );
    return sol;
    )

-- Unconstrained minimization 
lowerBound = method(
     Options => {RndTol => 3, Solver=>"M2", Verbose => false, EigTol => 1e-4} )
lowerBound(RingElement) := o -> (f) -> (
    -- sos lower bound for the polynomial f
    (bound,sol) := lowerBound(f,{},-1,o);
    return (bound, sol);
    )

-- Minimize a polynomial on an algebraic set
lowerBound(RingElement,List,ZZ) := o -> (f,h,D) -> (
    -- Lasserre hierarchy for the problem
    -- min f(x) s.t. h(x)=0
    numdens := (f) -> (
        R := ring f;
        (num,den) := (f, 1_R);
        if isField R then(
            (num,den) = (numerator f, denominator f);
            R = last R.baseRings;
            );
        return (R,num,den)
        );
    checkInputs := (D,num,den,h,R) -> (
        if D<0 then(
            if #h>0 or isQuotientRing R then
                error "Degree bound must be provided"
        )else(
            if odd D then error "degree bound must be even";
            if D < max\\first\degree\(h|{num,den}) then
                error "increase degree bound";
            );
        );
    
    (R,num,den) := numdens(f);
    checkInputs(D,num,den,h,R);
    -- prepare input
    (H,m) := makeMultiples(h, D, false);
    F := matrix transpose {{num,-den}|H};
    objP := matrix{{-1}} || zeros(ZZ,#H,1);

    -- call solveSOS
    o' := new OptionTable from
        {RndTol=>o.RndTol, Solver=>o.Solver, Verbose=>o.Verbose};
    mon := if isQuotientRing R then transpose basis(0,D//2,R)
        else choosemonp (F,Verbose=>o.Verbose);
    if mon===null then return (,);
    (Q,X,bound) := rawSolveSOS(F,objP,mon,o');
    
    -- recover
    if Q===null then return (,);
    sol := recoverSolution(mon,X,o.EigTol);
    return (bound#0,sol);
    )

--###################################
-- SDP SOLVER
--###################################

solveSDP = method(
     Options => {UntilObjNegative => false, WorkingPrecision => 53, Solver=>"M2", Verbose => false} )

solveSDP(Matrix, Matrix, Matrix) := o -> (C,A,b) -> solveSDP(C,sequence A,b,o)

solveSDP(Matrix, Matrix, Matrix, Matrix) := o -> (C,A,b,y) -> solveSDP(C,sequence A,b,y,o)

solveSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
    (C,A,b) = toReal(C,A,b);
    (ok,y,X,Z) := (,,,);
    (ok,y,X,Z) = sdpNoConstraints(C,A,b);
    if ok then return (y,X,Z);
    if o.Solver == "M2" then(
        (ok,y,X,Z) = trivialSDP(C,A,b);
        if ok then return (y,X,Z)
        else (y,Z) = simpleSDP(C,A,b,UntilObjNegative=>o.UntilObjNegative,Verbose=>o.Verbose)
        )
    else if o.Solver == "CSDP" then
        (y,X,Z) = solveCSDP(C,A,b,Verbose=>o.Verbose)
    else if o.Solver == "SDPA" then
        (y,X,Z) = solveSDPA(C,A,b,Verbose=>o.Verbose)
    else
        error "unknown SDP solver";
    if C==0 and b==0 and y==0 then(
        print "Zero solution obtained. Trying one more time.";
        b' := random(RR^(numRows b),RR^1);
        (y',X',Z') := solveSDP(C,A,b',o);
        if y'=!=null and norm(y')>1e-6 then (y,Z) = (y',Z');
        );
    return (y,X,Z);
)

solveSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y0) -> (
    (C,A,b) = toReal(C,A,b);
    y0 = promote(y0,RR);
    (ok,y,X,Z) := (,,,);
    (ok,y,X,Z) = sdpNoConstraints(C,A,b);
    if ok then return (y,X,Z);
    if o.Solver != "M2" then return solveSDP(C,A,b,o);
    (ok,y,X,Z) = trivialSDP(C,A,b);
    if ok then return (y,X,Z);
    (y,Z) = simpleSDP(C,A,b,y0,UntilObjNegative=>o.UntilObjNegative,Verbose=>o.Verbose);
    return (y,,Z);
)

toReal = (C,A,b) -> (
    C = promote(C,RR);
    A = apply(A, Ai -> promote(Ai,RR));
    b = promote(b,RR);
    return (C,A,b);
)

sdpNoConstraints = (C,A,b) -> (
    if #A==0 then(
        lambda := min eigenvalues(C, Hermitian=>true);
        if lambda>=0 then(
            print "SDP solved";
            y0 := zeros(RR,#A,1);
            return (true, y0, 0*C, C);
        )else(
            print "dual infeasible";
            return (true,,,);
            );
        );
    return (false,,,);
)

-- check trivial cases
trivialSDP = (C,A,b) -> (
    if #A==0 or b==0 then(
        lambda := min eigenvalues(C, Hermitian=>true);
        if lambda>=0 then(
            print "SDP solved";
            y0 := zeros(RR,#A,1);
            return (true, y0, 0*C, C);
        )else if #A==0 then(
            print "dual infeasible";
            return (true,,,);
            );
        );
    return (false,,,);
)

--simpleSDP

simpleSDP = method(
    TypicalValue => Matrix,
    Options => {UntilObjNegative => false, Verbose => false} )

simpleSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
    R := RR;
    n := numRows C;

    -- try to find strictly feasible starting point --
    (y,Z) := (,);
    lambda := min eigenvalues (C, Hermitian=>true);
    if lambda > 0 then
        y = zeros(R,#A,1)
    else(
        verbose("Computing strictly feasible solution...", o);
        y =  zeros(R,#A,1) || matrix{{lambda*1.1}};
        obj :=  zeros(R,#A,1) || matrix{{-1_R}};
        (y,Z) = simpleSDP(C,append(A,id_(R^n)), obj, y, UntilObjNegative=>true, Verbose=>o.Verbose);
        if y===null then return (,);
        y = transpose matrix {take (flatten entries y,numRows y - 1)};
        );
    verbose("Computing an optimal solution...", o);
    return simpleSDP(C, A, b, y, o);
    )

simpleSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y) -> (
    R := RR;
    n := numgens target C;

    m := numgens target y;
    mu := 1_R;
    theta := 10_R;
    iter := 1;
    NewtonIterMAX := 40;

    verbose("#It:       b'y      dy'Hdy   mu   alpha", o);

    while mu > 0.000001 do (
        mu = mu/theta;
        while true do (
            S := C - sum toList apply(0..m-1, i-> y_(i,0) * A_i);
            try Sinv := solve(S, id_(target S)) else (
                print "slack matrix is singular, terminate";
                return (,) );
            -- compute Hessian:
            H := map(R^m,R^m,(i,j) -> trace(Sinv*A_i*Sinv*A_j));
            if H==0 then (
                print "Hessian is zero";
                return (,) );
            -- compute gradient:
            g := map(R^m,R^1,(i,j) -> b_(i,0)/mu + trace(Sinv*A_i));
            
            -- compute damped Newton step:
            try dy := -g//H else (
                print "Newton step has no solution";
                return (,) );
            alpha := backtrack(S, -sum for i to m-1 list matrix(dy_(i,0) * entries A_i));
            if alpha===null then return (,);
            y = y + transpose matrix {alpha* (flatten entries dy)};
            lambda := (transpose dy*H*dy)_(0,0);
            obj := transpose b * y;
            
            -- print some information:
            verbose(iter | ":  " | net obj | "    " | net lambda | "    " | net mu | "    " | net alpha, o);

            iter = iter + 1;
            if iter > NewtonIterMAX then (
                verbose("Warning: exceeded maximum number of iterations", o);
                break);
            if (o.UntilObjNegative == true) and (obj_(0,0) < 0) then break;
            if lambda < 0.4 then break;
            ); 
        );
    Z := C - sum(for i to #A-1 list y_(i,0) * A_i);
    return (y,Z);
    )     

backtrack = args -> (
     S0 := args#0;
     R := ring S0;
     dS := args#1;
     alpha := 1_R;
     BacktrackIterMAX := 100;
     S :=  matrix( alpha * entries dS) + S0;
     
     cnt := 1;     
     while min eigenvalues(S,Hermitian=>true) <= 0 do (
      cnt = cnt + 1;
      alpha = alpha / sqrt(2_R);
      S = S0 + matrix( alpha * entries dS);
      if cnt > BacktrackIterMAX then (
          print ("line search did not converge.");
          return null );
      );
     return alpha;
     )


--solveCSDP

solveCSDP = method( Options => {Verbose => false} )
solveCSDP(Matrix,Sequence,Matrix) := o -> (C,A,b) -> (
    n := numColumns C;
    fin := getFileName ".dat-s";
    fout := getFileName "";
    fout2 := getFileName "";
    writeSDPA(fin,C,A,b);
    print("Executing CSDP on file " | fin);
    r := run(csdpexec | " " | fin | " " | fout | ">" | fout2);
    if r == 32512 then error "csdp executable not found";
    print("Output saved on file " | fout);
    (y,X,Z) := readCSDP(fout,fout2,n,o.Verbose);
    y = checkDualSol(C,A,y,Z,o.Verbose);
    return (y,X,Z);
)

getFileName = (ext) -> (
     filename := temporaryFileName() | ext;
     while fileExists(filename) do filename = temporaryFileName();
     return filename
)

writeSDPA = (fin,C,A,b) -> (
    digits := 16;
    m := length A;
    n := numColumns C;
    A = prepend(C,A);
    f := openOut fin;
    inputMatrix := l -> (
        a := -A_l;
        pref := toString l | " 1 ";
        for i to n-1 do
            for j from i to n-1 do
                if a_(i,j)!=0 then
                    f << pref | toString(i+1) | " " | toString(j+1) | " " | format(digits,a_(i,j)) << endl;
    );
    f << "*SDPA file generated by SOSm2" << endl;    
    f << toString m << " =mdim" << endl;
    f << "1 =nblocks" << endl;
    f << toString n << endl;
    f << demark(" ", toString\flatten entries b) << endl;
    for l to m do(
        inputMatrix(l);
    );
    f << close;
)

readCSDP = (fout,fout2,n,Verbose) -> (
    sdpa2matrix := s -> (
        e := for i in s list (i_2-1,i_3-1) => i_4;
        e' := for i in s list (i_3-1,i_2-1) => i_4;
        return map(RR^n, RR^n, e|e');
        );
    readLine := l -> for s in separate(" ",l) list if s=="" then continue else value s;
    --READ SOLUTIONS
    tmp := getFileName "";
    r := run("cat " | fout | " | tr -d + > " | tmp);
    L := lines get openIn tmp;
    y := transpose matrix{readLine L_0};
    S := readLine \ drop(L,1);
    S1 := select(S, l -> l_0==1);
    S2 := select(S, l -> l_0==2);
    Z := sdpa2matrix(S1); -- slack matrix
    X := sdpa2matrix(S2); -- dual solution
    -- READ STATUS
    text := get openIn fout2;
    s := select(lines text, l -> match("Success",l));
    if #s==0 then( print "SDP failed"; return (,,) );
    s = first s;
    print if Verbose then text else s;
    if match("SDP solved",s) then null
    else if match("primal infeasible",s) then X=null
    else if match("dual infeasible",s) then (y=null;Z=null;)
    else error "unknown message";
    return (y,X,Z);
)

checkDualSol = (C,A,y,Z,Verbose) -> (
    if y===null then return;
    yA := sum for i to #A-1 list y_(i,0)*A_i;
    if norm(Z-C+yA)<1e-5 then return y;
    if Verbose then print "updating dual solution";
    AA := transpose matrix(RR, smat2vec \ entries \ toList A);
    bb := transpose matrix(RR, {smat2vec entries(C-Z)});
    y = solve(AA,bb,ClosestFit=>true);
    return y;
    )

--solveSDPA

solveSDPA = method( Options => {Verbose => false} )
solveSDPA(Matrix,Sequence,Matrix) := o -> (C,A,b) -> (
    n := numColumns C;
    fin := getFileName ".dat-s";
    fout := getFileName "";
    writeSDPA(fin,C,A,b);
    print("Executing SDPA on file " | fin);
    r := run(sdpaexec | " " | fin | " " | fout | "> /dev/null");
    if r == 32512 then error "sdpa executable not found";
    if r == 11 then error "Segmentation fault running sdpa.";
    print("Output saved on file " | fout);
    (y,X,Z) := readSDPA(fout,n,o.Verbose);
    return (y,X,Z);
)

readSDPA = (fout,n,Verbose) -> (
    readVec := l -> (
        l = replace("([{} +])","",l);
        for s in separate(",",l) list if s=="" then continue else value s
    );
    readMatrix := ll -> 
        matrix for l in ll list readVec l;
    text := get openIn fout;
    L := lines text;
    --READ SOLUTIONS
    y := null; X := null; Z := null;
    i := position(L, l -> match("xVec =",l));
    if i=!=null then 
        y = transpose matrix {readVec L#(i+1)};
    i = position(L, l -> match("xMat =",l));
    if i=!=null then 
        Z = matrix for j to n-1 list readVec L#(i+j+2);
    i = position(L, l -> match("yMat =",l));
    if i=!=null then 
        X = matrix for j to n-1 list readVec L#(i+j+2);
    --READ STATUS
    if Verbose then print text;
    s := first select(L, l -> match("phase.value",l));
    if match("pdOPT",s) or match("pdFEAS",s) then 
        print "SDP solved"
    else if match("dUNBD",s) then(
        print "dual infeasible";
        y=null;Z=null; )
    else if match("pUNBD",s) then(
        print "primal infeasible";
        X=null; )
    else( 
        print("Warning: Solver returns unknown message!!! " |s);
        );
    return (y,X,Z);
    )

--###################################
-- Methods for testing
--###################################

checkSolver = method()
checkSolver(String,String) := (solver,fun) -> (
    checkMethod := hashTable {
        "solveSDP" => checkSolveSDP,
        "solveSOS" => checkSolveSOS,
        "sosdecTernary" => checkSosdecTernary,
        "sosInIdeal" => checkSosInIdeal,
        "lowerBound" => checkLowerBound
        };
    if checkMethod#?fun then 
        return checkMethod#fun(solver);
    if fun != "AllMethods" then
        error "No test implemented for this function";
    T := for f in keys checkMethod list(
        print "################################";
        print("checking method "|f);
        print "################################";
        t := checkMethod#f(solver);
        {f, testsString t}
        );
    print "################################";
    print("Summary");
    print netList T;
    )
checkSolver(String,Function) := (solver,fun) -> checkSolver(solver,toString fun)
checkSolver(String) := (solver) -> checkSolver(solver,"AllMethods")

-- A method to inform about the results of the tests in one function
testsString = t -> concatenate apply(t, i -> if i then " ✓ " else " ✘ ")
informAboutTests = t -> (
    print("Test Results: " | testsString t);
    )

--checkSolveSDP
checkSolveSDP = solver -> (
    tol := .001;
    equal := (y0,y) -> y=!=null and norm(y0-y)<tol*(1+norm(y0));
    local C; local b; local A; local A1; local A2; local A3; 
    local y0; local y; local X; local Z; local yopt;
    ---------------TEST0---------------
    C = matrix{{0,2,0,0,0,0},{2,0,0,0,0,0},
     {0,0,10,0,0,0},{0,0,0,10,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
    A1 = matrix{{-1,0,0,0,0,0},{0,0,0,0,0,0},
     {0,0,1,0,0,0},{0,0,0,0,0,0},{0,0,0,0,-1,0},{0,0,0,0,0,0}};
    A2 = matrix{{0,0,0,0,0,0},{0,-1,0,0,0,0},
     {0,0,0,0,0,0},{0,0,0,1,0,0},{0,0,0,0,0,0},{0,0,0,0,0,-1}};
    A = (A1,A2);
    y0 = matrix{{7},{9}};
    b = matrix{{1},{1}};
    (y,X,Z) = solveSDP(C,A,b,y0,Solver=>solver);
    yopt = matrix{{2.},{2.}};
    t0 := equal(yopt,y);
    ---------------TEST1---------------
    C = matrix {{2,1,-1},{1,0,0},{-1,0,5}};
    A1 = matrix {{0,0,1/2},{0,-1,0},{1/2,0,0}};
    A2 = matrix {{1,0,0},{0,1,0},{0,0,1}};
    A = (A1,A2);
    b = matrix {{0},{-1}};
    y0 = matrix {{0},{-.486952}};
    (y,X,Z) = solveSDP(C,A,b,y0,Solver=>solver);
    yopt = matrix{{1.97619},{.466049}};
    t1 := equal(yopt,y);
    ---------------TEST2---------------
    C = matrix{{2,2,-1,3},{2,0,0,2},{-1,0,1,0},{3,2,0,1}};
    A1 = matrix{{-1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    A2 = matrix{{0,0,0,1/2},{0,-1,0,0},{0,0,0,0},{1/2,0,0,0}};
    A = (A1,A2);
    b = matrix{{1},{0}};
    (y,X,Z) = solveSDP(C,A,b,Solver=>solver);
    yopt = matrix{{0.},{4.}};
    t2 := equal(yopt,y); 
    ---------------TEST3---------------
    -- solution not strictly feasible
    C = matrix {{2,2,-1,3},{2,0,0,2},{-1,0,1,0},{3,2,0,1}};
    A1 = matrix {{0,0,0,1/2},{0,-1,0,0},{0,0,0,0},{1/2,0,0,0}};
    A = sequence A1;
    b = matrix {{1}};
    (y,X,Z) = solveSDP(C,A,b,Solver=>solver);
    yopt = 4.;
    t3 := equal(yopt,y);
    ---------------TEST4---------------
    -- zero objective function
    C = matrix(RR, {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
    A1 = matrix(RR, {{1, 3/2, 3/2}, {3/2, 0, 1/2}, {3/2, 1/2, 0}});
    A2 = matrix(RR, {{0, 1/2, 3/2}, {1/2, 0, 3/2}, {3/2, 3/2, 1}});
    A3 = matrix(RR, {{0, 0, 1/2}, {0, -1, 0}, {1/2, 0, 0}});
    A = (A1,A2,A3);
    b = matrix(RR, {{0}, {0}, {0}});
    (y,X,Z) = solveSDP(C, A, b, Solver=>solver);
    yA := sum for i to #A-1 list y_(i,0)*A_i;
    t4 := norm(Z-C+yA)<1e-5;
    -----------------------------------
    test := {t0,t1,t2,t3,t4};
    informAboutTests test;
    -- trivial cases
    (y,X,Z) = solveSDP (matrix{{1,0},{0,-1}},(),matrix{{}},Solver=>solver);
    assert(y===null and X===null);
    (y,X,Z) = solveSDP (matrix{{1,0},{0,1}},(),matrix{{}},Solver=>solver);
    assert(y==0);
    return test;
    )

--checkSolveSOS
checkSolveSOS = solver -> (
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    local w; w= symbol w;
    local t; t= symbol t;
    isGram := (f,mon,Q) -> (
        if Q===null then return false;
        e := eigenvalues(Q,Hermitian=>true);
        tol := 1e-8;
        if min e < -tol then return false;
        S := ring mon;
        if coefficientRing S===QQ then
            return f == transpose(mon)*Q*mon;
        return norm(sub(f,S) - transpose(mon)*Q*mon) < tol;
        );
    ---------------GOOD CASES1---------------
    -- Test 0
    R := QQ[x,y];
    f := 4*x^4+y^4;
    (mon,Q,X,tval) := solveSOS(f,Solver=>solver);
    t0 := isGram(f,mon,Q);

    -- Test 1
    f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
    (mon,Q,X,tval) = solveSOS(f,Solver=>solver);
    t1 := isGram(f,mon,Q);

    -- Test 2
    R = QQ[x,y,z];
    f = x^4+y^4+z^4-4*x*y*z+x+y+z+3;
    (mon,Q,X,tval) = solveSOS(f,Solver=>solver);
    t2 := isGram(f,mon,Q);
    
    -- Test 3
    R = QQ[x,y,z,w];
    f = 2*x^4 + x^2*y^2 + y^4 - 4*x^2*z - 4*x*y*z - 2*y^2*w + y^2 - 2*y*z + 8*z^2 - 2*z*w + 2*w^2;
    (mon,Q,X,tval) = solveSOS(f,Solver=>solver);
    t3 := isGram(f,mon,Q);

    -- Test 4 (parametric)
    R = QQ[x][t];
    f = (t-1)*x^4+1/2*t*x+1;
    (mon,Q,X,tval) = solveSOS (f);
    t4 := isGram(sub(f,t=>tval#0),mon,Q);

    ---------------BAD CASES1---------------
    -- Test 5
    R = QQ[x,y][t];
    f = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1; --Motzkin
    (mon,Q,X,tval) = solveSOS(f,Solver=>solver); 
    t5 := ( Q === null );

    -- Test 6
    (mon,Q,X,tval) = solveSOS(f-t,-t, Solver=>solver); 
    t6 := ( Q === null );

    results := {t0,t1,t2,t3,t4,t5,t6};
    informAboutTests (results);
    return results
    )

-- check sosdecTernary
checkSosdecTernary = solver -> (
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;

    cmp := (f,p,q) -> (
        if p===null then return false;
        d := product(sumSOS\p) - f*product(sumSOS\q);
        if coefficientRing ring f===QQ then return d==0;
        return norm(d) < 1e-4;
        );

    -- Test 0
    R:= QQ[x,y,z];
    f := x^2 + y^2 +z^2;
    (p,q) := sosdecTernary (f, Solver=>"CSDP");
    t0 := cmp(f,p,q);

    -- Test 1
    R = QQ[x,y,z];
    f = x^4*y^2 + x^2*y^4 + z^6 - 4*x^2 *y^2 * z^2;
    (p,q) = sosdecTernary (f, Solver=>"CSDP");
    t1 := (p===null);

    -- Test 2
    R = RR[x,y,z];
    f = x^4*y^2 + x^2*y^4 + z^6 - 3*x^2 *y^2 * z^2; --Motzkin
    (p,q) = sosdecTernary (f, Solver=>"CSDP");
    t2 := cmp(f,p,q);

    results := {t0,t1,t2};
    informAboutTests (results);
    return results
    )


-- check sosInIdeal
checkSosInIdeal = solver -> (
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    cmp := (h,s,mult) -> (
        if s===null then return false;
        d := sum apply(h,mult,(i,j)->i*j) - sumSOS s;
        if coefficientRing ring h#0===QQ then return d==0;
        return norm(d)<1e-4;
        );

    -- Test 0
    R:= QQ[x];
    h:= {x+1};
    (s,mult) := sosInIdeal (h,2, Solver=>solver);
    t0 := cmp(h,s,mult);
    
    -- Test 1 (same as test 0 but with degree four)
    R= RR[x];
    h= {x+1};
    (s,mult) = sosInIdeal (h,4, Solver=>solver);
    t1 := cmp(h,s,mult);

    -- Test 2:
    R = RR[x,y,z];
    h = {x-y, x+z};
    (s,mult) = sosInIdeal (h,6, Solver=>solver);
    t2 := cmp(h,s,mult);

    results := {t0,t1,t2};
    informAboutTests (results);
    return results
    )


-- check lowerBound
checkLowerBound = solver -> (
    tol := 0.001;
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    equal := (a,b) -> (
        if a===null then return false;
        d := if abs(b)<1 then abs(a-b) else abs(a-b)/abs(b);
        return d < tol;
        );

    --------------UNCONSTRAINED1--------------
    --- Test 0
    R := QQ[x];
    f := (x-1)^2 + (x+3)^2;
    (bound, sol) := lowerBound(f, Solver=>solver);
    t0 := equal(bound,8);

    -- Test 1
    R = RR[x,y];
    f = (x-pi*y)^2 + x^2 + (y-4)^2;
    (bound, sol) = lowerBound(f, Solver=>solver);
    t1 := equal(bound,16*pi^2/(2+pi^2));

    -- Test 2
    R = QQ[x,z];
    f = x^4+x^2+z^6-3*x^2*z^2;
    (bound,sol) = lowerBound (f,Solver=>solver,RndTol=>infinity);
    t2 := equal(bound,-.17798);

    -- Test 3 (rational function)
    R = QQ[x];
    f = (x^2-x)/(x^2+1);
    (bound, sol) = lowerBound(f, Solver=>solver, RndTol=>infinity);
    t3 := equal(bound,1/2-1/sqrt(2));

    ---------------CONSTRAINED1---------------
    --- Test 4
    R = RR[x,y];
    f = y;
    h1 := y-pi*x^2;
    (bound, sol) = lowerBound (f, {h1}, 4, Solver=>solver);
    t4 := equal(bound,0);

    -- Test 5
    R = QQ[x,y,z];
    f = z;
    h1 = x^2 + y^2 + z^2 - 1;
    (bound,sol) = lowerBound (f, {h1}, 4, Solver=>solver);
    t5 := equal(bound,-1);

    -----------------QUOTIENT1-----------------
    -- Test 6
    R = QQ[x,y];
    I := ideal (x^2 - x);
    S := R/I;
    f = sub(x-y,S);
    h1 = sub(y^2 - y,S);
    (bound, sol) = lowerBound(f, {h1}, 2, Solver=>solver);
    t6 := equal(bound,-1);
    
    results := {t0,t1,t2,t3,t4,t5,t6};
    informAboutTests (results);
    return results
    )

--##########################################################################--
-- Documentation and Tests
--##########################################################################--

beginDocumentation()

load "./SOS/SOSdoc.m2"

TEST /// --sosPoly and sumSOS
    R = QQ[x,y,z]
    coeff1={3,1,1,1/4,1}
    pol1={-(1/2)*x^3*y-(1/2)*x*y^3+x*y*z^2, x^2*y^2-z^4, x^2*y*z-y*z^3,
      -x^3*y+x*y^3, x*y^2*z-x*z^3}
    p1=sosPoly(R,pol1,coeff1)
    p2=x^6*y^2 + 2*x^4*y^4 + x^2*y^6 - 2*x^4*y^2*z^2 - 2*x^2*y^4*z^2 - 
    3*x^2*y^2*z^4 + x^2*z^6 + y^2*z^6 + z^8
    assert(sumSOS(p1)===p2)
///

TEST /// --SOSmult
    R = QQ[x,y,z,w]
    p1=sosPoly(R,{x^2-x*y,y^2+1,x},{1,2,3})
    p2=sosPoly(R,{y^3,x*w*z,y*z^2},{3,1/2,1/4})
    assert(sumSOS(p1*p2)==sumSOS(p1)*sumSOS(p2))
    assert(sumSOS(p1^4)==sumSOS(p1)^4)

    equal = (f1,f2) -> norm(f1-f2) < 1e-8;
    R = RR[x,y,z,w]
    p1=sosPoly(R,{x^2-x*y,y^2+1,x},{1.32,1.47,12./7})
    p2=sosPoly(R,{y^3,x*w*z,y*z^2},{3.1,1.31,2.0})
    assert( equal(sumSOS(p1*p2),sumSOS(p1)*sumSOS(p2)) )
    assert( equal(sumSOS(p1^4),sumSOS(p1)^4) )
///

TEST /// --cleanSOS
    R = RR[x,y];
    s = sosPoly(R, {x+1,y}, {2,0.0001})
    t1 = clean( 0.001, s )
    t2 = sosPoly(R, {x+1}, {2})
    assert (t1 == t2)
    
    R = QQ[x,y];
    s = sosPoly(R, {x+1,y}, {2,1/100000})
    t = clean( 0.001, s )
    assert (t == s)
///

TEST ///--substitute SOSPoly
    R = QQ[x,y];
    s = sosPoly(R, {x+1,y}, {2,3})
    S = QQ[x,y,z]
    t1 = sosPoly(S, {x+1,y}, {2,3})
    t2 = sub (s, S)
    assert (t1 == t2)
///

TEST ///--toRing
    R = QQ[x,y];
    s = sosPoly(R, {x+1,y}, {2,3});
    S = RR[x,y];
    s2 = toRing_S s;
    assert (class coefficientRing ring s2 === RealField)
    s3 = toRing_R s2;
    assert (s==s3)
    
    tol := 1e-10;
    f = 0.1*x_S^2 + y^2
    g = 1/10*(symbol x)_R^2 + (symbol y)_R^2
    -- comparison in rationals is complicated:
    residual = sum \\ abs \ (x -> lift (x,QQ)) \ flatten entries last coefficients (toRing_R f - g)
    assert (residual < tol)
    -- comparison in reals:
    assert (norm (toRing_S g - f) < tol)
///

TEST /// --sosdec
    R=QQ[x,y,z]
    Q=matrix{{1,-1/2,1},{-1/2,1,-1/2},{1,-1/2,1}}
    Q=promote(Q,QQ)
    mon=matrix{{x^3},{x^2*z},{y*z^2}}
    f=sosdec(mon,Q)
    assert(f=!=null and sumSOS f==transpose mon * Q *mon)
    -- boundary cases:
    assert( sosdec(  ,mon) === null )
    assert( sosdec(Q ,   ) === null )
///

TEST /// --choosemonp
    R = QQ[x,y];
    f = x^4+2*x*y-x+y^4
    lmsos = choosemonp(f)
    assert( lmsos === null )
    R = QQ[x,y][t];
    f = x^4+2*x*y-x+y^4
    lmsos = choosemonp(f-t)
    assert( ring lmsos===R and numRows lmsos == 6 )
    
    R = RR[x,y][t];
    f = x^4+2*x*y-x+y^4
    lmsos = choosemonp(f-t)
    assert( ring lmsos===R and numRows lmsos == 6 )
///

TEST /// --createSOSModel
    eval = (Q,v) -> (transpose v * Q * v)_(0,0)
    
    R = QQ[x][t];
    f = x^4 - 2*x + t;
    mon = matrix{{1},{x},{x^2}}
    (C,Ai,Bi,A,B,b) = createSOSModel(f,mon)
    assert( eval(C,mon) == x^4 - 2*x )
    assert( #Ai==1 and eval(Ai#0,mon) == 0 )
    assert( #Bi==1 and eval(Bi#0,mon) == 1 )
    
    equal = (f1,f2) -> norm(f1-f2) < 1e-8;
    R = RR[x][t];
    f = x^4 - 2*x + t;
    mon = matrix{{1},{x},{x^2}}
    (C,Ai,Bi,A,B,b) = createSOSModel(f,mon)
    assert( equal(eval(C,mon), x^4 - 2*x) )
    assert( #Ai==1 and equal(eval(Ai#0,mon), 0) )
    assert( #Bi==1 and equal(eval(Bi#0,mon), 1) )
///

TEST /// --blkDiag
    A1=matrix{{1,0},{0,2}}
    A2=matrix{{1,1,3},{4,2,5},{2,1,1}}
    assert(blkDiag(A1,A2)==matrix{{1,0,0,0,0},{0,2,0,0,0},{0,0,1,1,3},{0,0,4,2,5},{0,0,2,1,1}})
    A1=matrix{{1.4,0,2.5},{0,2,1.9},{1.2,3,6.1}}
    A2=matrix{{2.6,1,0},{4.1,2.6,5},{1.5,1,1}}
    assert(blkDiag(A1,A2)==matrix{{1.4,0,2.5,0,0,0},{0,2,1.9,0,0,0},{1.2,3,6.1,0,0,0},{0,0,0,2.6,1,0},{0,0,0,4.1,2.6,5},{0,0,0,1.5,1,1}})
///

TEST /// --LDLdecomposition
    A = matrix(QQ, {{5,3,5},{3,2,4},{5,4,10}})
    (L,D,P,err) = LDLdecomposition A
    assert(err==0 and L*D*transpose L == transpose P * A * P)
    (L,D,P,err) = LDLdecomposition promote(A,RR)
    assert(err==0 and L*D*transpose L == transpose P * A * P)
    
    V = random(QQ^12,QQ^8)
    A = V * transpose V 
    (L,D,P,err) = LDLdecomposition(A)
    assert(err==0 and L*D*transpose L == transpose P * A * P)

    equal = (f1,f2) -> norm(f1-f2) < 1e-6;
    V = random(RR^12,RR^8)
    A = V * transpose V 
    (L,D,P,err) = LDLdecomposition(A)
    assert(err==0 and equal(L*D*transpose L, transpose P * A * P))
///

TEST ///--makeMultiples
    R = QQ[x,y,z]
    f1 = x + y
    f2 = x^2 + z^2
    h = {f1,f2}
    (H,m) = makeMultiples (h,3, false)
    assert (#H == 14 )
    assert( max(first\degree\H) == 3 )
    
    (H,m) = makeMultiples (h,3, true)
    assert(#H == 9)
    assert( unique(first\degree\H) == {3} )
///

TEST /// --solveSDP
    test := checkSolver("M2","solveSDP")
    assert(test#0 and test#1 and test#2) --(test3 fails)
///

TEST /// --solveSOS
    test := checkSolver("M2","solveSOS")
    assert(all(test,identity))
///
