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
    PackageImports => {"SimpleDoc","FourierMotzkin"},
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
    "checkSolveSDP",
    "genericCombination",
    "cleanSOS",
    "sospolyIdeal",
    "lowerBound",
    "lasserreHierarchy",
    "checkLasserreHierarchy",
--debugging
    "createSOSModel",
    "choosemonp",
    "project2linspace",
    "roundPSDmatrix",
    "toRing",
--Method options
    "RndTol",
    "UntilObjNegative",
    "WorkingPrecision",
    "Solver",
    "TraceObj",
    "EigTol",
    "sos"
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

sub(SOSPoly,Ring) := (S,R) -> 
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
    if R =!= ring S' then error "different rings";
    return sosPoly(R,S#gens|S'#gens, S#coefficients|S'#coefficients);
    )

sumSOS = method()

sumSOS (List, List) := (g,d) -> sum for i to #g-1 list g_i^2 * d_i

sumSOS SOSPoly := a -> sum for i to #(a#gens)-1 list a#gens_i^2 * a#coefficients_i

cleanSOS = method()
cleanSOS(SOSPoly,Number) := (s,tol) -> (
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

sosdec = (Q,mon) -> (
    if mon===null or Q===null then return (,);
    (L,D,P,err) := LDLdecomposition(Q);
    if err != 0 then error ("Gram Matrix is not positive semidefinite");
    n := numRows Q;
    g := toList flatten entries (transpose mon * transpose inverse P * L);
    d := apply(toList(0..n-1),i->D_(i,i));
    idx := positions (d, i->i!=0);
    d = d_idx;
    g = g_idx;
    return sosPoly(ring mon,g,d);
    )

-- This is the main method to decompose a polynomial as a 
-- sum of squares using an SDP solver.
solveSOS = method(
     Options => {RndTol => -3, Solver=>"M2", Verbose => false, TraceObj => false} )
 
solveSOS(RingElement,List,RingElement,List) := o -> (f,p,objFcn,bounds) -> (
    -- f is a polynomial to decompose
    -- p is a list of variables of f that are interpreted as parameters
    -- objFcn is a linear objective function for the SDP solver
    -- bounds can be empty or contain upper and lower bounds for the parameters
    if first degree objFcn > 1 then error("Only linear objective function allowed.");
    parBounded := false;
    if #bounds==2 then (
        lB := promote(bounds#0,QQ);
        uB := promote(bounds#1,QQ);
        parBounded = true;
    )else if #bounds!=0 then 
        error "expected a list with two elements";

    kk := coefficientRing ring f;
    if kk=!=QQ then (
        S := changeRingField(QQ,ring f);
        f = toRing_S f;
        p = toRing_S \ p;
        );
         
    -- build SOS model --     
    (C,Ai,Bi,A,B,b,mon,GramIndex,LinSpaceIndex) := createSOSModel(f,p,Verbose=>o.Verbose);
    if #mon==0 then return (,mon,,);

    ndim := numRows C;
    mdim := #Ai;

    if #p!=0 and objFcn!=0 then (
        -- compute an optimal solution --
        verbose("Solving SOS optimization problem...", o);
        (objMon,objCoef) := coefficients objFcn;
        objP := transpose matrix {apply(p,
            i->sub(-coefficient(i,objFcn),QQ))};
        obj := objP || map(QQ^#Ai,QQ^1,i->0);
        
        if parBounded then (
            C = blkDiag(C,diagonalMatrix(-lB),diagonalMatrix(uB));
            Ai = apply(0..#Ai-1,i -> blkDiag(Ai_i, map(QQ^(2*#p), QQ^(2*#p), j -> 0)));
            Bi = apply(0..#Bi-1,i -> blkDiag(Bi_i, 
                map(QQ^#p,QQ^#p, (j,k) -> if j==k and j==i then 1_QQ else 0_QQ),
                map(QQ^#p,QQ^#p, (j,k) -> if j==k and j==i then -1_QQ else 0_QQ)));
        );
    )else if o.TraceObj then (
        -- compute an optimal solution --
        verbose("Solving SOS optimization problem...", o);
        obj = map(RR^(#Bi),RR^1,(i,j)-> -trace Bi#i) || map(RR^(#Ai),RR^1,(i,j)-> -trace Ai#i);
    )else (
        -- compute a feasible solution --
        verbose( "Solving SOS feasibility problem...", o);
        obj = map(RR^(#Ai+#Bi),RR^1,i->0);
    );
    (my,X,Q) := solveSDP(C, Bi | Ai, obj, Solver=>o.Solver, Verbose=>o.Verbose);
    if Q===null then return (Q,mon,X,);
    y := -my;
    if parBounded then Q = Q^{0..ndim-1}_{0..ndim-1};
    pvec0 := flatten entries y^(toList(0..#p-1));

    if o.RndTol==infinity then
        return (Q,changePolyField(RR,mon),X,pvec0);
    if kk=!=QQ then
        return (Q,changePolyField(kk,mon),X,pvec0);

    -- rational rounding --
    (ok,Qp,pVec) := roundSolution(y,Q,A,B,b,GramIndex,LinSpaceIndex,o.RndTol);
    if ok then return (Qp,mon,X,pVec);
    print "rounding failed, returning real solution";
    return (Q,changePolyField(RR,mon),X,pvec0);
    )

-- Variants of solveSOS with fewer arguments
solveSOS(RingElement,List,RingElement) := o -> (f,p,objFcn) -> 
    solveSOS(f,p,objFcn,{},o)
solveSOS(RingElement,List) := o -> (f,p) -> 
    solveSOS(f,p,0_(ring f),o)
solveSOS(RingElement) := o -> (f) -> 
    drop(solveSOS(f,{},o),-1)

changeRingField = (kk,R) -> kk(monoid[gens R])

changePolyField = (kk,f) -> toRing(changeRingField(kk,ring f), f)

toRing = (S,f) -> (
    -- maps f to ring S
    R := ring f;
    kk := coefficientRing R;
    phi := map(S,R);
    -- QQ => RR
    if kk===QQ then return phi(f);
    -- RR => QQ
    liftmon := (x) -> S_(first exponents x);
    (mon,coef) := coefficients f;
    mon = matrix {liftmon \ flatten entries mon};
    K := 2^32;
    coef = matrix(QQ, {for c in flatten entries coef list round(K*sub(c,RR))/K});
    f' := (mon*transpose coef)_(0,0);
    return f';
    )

roundSolution = {Verbose=>false} >> o -> (y,Q,A,B,b,GramIndex,LinSpaceIndex,RndTol) -> (
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

        (ok,Qp) := roundPSDmatrix(Qnum,A,bPar,d,GramIndex,LinSpaceIndex,Verbose=>o.Verbose);
        if ok then break else d = d + 1;
        );
    pVec = if np!=0 then flatten entries pVec else null;
    return (ok,Qp,pVec);
    )

createSOSModel = {Verbose=>false} >> o -> (f,p) -> (
    -- Degree and number of variables
    n := numgens ring f;
    d := (first degree f)//2;
    -- Get list of monomials for SOS decomposition
    (lmf,lm) := choosemonp (f,p,o);
    if #lm==0 then return (,,,,,,lm,,);
    
    Hm := hashTable apply(lm, toList(1..#lm), identity);
    HHm := combine(Hm,Hm,times,(j,k)-> if j>=k then (j,k) else () , join );
    HHHm := applyValues(HHm,k->pack(2,k));
    -- This is a hash table that maps monomials into pairs of indices
    -- Writes the matrix, in sparse form
    ndim := #lm; -- binomial(n+d,d);
    mdim := #HHHm; --binomial(n+2*d,2*d);
    
    -- A hash table with the coefficients (should be an easier way)
    cf := coefficients f ;
    Hf := hashTable transpose {flatten entries cf#0, flatten entries cf#1};
    K := keys HHHm;
    
    -- Linear constraints: b
    b := transpose matrix(QQ,{apply (K, k-> (if Hf#?k then substitute(Hf#k,QQ) else 0))});
    
    -- Linear constraints: Ai, Bi
    Ah := new MutableHashTable;
    Bh := new MutableHashTable;
    LinSpaceDim := floor(ndim^2/2+ndim/2);
    LinSpaceIndex := hashTable apply (flatten values HHHm, toList(0..LinSpaceDim-1),identity);
    GramIndex := applyPairs (LinSpaceIndex, (i,j)->(j,i));
    for k from 0 to #K-1 do (
	-- Set constraints for monomial K#k 
	PairsEntries := toList HHHm#(K#k) ;
      	scan(PairsEntries, p -> (
                if p_0 == p_1 then Ah#(k,LinSpaceIndex#p)=1_QQ else Ah#(k,LinSpaceIndex#p)=2_QQ;)
	    );
    -- Consider search-parameters:
    for i from 0 to #p-1 do (
	mp := K#k*p_i;
	if Hf#?mp then Bh#(k,i) = -leadCoefficient Hf#mp;
	);
    );

    A := map(QQ^#K,QQ^(LinSpaceDim),(i,j) -> if Ah#?(i,j) then Ah#(i,j) else 0);
    B := map(QQ^#K,QQ^#p,(i,j) -> if Bh#?(i,j) then Bh#(i,j) else 0);

    -- compute the C matrix
    c := b//A;
    C := map(QQ^ndim,QQ^ndim, (i,j) -> if i>=j then c_(LinSpaceIndex#{i+1,j+1},0)
      	else c_(LinSpaceIndex#{j+1,i+1},0));
    -- compute the B_i matrices
    if #p!=0 then (
	bi := -B//A;
	Bi := apply(0..#p-1, k->
	    map(QQ^ndim,QQ^ndim, (i,j) -> if i>=j then bi_(LinSpaceIndex#{i+1,j+1},k)
		else bi_(LinSpaceIndex#{j+1,i+1},k)));
      	) else Bi = ();
    -- compute the A_i matrices     
    v := - generators kernel A;

    Ai := apply(0..(rank v) - 1,k ->
	map(QQ^ndim,QQ^ndim, (i,j) -> if i>=j then v_(LinSpaceIndex#{i+1,j+1},k) 
            else v_(LinSpaceIndex#{j+1,i+1},k))); 

    (C,Ai,Bi,A,B,b,transpose matrix {lm},GramIndex,LinSpaceIndex)
    )

choosemonp = {Verbose=>false} >> o -> (f,p) -> (
     filterVerts := (verts) -> (
         -- only consider those without parameters (this is a hack!)
         R := ring f;
         p0 := apply(p,i->i=>0);
         lm0 := set(apply(flatten entries (coefficients f)_0,i->(substitute(i,p0))));
         return select(verts, v -> member(R_v,lm0));
     );
     -- Get rid of parameters in polynomial:
     X := gens ring f;
     genpos := positions(X,i->not any(p,j->j==i));
     ringf := QQ(monoid[X_genpos]);
     n := #genpos;
     p1 := apply(p,i->i=>1);
     lmf := unique(apply(flatten entries (coefficients f)_0,i->(substitute(i,p1))));
     falt := sum lmf;
     
     -- Get exponent-points for Newton polytope:
     points := substitute(matrix (transpose exponents falt)_genpos,QQ);
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
     basVtrans := mingens kernel transpose basV;
     
     -- Compute Newton polytope:
     liftedpts := T*V || map (QQ^1,QQ^(size falt),i->1);
     dualpolytope := transpose substitute(first fourierMotzkin(liftedpts),QQ);
     argmin := L -> (m:= min L; set positions(L, l->l==m));
     idx := sum apply(entries(dualpolytope * liftedpts), i->argmin i);
     polytope := substitute(points_(toList idx),ZZ);
     oddverts := select(entries transpose polytope, i->any(i,odd));
     if #filterVerts(oddverts)>0 then(
         print("Newton polytope has odd vertices. Terminate.");
         return (lmf,{});
         );

     -- Get candidate points from basis of f:
     mon := flatten apply( toList(mindeg..maxdeg), k -> flatten entries basis(k, ringf));
     cp := apply (mon, i -> flatten exponents (i));
     verbose("#candidate points: " | #cp, o);

     -- Filter candidate points:
     -- Only the even ones within the box of degrees[mindegs:maxdegs]:
     cpf := select(cp,i-> all(i,even) and all(i-mindegs,j->j>=0) and all(maxdegs-i,j->j>=0)); 
     verbose("#points (even and within box of polynomial degrees): " | #cpf, o);
     -- Drop points that do not live on the subspace: 
     cpf2 := select(cpf,i-> matrix{i-shift}*basVtrans==0);
     verbose("#points in subspace of exponent-points: " | #cpf2, o);
     
     -- Find points within the polytope:
     lexponents := select(cpf2, i-> 
           max flatten entries (dualpolytope * ((T * transpose matrix {i-shift})||1)) <=0)/2;
     lmSOS := apply(lexponents, i-> product(n,j->(
        assert (denominator i#j==1);
         (ring f)_(genpos#j)^(numerator i#j)
         )));
     verbose("#points inside Newton polytope: " | #lmSOS, o);
     
     return (lmf,lmSOS);
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
     
roundPSDmatrix = {Verbose=>false} >> o -> (Q,A,b,d,GramIndex,LinSpaceIndex) -> (
     ndim := numRows Q;
     
     verbose("Rounding precision: " | d, o);
     Q0 := matrix pack (apply(flatten entries Q, i -> round(i*2^d)/2^d),ndim);
     x0 := transpose matrix {{apply(0..numgens source A-1, i -> Q0_(toSequence (GramIndex#i-{1,1})))}};
     t := timing (xp := project2linspace(A,b,x0););
     verbose("Time needed for projection: " | net t#0, o);
     Q = map(QQ^ndim,QQ^ndim, (i,j) -> if i>=j then xp_(LinSpaceIndex#{i+1,j+1},0) 
           else xp_(LinSpaceIndex#{j+1,i+1},0));

     t = timing((L,D,P,Qpsd) := LDLdecomposition(Q););
     verbose("Time needed for LDL decomposition: " | net t#0, o);
     if Qpsd == 0 then (true, Q) else (false,Q)
     )

LDLdecomposition = (A) -> (
     kk := ring A;
     if kk=!=QQ and kk=!=RR and not instance(kk,RealField) then 
        error "field must be QQ or RR";
     if transpose A != A then error("Matrix must be symmetric.");
      
     n := numRows A;
     Ah := new MutableHashTable; map (kk^n,kk^n,(i,j)->Ah#(i,j) = A_(i,j));
     v := new MutableList from toList apply(0..n-1,i->0_kk);
     d := new MutableList from toList apply(0..n-1,i->0_kk);
     piv := new MutableList from toList(0..n-1);
     err := 0;
     
     for k from 0 to n-1 do (
      q := maxPosition apply(k..n-1, i->Ah#(i,i)); q = q + k;
      -- Symmetric Matrix Permutation:
      tmp := piv#q; piv#q = piv#k; piv#k = tmp;
      scan(0..n-1, i-> (tmp := Ah#(i,q); Ah#(i,q) = Ah#(i,k); Ah#(i,k) = tmp;));
      scan(0..n-1, i-> (tmp := Ah#(q,i); Ah#(q,i) = Ah#(k,i); Ah#(k,i) = tmp;));
           
      --  positive semidefinite?
      if Ah#(k,k) < 0 then (err = k+1; break;);
      if (Ah#(k,k)==0) and (number(apply(0..n-1,i->Ah#(i,k)),f->f!=0)!=0) then (
           err = k+1; break;);
      
      -- Perform LDL factorization step:
      if Ah#(k,k) > 0 then (
                 scan(0..k-1, i -> v#i = Ah#(k,i)*Ah#(i,i));      
           Ah#(k,k) = Ah#(k,k) - sum apply(toList(0..k-1), i -> Ah#(k,i)*v#i);
           if Ah#(k,k) < 0 then (err = k+1; break;);
           if Ah#(k,k) > 0 then (
            scan(k+1..n-1, i ->
             (Ah#(i,k) = (Ah#(i,k)-sum apply(toList(0..k-1),j->Ah#(i,j)*v#j)) 
             / Ah#(k,k);))
                   );
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

dotProd = (a, b) -> sum apply (a, b, (i,j)->i*j);
genericCombination = (h, D, homog) -> (
    -- h is a list of polynomials
    -- D is a maximumd degree
    -- returns generic combination of the
    if #h==0 then error "list of polynomials is empty";
    if D < max\\first\degree\h then
        error "increase degree bound";
    R := ring h#0;
    -- compute monomials
    p := symbol p;
    mon := for i to #h-1 list (
        di := D - first degree h#i;
        b := if homog then basis(di,R) else basis(0,di,R);
        flatten entries b
        );
    -- ring of parameters
    pvars := for i to #h-1 list
        toList( p_(i,0)..p_(i,#(mon#i)-1) );
    S := newRing (R, Variables=> gens R|flatten pvars);
    pvars = for i to #h-1 list apply (pvars#i, m->S_m);
    -- polynomial multipliers
    g := for i to #h-1 list
        dotProd ( pvars#i , apply (mon#i, m -> sub (m, S)) );
    h = apply(h, hi -> sub(hi,S));
    F := dotProd(h,g);
    return (F,flatten pvars,g);
    )

sospolyIdeal = method(
     Options => {RndTol => -3, Solver=>"CSDP", Verbose => false} )
sospolyIdeal(List,ZZ) := o -> (h,D) -> (
    -- h is a list of polynomials
    -- D is a degree bound
    -- returns sos polynomial in <h>
    if odd D then error "D must be even";
    homog := all(h, isHomogeneous);
    (f,p,mult) := genericCombination(h, D, homog);
    (Q,mon,X,tval) := solveSOS (f, p, o);
    if Q==0 or (ring Q=!=QQ and norm Q<1e-6) then (
        print("no sos polynomial in degree "|D);
        return (null,null);
	);
    a := sosdec(Q,mon);
    kk := ring Q;
    S := kk(monoid[gens ring h#0]);
    a = sub(a,S);
    h = for hi in h list sub(hi,S);
    -- get multipliers
    T := kk(monoid[gens ring f]);
    dic := for i to #p-1 list sub(p#i,T) => tval#i;
    mult = for m in mult list sub(sub(sub(m,T),dic),S);
    -- another way (using gb)
    -- mult = flatten entries (sumSOS a // gens ideal(h));
    return (a,mult);
    )

sosdecTernary = method(
     Options => {RndTol => -3, Solver=>"CSDP", Verbose => false} )
sosdecTernary(RingElement) := o -> (f) -> (
    if numgens ring f =!= 3 then error "polynomial must involve 3 variables";
    if not isHomogeneous f then error "polynomial must be homogeneous";
    fi := f;
    S := {};
    di := first degree fi;
    while di > 4 do(
        (Si,mult) := sospolyIdeal({fi},2*di-4,o);
        if Si===null then return;
        fi = mult#0;
        di = first degree fi;
        S = append(S,Si);
        );
    (Q,mon,X) := solveSOS fi;
    if Q===null then return;
    Si = sosdec(Q,mon);
    S = append(S,Si);
    nums := for i to #S-1 list if odd i then continue else S#i;
    dens := for i to #S-1 list if even i then continue else S#i;
    return (nums,dens);
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

lowerBound = method(
     Options => {RndTol => -3, Solver=>"M2", Verbose => false, EigTol => 1e-4} )

lowerBound(RingElement,List) := o -> (f,pars) -> (
    -- sos lower bound for the polynomial f
    o' := new OptionTable from 
        {RndTol=>o.RndTol, Solver=>o.Solver, Verbose=>o.Verbose};
    R := ring f;
    tvar := local t;
    newR := newRing (R, Variables=> gens R| {tvar});
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
lowerBound(RingElement) := o -> (f) -> lowerBound(f,{},o)

-- Minimize a polynomial on an algebraic set using SDP relaxations:
lasserreHierarchy = method(
     Options => {RndTol => -3, Solver=>"M2", Verbose => false, EigTol => 1e-4} )
lasserreHierarchy(RingElement,List,ZZ) := o -> (f,h,D) -> (
    -- Lasserre hierarchy for the problem
    -- min f(x) s.t. h(x)=0
    if odd D then error "D must be even";
    if D < max\\first\degree\(h|{f}) then
        error "increase degree bound";
    if all(h|{f}, isHomogeneous) then
        error "problem is homogeneous";
    R := ring f;
    (H,p,mult) := genericCombination(h, D, false);
    S := ring H;
    f = sub(f,S);
    (bound,sol) := lowerBound(f+H, p, o);
    if bound===null then return (,);
    sol = for p in sol list sub(p#0,R)=>p#1;
    return (bound,sol);
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
        error "unknown algorithm";
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
            y0 := map(RR^#A,RR^1,i->0);
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
            y0 := map(RR^#A,RR^1,i->0);
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
        y = map(R^#A,R^1,(i,j)-> 0)
    else(
        verbose("Computing strictly feasible solution...", o); 
        y =  map(R^#A,R^1,i->0) || matrix{{lambda*1.1}};
        obj :=  map(R^#A,R^1,i->0) || matrix{{-1_R}};
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
            alpha := backtrack(S, -sum toList apply(0..m-1, i -> matrix(dy_(i,0) * entries A_i)));
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
    yA := sum for i to #A-1 list y_(i,0)*A_i;
    if norm(Z-C+yA)<1e-5 then return y;
    if Verbose then print "updating dual solution";
    AA := transpose matrix(RR, smat2vec \ toList A);
    bb := transpose matrix(RR, {smat2vec(C-Z)});
    y = solve(AA,bb,ClosestFit=>true);
    return y;
    )

smat2vec = A -> (
    n := numColumns A;
    v := for i to n-1 list
        for j from i to n-1 list A_(i,j);
    return flatten v;
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

-- A method to inform about the results of the tests in one function
informAboutTests = t -> (
    for i to #t-1 do if not t#i then print("test"|i|" failed");
    )

--checkSolveSDP

checkSolveSDP = (solver) -> (
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

-- check lasserreHierarchy
checkLasserreHierarchy = solver -> (
    tol := 0.001;
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    --- Test 0
    R := QQ[x,y,z];
    f := -z;
    h1 := x^2 + y^2 + z^2 - 1;
    (minb, sol) := lasserreHierarchy (f, {h1}, 4, Solver=>solver);
    t0 := (abs(1-minb) < tol);

    --- Test 1
    R = QQ[x,y];
    f = -y;
    h1 = y-x^2;
    (minb, sol) = lasserreHierarchy (f, {h1}, 4, Solver=>solver);
    t1 := (abs(minb) < tol);
    
    results := {t0,t1};
    informAboutTests (results);
    return results
    )

--##########################################################################--
-- Documentation and Tests
--##########################################################################--

beginDocumentation()

load "./SOS/SOSdoc.m2"

TEST /// --solveSOS (good cases)
    R = QQ[x,y];
    f = 4*x^4+y^4;
    (Q,mon,X) = solveSOS f
    a = sosdec(Q,mon)
    assert( f == sumSOS a )

    f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
    (Q,mon,X) = solveSOS f
    a = sosdec(Q,mon)
    assert( f == sumSOS a )

    R = QQ[x,y,z];
    f = x^4+y^4+z^4-4*x*y*z+x+y+z+3;
    (Q,mon,X) = solveSOS f
    a = sosdec(Q,mon)
    assert( f == sumSOS a )
    
    R = QQ[x,y,z,w];
    f = 2*x^4 + x^2*y^2 + y^4 - 4*x^2*z - 4*x*y*z - 2*y^2*w + y^2 - 2*y*z + 8*z^2 - 2*z*w + 2*w^2;
    (Q,mon,X) = solveSOS f
    a = sosdec(Q,mon)
    assert( f == sumSOS a )

    R = QQ[x,z,t];
    f = x^4+x^2+z^6-3*x^2*z^2-t;
    (Q,mon,X,tval) = solveSOS (f,{t},-t,RndTol=>12);
    assert( tval#0 == -729/4096 )
    assert( sub(f,t=>tval#0) == transpose(mon)*Q*mon )
///

TEST /// --solveSOS (bad cases)
    R = QQ[x,y,t];
    f = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1 --Motzkin
    (Q,mon,X) = solveSOS(f); 
    assert( Q === null )
    (Q,mon,X,tval) = solveSOS(f-t,{t},-t); 
    assert( Q === null )
///

TEST /// --choosemonp
    R = QQ[x,y,t];
    f = x^4+2*x*y-x+y^4
    (lmf,lmsos) = choosemonp(f,{}, Verbose=>true)
    assert( #lmsos == 0 )
    (lmf,lmsos) = choosemonp(f-t,{t}, Verbose=>true)
    assert( #lmsos == 6 )
///
        
TEST /// --LDLdecomposition
--  Simple example
    A = matrix(QQ, {{5,3,5},{3,2,4},{5,4,10}})
    (L,D,P,err) = LDLdecomposition A
    assert(L*D*transpose L == transpose P * A * P)
    (L,D,P,err) = LDLdecomposition promote(A,RR)
    assert(L*D*transpose L == transpose P * A * P)
    
--  Random low-rank matrix
    V = random(QQ^12,QQ^8)
    A = V * transpose V 
    (L,D,P,err) = LDLdecomposition(A)
    assert(L*D*transpose L == transpose P * A * P)
///

TEST /// --solveSDP
    test := checkSolveSDP("M2")
    assert(test#0 and test#1 and test#2) --(test3 fails)
///

TEST /// --lasserreHierarchy
    test := checkLasserreHierarchy ("M2")
    -- lasserreHierarchy fails with the default solver
    -- assert (test#0 and test#1)
///
