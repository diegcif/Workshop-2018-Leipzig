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
--Methods/Functions
    "solveSOS",
    "sosdec",
    "sumSOS",
    "blkDiag",
    "LDLdecomposition",
    "solveSDP",
    "checkSolver",
--debugging
    "createSOSModel",
    "choosemonp",
    "project2linspace",
    "getRationalSOS",
    "changeField",
--Method options
    "rndTol",
    "untilObjNegative",
    "workingPrecision",
    "Solver"
}

--##########################################################################--
-- GLOBAL VARIABLES 
--##########################################################################--

csdpexec=((options SOS).Configuration)#"CSDPexec"
sdpaexec=((options SOS).Configuration)#"SDPAexec"

--##########################################################################--
-- METHODS
--##########################################################################--

verbose = (s,o) -> if o.Verbose then print s

--###################################
-- solveSOS
--###################################

sumSOS = (g,d) -> sum for i to #g-1 list g_i^2 * d_i

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
     return (g,d);
     )     

solveSOS = method(
     Options => {rndTol => -3, Solver=>"M2", Verbose => false} )
solveSOS(RingElement,List,RingElement,List) := o -> (f,p,objFcn,bounds) -> (
    if coefficientRing ring f =!= QQ then error "polynomial must be over QQ";
    if first degree objFcn > 1 then error("Only linear objective function allowed.");
    parBounded := false;
    if #bounds==2 then (
        lB := promote(bounds#0,QQ);
        uB := promote(bounds#1,QQ);
        parBounded = true;
    )else if #bounds!=0 then 
        error "expected a list with two elements";
         
    -- build SOS model --     
    (C,Ai,Bi,A,B,b,mon,GramIndex,LinSpaceIndex) := createSOSModel(f,p,Verbose=>o.Verbose);
    if #mon==0 then return (false,,mon,);

    ndim := numRows C;
    mdim := #Ai;
    
    if #p!=0 and objFcn!=0 then (
        -- compute an optimal solution --
        verbose("Solving SOS optimization problem...", o);
        objFcnCF := coefficients objFcn;
        objFcnHT := hashTable transpose {flatten entries objFcnCF_0, flatten entries objFcnCF_1};
        obj := transpose matrix {apply(p,i->substitute(-objFcnHT#i,QQ))} || map(QQ^#Ai,QQ^1,i->0);
        
        if parBounded then (
            C = blkDiag(C,diagonalMatrix(-lB),diagonalMatrix(uB));
            Ai = apply(0..#Ai-1,i -> blkDiag(Ai_i, map(QQ^(2*#p), QQ^(2*#p), j -> 0)));
            Bi = apply(0..#Bi-1,i -> blkDiag(Bi_i, 
                map(QQ^#p,QQ^#p, (j,k) -> if j==k and j==i then 1_QQ else 0_QQ),
                map(QQ^#p,QQ^#p, (j,k) -> if j==k and j==i then -1_QQ else 0_QQ)));
        );
    )else (
        -- compute a feasible solution --
        verbose( "Solving SOS feasibility problem...", o);
        obj = map(RR^(#Ai+#Bi),RR^1,i->0);
    );
    (my,X,Q) := solveSDP(C, Bi | Ai, obj, Solver=>o.Solver, Verbose=>o.Verbose);
    if my===null then return (false,,mon,);
    y := -my;
    if parBounded then Q = Q^{0..ndim-1}_{0..ndim-1};
    pvec0 := y^(toList(0..#p-1));

    -- rational rounding --
    (ok,Qp,pVec) := roundSolution(y,Q,A,B,b,GramIndex,LinSpaceIndex,o.rndTol);
    if not ok then return (ok,Q,changeField(mon,RR),pvec0);
    return (ok,Qp,mon,pVec);
    )
solveSOS(RingElement,List,RingElement) := o -> (f,p,objFcn) -> 
    solveSOS(f,p,objFcn,{},o)
solveSOS(RingElement,List) := o -> (f,p) -> 
    solveSOS(f,p,0_(ring f),o)
solveSOS(RingElement) := o -> (f) -> 
    drop(solveSOS(f,{},o),-1)

roundSolution = {Verbose=>false} >> o -> (y,Q,A,B,b,GramIndex,LinSpaceIndex,rndTol) -> (
    -- round and project --
    Qnum := matrix applyTable(entries Q, a -> round(a*2^52)/2^52);

    dhi := 52;
    d := rndTol;
    np := numColumns B;
    
    while (d < dhi) do (
        verbose("rounding step #" | d, o);
        if np!=0 then (
           pVec := map(QQ^np,QQ^1,(i,j) -> round(y_(i,0)*2^d)/2^d);
           bPar := b - B*pVec;
           ) 
        else bPar= b;

        (ok,Qp) := getRationalSOS(Qnum,A,bPar,d,GramIndex,LinSpaceIndex,Verbose=>o.Verbose);
        if ok then break else d = d + 1;
        );
    pVec = if np!=0 then flatten entries pVec else null;
    return (ok,Qp,pVec);
    )

changeField = (f,kk) -> (
    R := ring f;
    X := gens R;
    S := kk(monoid[X]);
    phi := map(S,R);
    return phi(f);
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
     b := transpose matrix{apply (K, k-> (if Hf#?k then substitute(Hf#k,QQ) else 0))};

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
     if #oddverts>0 then(
         verbose("Newton polytope has odd vertices. Terminate.", o);
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
     
getRationalSOS = {Verbose=>false} >> o -> (Q,A,b,d,GramIndex,LinSpaceIndex) -> (
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
-- SDP SOLVER
--###################################

solveSDP = method(
     Options => {untilObjNegative => false, workingPrecision => 53, Solver=>"M2", Verbose => false} )

solveSDP(Matrix, Matrix, Matrix) := o -> (C,A,b) -> solveSDP(C,sequence A,b,o)

solveSDP(Matrix, Matrix, Matrix, Matrix) := o -> (C,A,b,y) -> solveSDP(C,sequence A,b,y,o)

solveSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
    (ok,y,X,Z) := trivialSDP(C,A,b);
    if ok then return (y,X,Z);
    if o.Solver == "M2" then
        (y,Z) = simpleSDP(C,A,b,untilObjNegative=>o.untilObjNegative,Verbose=>o.Verbose)
    else if o.Solver == "CSDP" then
        (y,X,Z) = solveCSDP(C,A,b,Verbose=>o.Verbose)
    else if o.Solver == "SDPA" then
        (y,X,Z) = solveSDPA(C,A,b,Verbose=>o.Verbose)
    else
        error "unknown algorithm";
    return (y,X,Z);
)

solveSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y0) -> (
    (ok,y,X,Z) := trivialSDP(C,A,b);
    if ok then return (y,X);
    if o.Solver != "M2" then return solveSDP(C,A,b,o);
    (y,Z) = simpleSDP(C,A,b,y0,untilObjNegative=>o.untilObjNegative,Verbose=>o.Verbose);
    return (y,,Z);
)

-- check trivial cases
trivialSDP = (C,A,b) -> (
    C = promote (C,RR);
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
    Options => {untilObjNegative => false, Verbose => false} )

simpleSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
    R := RR;
    C = promote(C,R);
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
        (y,Z) = simpleSDP(C,append(A,id_(R^n)), obj, y, untilObjNegative=>true, Verbose=>o.Verbose);
        if y===null then return (,);
        y = transpose matrix {take (flatten entries y,numRows y - 1)};   
        );
    verbose("Computing an optimal solution...", o);
    return simpleSDP(C, A, b, y, o);
    )

simpleSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y) -> (
    R := RR;
    C = promote(C,R);
    A = apply(0..#A-1, i -> promote(A_i,R));
    b = promote(b,R);
    n := numgens target C;                                    
    y = promote(y,R);

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
            Sinv :=  solve(S, id_(target S));
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
            if (o.untilObjNegative == true) and (obj_(0,0) < 0) then break;
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
        a := -sub(A_l,RR);
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
    s := first select(lines text, l -> match("Success",l));
    print if Verbose then text else s;
    if match("SDP solved",s) then null
    else if match("primal infeasible",s) then X=null
    else if match("dual infeasible",s) then (y=null;Z=null;)
    else error "unknown message";
    return (y,X,Z);
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
    else error "unknown message";
    return (y,X,Z);
    )

--checkSolver

checkSolver = (solver) -> (
    tol := .001;
    equal := (y0,y) -> y=!=null and norm(y0-y)<tol*(1+norm(y0));
    ---------------TEST0---------------
    C := matrix{{0,2,0,0,0,0},{2,0,0,0,0,0},
     {0,0,10,0,0,0},{0,0,0,10,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
    A1 := matrix{{-1,0,0,0,0,0},{0,0,0,0,0,0},
     {0,0,1,0,0,0},{0,0,0,0,0,0},{0,0,0,0,-1,0},{0,0,0,0,0,0}};
    A2 := matrix{{0,0,0,0,0,0},{0,-1,0,0,0,0},
     {0,0,0,0,0,0},{0,0,0,1,0,0},{0,0,0,0,0,0},{0,0,0,0,0,-1}};
    A := (A1,A2);
    y0 := matrix{{7},{9}};
    b := matrix{{1},{1}};
    (y,X,Z) := solveSDP(C,A,b,y0,Solver=>solver);
    yopt := matrix{{2.},{2.}};
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
    -----------------------------------
    test := {t0,t1,t2,t3};
    for i to 3 do
        if not test#i then print("test"|i|" failed");
    return test;
)

--##########################################################################--
-- Documentation and Tests
--##########################################################################--

beginDocumentation()

load "./SOS/SOSdoc.m2"

TEST /// --good cases
    R = QQ[x,y];
    f = 4*x^4+y^4;
    (ok,Q,mon) = solveSOS f
    (g,d) = sosdec(Q,mon)
    assert( f == sumSOS(g,d) )

    f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
    (ok,Q,mon) = solveSOS f
    (g,d) = sosdec(Q,mon)
    assert( f == sumSOS(g,d) )

    R = QQ[x,y,z];
    f = x^4+y^4+z^4-4*x*y*z+x+y+z+3;
    (ok,Q,mon) = solveSOS f
    (g,d) = sosdec(Q,mon)
    assert( f == sumSOS(g,d) )
    
    R = QQ[x,y,z,w];
    f = 2*x^4 + x^2*y^2 + y^4 - 4*x^2*z - 4*x*y*z - 2*y^2*w + y^2 - 2*y*z + 8*z^2 - 2*z*w + 2*w^2;
    (ok,Q,mon) = solveSOS f
    (g,d) = sosdec(Q,mon)
    assert( f == sumSOS (g,d) )

    R = QQ[x,z,t];
    f = x^4+x^2+z^6-3*x^2*z^2-t;
    (ok,Q,mon,tval) = solveSOS (f,{t},-t,rndTol=>12);
    assert( tval#0 == -729/4096 )
    assert( sub(f,t=>tval#0) == transpose(mon)*Q*mon )
///

TEST /// --bad cases
    R = QQ[x,y,t];
    f = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1 --Motzkin
    (ok,Q,mon) = solveSOS(f); 
    assert( ok == false )
    (ok,Q,mon,tval) = solveSOS(f-t,{t},-t); 
    assert( ok == false )
///

TEST /// --Newton polytope
    R = QQ[x,y,t];
    f = x^4+2*x*y-x+y^4
    (lmf,lmsos) = choosemonp(f,{}, Verbose=>true)
    assert( #lmsos == 0 )
    (lmf,lmsos) = choosemonp(f-t,{t}, Verbose=>true)
    assert( #lmsos == 6 )
///
        
TEST /// --LDL
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
    test := checkSolver("M2")
    assert(test#0 and test#1 and test#2)
///
