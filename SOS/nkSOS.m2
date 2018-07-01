--
restart
needsPackage ("SOS", Configuration=>{"CSDPexec"=>"/Users/kaihnsa/Documents/Workshop-2018-Leipzig-SOS/SOS/CSDP/csdp"})
R=QQ[x,y,z]
sospolyMultiply = method()
    
sospolyMultiply(SOSPoly,SOSPoly):= (g1,g2)-> (
    if g1#ring =!= g2#ring then error "different rings";
    q1=for i from 0 to #g1#gens-1 list(
        for j from 0 to #g2#gens-1 list g1#gens_i*g2#gens_j);
    q2=for i from 0 to #g1#coefficients-1 list(
        for j from 0 to #g2#coefficients-1 list g1#coefficients_i*g2#coefficients_j);
    return sosPoly(g1#ring, flatten(q1),flatten(q2));
    )


--EXAMPLE
R = QQ[x,y,z]
h = x^2 + y^2 + z^2
-- the Motzkin polynomial + Co
f1 = x^4*y^2 + x^2*y^4 + z^6 - 3*x^2 *y^2 * z^2 --Motzkin
(Q1,mon1,X1) = solveSOS (f1*h, Solver=>"CSDP")
g1 = sosdec (Q1,mon1)

--Robinson
f2 = x^6 + y^6 + z^6 - (x^4*y^2 + x^2*y^4 + x^4*z^2 + x^2*z^4 + y^4*z^2 + y^2*z^4) + 3*x^2*y^2*z^2 --Robinson
(Q2,mon2,X2) = solveSOS (f2*h, Solver=>"CSDP")
g2 = sosdec (Q2,mon2)

p1=sospolyMultiply(g1,g2)

sumSOS(p1)==sumSOS(g1)*sumSOS(g2)

----xxxx-----

sospolyPower = method(
      )
sospolyPower(SOSPoly,ZZ) := (p1,D)->(
    if D==2 then return sosPoly({sumSOS(p1)},{1});
    if D<=0 then error "power should be a positive integer.";
    if odd D then error "power should be an even integer.";
    d=D//2; 
    sq=sumSOS(p1);
    p2=(sq)^d;
    return sosPoly({p2},{1});
    )
sospolyPower(g1,4)

TEST /// --SOSmult
    R = QQ[x,y,z]
    coeff1={1,2,3,4}
    pol1={x^2,x*y,y^2,x}
    p1=sosPoly(pol1,coeff1)
    S=QQ[x,y,z,w]
    coeff2={3,1/2,1/4}
    pol2={y^3,x*w*z,y*z^2}
    p2=sosPoly(S,pol2,coeff2)
    p4=sub(p1,S)
    p3=sospolyMultiply(p4,p2)
    assert(sumSOS(p3)===sumSOS(p4)*sumSOS(p2))
///    

TEST /// --SOSmult2
    R = QQ[x,y,z,w]
    coeff1={1,2,3,4}
    pol1={x^2,x*y,y^2,x}
    p1=sosPoly(pol1,coeff1)
    coeff2={3,1/2,1/4}
    pol2={y^3,x*w*z,y*z^2}
    p2=sosPoly(pol2,coeff2)
    p3=sospolyMultiply(p1,p2)
    assert(sumSOS(p3)===sumSOS(p1)*sumSOS(p2))
///    
  
  


TEST /// --SOSmult4
    R = QQ[x,y,z]
    coeff1={1,2,3,4}
    pol1={x^2,x*y,y^2,x}
    p1=sosPoly(pol1,coeff1)
    S=RR[x,y,z,w]
    coeff2={3.1,1.31,2.0}
    pol2={y^3,x*w*z,y*z^2}
    p2=sosPoly(S,pol2,coeff2)
    p4=sub(p1,S)
    p3=sospolyMultiply(p4,p2)
    assert(sumSOS(p3)===sumSOS(p4)*sumSOS(p2))
///

TEST /// --sub(SOS)1
    R = QQ[x,y,z]
    coeff1={1,2,3,4}
    pol1={x^2,x*y,y^2,x}
    p1=sosPoly(pol1,coeff1)
    S=QQ[x,y,z,w]
    p4=sub(p1,S)
    assert(p4#ring===S)
///

TEST /// --sub(SOS)2
    R = QQ[x,y,z]
    coeff1={1,2,3,4}
    pol1={x^2,x*y,y^2,x}
    p1=sosPoly(pol1,coeff1)
    p1#ring===R
    S=RR[x,y]
    p4=sub(p1,S)
    assert(p4#ring===S)
///  

TEST /// --sospoly and sum of squares
    R = QQ[x,y,z]
    coeff1={3,1,1,1/4,1}
    pol1={-(1/2)*x^3*y-(1/2)*x*y^3+x*y*z^2, x^2*y^2-z^4, x^2*y*z-y*z^3,
      -x^3*y+x*y^3, x*y^2*z-x*z^3}
    p1=sosPoly(R,pol1,coeff1)
    p2=x^6*y^2 + 2*x^4*y^4 + x^2*y^6 - 2*x^4*y^2*z^2 - 2*x^2*y^4*z^2 - 
    3*x^2*y^2*z^4 + x^2*z^6 + y^2*z^6 + z^8
    assert(sumSOS(p1)===p2)
///

TEST /// --sospoly and sum of squares
    R = QQ[x,y,z]
    coeff1={1,1,1,3/4,1}
    pol1={x^4-(1/2)*x^2*y^2-(1/2)*y^4-(1/2)*x^2*z^2+y^2*z^2-(1/2)*z^4,
       x^3*y-x*y^3, x^3*z-x*z^3, x^2*y^2-y^4-x^2*z^2+z^4, y^3*z-y*z^3}
    p1=sosPoly(R,pol1,coeff1)
    p2=x^8 - 2*x^4*y^4 + y^8 + x^4*y^2*z^2 + x^2*y^4*z^2 - 2*x^4*z^4 + x^2*y^2*z^4 - 
    2*y^4*z^4 + z^8
    assert(sumSOS(p1)===p2)
///

TEST /// --sospolypower
    R = QQ[x,y,z]
    coeff1={1,1,1,3/4,1}
    pol1={x^4-(1/2)*x^2*y^2-(1/2)*y^4-(1/2)*x^2*z^2+y^2*z^2-(1/2)*z^4,
       x^3*y-x*y^3, x^3*z-x*z^3, x^2*y^2-y^4-x^2*z^2+z^4, y^3*z-y*z^3}
    p1=sosPoly(R,pol1,coeff1)
    p2=sospolyPower(p1,2)
    q1=sumSOS(p1)
    q2=sumSOS(p2)
    assert(q2==q1^2)
///

TEST /// --sospolypower
    R = QQ[x,y,z]
    coeff1={1,1,1,3/4,1}
    pol1={x^4-(1/2)*x^2*y^2-(1/2)*y^4-(1/2)*x^2*z^2+y^2*z^2-(1/2)*z^4,
       x^3*y-x*y^3, x^3*z-x*z^3, x^2*y^2-y^4-x^2*z^2+z^4, y^3*z-y*z^3}
    p1=sosPoly(R,pol1,coeff1)
    p2=sospolyPower(p1,6)
    q1=sumSOS(p1)
    q2=sumSOS(p2)
    assert(q2==q1^6)
///

TEST /// --sospolypower
    R = QQ[x,y,z]
    coeff1={3,1,1,1/4,1}
    pol1={-(1/2)*x^3*y-(1/2)*x*y^3+x*y*z^2, x^2*y^2-z^4, x^2*y*z-y*z^3,
      -x^3*y+x*y^3, x*y^2*z-x*z^3}
    p1=sosPoly(R,pol1,coeff1)
    p2=sospolyPower(p1,8)
    q1=sumSOS(p1)
    q2=sumSOS(p2)
    assert(q2==q1^8)
///

TEST /// --cleanSOS
    R = QQ[x,y,z]
    coeff1={1/2,0,0,4,6,8}
    pol1={x^2,y^2,x*z,x*y,z^3,y^2*z}
    p1=sosPoly(R,pol1,coeff1)
    p2=cleanSOS(p1,0)
    q1=sumSOS(p1)
    q2=sumSOS(p2)
    assert(q2==q1)
///

TEST /// --toRing
    R=RR[x,y,z]
    f=x^2+2.21*x*y-4*y*z
    S=QQ[x,y,z]
    f'=toRing(S,f)
    t =toRing(R,f')
    assert(norm(f-t)<1e-10)
    assert(class f'===S)
///


TEST /// --toRing
    R=QQ[x,y,z]
    f=x^2+2/9*x*y-4/21*y*z
    S=RR[x,y,z]
    f'=toRing(S,f)
    t =toRing(R,f')
    (mon,coeff)=coefficients(f-t)
    t1=flatten entries coeff
    t2=sum(for i in t1 list i^2)
    assert(1/10^(18)-t2>0)
    assert(class f'===S)
///

TEST /// --roundPSDmatrix
    Q=matrix{{2.01,0,0},{0,1.1,0},{0,0,2}}
    A=matrix{{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0}}
    b=matrix{{2},{1},{2}}
    d=100
    Gramind=hashTable {0 => {1,1},3 => {2,1},1 => {2,2},5=>{3,1},4=>{3,2},2=>{3,3}}
    LinIndex =applyPairs(Gramind,(i,j)->(j,i))
    (boolv,Qpsd)=roundPSDmatrix(Q,A,b,d,Gramind,LinIndex)
    e=eigenvalues Qpsd
    boolv2=if all(e, i -> i > 0) then 0 else 1
    assert((boolv,boolv2)==(true,0))
///

TEST /// --roundPSDmatrix2
    Q=matrix{{2.01,0,1},{0,1.1,0},{1,0,2}}
    eigenvalues Q
    A=matrix{{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,0,1,0}}
    b=matrix{{2},{1},{2},{1}}
    d=10
    Gramind=hashTable {0 => {1,1},3 => {2,1},1 => {2,2},4=>{3,1},5=>{3,2},2=>{3,3}}
    LinIndex =applyPairs(Gramind,(i,j)->(j,i))
    (boolv,Qpsd)=roundPSDmatrix(Q,A,b,d,Gramind,LinIndex)
    e=eigenvalues Qpsd
    boolv2=if all(e, i -> i > 0) then 0 else 1
    assert((boolv,boolv2)==(true,0))
///

TEST /// --blkDiag1
    A1=matrix{{1,0},{0,2}}
    A2=matrix{{1,1,3},{4,2,5},{2,1,1}}
    assert(blkDiag(A1,A2)==matrix{{1,0,0,0,0},{0,2,0,0,0},{0,0,1,1,3},{0,0,4,2,5},{0,0,2,1,1}})
///

TEST /// --blkDiag2
    A1=matrix{{1.4,0,2.5},{0,2,1.9},{1.2,3,6.1}}
    A2=matrix{{2.6,1,0},{4.1,2.6,5},{1.5,1,1}}
    assert(blkDiag(A1,A2)==matrix{{1.4,0,2.5,0,0,0},{0,2,1.9,0,0,0},{1.2,3,6.1,0,0,0},{0,0,0,2.6,1,0},{0,0,0,4.1,2.6,5},{0,0,0,1.5,1,1}})
///

TEST /// --sosdec
    R=QQ[x,y,z]
    Q=matrix{{1,0,0},{0,1,0},{0,0,1}}
    Q=promote(Q,QQ)
    mon=matrix{{x^2},{x*z},{y^2}}
    f=sosdec(Q,mon)
    assert(sumSOS f==transpose mon * Q *mon)
///

TEST /// --sosdec
    R=QQ[x,y,z]
    Q=matrix{{1,-1/2,1},{-1/2,1,-1/2},{1,-1/2,1}}
    --eigenvalues Q
    Q=promote(Q,QQ)
    mon=matrix{{x^3},{x^2*z},{y*z^2}}
    f=sosdec(Q,mon)
    assert(sumSOS f==transpose mon * Q *mon)
///

TEST /// --roundPSDmatrix3 (doesn't work for d=1, bug in LDL)
    Q=matrix{{1,-0.75,-0.75},{-0.75,1,0.99},{-0.75,0.99,1}}
    eigenvalues Q
    A=matrix{{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,0,1,0}}
    b=matrix{{1},{1},{1},{1}}
    d=10
    Gramind=hashTable {0 => {1,1},3 => {2,1},1 => {2,2},4=>{3,1},5=>{3,2},2=>{3,3}}
    LinIndex =applyPairs(Gramind,(i,j)->(j,i))
    (boolv,Qpsd)=roundPSDmatrix(Q,A,b,d,Gramind,LinIndex)
    e=eigenvalues Qpsd
    boolv2=if all(e, i -> i > 0) then 0 else 1
    assert((boolv)==(false))
///

TEST /// --roundPSDmatrix4 
    Q=matrix{{2,-1,-0.75},{-1,2,-1.1},{-0.75,-1.1,2}}
    eigenvalues Q
    A=matrix{{1,0,0,3,0,0},{2,1,0,0,5,0},{0,4,1,0,0,-3},{6,0,-5,0,1,0}}
    b=matrix{{1},{4},{-3},{1}}
    d=100
    Gramind=hashTable {0 => {1,1},3 => {2,1},1 => {2,2},4=>{3,1},5=>{3,2},2=>{3,3}}
    LinIndex =applyPairs(Gramind,(i,j)->(j,i))
    (boolv,Qpsd)=roundPSDmatrix(Q,A,b,d,Gramind,LinIndex)
    assert((boolv)==(false))
///


