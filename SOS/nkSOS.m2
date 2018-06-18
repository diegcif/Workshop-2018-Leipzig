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
    sqpol1=sumSOS(p1);
    p2=sqpol1^(D-2);
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
    p4#ring===S
///

TEST /// --sub(SOS)2
    R = QQ[x,y,z]
    coeff1={1,2,3,4}
    pol1={x^2,x*y,y^2,x}
    p1=sosPoly(pol1,coeff1)
    p1#ring===R
    S=RR[x,y]
    p4=sub(p1,S)
    p4#ring===S
///    
