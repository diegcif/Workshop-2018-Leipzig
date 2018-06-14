--
restart
needsPackage ("SOS", Configuration=>{"CSDPexec"=>"/Users/kaihnsa/Documents/Workshop-2018-Leipzig-SOS/SOS/CSDP/csdp"})
R=QQ[x,y,z]
sospolyMultiply = method()
    
sospolyMultiply(SOSPoly,SOSPoly):= (g1,g2)-> (
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




sospolyPower = method(
      )
sospolyPower(SOSPoly,ZZ) := (p1,D)->(
    if D<0 then error "power should be a positive integer.";
    sqpol1=apply(p1#gens,i->i^2);
    p2=sum(for i from 0 to #p1#coefficients-1 list p1#coefficients_i*sqpol1_i);
    (Q1,mon1,X1) = solveSOS (p2^D, Solver=>"CSDP");
    g1 = sosdec (Q1,mon1);
    return(g1);
    )
sospolyPower(g2,2)




