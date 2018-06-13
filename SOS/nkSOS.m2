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



--p1=sospolyMultiply(g1,g2)




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




