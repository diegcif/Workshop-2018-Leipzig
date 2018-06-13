--
restart
needsPackage ("SOS", Configuration=>{"CSDPexec"=>"/Users/kaihnsa/Documents/Workshop-2018-Leipzig-SOS/SOS/CSDP/csdp"})
R=QQ[x,y,z]
sospolyMultiply = method()
    
sospolyMultiply(SOSPoly,SOSPoly):= (p1,p2)-> (
    sqpol1=apply(p1#gens,i->i^2);
    q1=sum(for i from 0 to #p1#coefficients-1 list p1#coefficients_i*sqpol1_i);
    sqpol2=apply(p2#gens,i->i^2);
    q2=sum(for i from 0 to #p2#coefficients-1 list p2#coefficients_i*sqpol2_i);
    return(q1*q2);
    )


--p1=sospolyMultiply(g1,g2)




sospolyPower = method(
      )
sospolyPower(SOSPoly,ZZ) := (p1,D)->(
    if D<0 then error "power should be a positive number.";
    sqpol1=apply(p1#gens,i->i^2);
    p2=sum(for i from 0 to #p1#coefficients-1 list p1#coefficients_i*sqpol1_i);
    return(p2^D);
    )
--sospolyPower(g2,1)




