--
restart
needsPackage ("SOS", Configuration=>{"CSDPexec"=>"/Users/kaihnsa/Documents/Workshop-2018-Leipzig-SOS/SOS/CSDP/csdp"})
R=QQ[x,y,z]
sospolyMultiply = method(
     Options => {RndTol => -3, Solver=>"CSDP", Verbose => false} )
sospolyMultiply(List,List,List,List):= o -> (coeff1,pol1,coeff2,pol2)-> (
    if #coeff1!=#pol1 then error "length of lists are not equal";
    if #coeff2!=#pol2 then error "length of lists are not equal";
    sqpol1=apply(pol1,i->i^2);
    p1=sum(for i from 0 to #coeff1-1 list coeff1_i*sqpol1_i);
    sqpol2=apply(pol2,i->i^2);
    p2=sum(for i from 0 to #coeff2-1 list coeff2_i*sqpol2_i);
    return(p1*p2);
    )

coeff1={1,2,3,4} 
pol1={x^2,y^2,x*y,x*z}
coeff2={2,3,4}
pol2={z*y,z^2,z^2+x^2}
p1=sospolyMultiply(coeff1,pol1,coeff2,pol2)
p2=sub(p1,R)

(Q2,mon2,X2) = solveSOS(p1, Solver=>"CSDP")
g2 = sosdec (Q2,mon2)
p1

sospolyPower = method(
     Options => {RndTol => -3, Solver=>"CSDP", Verbose => false} )
sospolyPower(List,List,ZZ) :=o -> (coeff,pol,D)->(
    if D<0 then error "power should be a positive number.";
    sqpol1=apply(pol1,i->i^2);
    p1=sum(for i from 0 to #coeff1-1 list coeff1_i*sqpol1_i);
    return(p1^D);
    )

sospolyPower(coeff1,pol1,4)

