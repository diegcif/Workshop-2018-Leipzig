needsPackage( "SOS", Configuration=>{"CSDPexec"=>"CSDP/csdp"} )


R=QQ[x];
--f = x^2+y^2+1;
f = (x-1)^2 + (x+3)^2;
r = lowerBound(f, {}, RndTol=>12, Solver => "CSDP");
print(r);

-- minimize f(x) subject to h(x) = 0
R=QQ[x,y,z];
f = 10 - x^2 - y
h = {x^2+y^2+z^2-1}
r = lasserreHierarchy(f,h,2, RndTol=>12, Solver => "CSDP")
print(r);

R=QQ[x,y];
f = x^2 + y^2
h = {x^2*y-1}
r = lasserreHierarchy(f,h,4, RndTol=>12, Solver => "CSDP")
print(r);
