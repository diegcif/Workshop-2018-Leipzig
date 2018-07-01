needsPackage ("SOS" , Configuration => { "CSDPexec" => "../CSDP/csdp"})
R = QQ[x,y]
f = (x-y)^2 + x^2 + (y-4)^2;
(bound, sol) = lowerBound(f, Solver=>"CSDP");
-- the minimum is 16/3.
result = (abs (16/3-bound) < 0.001)

