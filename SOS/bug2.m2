needsPackage ("SOS" , Configuration => { "CSDPexec" => "CSDP/csdp", "SDPAexec" => "./sdpa"})
R = QQ[x];
L = sosInIdeal({x+1},2, Solver=>"CSDP")
print L

