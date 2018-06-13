restart

needsPackage ("SOS", Configuration=>{"CSDPexec"=>"/Users/kaihnsa/Documents/Workshop-2018-Leipzig-SOS/SOS/CSDP/csdp"})
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
(p2,q2) = sosdecTernary(f2)
f3=x*y*z+x^4
f={f1,f2,f3}
p=sospolyIdeal(f,8)

--Lax-Lax
restart
needsPackage ("SOS", Configuration=>{"CSDPexec"=>"/Users/kaihnsa/Documents/Workshop-2018-Leipzig-SOS/SOS/CSDP/csdp"})
R = QQ[a,b,c,d]
h = a^2 + b^2 + c^2 +d^2
f1=(a-b)*(a-c)*(a-d)*a+(b-a)*(b-c)*(b-d)*b+(c-a)*(c-b)*(c-d)*c+(d-a)*(d-b)*(d-c)*d+a*b*c*d
(Q1,mon1,X1) = solveSOS (f1*h, Solver=>"CSDP")
g1 = sosdec (Q1,mon1)

--Harris polynomial
R = QQ[x,y,z]
-- (a,b,c,d,e) = (1,1,1,1,1) -- works
--(a,b,c,d,e) = (1,1,-2,0,0)
--a -> 88, b -> -71, c -> -17, d -> 504, e -> -590
--Harris pol
(a,b,c,d,e) = (16,-36,20,57,-38)
f4 = a*( x^10 + y^10 + z^10)+ 
    b*( x^8* y^2 + x^2* y^8 + x^8* z^2 + x^2* z^8 + y^8* z^2 + y^2* z^8 ) +
    c*( x^6* y^4 + x^4* y^6 + x^6* z^4 + x^4* z^6 + y^6* z^4 + y^4* z^6 ) + 
    d*( x^6* y^2* z^2 + x^2* y^6* z^2 + x^2* y^2* z^6) +
    e*( x^4* y^4* z^2 + x^4* y^2* z^4 + x^2* y^4* z^4)
(Q4,mon4,X4) = solveSOS (f4*h, Solver=>"CSDP")
g4 = sosdec (Q4,mon4)
(p4,q4) = sosdecTernary (f4)
sumSOS(g4,d4) == f4*h


