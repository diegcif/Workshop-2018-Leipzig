restart
load "realzeros.m2"
-- K = QQ[c_1,c_2,c_3]
K = QQ
x = symbol x
R = K[x_1..x_9]

k_1 = 1
k_2 = 1
k_3 = 1
k_4 = 1
k_5 = 1
k_6 = 2
k_7 = 1
k_8 = 1
k_9 = 1
k_10 = 1
k_11 = 1
k_12 = 1

f1 = -k_1*x_1*x_2+k_2*x_3+k_12*x_9
f2 = -k_1*x_1*x_2+(k_2+k_3)*x_3-k_4*x_2*x_4+(k_5+k_6)*x_5
f3 = k_1*x_1*x_2-(k_2+k_3)*x_3
f4 = k_3*x_3-k_4*x_2*x_4+k_5*x_5+k_9*x_8-k_10*x_4*x_7+k_11*x_9
f5 = k_4*x_2*x_4-(k_5+k_6)*x_5
f6 = k_6*x_5-k_7*x_6*x_7+k_8*x_8
f7 = -k_7*x_6*x_7+(k_8+k_9)*x_8-k_10*x_4*x_7+(k_11+k_12)*x_9
f8 = k_7*x_6*x_7-(k_8+k_9)*x_8
f9 = k_10*x_4*x_7-(k_11+k_12)*x_9
(cc1,cc2,cc3) = (1,1,1)
cons1 = x_2 + x_3 + x_5 - cc1
cons2 = x_7 + x_8 + x_9  -cc2
cons3 = x_1 + x_3 + x_4 + x_5 + x_6 + x_8 + x_9 - cc3

I = ideal(f1,f2,f3,f4,f5,f6,f7,f8,f9,cons1,cons2,cons3)
inputpolys = I_*
pts = realZeros(inputpolys, 8, CleanTol=>1e-3, ResTol=>1e-2)
