# khai bao bien --------------
muichieu <- 3
a <- 6378137.0
f <- 1/298.257223563
b <- a*(1-f)
E0 <- 500000
N0 <- 0
# Cac thong so chuyen theo thong tu so 973/2001/TT-TCĐC
#translation parallel to X, Y, Z (unit: metres)
DX <- 191.9044143 
DY <- 39.30318279
DZ <- 111.45032835

#rotation about X, Y, Z (unit: seconds of arc)
X_rot <- -0.00928836 
Y_rot <- 0.01975479
Z_rot <- -0.00427372

#scale change (unit: parts per million)
S <- 0.252906278

PHI0 <- 0