# GROUP1 : Lat long H to XYZ WGS84 ----------------------------

Lat_Long_H_to_X <- function(PHI, LAM, H, a, b){
#Convert angle measures to radians
    RadPHI = PHI * (pi / 180)
    RadLAM = LAM * (pi / 180)

#Compute eccentricity squared and nu
  e2 = ((a ^ 2) - (b ^ 2)) / (a ^ 2)
  V = a / (sqrt(1 - (e2 * ((sin(RadPHI)) ^ 2))))
  
#Compute X
    kq <- (V + H) * (cos(RadPHI)) * (cos(RadLAM))
    return(kq)
}


Lat_Long_H_to_Y <- function(PHI, LAM, H, a, b){
  #Convert angle measures to radians
  RadPHI = PHI * (pi / 180)
  RadLAM = LAM * (pi / 180)
  
  #Compute eccentricity squared and nu
  e2 = ((a ^ 2) - (b ^ 2)) / (a ^ 2)
  V = a / (sqrt(1 - (e2 * ((sin(RadPHI)) ^ 2))))
  
  #Compute Y
  kq <- (V + H) * (cos(RadPHI)) * (sin(RadLAM))
  return(kq)
}


Lat_H_to_Z <- function(PHI, H, a, b){
  #Convert angle measures to radians
  RadPHI = PHI * (pi / 180)
  
  #Compute eccentricity squared and nu
  e2 = ((a ^ 2) - (b ^ 2)) / (a ^ 2)
  V = a / (sqrt(1 - (e2 * ((sin(RadPHI)) ^ 2))))
  
  #Compute Y
  kq <- ((V * (1 - e2)) + H) * (sin(RadPHI))
  return(kq)
}



#GROUP 2: XYZ to X1Y1Z1 VN2000-------------

Helmert_X <- function(X, Y, Z, DX, Y_rot, Z_rot, S){
  #Convert rotations to radians and ppm scale to a factor
  sfactor = S * 0.000001
  RadY_Rot = (Y_rot / 3600) * (pi / 180)
  RadZ_Rot = (Z_rot / 3600) * (pi / 180)
  
  #Compute transformed X coord
  Helmert_X = X + (X * sfactor) - (Y * RadZ_Rot) + (Z * RadY_Rot) + DX
  return(Helmert_X)
}

Helmert_Y <- function(X, Y, Z, DX, Y_rot, Z_rot, S){
  #Convert rotations to radians and ppm scale to a factor
  sfactor = S * 0.000001
  RadX_Rot = (X_rot / 3600) * (pi / 180)
  RadZ_Rot = (Z_rot / 3600) * (pi / 180)
  
  #Compute transformed X coord
  Helmert_Y = (X * RadZ_Rot) + Y + (Y * sfactor) - (Z * RadX_Rot) + DY
  return(Helmert_Y)
}


Helmert_Z <- function(X, Y, Z, DX, Y_rot, Z_rot, S){
  #Convert rotations to radians and ppm scale to a factor
  sfactor = S * 0.000001
  RadX_Rot = (X_rot / 3600) * (pi / 180)
  RadY_Rot = (Y_rot / 3600) * (pi / 180)
  
  #Compute transformed X coord
  Helmert_Z = (-1 * X * RadY_Rot) + (Y * RadX_Rot) + Z + (Z * sfactor) + DZ
  return(Helmert_Z)
}






#GROUP 3 X1Y1Z1 to Lat Long H VN2000 -------------
Iterate_XYZ_to_Lat <- function(a, e2, PHI1, Z, RootXYSqr){
  V = a / (sqrt(1 - (e2 * ((sin(PHI1)) ^ 2))))
  PHI2 = atan((Z + (e2 * V * (sin(PHI1)))) / RootXYSqr)
  
  while(abs(PHI1 - PHI2) > 0.000000001){
    PHI1 = PHI2
    V = a / (sqrt(1 - (e2 * ((sin(PHI1)) ^ 2))))
    PHI2 = atan((Z + (e2 * V * (sin(PHI1)))) / RootXYSqr)
  }
  return(PHI2)
}






XYZ_to_Lat <- function(X,Y,Z,a,b){
  RootXYSqr = sqrt((X ^ 2) + (Y ^ 2))
  e2 = ((a ^ 2) - (b ^ 2)) / (a ^ 2)
  PHI1 = atan(Z / (RootXYSqr * (1 - e2)))
  
  PHI = Iterate_XYZ_to_Lat(a, e2, PHI1, Z, RootXYSqr)
  XYZ_to_Lat = PHI * (180 / pi)
  return(XYZ_to_Lat)
  
}


XYZ_to_Long <- function(X,Y){
  if (X >= 0){
    #longitude is in the W90 thru 0 to E90 hemisphere
    XYZ_to_Long = (atan()) * (180 / pi)
  }else if(X < 0 & Y >= 0){
    #longitude is in the E90 to E180 quadrant
    XYZ_to_Long = ((atan(Y / X)) * (180 / pi)) + 180
  }else if(X < 0 & Y < 0){
    #longitude is in the E180 to W90 quadrant
    XYZ_to_Long = ((atan(Y / X)) * (180 / pi)) - 180
  }
  return(XYZ_to_Long)
}


XYZ_to_H <- function(X,Y,Z,a,b){
  #Compute PHI (Dec Degrees) first
  PHI = XYZ_to_Lat(X, Y, Z, a, b)

  #Convert PHI radians
  RadPHI = PHI * (pi / 180)
  
  #Compute H
  RootXYSqr = sqrt((X ^ 2) + (Y ^ 2))
  e2 = ((a ^ 2) - (b ^ 2)) / (a ^ 2)
  V = a / (sqrt(1 - (e2 * ((sin(RadPHI)) ^ 2))))
  H = (RootXYSqr / cos(RadPHI)) - V
    
  return(H)
  
}





# GROUP 4 : Lat Long H to E N----------
Lat_Long_to_East <- function(PHI, LAM, a, b, E0, f0, PHI0, LAM0){
  
  RadPHI = PHI * (pi / 180)
  RadLAM = LAM * (pi / 180)
  RadPHI0 = PHI0 * (pi / 180)
  RadLAM0 = LAM0 * (pi / 180)
  
  af0 = a * f0
  bf0 = b * f0
  e2 = ((af0 ^ 2) - (bf0 ^ 2)) / (af0 ^ 2)
  n = (af0 - bf0) / (af0 + bf0)
  nu = af0 / (sqrt(1 - (e2 * ((sin(RadPHI)) ^ 2))))
  rho = (nu * (1 - e2)) / (1 - (e2 * (sin(RadPHI)) ^ 2))
  eta2 = (nu / rho) - 1
  p = RadLAM - RadLAM0
  
  IV = nu * (cos(RadPHI))
  V = (nu / 6) * ((cos(RadPHI)) ^ 3) * ((nu / rho) - ((tan(RadPHI) ^ 2)))
  VI = (nu / 120) * ((cos(RadPHI)) ^ 5) * (5 - (18 * ((tan(RadPHI)) ^ 2)) + ((tan(RadPHI)) ^ 4) + (14 * eta2) - (58 * ((tan(RadPHI)) ^ 2) * eta2))
  
  Lat_Long_to_East = E0 + (p * IV) + ((p ^ 3) * V) + ((p ^ 5) * VI)
  return(Lat_Long_to_East)
}

Marc <- function(bf0, n, PHI0, PHI){
  Marc = bf0 * (((1 + n + ((5 / 4) * (n ^ 2)) + ((5 / 4) * (n ^ 3))) * (PHI - PHI0))-
                  (((3 * n) + (3 * (n ^ 2)) + ((21 / 8) * (n ^ 3))) * (sin(PHI - PHI0)) * (cos(PHI + PHI0)))+
                  ((((15 / 8) * (n ^ 2)) + ((15 / 8) * (n ^ 3))) * (sin(2 * (PHI - PHI0))) * (cos(2 * (PHI + PHI0))))-
                  (((35 / 24) * (n ^ 3)) * (sin(3 * (PHI - PHI0))) * (cos(3 * (PHI + PHI0)))))
  return(Marc)
}



Lat_Long_to_North <- function(PHI, LAM, a, b, E0, N0, f0, PHI0, LAM0){

  RadPHI = PHI * (pi / 180)
  RadLAM = LAM * (pi / 180)
  RadPHI0 = PHI0 * (pi / 180)
  RadLAM0 = LAM0 * (pi / 180)
  
  af0 = a * f0
  bf0 = b * f0
  e2 = ((af0 ^ 2) - (bf0 ^ 2)) / (af0 ^ 2)
  n = (af0 - bf0) / (af0 + bf0)
  nu = af0 / (sqrt(1 - (e2 * ((sin(RadPHI)) ^ 2))))
  rho = (nu * (1 - e2)) / (1 - (e2 * (sin(RadPHI)) ^ 2))
  eta2 = (nu / rho) - 1
  p = RadLAM - RadLAM0
  M = Marc(bf0, n, RadPHI0, RadPHI)
  
  I = M + N0
  II = (nu / 2) * (sin(RadPHI)) * (cos(RadPHI))
  III = ((nu / 24) * (sin(RadPHI)) * ((cos(RadPHI)) ^ 3)) * (5 - ((tan(RadPHI)) ^ 2) + (9 * eta2))
  IIIA = ((nu / 720) * (sin(RadPHI)) * ((cos(RadPHI)) ^ 5)) * (61 - (58 * ((tan(RadPHI)) ^ 2)) + ((tan(RadPHI)) ^ 4))
  
  Lat_Long_to_North = I + ((p ^ 2) * II) + ((p ^ 4) * III) + ((p ^ 6) * IIIA)
  return(Lat_Long_to_North)
}


# main funtion ----------------------------------------

wgs84_to_vn2000 <- function(long,lat,alt){
  source("R/parameter.R")
  f0 <<- find_f0(muichieu)
  # Lat long H to XYZ WGS84 ---------------------
  X <- Lat_Long_H_to_X(PHI = lat, LAM = long, H=alt, a=a, b=b)
  Y <- Lat_Long_H_to_Y(PHI = lat, LAM = long, H=alt, a=a, b=b)
  Z <- Lat_H_to_Z(PHI = lat, H = alt, a=a, b=b)
  
  #XYZ to X1Y1Z1 VN2000-------------
  X1 = Helmert_X(X, Y, Z, DX, Y_rot, Z_rot, S)
  Y1 = Helmert_Y(X, Y, Z, DY, X_rot, Z_rot, S)
  Z1 = Helmert_Z(X, Y, Z, DZ, X_rot, Y_rot, S)
  
  #X1Y1Z1 to Lat Long H VN2000
  B1 = XYZ_to_Lat(X1, Y1, Z1, a, b)
  L1 = XYZ_to_Long(X1, Y1)
  H1 = XYZ_to_H(X1, Y1, Z1, a, b)
  
  
  
  #Lat Long H to E N
  VN2000X = Lat_Long_to_East(B1, L1, a, b, E0, f0, PHI0, LAM0)
  VN2000Y = Lat_Long_to_North(B1, L1, a, b, E0, N0, f0, PHI0, LAM0)
  VN2000H = H1
  
  vn2000 <-c(VN2000X,VN2000Y,VN2000H)
  names(vn2000) <- c("x","y","z")
  return(vn2000)
}


find_f0 <- function(muichieu){
  if(muichieu == 3){
    f0 <- 0.999600000000
  }else if(muichieu == 6){
    f0 <- 0.999900000000
  }else{
    stop("muichieu must be 3 or 6")
  }
  return(f0)
}


