rm(list=ls())

province <- read.csv("data/province.csv")
#tinh Tay Ninh
LAM0 <- 105.50
# de tim thong so cua cac tinh thanh khac xem trong df "province"

lat <- 11.58976
long <- 105.88743
alt <- 0

wgs84_to_vn2000(long = long, lat=lat, alt = 0)

