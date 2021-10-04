library("ETAS")

gregion <- list(lat = c(26, 25, 29, 38, 35), long = c(52, 59, 58, 45, 43))
iran.cat <- catalog(iran.quakes, study.start = "1991/01/01", study.end = "2011/01/01", region.poly = gregion, mag.threshold = 4.5)
, roundoff = FALSE)

param01 <- c(0.5, 0.2, 0.05, 2.7, 1.2, 0.02, 2.3, 0.03)
iran.fit <- etas(iran.cat, param0 = param01, cxxcode=FALSE)