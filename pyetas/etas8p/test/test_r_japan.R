
library("ETAS")
#japan.cat <- catalog(japan.quakes)


jpoly <- list(long = c(134.0, 137.9, 143.1, 144.9, 147.8, 137.8, 137.4, 135.1, 130.6), lat = c(31.9, 33.0, 33.2, 35.2, 41.3, 44.2, 40.2, 38.0, 35.4))
japan.cat <- catalog(japan.quakes, study.start = "1953-05-26", study.end = "1990-01-08", region.poly = jpoly, mag.threshold = 4.5)




param0=NULL

param0 <- c(0.592844590, 0.204288231, 0.022692883, 1.495169224, 1.109752319, 0.001175925, 1.860044210, 1.041549634)

bwd = NULL
nnp = 5
bwm = 0.05
verbose = TRUE
plot.it = FALSE
ndiv = 1000
no.itr = 11
rel.tol=1e-03
eps = 1e-06
cxxcode = FALSE
nthreads = 1

  ptm <- proc.time()
  spatstat::verifyclass(object, "catalog")
  revents <-  object$revents
  rpoly <-    object$rpoly
  rtperiod <- object$rtperiod
  m0 <- object$mag.threshold
  win <- object$region.win

  if (nthreads > parallel::detectCores())
  {
    stop(paste("nthreads can not be greater than", parallel::detectCores(),
               "on this machine!"))
  }

  # initial prameter values
  if (is.null(param0))
  {
    mu0 <- nrow(revents)/(4 * diff(rtperiod) * spatstat::area.owin(win))
    param0 <- c(mu=mu0, A=0.01, c=0.01, alpha=1, p=1.3, D=0.01, q=2,
                gamma=1)
    if (object$dist.unit == "km")
      param0["D"] <- 111^2 * param0["D"]
    if (verbose)
    {
      cat("using non-informative initial parameter values:\n")
      print(param0)
    }
    warning("the algorithm is very sensitive to the choice of starting point")
  }
  # bandwidths for smoothness and integration
  if (is.null(bwd))
  {
    if (object$dist.unit == "km")
      bwm <-  6371.3 * pi / 180 * bwm
    rbwd <- spatstat::nndist.default(revents[, 2], revents[, 3], k=nnp)
    rbwd <- pmax(rbwd, bwm)
  }
  else
  {
    stopifnot(is.numeric(bwd), length(bwd) != nrow(revents))
    rbwd <- bwd
  }

  # check initial values for the model parameters
  if (!is.numeric(param0) || length(param0) != 8 || any(param0 < 0))
    stop("param0 must be a numeric vector of length 8 with positive components")

  param1 <- param0
  thetar <- asd <- matrix(NA, nrow=no.itr, ncol=8)
  par.names <- c("mu", "A", "c", "alpha", "p", "D", "q", "gamma")
  names(param1) <- colnames(thetar) <- colnames(asd) <- par.names
  loglikfv <- numeric(no.itr)
  rownames(thetar) <- rownames(asd) <- names(loglikfv) <- paste("iteration", 1:no.itr)
  ihess <- diag(8)
  bk <- numeric(nrow(revents))

  for (itr in 1:no.itr)
  {
    cat("declustering:\n")
    bkgpbar <- utils::txtProgressBar(min=0, max=no.itr + 1 - itr, style=3)
    for (l in 1:(no.itr + 1 - itr))
    {
      bkg <- decluster(param1, rbwd, revents, rpoly, rtperiod, ndiv, cxxcode)
      revents <- bkg$revents
      utils::setTxtProgressBar(bkgpbar, l)
	cat(bkg$integ0)
    }
    close(bkgpbar)
    integ0 <- bkg$integ0
    dbk <- bk - revents[, 6]
    bk <- revents[, 6]
    pb <- revents[, 7]
    if (verbose)
    {
      cat("iteration: ", itr, "\n")
      cat("======================================================\n")
      cat("background seismicity rate:\n")
      print(summary(bk))
      cat("probability of being a background event:\n")
      print(summary(pb))
      cat("integral of background seismicity rate: ", integ0, "\n")
      cat("======================================================\n")
    }
