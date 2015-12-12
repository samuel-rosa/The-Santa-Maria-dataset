# USER DEFINED FUNCTIONS #######################################################
rboxcox <- 
  function (n, lambda, lambda2 = NULL, mean = 0, sd = 1, verbose = TRUE) {
    if (is.null(lambda2) || is.na(lambda2))
      lambda2 <- 0
    xn <- rnorm(n = n, mean = mean, sd = sd)
    if (isTRUE(all.equal(unname(lambda), 0))) 
      xbc <- exp(xn)
    else {
      xbc <- rep(NA, n)
      ind <- xn < -1/lambda
      sum.ind <- sum(ind)
      if (sum.ind > 0)
        if (verbose) {
          cat(paste("rboxcox: WARNING ", sum.ind, "values truncated to 0\n"))
        }
      xn[ind] <- -1/lambda
      xbc <- ((xn * lambda) + 1)^(1/lambda)
    }
    return(xbc - lambda2)
  }

negative <- 
  function (x) {
    x < 0
  }
anyZero <-
  function (x) {
    any(x == 0)
  }
anyNegative <- 
  function (x) {
    any(x < 0)
  }
whichNegative <-
  function (x) {
    which(x < 0)
  }
# This function is similar to stats::terms
trend.terms <- 
  function (x) {
    cl <- class(x)
    if (all(cl == c("likGRF","variomodel"))) {
      res <- all.vars(x$trend)
    }
    return (res)
  }
# This function is similar to stats::model.frame
trend.matrix <-
  function (x) {
    cl <- class(x)
    if (all(cl == c("likGRF","variomodel"))) {
      res <- x$trend.matrix[, -1]
      colnames(res) <- trend.terms(x)
    }
    return (res)
  }
# Automatically prepare initial covariance parameters for likfit
likfit_ini.cov.pars <- 
  function (obj, z) {
    v <- var(na.omit(obj@data[, z]))
    v <- c(0.5 * v, 0.75 * v, v)
    r <- max(dist(obj@coords)) / 3
    r <- c(0.1 * r, 0.2 * r, 0.3 * r)
    res <- as.matrix(expand.grid(v, r))
    return (res)
  }
likfit_nugget <-
  function (obj, z) {
    n <- var(na.omit(obj@data[, z]))
    n <- c(0, 0.25 * n, 0.5 * n)
    return (n)
  }
fit_reml <-
  function (obj, z, nugget, lik.method = "REML", ...) {
    
    # Nugget
    if (missing(nugget)) {
      nugget <- likfit_nugget(obj, z)
      fix.nugget <- FALSE
    } else {
      fix.nugget <- TRUE
    }
    
    # Fit REML
    res <- geoR::likfit(geodata = as.geodata(obj, data.col = z),
                        ini.cov.pars = likfit_ini.cov.pars(obj, z),
                        nugget = nugget, fix.nugget = fix.nugget,
                        lik.method = lik.method, ...)
    
    # Output
    return (res)
  }
# spatially correlated variance ################################################
vgmSCV <- 
  function (x, digits = 4) {
    res <- x$sigmasq / (x$sigmasq + x$nugget)
    res <- round(res, digits = digits)
    return (res)
  }
# compute the nadir and utopia points for ACDC #################################
# nadirACDC <-
#   function (osc = list(CORR, DIST), candi, covars) {
#     
#     # Convert numeric covariates into factor covariates
#     if (pedometrics::anyFactor(covars) && !pedometrics::allFactor(covars)) {
#       id <- which(!sapply(covars, is.factor))
#       message(paste("converting ", length(id), 
#                     " numeric covariates into factor covariates...", sep = ""))
#       covars[, id] <- 
#         pedometrics::stratify(x = covars[, id], n = nrow(osc$CORR))
#     }
#     
#     # Compute objective function values
#     obj_dist <- sapply(osc, objDIST, candi = candi, covars = covars)
#     obj_corr <- sapply(osc, objCORR, candi = candi, covars = covars)
#     
#     # Prepare output
#     res <- data.frame(CORR = obj_corr, DIST = obj_dist)
#     return (res)
#   }
# compute the nadir and utopia points for SPAN #################################
# nadirSPAN <-
#   function (osc = list(CORR, DIST, PPL, MSSD), candi, covars, x.max, y.max) {
#     
#     # Convert numeric covariates to factor covariates
#     if (pedometrics::anyFactor(covars) && !pedometrics::allFactor(covars)) {
#       id <- which(!sapply(covars, is.factor))
#       covars[, id] <- pedometrics::stratify(x = covars[, id], n = nrow(osc$CORR))
#     }
#     
#     # Compute objective function values
#     obj_dist <- sapply(osc, objDIST, candi = candi, covars = covars)
#     obj_corr <- sapply(osc, objCORR, candi = candi, covars = covars)
#     obj_ppl <- sapply(osc, objPPL, candi = candi, x.max = x.max, y.max = y.max)
#     obj_mssd <- sapply(osc, objMSSD, candi = candi)
#     
#     # Prepare output
#     res <- data.frame(DIST = obj_dist, CORR = obj_corr, PPL = obj_ppl,
#                       MSSD = obj_mssd)
#     return (res)
#   }
# Plot optimized sample configuration
# plot.OSC <-
#   function (osc, which = 1:2, boundary) {
#     
#     par0 <- par()
#     on.exit(suppressWarnings(par(par0)))
#     if (all(which == 1:2)) {
#       graphics::par(mfrow = c(1, 2))
#     }
#     
#     # Plot the energy states
#     if (all(which == 1:2) || which == 1) {
#       k <- attr(osc, "iterations")
#       a <- objSPSANN(OSC = osc, n = k + 1)
#       
#       # Multi-objective optimization problem
#       if (is.data.frame(a)) {
#         l <- colnames(a)
#         n <- ncol(a)
#         graphics::plot(1, type = 'n', xlim = c(0, k), 
#                        ylim = c(0, max(sapply(a, max)) * 1.1), 
#                        xlab = "iteration", ylab = "energy state")
#         graphics::legend("topright", legend = l, lwd = 1, lty = 1:n)
#         for(i in 1:ncol(a)) {
#           graphics::lines(a[, i] ~ seq(1, k + 1), type = "l", lty = i)
#         }
#         
#         # Single-objective optimization problem
#       } else {
#         graphics::plot(a ~ c(0:k), type = "l", xlab = "iteration", 
#                        ylab = "energy state")
#       }
#     }
#     
#     # Plot optimized sample configuration
#     if (which == 1:2 || which == 2) {
#       if (!missing(boundary)) {
#         bb <- sp::bbox(boundary)
#         if (class(boundary) == "SpatialPoints") {
#           sp::plot(boundary, pch = 20, cex = 0.1)
#         } else {
#           sp::plot(boundary)
#         }
#         graphics::points(osc[, 2], osc[, 3], pch = 20, cex = 0.5)  
#       } else {
#         plot(osc[, 2:3], pch = 20, cex = 0.5)  
#       }
#     }
#   }
# Prepare argument 'newdata' for simulations
newdata4sim <-
  function (dnos_raster, covars, model, trend = "lmm") {
    if (trend == "lmm") {
      newdata <- cbind(dnos_raster, covars[, trend.terms(model)])
    } else {
      newdata <- data.frame(dnos_raster)
    }
    gridded(newdata) <- ~ x + y
    return (newdata)
  }

# prepare ini.cov.pars
vgmICP <- 
  function (z, coords, lags = pedometrics::vgmLags(coords)[-1], 
            method = "asr", min.npairs = 30, cov.model = "RMmatern", nu = 0.5,
            scv = 0.666, 
            expand = FALSE, by = list(tausq = 100, sigmasq = 100, phi = 100),
            ...) {
    # Compute variogram
    v <- georob::sample.variogram(
      response = z, locations = coords, lag.class.def = lags)
    
    # Merge lag-distance classes that have too few point-pairs
    if (any(v$npairs < min.npairs)) {
      idx <- which(v$npairs < min.npairs)
      i <- 1
      while (length(idx) >= 1) {
        lags <- lags[-idx[i]]
        v <- georob::sample.variogram(
          response = z, locations = coords, lag.class.def = lags)
        idx <- which(v$npairs < min.npairs)
        i <- i + 1
      }
    }
    
    # RANGE
    # In general, the initial guess for the range (scale) parameter is made 
    # based on the lag-distance classes and on the dimensions of the study area.
    # The most commom rule is to compute the initial range as half the 
    # maximum distance up to which lag-distance classes have been defined 
    # (Jian et al., 1996; Larrondo et al., 2003; Desassis & Renard, 2012). 
    # Others set the initial range to a proportion of the diagonal of the study 
    # area, say 0.1, 0.35 or 0.5 (Hiemstra et al., 2009).
    # I think that this rather arbitrary and, possibly, suboptimal because the 
    # lag-distance classes usually are defined by some authomatic proceedure
    # implemented in the software being used, which does not account for the
    # features of the data that is being analysed.
    # I propose using the estimate of the variance and semivariance of the data
    # to make an initial guess for the range parameter. The variance is used
    # here because it is the initial guess of the total sill (see bellow).
    # I start computing the absolute difference between the semivariogram and 
    # the variance in each lag-distance classe (except for the first). Then, I
    # record the index of the lag-distance class where the semivariance is 
    # closest to the variance. The separation distance at the centre of this 
    # lag-distance class is used as the initial guess for the range parameter.
    range <- switch(
      method, 
      a = { # Samuel-Rosa
        v$lag.dist[which.min(abs(v$gamma[-c(1, length(v$gamma))] - var(z))) + 1] 
      }, 
      b = { # JianEtAl1996
        lags[length(lags)] * 0.5 
      },
      c = { # HiemstraEtAl2009
        sp::spDists(t(sp::bbox(SpatialPoints(coords))))[1, 2] * 0.1
      },
      d = { # DesassisEtAl2012
        lags[length(lags)] * 0.5
      },
      e = { # LarrondoEtAl2003
        lags[length(lags)]
      }
    )
    # Correct the initial guess of the range parameter for exponential...
    if (cov.model %in% c("exponential", "Exponential", "exp", "Exp")) {
      range <- range / 3
    }
    # ... And Gaussian covariance models
    if (cov.model %in% c("gaussian", "Gaussian", "gauss", "Gauss")) {
      range <- range / sqrt(3)
    }
    if (cov.model %in% "RMmatern") {
      if (nu == 0.5) { range <- range / 3 }
      if (nu == 2.0) { range <- range / sqrt(3) }
    }
    # NUGGET
    # The initial guess for the nugget variance is commonly made using one of
    # the following rules:
    # 1) use the minimum semivariance in the sample variogram (Hiemstra et al.,
    #    2009)
    # 2) set the initial nugget value to zero (Larrondo et al., 2003)
    # We can also find rules that take into account the difference in the
    # semivariance between the first and second lag-distance classes ponderated
    # by the difference in the size of these lag-distance classes (Jian et al.,
    # 1996). The resulting initial guess for the nugget variance is always lower
    # than the minimum semivariance value.
    nugget <- switch(
      method,
      a = { # Samuel-Rosa
        
        # Is the smallest semivariance in lags others than the first?
        if (which.min(v$gamma) > 1) {
          
          # Is the semivariance of the fist lag larger than that of the second
          # lag? Or is the semivariance of the fist lag larger than the variance of
          # z? Or is the semivariance of the fist and second lags larger than that
          # of the third? This looks like a messy sample variogram.
          if (v$gamma[1] >= v$gamma[2] || v$gamma[1] >= var(z) || 
              which.min(v$gamma[1:3]) == 3) {
            
            # Well... Let us imply use the minimum gamma value.
            min(v$gamma)
            
          } else {
            
            # In the worst case we resort to the user-estimate of the amount of 
            # spatially correlated variance
            if (scv >= 0 && scv <= 1) {
              var(z) * (1 - scv)
            } else {
              stop("'scv' must be in the 0:1 interval")
            }
          }
          
        } else {
          # Is the difference in the semivariance of the first and second 
          # lags larger than the difference between the semivariance in the 
          # second lag and the variance of z? If TRUE, it may be that the 
          # semivariance in the first lag-distance class is spuriously low 
          # (or high).
          if (diff(v$gamma[1:2]) > abs(diff(c(v$gamma[2], var(z))))) {
            
            # Let us check if the semivariance increases monotonicaly from the
            # first to the third lag-distance class. If not, then it means that 
            # it will be hard to estimate the shape of the variogram close to 
            # the origin. It
            # may be better to use a conservative initial guess, such as the mean 
            # of the semivariance at the first lag and the minimum semivariance 
            # among the other lags (up to the third?).
            if (!identical(order(v$gamma[1:3]), 1:3)) {
              mean(c(v$gamma[1], min(v$gamma[-1])))
              # mean(c(v$gamma[1], min(v$gamma[2:3])))
              
            } else {
              if (scv >= 0 && scv <= 1) {
                var(z) * (1 - scv)
              } else {
                stop("'scv' must be in the 0:1 interval")
              }
            }
          } else {
            v$gamma[1] - diff(c(0, v$gamma[1])) * 0.5
            # if (identical(order(v$gamma[1:3]), 1:3)) {
            # if (scv >= 0 && scv <= 1) {
            # var(z) * (1 - scv)
            # } else {
            # stop("'scv' must be in the 0:1 interval")
            # }
            # } else {
            # v$gamma[1] - diff(c(0, v$gamma[1])) * 0.5
            # }
          }
        }
      },
      b = { # JianEtAl1996
        max(0, c(v$gamma[1] - (lags[1] / diff(lags[1:2])) * diff(v$gamma[1:2])))
      },
      c = { # HiemstraEtAl2009
        min(v$gamma)
      },
      d = { # DesassisEtAl2012
        1e-12
      },
      e = { # LarrondoEtAl2003
        1e-12
      }
    )
    
    # SILL
    # The initial guess for the sill should be the easiest among all
    # parameters needed to fit a covariance model. Several rules are used in the
    # literature:
    # 1) the average of the maximum semivariance and the median semivariance in 
    #    the sample variogram (Hiemstra et al., 2009);
    # 2) the average of all the experimental points (Larrondo et al., 2003);
    # 3) the average  semivariance of the three last lag-distance classes
    #    (Jian et al., 1996)
    # 4) the total variance (Desassis & Renard, 2012).
    # I define the partial sill as the difference between the variance
    # of the data minus the nugget variance.
    sill <- switch(
      method,
      a = { # Samuel-Rosa
        # The average of the variance of the data and gamma at the two last 
        # lag-distance classes.
        mean(c(var(z), v$gamma[c(length(v$gamma), length(v$gamma) - 1)]))
      },
      b = { # JianEtAl1996
        mean(v$gamma[seq(length(v$gamma) - 2, length(v$gamma))])
      },
      c = { # HiemstraEtAl2009
        mean(c(max(v$gamma), median(v$gamma)))
      },
      d = { # DesassisEtAl2012
        var(z)
      },
      e = { # LarrondoEtAl2003
        mean(v$gamma)
      }
    )
    p_sill <- sill - nugget
    
    # Prepare output
    if (expand) { # output for geoR-package
      range <- c(range - by$range, range, range + by$range)
      if (range[1] <= 0) range <- range + abs(range[1]) * 1.2
      
      p_sill <- c(p_sill - by$p_sill, p_sill, p_sill + by$p_sill)
      if (p_sill[1] <= 0) p_sill <- p_sill + abs(p_sill[1]) * 1.1
      
      nugget <- c(nugget - by$nugget, nugget, nugget + by$nugget)
      if (nugget[1] <= 0) nugget <- nugget + abs(nugget[1]) * 1.1
      
      cov.pars <- as.matrix(expand.grid(p_sill = p_sill, range = range))
      res <- list(cov.pars = cov.pars, nugget = nugget)
      
    } else { # output for other routines
      res <- list(range = range, p_sill = p_sill, nugget = nugget)
    }
    
    return (res)
  }

# Plot several optimized sample configurations #################################
plotManyOSC <-
  function (osc, by = "jitter") {
    
    # Compute plot dimensions
    ylim <- sapply(1:length(osc), 
                   function (i) range(osc[[i]]@objective$energy$obj))
    
    if (by == "jitter") {
      xlim <- sapply(1:length(osc), 
                     function (i) osc[[i]]@spsann$chains$used * 
                       osc[[i]]@spsann$chains$length *
                       nrow(osc[[i]]@points))
    } else {
      xlim <- sapply(1:length(osc), function (i) osc[[i]]@spsann$chains$used)
    }
    
    # Prepare plotting area
    xlab <- ifelse(by == "jitter", "jitter", "chain")
    plot(x = "", y = "", type = "n", ylim = range(ylim), xlim = c(1, max(xlim)),
         main = osc[[1]]@objective$name, xlab = xlab, ylab = "energy state")
    
    # Prepare legend
    col <- gray(seq(0, 0.5, length.out = length(osc)))
    leg <- lapply(1:length(osc), function (i) nrow(osc[[i]]@points))
    legend("topright", legend = unlist(leg), lty = 1:length(osc), col = col)
    
    # Plot lines
    for (i in 1:length(osc)) {
      if (by == "jitter") {
        lines(osc[[i]]@objective$energy$obj, lty = i, col = col[i])
      } else {
        if (by == "chain") {
          id <- seq(osc[[i]]@spsann$chains$length * nrow(osc[[i]]@points),
                    osc[[i]]@spsann$chains$length * nrow(osc[[i]]@points) *
                      osc[[i]]@spsann$chains$used,
                    length.out = osc[[i]]@spsann$chains$used)
          id <- c(1, id + 1)
          lines(osc[[i]]@objective$energy$obj[id], lty = i, col = col[i])
        }
      }
      
    }
  }
