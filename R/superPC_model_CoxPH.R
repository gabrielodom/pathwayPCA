#' Train Cox Proportional Hazards model for supervised PCA
#'
#' @description Main and utility functions for training the Cox PH model.
#'
#' @param x A "tall" pathway data frame (\eqn{p \times n}).
#' @param y A response vector of follow-up / event times.
#' @param censoring.status A censoring vector.
#' @param s0.perc A stabilization parameter. This is an optional argument to
#'   each of the functions called internally. Defaults to \code{NULL}.
#'
#' @return A list containing:
#' \itemize{
#'   \item{\code{tt} : }{The scaled p-dimensional score vector: each value has
#'      been divided by the respective standard deviation plus the \code{fudge}
#'      value.}
#'   \item{\code{numer} : }{The original p-dimensional score vector. From the
#'      internal \code{.coxscor} function.}
#'   \item{\code{sd} : }{The standard deviations of the scores. From the
#'      internal \code{.coxvar} function.}
#'   \item{\code{fudge} : }{A regularization scalar added to the standard
#'      deviation. If \code{s0.perc} is supplied,
#'      \code{fudge = quantile(sd, s0.perc)}.}
#' }
#'
#' @details See \url{https://web.stanford.edu/~hastie/Papers/spca_JASA.pdf},
#'    Section 5, for a description of Supervised PCA applied to survival data.
#'    The internal utility functions defined in this file (\code{.coxscor},
#'    \code{.coxvar}, and \code{.coxstuff}) are not called anywhere else, other
#'    than in the \code{coxTrain_fun} function itself. Therefore, we do not
#'    document these functions.
#'
#'    NOTE: No missing values allowed.
#'
#' @keywords internal
#'
#' @importFrom stats median
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use SuperPCA_pVals() instead

coxTrain_fun <- function(x, y, censoring.status, s0.perc = NULL){
  # browser()

  junk <- .coxscor(x, y, censoring.status)
  scor <- junk$scor
  coxVar_num <- .coxvar(
    x, y, censoring.status,
    coxstuff.obj = junk$coxstuff.obj
  )
  # The .coxvar() function can return practical 0s that are numerically negative
  coxVar_num[coxVar_num < 0] <- 0
  sd <- sqrt(coxVar_num)

  if(is.null(s0.perc)) fudge <- median(sd)
  if(!is.null(s0.perc)){

    if(0 <= s0.perc && s0.perc <= 1){
      fudge <- quantile(sd,s0.perc)
    } else {

      fudge <- 0
      warning("s0.perc parameter must in [0,1].")

    }

  }

  tt <- scor / (sd + fudge)

  return(list(tt = tt, numer = scor, sd = sd, fudge = fudge))

}


.coxscor <- function(X1, Y1, IC1, Offset1 = rep(0, length(Y1))){
  # browser()

  # computes cox scor function for rows of nx by n matrix  X1
  # first put everything in time order
  n <- length(Y1)
  nx <- nrow(X1)
  yy <- Y1 + (IC1 == 0) * (1e-05)
  otag <- order(yy)
  Y1 <- Y1[otag]
  IC1 <- IC1[otag]
  X1 <- X1[, otag, drop = FALSE]
  #compute  unique failure times, d=# of deaths at each failure time,
  #dd= expanded version of d to length n, s=sum of covariates at each
  # failure time, nn=#obs in each risk set, nno=sum(exp(Offset1)) at each failure time
  Offset1 <- Offset1[otag]
  a <- .coxstuff(X1, Y1, IC1, Offset1)
  nf <- a$nf
  fail.times <- a$fail.times
  s <- a$s
  d <- a$d
  dd <- a$dd
  nn <- a$nn
  nno <- a$nno
  w <- rep(0, nx)

  for(i in (1:nf)){

    w <- w + s[, i]
    oo <- (1:n)[Y1 >= fail.times[i]]
    r <- rowSums(X1[, oo, drop = F] * exp(Offset1[oo]))
    w <- w - (d[i] / nno[i]) * r

  }

  return(list(scor = w, coxstuff.obj = a))
}


.coxvar <- function(X2, Y2, IC2,
                   Offset2 = rep(0, length(Y2)),
                   coxstuff.obj = NULL){

  # computes information elements (var) for cox
  # X2 is nx by n matrix of expression  values
  nx <- nrow(X2)
  n <- length(Y2)
  yy <- Y2 + (IC2 == 0) * (1e-06)
  otag <- order(yy)
  Y2 <- Y2[otag]
  IC2 <- IC2[otag]
  X2 <- X2[, otag, drop = FALSE]
  Offset2 <- Offset2[otag]

  if(is.null(coxstuff.obj)) {
    coxstuff.obj <- .coxstuff(X2, Y2, IC2, Offset2)
  }

  nf <- coxstuff.obj$nf
  fail.times <- coxstuff.obj$fail.times
  s <- coxstuff.obj$s
  d <- coxstuff.obj$d
  dd <- coxstuff.obj$dd
  nn <- coxstuff.obj$nn
  nno <- coxstuff.obj$nno

  X2sq <- X2 ^ 2
  oo <- (1:n)[Y2 >= fail.times[1] ]
  sx <- (1 / nno[1]) * rowSums(X2[, oo] * exp(Offset2[oo]))
  s <- (1 / nno[1]) * rowSums(X2sq[, oo] * exp(Offset2[oo]))
  w <- d[1] * (s - sx * sx)

  for(i in 2:nf){

    oo <- (1:n)[Y2 >= fail.times[i - 1] & Y2 < fail.times[i]]
    sx <- (1 / nno[i]) *
      (nno[i - 1] * sx - rowSums(X2[, oo, drop = FALSE] * exp(Offset2[oo])))
    s <- (1 / nno[i]) *
      (nno[i - 1] * s - rowSums(X2sq[, oo, drop = FALSE] * exp(Offset2[oo])))
    w <- w + d[i] * (s - sx * sx)

  }

  return(w)

}

.coxstuff <- function(X3, Y3, IC3, Offset3 = rep(0, length(Y3))){

  fail.times <- unique(Y3[IC3 == 1])
  nf <- length(fail.times)
  n <- length(Y3)
  nn <- rep(0, nf)
  nno <- rep(0, nf)

  for(i in 1:nf){

    nn[i] <- sum(Y3 >= fail.times[i])
    nno[i] <- sum(exp(Offset3)[Y3 >= fail.times[i]])

  }

  s <- matrix(0, ncol = nf, nrow = nrow(X3))
  d <- rep(0, nf)

  #expand d out to a vector of length n
  for(i in 1:nf){

    # Try this:
    # d[i] <- sum((Y3 == fail.times[i]) & (IC3 == 1))
    # At each even time, we want to count the number of events
    o <- (1:n)[(Y3 == fail.times[i]) & (IC3 == 1)]
    d[i] <- length(o)

  }

  oo <- match(Y3, fail.times)
  oo[IC3 == 0] <- NA
  oo[is.na(oo)] <- max(oo[!is.na(oo)]) + 1
  s <- t(rowsum(t(X3), oo))

  if(ncol(s) > nf) s <- s[, -ncol(s)]

  dd <- rep(0, n)

  for(j in 1:nf){
    dd[(Y3 == fail.times[j]) & (IC3 == 1)] <- d[j]
  }

  return(list(fail.times = fail.times,
              s = s, d = d, dd = dd,
              nf = nf, nn = nn, nno = nno))

}
