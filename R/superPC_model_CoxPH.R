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
  sd <- sqrt(.coxvar(x, y, censoring.status,
                    coxstuff.obj = junk$coxstuff.obj))

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


.coxscor <- function(x, y, ic, offset = rep(0, length(y))){
  # browser()

  # computes cox scor function for rows of nx by n matrix  x
  # first put everything in time order
  n <- length(y)
  nx <- nrow(x)
  yy <- y + (ic == 0) * (1e-05)
  otag <- order(yy)
  y <- y[otag]
  ic <- ic[otag]
  x <- x[, otag, drop = FALSE]
  #compute  unique failure times, d=# of deaths at each failure time,
  #dd= expanded version of d to length n, s=sum of covariates at each
  # failure time, nn=#obs in each risk set, nno=sum(exp(offset)) at each failure time
  offset <- offset[otag]
  a <- .coxstuff(x, y, ic, offset = offset)
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
    oo <- (1:n)[y >= fail.times[i]]
    r <- rowSums(x[, oo, drop = F] * exp(offset[oo]))
    w <- w - (d[i] / nno[i]) * r

  }

  return(list(scor = w, coxstuff.obj = a))
}


.coxvar <- function(x, y, ic,
                   offset = rep(0, length(y)),
                   coxstuff.obj = NULL){

  # computes information elements (var) for cox
  # x is nx by n matrix of expression  values
  nx <- nrow(x)
  n <- length(y)
  yy <- y + (ic == 0) * (1e-06)
  otag <- order(yy)
  y <- y[otag]
  ic <- ic[otag]
  x <- x[, otag, drop = FALSE]
  offset <- offset[otag]

  if(is.null(coxstuff.obj)) {
    coxstuff.obj <- .coxstuff(x, y, ic, offset = offset)
  }

  nf <- coxstuff.obj$nf
  fail.times <- coxstuff.obj$fail.times
  s <- coxstuff.obj$s
  d <- coxstuff.obj$d
  dd <- coxstuff.obj$dd
  nn <- coxstuff.obj$nn
  nno <- coxstuff.obj$nno

  x2 <- x ^ 2
  oo <- (1:n)[y >= fail.times[1] ]
  sx <- (1 / nno[1]) * rowSums(x[, oo] * exp(offset[oo]))
  s <- (1 / nno[1]) * rowSums(x2[, oo] * exp(offset[oo]))
  w <- d[1] * (s - sx * sx)

  for(i in 2:nf){

    oo <- (1:n)[y >= fail.times[i - 1] & y < fail.times[i]]
    sx <- (1 / nno[i]) *
      (nno[i - 1] * sx - rowSums(x[, oo, drop = FALSE] * exp(offset[oo])))
    s <- (1 / nno[i]) *
      (nno[i - 1] * s - rowSums(x2[, oo, drop = FALSE] * exp(offset[oo])))
    w <- w + d[i] * (s - sx * sx)

  }

  return(w)

}

.coxstuff <- function(x, y, ic, offset = rep(0, length(y))){

  fail.times <- unique(y[ic == 1])
  nf <- length(fail.times)
  n <- length(y)
  nn <- rep(0, nf)
  nno <- rep(0, nf)

  for(i in 1:nf){

    nn[i] <- sum(y >= fail.times[i])
    nno[i] <- sum(exp(offset)[y >= fail.times[i]])

  }

  s <- matrix(0, ncol = nf, nrow = nrow(x))
  d <- rep(0, nf)

  #expand d out to a vector of length n
  for(i in 1:nf){

    # Try this:
    # d[i] <- sum((y == fail.times[i]) & (ic == 1))
    # At each even time, we want to count the number of events
    o <- (1:n)[(y == fail.times[i]) & (ic == 1)]
    d[i] <- length(o)

  }

  oo <- match(y, fail.times)
  oo[ic == 0] <- NA
  oo[is.na(oo)] <- max(oo[!is.na(oo)]) + 1
  s <- t(rowsum(t(x), oo))

  if(ncol(s) > nf) s <- s[, -ncol(s)]

  dd <- rep(0, n)

  for(j in 1:nf){
    dd[(y == fail.times[j]) & (ic == 1)] <- d[j]
  }

  return(list(fail.times = fail.times,
              s = s, d = d, dd = dd,
              nf = nf, nn = nn, nno = nno))

}
