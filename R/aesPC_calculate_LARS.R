#' Least Angle and Lasso Regression models
#'
#' @description These are all variants of Lasso, and provide the entire
#'   sequence of coefficients and fits, starting from zero, to the least
#'   squares fit.
#'
#' @param Sigma0 A Grammian / covariance matrix of pathway predictors
#' @param b0 An eigenvector of \code{Sigma0}.
#' @param n The sample size
#' @param type Option between "lar" and "lasso". Defaults to "lasso"
#' @param max.steps How many steps should the LAR or LASSO algorithms take?
#'   Defaults to 8 times the pathway dimension.
#' @param eps What should we consider to be numerically 0? Defaults to the
#'   machine's default error limit for doubles (\code{.Machine$double.eps}).
#' @param adaptive Ignore
#' @param para Ignore
#'
#' @return An object of class \code{"lars"}.
#'
#' @details LARS is described in detail in Efron, Hastie, Johnstone and
#'   Tibshirani (2002). With the "lasso" option, it computes the complete lasso
#'   solution simultaneously for ALL values of the shrinkage parameter in the
#'   same computational cost as a least squares fit. This function is adapted
#'   from the \code{lars()} function in the \code{lars} package to apply to
#'   covariance or Grammian pathway design matrices.
#'
#' @seealso \url{https://web.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}
#'
#' @importFrom lars updateR
#' @importFrom lars backsolvet
#' @importFrom lars downdateR
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() instead
lars.lsa <- function(Sigma0, b0, n,
                     type = c("lar","lasso"),
                     max.steps = NULL,
                     eps = .Machine$double.eps,
                     adaptive = TRUE,
                     para = NULL){

  type <- match.arg(type)
  TYPE <- switch(type, lasso = "LASSO", lar = "LAR")

  n1 <- dim(Sigma0)[1]

  Sigma <- Sigma0
  b <- b0

  if(adaptive) {
    # If for any reason b0 is a scalar, this will fail
    Sigma <- diag(abs(b)) %*% Sigma %*% diag(abs(b))
    b <- sign(b)
  }

  nm <- dim(Sigma)
  m <- nm[2]
  im <- inactive <- seq(m)

  Cvec <- drop(t(b) %*% Sigma)
  ssy <- sum(Cvec * b)

  if(is.null(max.steps)){
    max.steps <- 8 * m
  }

  beta <- matrix(0, nrow = max.steps + 1, ncol = m)
  Gamrat <- NULL
  arc.length <- NULL
  R2 <- 1
  RSS <- ssy
  first.in <- integer(m)
  active <- NULL
  actions <- as.list(seq(max.steps))
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  ignores <- NULL

  ###  The while() Loop  ###
  while((k < max.steps) & (length(active) < m)){

    action <- NULL
    k <- k + 1
    C <- Cvec[inactive]
    Cmax <- max(abs(C))

    if(!any(drops)){

      new <- abs(C) >= Cmax - eps
      C <- C[!new]
      new <- inactive[new]

      for(inew in new){
        # lars::updateR
        R <- updateR(Sigma[inew, inew], R,
                     drop(Sigma[inew, active]),
                     Gram = TRUE,
                     eps = eps)

        if(attr(R, "rank") == length(active)){

          # singularity; back out
          nR <- seq(length(active))
          R <- R[nR, nR, drop = FALSE]
          attr(R, "rank") <- length(active)
          ignores <- c(ignores, inew)
          action <- c(action, -inew)

          # END IF
        } else {

          if(first.in[inew] == 0){
            first.in[inew] <- k
          }

          active <- c(active, inew)
          Sign <- c(Sign, sign(Cvec[inew]))
          action <- c(action, inew)

          # END IF/ELSE
        }

        # END FOR
      }

      # END IF
    }	else {
      # any(drops) is TRUE
      action <- -dropid
      # END IF(!any(drops)) ELSE
    }

    # backsolve is from base, while backsolvet is from lars
    Gi1 <- backsolve(R, backsolvet(R, Sign))
    dropouts <- NULL
    A <- 1 / sqrt(sum(Gi1 * Sign))  # This can easily be NaN
    w <- A * Gi1

    if (length(active) >= m) {
      gamhat <- Cmax / A
    }	else {

      a <- drop(w %*% Sigma[active, -c(active, ignores), drop = FALSE])
      gam <- c((Cmax - C) / (A - a), (Cmax + C) / (A + a))
      gamhat <- min(gam[gam > eps], Cmax / A)

    }

    ###  Restart Inspection HERE  ###
    # browser()
    # # The issue is that A can be a vector of NaNs, and there's no escape

    if(type == "lasso"){

      dropid <- NULL
      b1 <- beta[k, active]
      z1 <- -b1 / w
      zmin <- min(z1[z1 > eps], gamhat)

      if (zmin < gamhat) {

        gamhat <- zmin
        drops <- z1 == zmin

      }	else {
        drops <- FALSE
      }

      # END IF
    }

    beta[k + 1, ] <- beta[k, ]
    beta[k + 1, active] <- beta[k + 1, active] + gamhat * w

    Cvec <- Cvec - gamhat * Sigma[, active, drop = FALSE] %*% w
    Gamrat <- c(Gamrat, gamhat / (Cmax / A))

    arc.length <- c(arc.length, gamhat)

    if (type == "lasso" && any(drops)) {

      dropid <- seq(drops)[drops]

      for(id in rev(dropid)){
        # lars::downdateR
        R <- downdateR(R, id)
      }

      dropid <- active[drops]
      beta[k + 1, dropid] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops]

      # END IF
    }

    actions[[k]] <- action
    inactive <- im[-c(active)]

    # END WHILE
  }

  beta <- beta[seq(k + 1), ]
  dff <- b - t(beta)
  RSS <- diag(t(dff) %*% Sigma %*% dff)

  if(adaptive){
    beta <- t(abs(b0) * t(beta))
  }

  # Model Information Criteria
  dof <- apply(abs(beta) > eps, 1, sum)
  BIC <- RSS + log(n) * dof
  AIC <- RSS + 2 * dof

  if(is.null(para)){

    tmp <- sort(BIC, ind = TRUE)
    beta.bic <- beta[tmp$ix[1], ]

  }	else {

    tmp1 <- dof == para
    beta.bic <- beta[tmp1, ]

  }

  ###  Return  ###
  list(beta.bic = beta.bic,
       BIC = BIC,
       beta = beta)
}
