#' Adaptive, Elastic-Net PCA
#'
#' @description A function to perform adaptive, elastic-net principal component
#'   analysis (AES-PCA).
#'
#' @param X A pathway design matrix: either the data matrix itself (when
#'   \code{type = "predictor"}) or the Gram matrix (when \code{type = "Gram"}).
#'   The data matrix should be $n x p$, where $n$ is the sample size and $p$ is
#'   the number of variables included in the pathway.
#' @param n The sample size. Needed when \code{X} is the Gram matrix (for
#'    computing BIC). Due to a coding error, it is a required argument even
#'    when \code{X} is not a Grammian. FIX THIS.
#' @param d The number of PCs to extract from the pathway. Defaults to 1.
#' @param lambda The ridge regression penalty. Defaults to 10 ^ -4.
#' @param type Is \code{X} a pathway design matrix or a Grammian? Defaults to
#'    both: \code{c("predictor", "Gram")}. FIX THIS TOO.
#' @param corr If \code{type = Gram}, is the matrix \code{X} actually a
#'    correlation matrix? Defaults to FALSE.
#' @param max.iter The maximum number of times an internal \code{while()} loop
#'    can make calls to the \code{lars.lsa()} function. Defaults to 10.
#' @param adaptive Internal argument of the \code{lars.lsa()} function. Defaults
#'    to TRUE.
#' @param para Internal argument of the \code{lars.lsa()} function. Defaults to
#'    NULL.
#'
#' @return What does the function return?
#'
#' @details A thorough explination of how the function works
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() instead
aespca <- function(X, n, d = 1,
                   lambda = 0.0001,
                   type = c("predictor", "Gram"),
                   corr = FALSE,
                   max.iter = 10,
                   adaptive = TRUE,
                   para = NULL){

  # browser()

  eps.conv <- 1e-3
  call <- match.call()
  type <- match.arg(type)      # Why?
  # method <- match.arg(method)  # Again, why?
  vn <- dimnames(X)[[2]]       # colnames()
  p <- dim(X)[2]               # ncol()

  # Calculate the regularized Grammian
  xtx <- switch(type,
                predictor = {
                  # This overwrites the original data
                  X <- scale(X, center = TRUE, scale = TRUE)
                  xtx <- (t(X) %*% X + diag(lambda, p, p)) / (1 + lambda)
                },
                Gram = {
                  xtx <- X + diag(lambda,p,p)
                }
  )
  # Calculate the SVD of the Grammian. Why? it's more efficient to calculate
  #   SVD of the scaled data matrix (if n < p)
  svdobj <- svd(xtx)
  v <- svdobj$v
  # Matrix of the first d eigenvectors
  A <- as.matrix(v[, 1:d, drop = FALSE])

  # For each eigenvector, swap the signs of the vector elements if the first
  #   entry is negative. This has popped up all over the place, so it needs its
  #   own function.
  for(i in 1:d){
    A[, i] <- A[, i] * sign(A[1, i])
  }
  # ASIDE: this "sign flipping" thing changes the eigenvectors of xtx into the
  #   eigenvectors of scale(X). What about this: instead of calculating the
  #   Grammian, regularising it, taking the SVD of it, and extracting the first
  #   d eigenvectors, why don't we just extract the first d singular vectors
  #   directly from the scaled data matrix itself? The regularisation effect
  #   only inflates the singular values anyway, and it has no effect on the
  #   singular vectors in any way. Moreover, we don't even call for the eigen-
  #   values anywhere in this code, so this whole process is supurfluous
  # EDIT: It doesn't even matter, because the damned lars code requires an
  #   estimate of the covariance matrix. Ugh.
  # EDIT 2: The lars::lars() function *can* take in the full data, instead of
  #   just a Grammian. As an enhancement, we should either update our copy of
  #   the lars() function, or make a call to the exported lars() function. We
  #   will worry about getting the rest of the package running before we spend
  #   times optimising.

  ###  Initial Values for while() Loop  ###
  # Initialize matrices of eigenvectors
  B <- A; temp <- B; A0 <- B

  # Calculate the Grammian (again?)
  xtx.1 <- xtx
  if(corr){
    xtx.1 <- xtx.1 * n
  }

  II <- diag(p) * n       # The diag() function breaks if p = 1
  k <- 0
  diff <- 1
  # tr.x <- tr(xtx)  # Not called anywhere.

  ###  The while() Loop  ###
  while((k < max.iter) & (diff > eps.conv)){
    # browser()

    k <- k + 1

    for(i in 1:d){

      # # One switch? Shouldn't this be an ifelse or if()?
      # xtx1 <- switch(method,
      #                SPCA = {
      #                  xtx1 <- xtx.1
      #                })
      xtx1 <- xtx.1
      lfit <- lars.lsa(Sigma0 = xtx1,
                       b0 = A[, i],
                       n = n,
                       adaptive = adaptive,
                       type = "lasso",
                       para = para[i])
      B[, i] <- lfit$beta.bic

      # END FOR
    }

    # Part of the internal calls of the normalize() function. See comments below
    normB <- sqrt(apply(B ^ 2, 2, sum))
    normB[normB == 0] <- 1
    B2 <- t(t(B) / normB)
    diff <- max(abs(B2 - temp))
    temp <- B
    A <- xtx %*% B
    z <- svd(A)
    A <- (z$u) %*% t(z$v)

    # END WHILE
  }

  B <- normalize(B, d)
  dimnames(B) <- list(vn, paste("PC", 1:d, sep = "")) # use paste0() instead
  dimnames(A0) <- list(vn, paste("PC", 1:d, sep = ""))
  # Initialise an empty score array. Why not initialise oldscore here too?
  score <- array(0, dim = c(dim(X)[1], d))

  for(i in 1:d){
    score[, i] <- X %*% B[, i]
  }

  oldscore <- array(0, dim = c(dim(X)[1], d))
  for(i in 1:d){
    oldscore[, i] <- X %*% A0[, i]
  }

  ###  Return  ###
  obj <- list(loadings = B,
              B0 = A0,
              score = score,
              oldscore = oldscore)
  obj

}
