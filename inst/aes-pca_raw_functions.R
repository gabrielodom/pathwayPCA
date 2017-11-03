######  The AES_PCA Function  ######

#==============================================================================
# Input: x - either the data matrix (when type = "predictor") or the Gram
#            matrix (when type = "Gram"). The data matrix should be n * p,
#            where n is the sample size and p is the number of variables.
#        n - sample size, needed especially when x is the Gram matrix for
#            computing BIC
#        d - the number of PCs to retain. Defaults to 1.
#        type - . Defaults to "predictor" and "Gram".
#        max.iter - . Defaults to 10.
#        method - . Defaults to "SPCA".
#        adaptive - whether to use adaptive penalty. Defaults to TRUE.
#        corr - whether to use correlation matrix. Defaults to FALSE.
#        para - a vector of integers specifying how many how zero coefficients
#               should be retained for each loading, inherited from function
#               elasticnet::spca(). Defaults to NULL.
#        lambda - ridge regression penalty. Defaults to 10 ^ -4.
# Output: loadings - the estimated loadings
#         B0 - the estimated loadings of PCA
#         score - adaptive elastic net PCA score
#         oldscore - standard PCA score
#==============================================================================

aes.pca <- function(x, n, d = 1,
                    type = c("predictor", "Gram"),
                    max.iter = 10,
                    method = "SPCA",
                    adaptive = TRUE,
                    para = NULL,
                    corr = FALSE,
                    seed = NULL,
                    lambda = 0.0001){

  eps.conv <- 1e-3
  call <- match.call()
  type <- match.arg(type)      # Why?
  method <- match.arg(method)  # Again, why?
  vn <- dimnames(x)[[2]]       # colnames()
  p <- dim(x)[2]               # ncol()

  # Calculate the regularized Grammian
  xtx <- switch(type,
                predictor = {
                  # This overwrites the original data
                  x <- scale(x, center = TRUE, scale = TRUE)
                  xtx <- (t(x) %*% x + diag(lambda, p, p)) / (1 + lambda)
                },
                Gram = {
                  xtx <- x + diag(lamda,p,p)
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
  #   eigenvectors of scale(x). What about this: instead of calculating the
  #   Grammian, regularising it, taking the SVD of it, and extracting the first
  #   d eigenvectors, why don't we just extract the first d singular vectors
  #   directly from the scaled data matrix itself? The regularisation effect
  #   only inflates the singular values anyway, and it has no effect on the
  #   singular vectors in any way. Moreover, we don't even call for the eigen-
  #   values anywhere in this code, so this whole process is supurfluous

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
  tr.x <- tr(xtx)

  ###  The while() Loop  ###
  while((k < max.iter) & (diff > eps.conv)){

    k <- k + 1

    for(i in 1:d){

      # One switch? Shouldn't this be an ifelse or if()?
      xtx1 <- switch(method,
                     SPCA = {
                       xtx1 <- xtx.1
                     })
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
  score <- array(0, dim = c(dim(x)[1], d))

  for(i in 1:d){
    score[, i] <- x %*% B[, i]
  }

  oldscore <- array(0, dim = c(dim(x)[1], d))
  for(i in 1:d){
    oldscore[, i] <- x %*% A0[, i]
  }

  ###  Return  ###
  obj <- list(loadings = B,
              B0 = A0,
              score = score,
              oldscore = oldscore)
  obj

}



######  Utility Functions  ######
# The matrix trace function
tr <- function(x) sum(diag(x))

# The matrix square root function
rootmatrix <- function(x) {
  x.eigen <- eigen(x)
  d <- x.eigen$values
  # Zero-out any negative eigenvalues
  d <- (d + abs(d))/2
  v <- x.eigen$vectors
  return(v %*% diag(sqrt(d)) %*% t(v))
}

# A function that norms a matrix, but I don't understand all of it. I met with
#   James and Steven on 26 September and neither of them understood the sign
#   reversal in the last line of the internal for() loop. Before this function
#   can be extracted and documented, I need to understand why it's doing what
#   it's doing.
normalize <- function(B,d)
{
  # browser()
  # Square root of the column sums of B ^ 2 (the entries squared). Why though?
  normB <- sqrt(apply(B^2, 2, sum))
  # Replace any 0 values with 1. Once again, why? The domain is [0, infty), so
  #   shouldn't we replace 0s with the minimum of the non-zero values? Why 1?
  normB[normB == 0] <- 1
  # Divide each column of B by its "norm" value
  B <- t(t(B)/normB)
  # Replace any really small "normed" B values with 0
  B[abs(B)<1e-6] <- 0
  # This loop moves over columns: it swaps the signs of each value in a column
  #   if the first value in that column is < 0. I still don't get what the point
  #   is of what we're doing.
  for (i in 1:d) {
    # browser()
    # Check for non-"zero" values in the ith column
    tmp <- abs(B[,i])>1e-6;
    # If we have any non-"zero" entries in the ith column, then
    if(sum(tmp)!=0) {
      # we find the indices of these non-small values and
      tmp1 <- B[tmp,i]
      # reverse their sign if the sign of the first non-small entry in B[,i] is
      #   negative. What is going on?? This is to negate the sign swap from the
      #   Grammian matrix. Now that we have changed to the svd of the data, we
      #   should be able to remove it.
      B[,i] <- sign(tmp1[1])*B[,i]
    }
  }
  B
}



######  The LAR Function  ######

lars.lsa <- function(Sigma0, b0, n,
                     adaptive = TRUE,
                     type = c("lar","lasso"),
                     eps = .Machine$double.eps,
                     max.steps,
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

  if(missing(max.steps)){
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
        R <- lars::updateR(Sigma[inew, inew], R,
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

    # backsolve is from base, while backsolvet is from lars::
    Gi1 <- lars::backsolvet(R, lars::backsolvet(R, Sign))
    dropouts <- NULL
    A <- 1 / sqrt(sum(Gi1 * Sign))
    w <- A * Gi1

    if (length(active) >= m) {
      gamhat <- Cmax / A
    }	else {

      a <- drop(w %*% Sigma[active, -c(active, ignores), drop = FALSE])
      gam <- c((Cmax - C) / (A - a), (Cmax + C) / (A + a))
      gamhat <- min(gam[gam > eps], Cmax / A)

    }

    ###  Restart Inspection HERE  ###
    browser()

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
