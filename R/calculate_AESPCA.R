#' Adaptive, Elastic-Net PCA
#'
#' @description A function to perform adaptive, elastic-net principal component
#'   analysis (AES-PCA).
#'
#' @param x
#' @param d
#' @param max.iter
#' @param adaptive
#' @param para
#' @param seed
#' @param n
#' @param type
#' @param method
#' @param corr
#' @param lambda
#'
#' @return What does the function return?
#'
#' @details A thorough explination of how the function works
#'
#' @export
#'
#' @examples
#'    NULL
aes.pca <- function(x, n, d = 1,
                    type = c("predictor", "Gram"),
                    max.iter = 10,
                    method = "SPCA",
                    adaptive = TRUE,
                    para = NULL,
                    corr = FALSE,
                    seed = NULL,
                    lambda = 0.0001){

  # browser()

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
  # EDIT: It doesn't even matter, because the damned lars code requires an
  #   estimate of the covariance matrix. Ugh.

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
    # browser()

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
