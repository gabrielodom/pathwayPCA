#' Adaptive, elastic-net, sparse principal component analysis
#'
#' @description A function to perform adaptive, elastic-net, sparse principal
#'   component analysis (AES-PCA).
#'
#' @param X A pathway design matrix: the data matrix should be \eqn{n \times p},
#'    where \eqn{n} is the sample size and \eqn{p} is the number of variables
#'    included in the pathway.
#' @param d The number of principal components (PCs) to extract from the
#'    pathway. Defaults to 1.
#' @param max.iter The maximum number of times an internal \code{while()} loop
#'    can make calls to the \code{lars.lsa()} function. Defaults to 10.
#' @param eps.conv A numerical convergence threshold for the same \code{while()}
#'    loop. Defaults to 0.001.
#' @param adaptive Internal argument of the \code{lars.lsa()} function. Defaults
#'    to TRUE.
#' @param para Internal argument of the \code{lars.lsa()} function. Defaults to
#'    NULL.
#'
#' @return A list of four elements containing the loadings and projected
#'   predictors:
#'   \itemize{
#'     \item{\code{aesLoad} : }{A \eqn{d \times p} projection matrix of the
#'        \eqn{d} AES-PCs.}
#'     \item{\code{oldLoad} : }{A \eqn{d \times p} projection matrix of the
#'        \eqn{d} PCs from the singular value decomposition (SVD).}
#'     \item{\code{aesScore} : }{An \eqn{n \times d} predictor matrix: the
#'        original \eqn{n} observations loaded onto the \eqn{d} AES-PCs.}
#'     \item{\code{oldScore} : }{An \eqn{n \times d} predictor matrix: the
#'        original \eqn{n} observations loaded onto the \eqn{d} SVD-PCs.}
#'   }
#'
#' @details This function calculates the loadings and reduced-dimension
#'    predictor matrix using both the Singular Value Decomposition and AES-PCA
#'    Decomposition (as described in Efron et al (2003)) of the data matrix.
#'
#'    See \url{https://web.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}.
#'
#'    For potential enhancement details, see the comment in the "Details"
#'    section of \code{\link{normalize}}.
#'
#' @seealso \code{\link{normalize}}; \code{\link{lars.lsa}};
#'    \code{\link{ExtractAESPCs}}; \code{\link{AESPCA_pVals}}
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Call this function through AESPCA_pVals() instead.
aespca <- function(X,
                   d = 1,
                   max.iter = 10,
                   eps.conv = 1e-3,
                   adaptive = TRUE,
                   para = NULL){

  # browser()

  call <- match.call()
  vn <- colnames(X)
  p <- ncol(X)
  n <- nrow(X)


  ###  SVD  ###
  X <- scale(X, center = TRUE, scale = TRUE)
  xtx <- t(X) %*% X
  svdGram <- svd(xtx)
  A <- svdGram$v[, 1:d, drop = FALSE]


  ###  LARS Setup  ###
  # For each eigenvector, swap the signs of the vector elements if the first
  #   entry is negative. This is Equation (2.4) of Efron et al (2003).
  for(i in 1:d){
    A[, i] <- A[, i] * sign(A[1, i])
  }
  # SEE THE COMMENTS IN THE "Details" SECTION OF ?normalize


  ###  Initial Values for while() Loop  ###
  temp <- B <- A0 <- A
  dimnames(A0) <- dimnames(B) <- list(vn, paste0("PC", 1:d))
  k <- 0; diff <- 1


  ###  The while() Loop  ###
  while((k < max.iter) & (diff > eps.conv)){
    # browser()

    k <- k + 1
    LARSerror <- rep(FALSE, d)

    for(i in 1:d){

      xtx1 <- xtx
      lfit <- tryCatch(
        {
          lars.lsa(
            Sigma0 = xtx1,
            b0 = A[, i],
            n = n,
            adaptive = adaptive,
            type = "lasso",
            para = para[i]
          )
        },
        error = function(e) NULL
      )

      if(is.null(lfit)){

        LARSerror[i] <- TRUE
        next

      } else {
        B[, i] <- lfit$beta.bic
      }

    }

    # Escape while() if lars.lsa() craps out
    if(any(LARSerror)){

      message("LARS algorithm encountered an error. Using SVD instead.")

      A <- svdGram$v[, 1:d, drop = FALSE]
      for(i in 1:d){
        A[, i] <- A[, i] * sign(A[1, i])
      }
      B <- A0 <- A
      break

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

  ###  Format Results  ###
  B <- normalize(B, d)

  oldscore <- score <- matrix(0, nrow = n, ncol = d)

  for(i in 1:d){

    score[, i] <- X %*% B[, i]
    oldscore[, i] <- X %*% A0[, i]

  }


  ###  Return  ###
  obj <- list(
    aesLoad   = t(B),      # d x p
    oldLoad   = t(A0),     # d x p
    aesScore  = score,     # n x d
    oldScore  = oldscore   # n x d
  )

  obj

}
