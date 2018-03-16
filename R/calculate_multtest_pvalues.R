#' Adjust p-Values for Simple Multiple-Testing Procedures
#'
#' @description This is the \code{mt.rawp2adjp} function from the Bioconductor
#'    package \code{multtest}. We did not write this function. For more
#'    information, see
#'    \url{https://www.bioconductor.org/packages/3.7/bioc/manuals/multtest/man/multtest.pdf}
#'
#' @author Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}
#' @author Yongchao Ge, \url{yongchao.ge@@mssm.edu}
#' @author Houston Gilbert, \url{http://www.stat.berkeley.edu/~houston}
#'
#' @param rawp A vector of raw (unadjusted) p-values for each hypothesis under
#'    consideration. These could be nominal p-values, for example, from t-tables,
#'    or permutation p-values.
#' @param proc A vector of character strings containing the names of the
#'    multiple testing procedures for which adjusted p-values are to be
#'    computed. This vector should include any of the options listed in the
#'    "Details" Section. Adjusted p-values are computed for simple FWER- and
#'    FDR- controlling procedures based on a vector of raw (unadjusted) p-values.
#' @param alpha A nominal type I error rate, or a vector of error rates, used
#'    for estimating the number of true null hypotheses in the two-stage
#'    Benjamini & Hochberg procedure (\code{"TSBH"}). Default is 0.05.
#' @param na.rm An option for handling NA values in a list of raw p-values. If
#'    \code{FALSE}, the number of hypotheses considered is the length of the
#'    vector of raw p-values. Otherwise, if \code{TRUE}, the number of
#'    hypotheses is the number of raw p-values which were not \code{NA}s.
#' @param as.multtest.out Should the output match the output from the
#'    \code{mt.rawp2adjp} function? If not, the output will match the input (a
#'    vector). Defaults to \code{FALSE}.
#'
#' @details This function computes adjusted p-values for simple multiple testing
#'    procedures from a vector of raw (unadjusted) p-values. The procedures
#'    include the Bonferroni, Holm (1979), Hochberg (1988), and Sidak procedures
#'    for strong control of the family-wise Type I error rate (FWER), and the
#'    Benjamini & Hochberg (1995) and Benjamini & Yekutieli (2001) procedures
#'    for (strong) control of the false discovery rate (FDR). The less
#'    conservative adaptive Benjamini & Hochberg (2000) and two-stage Benjamini
#'    & Hochberg (2006) FDR-controlling procedures are also included.
#'
#'    The \code{proc} options are
#'    \itemize{
#'      \item{\code{"Bonferroni"} : }{Bonferroni single-step adjusted p-values
#'         for strong control of the FWER.}
#'      \item{\code{"Holm"} : }{Holm (1979) step-down adjusted p-values for
#'         strong control of the FWER.}
#'      \item{\code{"Hochberg"} : }{Hochberg (1988) step-up adjusted p-values
#'         for strong control of the FWER (for raw (unadjusted) p-values
#'         satisfying the Simes inequality).}
#'      \item{\code{"SidakSS"} : }{Sidak single-step adjusted p-values for
#'         strong control of the FWER (for positive orthant dependent test
#'         statistics).}
#'      \item{\code{"SidakSD"} : }{Sidak step-down adjusted p-values for strong
#'         control of the FWER (for positive orthant dependent test statistics).}
#'      \item{\code{"BH"} : }{Adjusted p-values for the Benjamini & Hochberg
#'         (1995) step-up FDR-controlling procedure (independent and positive
#'         regression dependent test statistics).}
#'      \item{\code{"BY"} : }{Adjusted p-values for the Benjamini & Yekutieli
#'         (2001) step-up FDR-controlling procedure (general dependency
#'         structures).}
#'      \item{\code{"ABH"} : }{Adjusted p-values for the adaptive Benjamini &
#'         Hochberg (2000) step-up FDR-controlling procedure. This method
#'         ammends the original step-up procedure using an estimate of the
#'         number of true null hypotheses obtained from p-values.}
#'      \item{\code{"TSBH"} : }{Adjusted p-values for the two-stage Benjamini &
#'         Hochberg (2006) step-up FDR-controlling procedure. This method
#'         ammends the original step-up procedure using an estimate of the
#'         number of true null hypotheses obtained from a first-pass application
#'         of \code{"BH"}. The adjusted p-values are a-dependent, therefore
#'         alpha must be set in the function arguments when using this
#'         procedure.}
#'    }
#'
#'
#' @return A vector of the same length and order as \code{rawp}, unless the user
#'    specifies that the output should match the output from the \code{multtest}
#'    package. In that case, the use should specify \code{as.multtest.out = TRUE}
#'    and this function will return output identical to that of the
#'    \code{mt.rawp2adjp} function from package \code{multtest}. That output is
#'    as follows:
#'    \itemize{
#'      \item{\code{adjp} : }{A matrix of adjusted p-values, with rows
#'         corresponding to hypotheses and columns to multiple testing
#'         procedures. Hypotheses are sorted in increasing order of their raw
#'         (unadjusted) p-values.}
#'      \item{\code{index} : }{A vector of row indices, between 1 and
#'         \code{length(rawp)}, where rows are sorted according to their raw
#'         (unadjusted) p-values. To obtain the adjusted p-values in the
#'         original data order, use \code{adjp\[order(index),\]}.}
#'      \item{\code{h0.ABH} : }{The estimate of the number of true null
#'         hypotheses as proposed by Benjamini & Hochberg (2000) used when
#'         computing adjusted p-values for the \code{"ABH"} procedure (see
#'         Dudoit et al., 2007).}
#'      \item{\code{h0.TSBH} : }{The estimate (or vector of estimates) of the
#'         number of true null hypotheses as proposed by Benjamini et al. (2006)
#'         when computing adjusted p-values for the \code{"TSBH"} procedure
#'         (see Dudoit et al., 2007).}
#'    }
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() or superPCA_pVals() instead
adjustRaw_pVals <- function (rawp,
                             proc = c("Bonferroni",
                                      "Holm",
                                      "Hochberg",
                                      "SidakSS",
                                      "SidakSD",
                                      "BH",
                                      "BY",
                                      "ABH",
                                      "TSBH"),
                             alpha = 0.05,
                             na.rm = FALSE,
                             as.multtest.out = FALSE){

  ###  Proc Setup  ###
  m <- length(rawp)

  if(na.rm){
    mgood <- sum(!is.na(rawp))
  } else {
    mgood <- m
  }

  n <- length(proc)
  a <- length(alpha)
  index <- order(rawp)
  h0.ABH <- NULL
  h0.TSBH <- NULL
  spval <- rawp[index]
  adjp <- matrix(0, m, n + 1)
  dimnames(adjp) <- list(NULL, c("rawp", proc))
  adjp[, 1] <- spval

  ###  Bonferroni  ###
  if(is.element("Bonferroni", proc)){

    tmp <- mgood * spval
    tmp[tmp > 1] <- 1
    adjp[, "Bonferroni"] <- tmp

  }

  ###  Holm  ###
  if(is.element("Holm", proc)){

    tmp <- spval
    tmp[1] <- min(mgood * spval[1], 1)

    for(i in 2:m){
      tmp[i] <- max(tmp[i - 1], min((mgood - i + 1) * spval[i], 1))
    }

    adjp[, "Holm"] <- tmp

  }

  ###  Hochberg  ###
  if(is.element("Hochberg", proc)){

    tmp <- spval

    for(i in (m - 1):1){

      tmp[i] <- min(tmp[i + 1],
                    min((mgood - i + 1) * spval[i], 1, na.rm = TRUE),
                    na.rm = TRUE)

      if(is.na(spval[i])){
        tmp[i] <- NA
      }

    }

    adjp[, "Hochberg"] <- tmp

  }

  ###  Sidak Single-Step  ###
  if(is.element("SidakSS", proc)){
    adjp[, "SidakSS"] <- 1 - (1 - spval) ^ mgood
  }

  ###  Sidak Step-Down  ###
  if(is.element("SidakSD", proc)){

    tmp <- spval
    tmp[1] <- 1 - (1 - spval[1]) ^ mgood

    for(i in 2:m){
      tmp[i] <- max(tmp[i - 1], 1 - (1 - spval[i]) ^ (mgood - i + 1))
    }

    adjp[, "SidakSD"] <- tmp

  }

  ###  Benjamini & Hochberg  ###
  if(is.element("BH", proc)){

    tmp <- spval

    for(i in (m - 1):1){

      tmp[i] <- min(tmp[i + 1],
                    min((mgood / i) * spval[i], 1, na.rm = TRUE),
                    na.rm = TRUE)

      if(is.na(spval[i])){
        tmp[i] <- NA
      }

    }

    adjp[, "BH"] <- tmp

  }

  ###  Benjamini & Yekutieli  ###
  if(is.element("BY", proc)){

    tmp <- spval
    a <- sum(1 / (1:mgood))
    tmp[m] <- min(a * spval[m], 1)

    for(i in (m - 1):1){

      tmp[i] <- min(tmp[i + 1],
                    min((mgood * a/i) * spval[i], 1, na.rm = TRUE),
                    na.rm = TRUE)

      if(is.na(spval[i])){
        tmp[i] <- NA
      }

    }

    adjp[, "BY"] <- tmp

  }

  ###  Adaptive Benjamini & Hochberg  ###
  if(is.element("ABH", proc)){

    tmp <- spval
    h0.m <- rep(0, mgood)

    for(k in 1:mgood){
      h0.m[k] <- (mgood + 1 - k) / (1 - spval[k])
    }

    grab <- min(which(diff(h0.m, na.rm = TRUE) > 0), na.rm = TRUE)
    h0.ABH <- ceiling(min(h0.m[grab], mgood))

    for(i in (m - 1):1){

      tmp[i] <- min(tmp[i + 1],
                    min((mgood / i) * spval[i], 1, na.rm = TRUE),
                    na.rm = TRUE)

      if(is.na(spval[i])){
        tmp[i] <- NA
      }

    }

    adjp[, "ABH"] <- tmp * h0.ABH / mgood

  }

  ###  Two-Stage Benjamini & Hochberg  ###
  if(is.element("TSBH", proc)){

    TS.spot <- which(proc == "TSBH")
    TSBHs <- paste("TSBH", alpha, sep = "_")
    newprocs <- append(proc, TSBHs, after = TS.spot)
    newprocs <- newprocs[newprocs != "TSBH"]
    adjp <- matrix(0, m, n + a)
    dimnames(adjp) <- list(NULL, c("rawp", newprocs))
    adjp[, 1] <- spval
    tmp <- spval

    for(i in (m - 1):1){

      tmp[i] <- min(tmp[i + 1],
                    min((mgood / i) * spval[i], 1, na.rm = TRUE),
                    na.rm = TRUE)

      if(is.na(spval[i])){
        tmp[i] <- NA
      }

    }

    h0.TSBH <- rep(0, length(alpha))
    names(h0.TSBH) <- paste("h0.TSBH", alpha, sep = "_")

    for(i in 1:length(alpha)){

      h0.TSBH[i] <- mgood - sum(tmp < alpha[i] / (1 + alpha[i]), na.rm = TRUE)
      adjp[, TS.spot + i] <- tmp * h0.TSBH[i] / mgood

    }

  }

  ###  Return  ###
  # The multtest package has the mt.rawp2adjp() function return
  if(as.multtest.out){

    list(adjp = adjp,
         index = index,
         h0.ABH = h0.ABH[1],
         h0.TSBH = h0.TSBH[1:length(alpha)])

  } else {
    adjp[order(index), ]
  }


}
