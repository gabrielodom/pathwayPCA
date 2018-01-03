###  Extreme Distribution Adjusted p-Values  ###
# Migrated from Test_supervisedPCA.R



######  Data Load  ############################################################
data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_Genesets4240_ls")
geneset <- supervised_Genesets4240_ls
data("tScores4240_mat")
data("tControl4240_mat")
rm(supervised_Tumors_df, supervised_Genesets4240_ls)



######  Largest and Smallest Absolute t-Scores  ###############################
absMax <- function(vec){
  vec[which.max(abs(vec))]
}

tScore_max <- apply(tScores4240_mat, MARGIN = 1, FUN = absMax)
tControl_max <- apply(tControl4240_mat, MARGIN = 1, FUN = absMax)



######  Setup for Likelihood Optimization  ####################################

###  Pathway Cardinality  ###
# Don't know what this does
calc_anbn <- function(length_vec){

  logn <- log(length_vec)
  an1 <- sqrt(2 * logn)
  top <- log(4 * pi) + log(logn)
  bottom <- 2 * logn
  bn1 <- an1 * (1 - 0.5 * top / bottom)

  list(an = an1, bn = bn1)

}

# Test
pathwaylength_vec <- unlist(geneset$setsize)
abn_ls <- calc_anbn(pathwaylength_vec)


###  The Likelihood Function  ###
# This function is equation (4) Chen et al (2008). The formulation here is
#   different from the paper in the following:
#   1. They define zi = (t - mu_i) / sigma_i
#   2. They divide p and (1 - p) by sigma 1 and 2, respectively.
# This likelihood is a mixture of two Gumbel Extreme Value pdfs, with mixing
#   proportion p. The values u and s are placeholders for the mean and standard
#   deviation, respectively, but they are not used as they should be: the s
#   values are used more as precision (being multiplied instead of used as a
#   divisor).
# COMPUTATIONAL NOTE: the "L-BFGS-B" optim routine requires a finite function,
#   so we put 10 ^ 200 in the ifelse() function instead of "Inf". As we are
#   attempting to minimise L, this maximum machine value is effectively Infinity
gumbelMixture <- function(p_vec, maxt_vec, an_vec, bn_vec){

  z1 <- (maxt_vec - bn_vec - p_vec["u1"]) * an_vec * p_vec["s1"]
  z2 <- (maxt_vec + bn_vec + p_vec["u2"]) * an_vec * p_vec["s2"]
  e  <- (p_vec["p"] * an_vec * p_vec["s1"]) *
    exp(-z1 - exp(-z1)) +
    ((1 - p_vec["p"]) * an_vec * p_vec["s2"]) *
    exp(z2 - exp(z2))

  ifelse(test = any(e <= 0), yes = 10 ^ 200, no = -sum(log(e)))

}

# Test at an initial value vector
p0 <- c(p = 0.5, u1 = 1, s1 = 0.5, u2 = 1, s2 = 0.5)
gumbelMixture(p_vec = p0,
        maxt_vec = tControl_max,
        an_vec = abn_ls$an,
        bn_vec = abn_ls$bn)


###  The Score Function  ###
gumbelMix_score <- function(p_vec, maxt_vec, an_vec, bn_vec){

  innard_formula <- as.formula(~ -log(
    (p * an * s1) *
      exp(-((x - bn - u1) * an * s1) - exp(-(x - bn - u1) * an * s1)) +
      ((1 - p) * an * s2) *
      exp(((x + bn + u2) * an * s2) - exp((x + bn + u2) * an * s2))
  )
  )

  lmix_fun <- deriv(innard_formula,
                    c("p", "u1", "s1", "u2", "s2"),
                    function(x, an, bn, p, u1, s1, u2, s2) NULL)

  gradient <- lmix_fun(x = maxt_vec,
                       an = an_vec, bn = bn_vec,
                       p = p_vec["p"],
                       u1 = p_vec["u1"], s1 = p_vec["s1"],
                       u2 = p_vec["u2"], s2 = p_vec["s2"])
  # Extract the gradient matrix from the gradient vector, then find the column
  #   sums of that matrix
  colSums(attr(gradient,"gradient"))

}

# Test at an initial value vector
gumbelMix_score(p_vec = p0,
                maxt_vec = tControl_max,
                an_vec = abn_ls$an,
                bn_vec = abn_ls$bn)


###  Constrained Optimization of the Likelihood  ###
# Put a constraint on the bounds of p to avoid a global minimum where p is not
#   an element of [0,1]. Otherwise, the solution will be p ~= -8k.
pOptim <- optim(par = p0,
                fn = gumbelMixture,
                gr = gumbelMix_score,
                maxt_vec = tControl_max,
                an_vec = abn_ls$an,
                bn_vec = abn_ls$bn,
                method = "L-BFGS-B",
                lower = c(0, -Inf, 0, -Inf, 0),
                upper = c(1, Inf, Inf, Inf, Inf))$par

pOptim



######  Do Things with the t-Scores  ##########################################
# What things? No idea.
# Clean, vectorised, and functional code that does some mystery calculation.
newP_fun <- function(tScore_vec, optimParams_vec, abCounts_ls){

  an_s1 <- abCounts_ls$an * optimParams_vec["s1"]
  arg_A <- -(tScore_vec - abCounts_ls$bn - optimParams_vec["u1"]) * an_s1
  A <- 1 - exp(-exp(arg_A)) # in [0, 1]

  an_s2 <- abCounts_ls$an * optimParams_vec["s2"]
  arg_B <- (tScore_vec + abCounts_ls$bn + optimParams_vec["u2"]) * an_s2
  B <- exp(-exp(arg_B)) # also in [0, 1]

  # These values will be between 0 and 1 if the "p" value is in [0, 1]
  tt1 <- optimParams_vec["p"] * A + (1 - optimParams_vec["p"]) * B
  tt2 <- 1 - tt1

  apply(cbind(tt1, tt2), MARGIN = 1, FUN = min)

}

# Because pOptim was calculated using the tControl_max vector, that is how we
#   adjust our tScores_max vector by the control data set.
newp <- newP_fun(tScore_vec = tScore_max,
                 optimParams_vec = pOptim,
                 abCounts_ls = abn_ls)


# These are the p-values per pathway as returned by the newP_fun() function
ntest <- data.frame(goterms = names(geneset$pathways),
                    setsize = pathwaylength_vec,
                    rawp = newp)
rownames(ntest) <- NULL

###  Adjust the p-Values  ###
bh <- multtest::mt.rawp2adjp(ntest$rawp, "BH")
adjustedP <- bh$adjp[order(bh$index), ]

ntest$FDR <- adjustedP[, 2]
ntest$terms <- unlist(geneset$TERMS)
ntest <- ntest[order(ntest$FDR, ntest$rawp), ]
spcaPathwayPvals_df <- ntest

devtools::use_data(spcaPathwayPvals_df)
# write.csv(ntest, "results.csv")
