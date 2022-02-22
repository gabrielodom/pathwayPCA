# Issue 65
# Gabriel Odom
# 2022-02-22 (it's TWOsday)

library(pathwayPCA)
data("colonSurv_df")
aespca(as.matrix(colonSurv_df[, 4:264]))

set.seed(8675309)
randCols_idx <- sample(4:ncol(colonSurv_df), size = 260, replace = FALSE)
test_mat <- as.matrix(colonSurv_df[, randCols_idx])
aespca(test_mat)
# 6 warnings

# Does this still happen if we scale the matrix first?
set.seed(8675309)
randCols_idx <- sample(4:ncol(colonSurv_df), size = 260, replace = FALSE)
testScaled_mat <- scale( as.matrix( colonSurv_df[, randCols_idx] ) )
aespca(testScaled_mat)
# Yes, and we're up to 22 warnings. Regardless, let's stick with this version 
#   because the data matrix will have been scaled by this point.


###  Dive into lars.lsa()  ###
set.seed(8675309)
randCols_idx <- sample(4:ncol(colonSurv_df), size = 260, replace = FALSE)
testScaled_mat <- scale( as.matrix( colonSurv_df[, randCols_idx] ) )

XtX <- t(testScaled_mat) %*% testScaled_mat
A_mat <- svd(XtX)$v[, 1, drop = FALSE]

lars.lsa(
  Sigma0 = XtX,
  b0 = A_mat[, 1] * sign(A_mat[1, 1]),
  n = ncol(testScaled_mat),
  type = "lasso"
)
# at k = 272, we finally see a different length for a and C. I think I have a
#   solution (adding the check for the length of "a" and "C"), but I want to now
#   try it on the full data. For the smaller example above, our "catch" happens
#   for iteration 272, 273, 274, 275, 276, and 277


testScaled2_mat <- scale( as.matrix( colonSurv_df[, -(1:3)] ) )
XtX2 <- t(testScaled2_mat) %*% testScaled2_mat
A2_mat <- svd(XtX2)$v[, 1, drop = FALSE]

test_lars <- lars.lsa(
  Sigma0 = XtX2,
  b0 = A2_mat[, 1] * sign(A2_mat[1, 1]),
  n = ncol(testScaled2_mat),
  type = "lasso"
)
# Now, our check happens for iteration 357 to 897.


# Compare the results
test2_lars <- lars.lsa(
  Sigma0 = XtX2,
  b0 = A2_mat[, 1] * sign(A2_mat[1, 1]),
  n = ncol(testScaled2_mat),
  type = "lasso"
) 

test_lars$beta.bic
test2_lars$beta.bic

all.equal(
  test_lars$beta.bic,
  test2_lars$beta.bic
)
# WHAT?? Ok, i'll take it.

dim(test_lars$beta)
dim(test2_lars$beta)

# It looks like the algorithm returns the same vector of best-fitting betas,
#   but takes another 200 iterations when we don't have the escape.


###  Back to aespca()  ###
# Let's now go back to the aespca() function and compare the results.
set.seed(8675309)
randCols_idx <- sample(4:ncol(colonSurv_df), size = 260, replace = FALSE)
test_mat <- as.matrix(colonSurv_df[, randCols_idx])

test_aespca  <- aespca(test_mat) # remember to cut the new if () statement
test2_aespca <- aespca(test_mat) # and we still get the 6 warnings, good.

all.equal(test_aespca$aesLoad, test2_aespca$aesLoad)
all.equal(test_aespca$oldLoad, test2_aespca$oldLoad)
all.equal(test_aespca$aesScore, test2_aespca$aesScore)
all.equal(test_aespca$oldScore, test2_aespca$oldScore)


###  Code comparison  ###
# # before
# if (length(active) >= m) {
#   gamhat <- Cmax / A
# }

# # after
# if ( length(C) != (m - length(dropCols_idx)) ) {
#   
#   gamhat <- Cmax / A
#   
# } else if (length(active) >= m) {
#   gamhat <- Cmax / A
# }