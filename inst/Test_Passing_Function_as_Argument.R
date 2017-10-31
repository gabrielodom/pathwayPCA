# Functions as Arguments

###  X  ###
X_mat <- matrix(rnorm(40, mean = 1), ncol = 2)
X_mat <- cbind(X_mat, 2 * rep(c(1, -1), each = 10), runif(20))



###  Response Vectors / Matrix  ###

# regression
regY_vec <- 2 +
  1.5 * X_mat[, 1] -
  3 * X_mat[, 2] +
  2 * X_mat[, 3] +
  rnorm(20, sd = 0.1)
plot(regY_vec)


# survival
library(survival)
survY_vec <- sapply(regY_vec, function(i){
  rnorm(1, mean = 78 + i, sd = 1.25)
})
plot(survY_vec)
censorProb <- 0.15
survDelta_vec <- sample(c(0, 1),
                        size = length(survY_vec),
                        replace = TRUE,
                        prob = c(censorProb, 1 - censorProb))
plot(survDelta_vec)
Y_surv <- Surv(time = survY_vec, event = survDelta_vec)


# logistic
binomY_vec <- round(1 / (1 + exp(-regY_vec + rnorm(length(regY_vec),
                                                   sd = 1))))
binomY_vec[6] <- 1
plot(binomY_vec)
binomY_vec <- as.factor(binomY_vec)
plot(binomY_vec)


# multinomial logistic
cuts <- quantile(regY_vec, c(0.3333, 0.6667))
jitt <- regY_vec + rnorm(20, sd = 0.5)
multinomY_vec <- sapply(regY_vec, function(x){
  if(x < cuts[1]){
    1
  } else if(x < cuts[2]){
    2
  } else {
    3
  }
})
plot(regY_vec)
plot(multinomY_vec)



###  The Modelling Function  ###
modelStuff <- function(formula = y ~ X, FUN = lm, ...){
  FUN(formula, ...)
}

# Least Squares
modelStuff(regY_vec ~ X_mat)
summary(modelStuff(regY_vec ~ X_mat))

# Survival
modelStuff(Y_surv ~ X_mat, FUN = coxph)
summary(modelStuff(Y_surv ~ X_mat, FUN = coxph))

# Logistic
# We are getting a "perfect fit" error, so go back to the binomial vector and
#   add a single mistake.
glm(binomY_vec ~ X_mat, family = binomial(link = "logit"))
modelStuff(binomY_vec ~ X_mat, FUN = glm, family = binomial)
summary(modelStuff(binomY_vec ~ X_mat, FUN = glm, family = binomial))

# Multinomial
library(nnet)
multinom(multinomY_vec ~ X_mat)
modelStuff(multinomY_vec ~ X_mat, FUN = multinom)
