p <- 500
n <- 50

x_mat <- matrix(rnorm(n * p), nrow = p, ncol = n)
x_df <- data.frame(x_mat)
time_int <- rpois(n, lambda = 365 * 2)
obs_logi <- sample(
  c(FALSE, TRUE),
  size = n,
  replace = TRUE,
  prob = c(0.2, 0.8)
)

coxTrain_fun(
  x = x_df,
  y = time_int,
  censoring.status = !obs_logi
)

olsTrain_fun(
  x = x_mat,
  y = time_int
)

glmTrain_fun(
  x = x_mat,
  y = obs_logi
)
