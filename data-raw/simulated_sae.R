## code to prepare `simulated_sae` dataset goes here

# Number of areas
M <- 100

# Reproducibility
set.seed(123)

# Spatial coordinates
lat <- runif(M, -11, 6)
lon <- runif(M, 95, 141)

# Auxiliary variables
x1 <- rbinom(M, 1, 0.9)
x2 <- runif(M, 1, 5)
x3 <- runif(M, -pi, pi)

# Model parameters
alpha <- 0.5
beta1 <- 1
beta2 <- 1

# Nonlinear component
g_x3 <- 2 * sin(x3)

# Spatial component
h_spatial <- 1.5 * sin(pi * lat / 180) * cos(pi * lon / 180)

# True area parameter
theta <- alpha + beta1 * x1 + beta2 * x2 + g_x3 + h_spatial

# ------------------------------------------
# Different sampling variances
# ------------------------------------------

# Different sample sizes for each area
n_i <- sample(10:50, M, replace = TRUE)

# Population sampling variance
sigma_e2 <- 1

# Sampling variance for each area
vardir <- sigma_e2 / n_i

# ------------------------------------------
# Direct estimator
# ------------------------------------------

y <- theta + rnorm(M, mean = 0, sd = sqrt(vardir))

# ------------------------------------------
# Population data
# ------------------------------------------

simulated_sae <- data.frame(
  area = sprintf("Area_%02d", 1:M),
  y = y,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  lat = lat,
  lon = lon,
  vardir = vardir,
  theta = theta
)

usethis::use_data(simulated_sae, overwrite = TRUE)
