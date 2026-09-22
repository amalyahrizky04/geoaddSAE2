test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("simulated data is available", {
  data(simulated_sae)

  expect_true(is.data.frame(simulated_sae))
  expect_equal(nrow(simulated_sae), 100)
})

test_that("geosae fits simulated data", {
  data(simulated_sae)

  fit <- geosae(
    data = simulated_sae,
    formula = y ~ x1,
    vardir = vardir,
    nonlinear = "x2",
    spatial = c("lat", "lon"),
    bootstrap = TRUE,
    B = 5,
    seed = 1
  )

  expect_s3_class(fit, "geosae")
  expect_equal(nrow(fit$estimation), 100)
})

test_that("spatial requires two variables", {
  data(simulated_sae)

  expect_error(
    geosae(
      data = simulated_sae,
      formula = y ~ x1 + x2,
      vardir = vardir,
      nonlinear = "x3",
      spatial = "lat",
      bootstrap = FALSE
    )
  )
})
