#' geoaddSAE2: Geoadditive Small Area Estimation for Area Level
#'
#' Fits area-level Geoadditive Small Area Estimation (SAE) models: a extension of the Fay-Herriot model that represents nonlinear covariate effects with and spatial variation, estimated by REML within a linear mixed model framework (via \pkg{mgcv}), with Mean Squared Error obtained by parametric bootstrap. The package also lets users compare the geoadditive model against the classical Fay-Herriot and Spatial Fay-Herriot models, both of which are provided by, and fitted with, the existing \pkg{sae} package.
#'
#' Main entry point: \code{\link{geosae}}.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(".inv_var")

#' Build a row-standardized spatial proximity matrix
#' @keywords internal
#' @noRd
build_proxmat <- function(coords, k = 4) {
  coords <- as.matrix(coords)

  if (!is.numeric(coords)) {
    stop("`coords` must contain numeric spatial coordinates.", call. = FALSE)
  }

  if (anyNA(coords)) {
    stop("`coords` must not contain missing values.", call. = FALSE)
  }

  n <- nrow(coords)

  if (n < 2) {
    stop("At least two areas are required to build a proximity matrix.", call. = FALSE)
  }

  if (length(k) != 1L || !is.numeric(k) || is.na(k) || k < 1) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }

  k <- as.integer(k)
  k <- min(k, n - 1L)
  d <- as.matrix(stats::dist(coords))
  W <- matrix(0, nrow = n, ncol = n)

  for (i in seq_len(n)) {
    d_i <- d[i, ]
    d_i[i] <- Inf
    nn <- order(d_i)[seq_len(k)]
    W[i, nn] <- 1
  }

  diag(W) <- 0
  rs <- rowSums(W)

  if (any(rs == 0)) {
    stop("At least one area has no spatial neighbours.", call. = FALSE)
  }

  # Row-standardization
  W <- W / rs

  # Final checks
  if (any(!is.finite(W))) {
    stop("The resulting proximity matrix contains non-finite values.", call. = FALSE)
  }

  W
}

#' Default number of knots (Ruppert's rule of thumb)
#' @keywords internal
#' @noRd
default_num_knots <- function(x, min_knots = 3, max_knots = 35) {
  nu <- length(unique(x[!is.na(x)]))
  k <- round(0.25 * nu)
  k <- max(min_knots, min(k, max_knots, nu - 1))
  as.integer(k)
}

#' @keywords internal
#' @noRd
.get_knots_for <- function(var_name, knots, default_k) {
  if (is.null(knots)) return(default_k)
  if (is.list(knots) || (!is.null(names(knots)) && any(nzchar(names(knots))))) {
    if (!is.null(knots[[var_name]])) return(as.integer(knots[[var_name]]))
    return(default_k)
  }
  if (length(knots) == 1) return(as.integer(knots))
  stop("`knots` must be NULL, a single number, or a named list/vector keyed by ",
       "nonlinear variable name.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.build_gam_formula <- function(response, linear_terms, nonlinear, nonlinear_k,
                               spatial, spatial_k) {
  rhs_parts <- character(0)

  if (length(linear_terms) > 0) {
    rhs_parts <- c(rhs_parts, linear_terms)
  }

  if (length(nonlinear) > 0) {
    smooth_terms <- vapply(nonlinear, function(v) {
      k <- nonlinear_k[[v]]
      sprintf('s(%s, bs = "ps", k = %d)', v, k)
    }, character(1))
    rhs_parts <- c(rhs_parts, smooth_terms)
  }

  if (!is.null(spatial)) {
    if (is.null(spatial_k)) {
      spatial_term <- sprintf('s(%s, %s, bs = "tp")', spatial[1], spatial[2])
    } else {
      spatial_term <- sprintf('s(%s, %s, bs = "tp", k = %d)', spatial[1], spatial[2], spatial_k)
    }
    rhs_parts <- c(rhs_parts, spatial_term)
  }

  if (length(rhs_parts) == 0) rhs_parts <- "1"

  stats::as.formula(paste(response, "~", paste(rhs_parts, collapse = " + ")))
}

#' @keywords internal
#' @noRd
.fmt <- function(x, digits = 4) {
  formatC(x, format = "f", digits = digits)
}

#' @keywords internal
#' @noRd
.fit_geoadditive_gam <- function(data, gam_formula, vardir_name, method) {
  data$.inv_var <- 1 / data[[vardir_name]]

  mgcv::gam(
    formula = gam_formula,
    data = data,
    weights = .inv_var,
    method = method,
    family = stats::gaussian(),
    scale = 1
  )
}

#' @keywords internal
#' @noRd
.bootstrap_mse_geoadditive <- function(fit, data, response, vardir_name,
                                       gam_formula, method, B = 100, seed = NULL,
                                       verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(data)
  Xp <- stats::predict(fit, type = "lpmatrix")
  coefs <- stats::coef(fit)
  vardir_vec <- data[[vardir_name]]

  smooth_list <- fit$smooth
  n_smooth <- length(smooth_list)
  sigma2_k <- if (n_smooth > 0) 1 / fit$sp else numeric(0)

  sq_err <- matrix(NA_real_, n, B)
  n_fail <- 0

  for (b in seq_len(B)) {
    coefs_star <- coefs
    if (n_smooth > 0) {
      for (k in seq_len(n_smooth)) {
        idx <- smooth_list[[k]]$first.para:smooth_list[[k]]$last.para
        coefs_star[idx] <- stats::rnorm(length(idx), mean = 0, sd = sqrt(sigma2_k[k]))
      }
    }
    theta_star <- as.vector(Xp %*% coefs_star)
    e_star <- stats::rnorm(n, mean = 0, sd = sqrt(vardir_vec))
    y_star <- theta_star + e_star

    data_star <- data
    data_star[[response]] <- y_star

    fit_star <- tryCatch(
      suppressWarnings(.fit_geoadditive_gam(data_star, gam_formula, vardir_name, method)),
      error = function(e) NULL
    )

    if (is.null(fit_star)) {
      n_fail <- n_fail + 1
      next
    }
    theta_hat_star <- stats::fitted(fit_star)
    sq_err[, b] <- (theta_hat_star - theta_star)^2
  }

  mse <- rowMeans(sq_err, na.rm = TRUE)
  list(mse = mse, n_fail = n_fail, B = B)
}

#' Fit an area-level Geoadditive Small Area Estimation model
#'
#' Fits area-level Geoadditive Small Area Estimation (SAE) models: a extension of the Fay-Herriot model that represents nonlinear covariate effects with and spatial variation, estimated by REML within a linear mixed model framework (via \pkg{mgcv}), with Mean Squared Error obtained by parametric bootstrap. The package also lets users compare the geoadditive model against the classical Fay-Herriot and Spatial Fay-Herriot models, both of which are provided by, and fitted with, the existing \pkg{sae} package.
#'
#' Optionally (\code{compare = TRUE}), the classical Fay-Herriot and Spatial Fay-Herriot models are also fitted -- using the existing, well-tested implementations in the \pkg{sae} package (\code{sae::mseFH}, \code{sae::eblupSFH}, \code{sae::mseSFH}) rather than re-implemented here -- so that all three models can be compared side by side.
#'
#' @param data A data frame containing the direct estimates, the known sampling variances, the covariates, and (if used) spatial coordinates.
#' @param formula A model formula \code{y ~ x1 + x2} giving the response (direct estimator) on the left-hand side and the *linear* covariates on the right-hand side. Use \code{y ~ 1} if there are no linear covariates.
#' @param vardir Name of the column in \code{data} holding the known sampling variances of the direct estimator.
#' @param nonlinear Character vector of covariate names to be modelled nonlinearly with a P-spline.
#' @param spatial Character vector of length 2 giving the names of the two spatial coordinate columns (e.g. \code{c("lat", "lon")}).
#' @param knots Optional. Controls the basis dimension (number of knots) of each nonlinear P-spline term. \code{NULL} (default) uses an automatic rule of thumb following Ruppert (2002).
#' @param spatial_k Optional basis dimension for the spatial thin-plate smooth. \code{NULL} uses \pkg{mgcv}'s default.
#' @param method Smoothing-parameter selection method passed to \code{mgcv::gam} (and to the \pkg{sae} comparison models when \code{compare = TRUE}). Default \code{"REML"}.
#' @param compare Logical. If \code{FALSE} (default), only the geoadditive model is fitted. If \code{TRUE}, comparison models are fitted and a comparison table (model | mse | rmse) is returned.
#' @param proxmat Optional spatial proximity matrix for the Spatial Fay-Herriot model. The matrix should have one row and one column for each area, with zero diagonal and row-standardized spatial weights. If \code{NULL} and \code{compare = TRUE} with \code{spatial} specified, the matrix is constructed automatically using k-nearest neighbours based on the spatial coordinates.
#' @param proxmat_k Number of nearest neighbours used to construct the automatic spatial proximity matrix when \code{proxmat = NULL}. Default is 4.
#' @param bootstrap Logical. Whether to estimate the MSE of the geoadditive predictor via parametric bootstrap (\code{TRUE}). If \code{FALSE}, an analytical model-based approximation using the standard errors of the fitted GAM is used. Default \code{TRUE}.
#' @param B Number of parametric bootstrap replicates. Default 100.
#' @param seed Optional integer seed for the bootstrap, for reproducibility.
#'
#' @returns An object of class \code{"geosae"}, a list containing the call, settings, estimation results, diagnostics, fixed-effect parameters, model comparisons (if requested), and underlying fitted model objects.
#'
#' @export
geosae <- function(data,
                   formula,
                   vardir,
                   nonlinear,
                   spatial,
                   knots       = NULL,
                   spatial_k   = NULL,
                   method      = c("REML"),
                   compare     = FALSE,
                   proxmat     = NULL,
                   proxmat_k   = 4,
                   bootstrap   = TRUE,
                   B           = 100,
                   seed        = NULL) {

  call <- match.call()
  method <- match.arg(method)

  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!inherits(formula, "formula"))stop("`formula` must be a formula, e.g. y ~ x1 + x2.", call. = FALSE)

  vardir_name <- deparse(substitute(vardir))

  if (length(vardir_name) != 1L || !vardir_name %in% names(data)) {
    stop("`vardir` must be the name of a single column present in `data`.", call. = FALSE)
  }

  if (!is.numeric(data[[vardir_name]]) || any(data[[vardir_name]] <= 0, na.rm = TRUE)) {
    stop("`vardir` must reference a positive numeric column (known sampling variances).", call. = FALSE)
  }

  response <- all.vars(formula)[1]
  if (!response %in% names(data)) {
    stop(sprintf("Response variable `%s` (left-hand side of `formula`) not found in `data`.", response),
         call. = FALSE)
  }

  linear_terms <- attr(stats::terms(formula), "term.labels")

  if (is.null(nonlinear) || length(nonlinear) == 0) {
    stop("`nonlinear` must contain at least one covariate for the Geoadditive SAE model.", call. = FALSE)
  }

  nonlinear <- as.character(nonlinear)

  missing_nl <- setdiff(nonlinear, names(data))
  if (length(missing_nl) > 0) {
    stop("`nonlinear` variable(s) not found in `data`: ", paste(missing_nl, collapse = ", "), call. = FALSE)
  }

  overlap <- intersect(nonlinear, linear_terms)
  if (length(overlap) > 0) {
    stop("Variable(s) listed in both `formula` and `nonlinear`: ", paste(overlap, collapse = ", "), ". Remove them from `formula` if they should be modelled nonlinearly.", call. = FALSE)
  }

  if (is.null(spatial) || length(spatial) != 2) {
    stop("`spatial` must be a character vector of length 2, e.g. c(\"lat\", \"lon\").", call. = FALSE)
  }

  spatial <- as.character(spatial)

  missing_sp <- setdiff(spatial, names(data))
  if (length(missing_sp) > 0) {
    stop("`spatial` variable(s) not found in `data`: ", paste(missing_sp, collapse = ", "), call. = FALSE)
  }

  if (!all(vapply(data[spatial], is.numeric, logical(1)))) {
    stop("Spatial coordinate variables must be numeric.", call. = FALSE)
  }

  nonlinear_k <- stats::setNames(vector("list", length(nonlinear)), nonlinear)

  if (length(nonlinear) > 0) {
    non_numeric <- nonlinear[
      !vapply(data[nonlinear], is.numeric, logical(1))
    ]
    if (length(non_numeric) > 0) {
      stop("`nonlinear` variables must be numeric: ", paste(non_numeric, collapse = ", "), call. = FALSE)
    }
  }

  for (v in nonlinear) {
    default_k <- default_num_knots(data[[v]])
    nonlinear_k[[v]] <- .get_knots_for(v, knots, default_k)
  }

  if (!is.numeric(B) || length(B) != 1L || is.na(B) ||
      B < 1 || B != as.integer(B)) {
    stop("`B` must be a positive integer.", call. = FALSE)
  }
  B <- as.integer(B)
  if (!is.numeric(proxmat_k) || length(proxmat_k) != 1L ||
      is.na(proxmat_k) || proxmat_k < 1 ||
      proxmat_k != as.integer(proxmat_k)) {
    stop("`proxmat_k` must be a positive integer.", call. = FALSE)
  }

  proxmat_k <- as.integer(proxmat_k)

  gam_formula <- .build_gam_formula(response, linear_terms, nonlinear, nonlinear_k, spatial, spatial_k)

  # Estimate Model
  fit <- .fit_geoadditive_gam(data, gam_formula, vardir_name, method)

  # Extract point estimate and analytical standard error
  pred <- stats::predict(fit, se.fit = TRUE)
  point_est <- as.numeric(pred$fit)

  mse_geo <- rep(NA_real_, nrow(data))
  boot_info <- NULL

  if (isTRUE(bootstrap)) {
    boot_info <- .bootstrap_mse_geoadditive(
      fit = fit, data = data, response = response, vardir_name = vardir_name,
      gam_formula = gam_formula, method = method, B = B, seed = seed
    )
    mse_geo <- boot_info$mse
  } else {
    mse_geo <- as.numeric(pred$se.fit)^2
  }

  area_id <- if ("area" %in% names(data)) {
    as.character(data$area)
  } else {
    as.character(seq_len(nrow(data)))
  }

  estimation <- data.frame(
    area            = area_id,
    direct_est      = data[[response]],
    direct_mse      = data[[vardir_name]],
    geoadditive_est = point_est,
    geoadditive_mse = mse_geo,
    stringsAsFactors = FALSE
  )
  rownames(estimation) <- NULL

  smry <- summary(fit)
  diagnostics <- list(
    method            = method,
    converged         = isTRUE(fit$converged),
    n_obs             = nrow(data),
    n_smooth_terms    = length(fit$smooth),
    edf_table         = if (length(fit$smooth) > 0) {
      data.frame(term = rownames(smry$s.table), edf = smry$s.table[, "edf"],
                 ref_df = smry$s.table[, "Ref.df"], p_value = smry$s.table[, "p-value"],
                 row.names = NULL, stringsAsFactors = FALSE)
    } else NULL,
    variance_components = if (length(fit$smooth) > 0) {
      vc_mat <- mgcv::gam.vcomp(fit, conf.lev = 0.95)$vc
      data.frame(
        term = rownames(vc_mat),
        std_dev = vc_mat[, "std.dev"],
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    } else NULL,
    r_sq_adj          = smry$r.sq,
    deviance_explained = smry$dev.expl,
    reml_score        = fit$gcv.ubre,
    bootstrap         = if (isTRUE(bootstrap)) list(B = B, n_failed = boot_info$n_fail) else NULL
  )

  pt <- smry$p.table
  parameters <- data.frame(
    term      = rownames(pt),
    estimate  = pt[, "Estimate"],
    std_error = pt[, "Std. Error"],
    p_value   = pt[, ncol(pt)],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  models <- list(geoadditive = fit)
  comparison <- NULL

  if (isTRUE(compare)) {
    rows <- list()

    # Fit Fay-Herriot
    fh_res <- .fit_fh(formula, vardir_name, data, method = method)
    models$fh <- fh_res
    rows$fh <- .comparison_row("Fay-Herriot", fh_res$estimates, fh_res$mse, data[[response]])

    estimation$fh_est <- fh_res$estimates
    estimation$fh_mse <- fh_res$mse

    if (is.null(proxmat)) proxmat <- build_proxmat(data[, spatial, drop = FALSE], k = proxmat_k)

    # Fit Spatial Fay-Herriot
    sfh_res <- .fit_sfh(formula, vardir_name, proxmat, data, method = method)
    models$sfh <- sfh_res
    rows$sfh <- .comparison_row("Spatial Fay-Herriot", sfh_res$estimates, sfh_res$mse, data[[response]])

    estimation$sfh_est <- sfh_res$estimates
    estimation$sfh_mse <- sfh_res$mse

    rows$geo <- .comparison_row("Geoadditive SAE", point_est, mse_geo, data[[response]])

    comparison <- do.call(rbind, rows)
    rownames(comparison) <- NULL
  }

  out <- list(
    call        = call,
    settings    = list(
      formula = formula, gam_formula = gam_formula, response = response,
      vardir = vardir_name, nonlinear = nonlinear, nonlinear_knots = nonlinear_k,
      spatial = spatial, spatial_k = spatial_k, method = method,
      compare = compare, bootstrap = bootstrap, B = B, seed = seed
    ),
    estimation  = estimation,
    diagnostics = diagnostics,
    parameters  = parameters,
    comparison  = comparison,
    models      = models
  )
  class(out) <- "geosae"
  out
}

#' Print method for geosae object
#'
#' @param x An object of class \code{geosae}
#' @param digits Number of decimal digits to print. Default 4.
#' @param ... Additional arguments passed to other methods.
#' @returns Invisibly returns the input object of class \code{"geosae"}.
#' The method prints the main model settings and area-level estimation
#' results to the console.
#' @export
print.geosae <- function(x, digits = 4, ...) {
  cat("=== Geoadditive Small Area Estimation Model ===\n\n")

  s <- x$settings
  if (!is.null(s$formula)) cat("Formula (linear) :", deparse(s$formula), "\n")
  cat("Nonlinear terms  :", paste(s$nonlinear, collapse = ", "), "\n")
  cat("Spatial terms    :", paste(s$spatial, collapse = ", "), "\n")
  cat("Number of Areas  :", nrow(x$estimation), "\n\n")

  print(x$estimation, digits)

  cat("\n* Note: Use summary() for model diagnostics or extract $estimation for full results.\n")
  invisible(x)
}

#' Summary method for geosae object
#'
#' @param object An object of class \code{geosae}
#' @param ... Additional arguments passed to other methods.
#' @returns Invisibly returns the input object of class \code{"geosae"}.
#' The method prints model diagnostics, fixed-effect parameter estimates,
#' smooth-term information, MSE and RMSE summaries, and model comparison
#' results when available.
#' @export
summary.geosae <- function(object, ...) {
  cat("=== Summary: Geoadditive Small Area Estimation Model ===\n\n")

  s <- object$settings

  if (!is.null(s$nonlinear) && length(s$nonlinear) > 0) {
    cat("Nonlinear Terms & Knots:\n")
    for (v in s$nonlinear) {
      k <- s$nonlinear_knots[[v]]
      if (!is.null(k)) { cat(sprintf("  - %s : k = %d\n", v, k))}
      else { cat(sprintf("  - %s\n", v))}
    }
  }

  d <- object$diagnostics

  cat("\nModel Diagnostics:\n")

  if (!is.null(d$converged)) {
    cat(sprintf("  Convergence     : %s\n", d$converged))
  }

  if (!is.null(d$deviance_explained)) {
    cat(sprintf("  Deviance Expl.  : %.2f%%\n", 100 * d$deviance_explained))
  }

  if (!is.null(d$r_sq_adj)) {
    cat(sprintf("  R-squared (adj) : %.4f\n", d$r_sq_adj))
  }

  cat("\nFixed-Effect Parameters:\n")

  if (!is.null(object$parameters) &&
      nrow(object$parameters) > 0) {
    print(object$parameters, row.names = FALSE)
  } else {
    cat("  No fixed-effect parameters available.\n")
  }

  if (!is.null(d$edf_table) &&
      nrow(d$edf_table) > 0) {
    cat("\nSmooth Terms:\n")
    print(d$edf_table, row.names = FALSE)
    cat("\n")
  }

  mse_type <- if (isTRUE(s$bootstrap)) { "Bootstrap" }
  else {"Analytical Model-Based"}

  cat(sprintf("MSE Summary across areas (%s):\n", mse_type))

  mse <- object$estimation$geoadditive_mse
  rmse <- sqrt(mse)

  err_summary <- data.frame(
    Min = c(min(mse, na.rm = TRUE), min(rmse, na.rm = TRUE)),
    Mean = c(mean(mse, na.rm = TRUE), mean(rmse, na.rm = TRUE)),
    Max = c(max(mse, na.rm = TRUE), max(rmse, na.rm = TRUE))
  )

  rownames(err_summary) <- c("MSE", "RMSE")

  print(err_summary)

  if (!is.null(object$comparison)) {
    cat("\nModel Comparison:\n")
    print(object$comparison, row.names = FALSE)
  }

  invisible(object)
}
