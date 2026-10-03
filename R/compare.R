#' @keywords internal
#' @noRd
.fit_fh <- function(formula, vardir_name, data, method = "REML") {
  data_sae <- as.data.frame(data)
  v_sym <- as.name(vardir_name)

  call_obj <- bquote(sae::mseFH(
    formula = .(formula),
    vardir  = .(v_sym),
    method  = .(method),
    data    = data_sae
  ))

  res <- eval(call_obj)

  if (isFALSE(res$est$fit$convergence)) {
    warning("The Fay-Herriot model did not converge.", call. = FALSE)
  }

  list(
    estimates = as.vector(res$est$eblup),
    mse       = as.vector(res$mse),
    fit       = res$est$fit
  )
}

#' @keywords internal
#' @noRd
.fit_sfh <- function(formula, vardir_name, proxmat, data, method = "REML") {
  data_sae <- as.data.frame(data)
  v_sym <- as.name(vardir_name)

  call_eblup <- bquote(sae::eblupSFH(
    formula = .(formula),
    vardir  = .(v_sym),
    proxmat = proxmat,
    method  = .(method),
    data    = data_sae
  ))
  res_eblup <- eval(call_eblup)

  call_mse <- bquote(sae::mseSFH(
    formula = .(formula),
    vardir  = .(v_sym),
    proxmat = proxmat,
    method  = .(method),
    data    = data_sae
  ))
  res_mse <- eval(call_mse)

  if (isFALSE(res_eblup$fit$convergence)) {
    warning("The Spatial Fay-Herriot model did not converge.", call. = FALSE)
  }

  list(
    estimates = as.vector(res_eblup$eblup),
    mse       = as.vector(res_mse$mse),
    fit       = res_eblup$fit
  )
}

# Summary row: mse = mean of the area-level MSEs;
# rmse = square root of that mean (consistent with the thesis tables).
#' @keywords internal
#' @noRd
.comparison_row <- function(model_name, est, mse, direct = NULL) {
  m <- mean(mse, na.rm = TRUE)
  data.frame(
    model = model_name,
    mse   = m,
    rmse  = sqrt(m),
    stringsAsFactors = FALSE
  )
}
