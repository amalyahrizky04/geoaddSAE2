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

  list(
    estimates = as.vector(res_eblup$eblup),
    mse       = as.vector(res_mse$mse),
    fit       = res_eblup$fit
  )
}

#' @keywords internal
#' @noRd
.comparison_row <- function(model_name, est, mse, direct) {
  data.frame(
    model = model_name,
    mse   = mean(mse, na.rm = TRUE),
    rmse  = mean(sqrt(mse), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
