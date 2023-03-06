
#' Calculate bivariate residuals in a poLCA model.
#' @param fit A poLCA model fit object
#' @param tol A tolerance (scalar double) to deal with zero counts
#' @param rescale_to_df In bivariate crosstables for *polytomous* variables,
#'  whether to divide the resulting 'chi-square' type value by its
#'  degrees of freedom. Defaults to true.
#' @importFrom stats as.formula
#' @importFrom stats xtabs
#' @importFrom utils combn
#' @export
bvr <- function(fit, tol = 1e-3, rescale_to_df = TRUE) {
  stopifnot(class(fit) == "poLCA")

  ov_names <- names(fit$predcell)[1:(ncol(fit$predcell) - 2)]
  ov_combn <- combn(ov_names, 2)

  get_bvr <- function(ov_pair) {
    form_obs <- as.formula(paste0("observed ~ ", ov_pair[1], " + ", ov_pair[2]))
    form_exp <- as.formula(paste0("expected ~ ", ov_pair[1], " + ", ov_pair[2]))

    counts_obs <- xtabs(form_obs, data = fit$predcell)
    counts_exp <- xtabs(form_exp, data = fit$predcell)
    counts_exp <- ifelse(counts_exp < tol, tol, counts_exp) # Prevent Inf/NaN

    bvr_df <- prod(dim(counts_exp) - 1)
    bvr_value <- sum((counts_obs - counts_exp)^2 / counts_exp)

    if(rescale_to_df) bvr_value <- bvr_value / bvr_df

    attr(bvr_value, "df") <- bvr_df

    bvr_value
  }

  bvr_pairs <- apply(ov_combn, 2, get_bvr)

  attr(bvr_pairs, "rescale_to_df") <- rescale_to_df
  attr(bvr_pairs, "class") <- "dist"
  attr(bvr_pairs, "Size") <- length(ov_names)
  attr(bvr_pairs, "Labels") <- ov_names
  attr(bvr_pairs, "Diag") <- FALSE
  attr(bvr_pairs, "Upper") <- FALSE

  bvr_pairs
}


#' Some convenience function for parametric bootstrapping
#' Given a formula, and other arguments to poLCA,
#' output a  function : data → statistics:
#'  the output function fits poLCA to data as requested and
#'  returns the statistics of interest.
#' @param f A formula to be used as input by poLCA
#' @param ... further arguments to poLCA, such as the number of classes `nclass`
#' @returns A function to be used as `statistic=` by `boot()`
#' @importFrom poLCA poLCA
#' @export
get_fitter_and_calculator <- function(f, ...) function(dat) {
  # Fit the model to the data:
  fit_boot <- poLCA(formula = f, data = dat, verbose = FALSE, ...)

  # Calculate statistics:
  as.vector(bvr(fit_boot)) # Just bivariate residuals for now
}

#' Given a poLCA fit object, output a function that generates data from it.
#' @param fit_polca A poLCA model fit object
#' @returns A function to be used as `ran.gen=` by `boot()`
#' @importFrom poLCA poLCA.predcell
#' @export
get_generator <- function(fit_polca) {
  stopifnot(class(fit_polca) == "poLCA")

  num_indicators <- ncol(fit_polca$y)
  num_obs <- fit_polca$Nobs

  function(dat, p) {
    # Generate all possible patterns (these are not in poLCA by default)
    df_pats <- fit_polca$y |> sapply(unique, simplify = FALSE) |> expand.grid()
    # Ask poLCA what it would have predicted for all patterns
    p_pred <- poLCA::poLCA.predcell(fit_polca, df_pats)[, 1]
    # Sample pattern id's
    idx_samp <- sample(1:nrow(df_pats), size = num_obs, replace = TRUE, prob = p_pred)
    # The sampled data is created by indexing the patterns data.frame:
    df_samp <- df_pats[idx_samp, ]

    df_samp
  }
}

