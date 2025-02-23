#' Prescriptive Disparity Analysis
#'
#' This function implements a presriptive analysis of group differences in an outcome, as proposed in Opacic Wei and Zhou (2025).
#' @param data_list A list of data frames, each representing an imputed dataset, or a single data frame for non-imputed data.
#' @param x A vector of characters corresponding to variable names.
#' @param a The treatment variable (a character string).
#' @param z The variable (a character string) for the conditional equalization intervention.
#' @param r The group indicator variable (a character string). A value of 1 should correspond to the disadvantaged group.
#' @param g The non-focal grouping variable for the affirmative action-type intervention (a character string).
#' @param y The outcome variable (a character string).
#' @param weight.var (Optional) The name of the survey weight variable (a character string).
#'    If not provided, all observations receive equal weight.
#' @param estimator  A character string specifying the estimator to use. Options are "ri" (regression imputation),
#'   "w" (weighting), "dr" (doubly robust estimation), or "dml" (double machine learning). Default is "ri".
#' @param K For DML estimation, the number of cross-fitting folds (default is 5).
#' @param B For parametric estimation, the number of bootstrap iterations (default is 250).
#'
#' @return A list containing 3 data frames with the results: interventional_disparity (Eq. 19 in Opacic Wei and Zhou 2025), gap_reduced (estimates of gaps reduced), and counterfactual_means (estimates of counterfactual means under the intervention).
#' @export
#'
#' @examples
#' \dontrun{
#'   # Prescriptive analysis using parametric (ri) estimator
#'   res_ri <- prescriptive(data_list = data, x = c("age", "female", "parinc"), a = "selective", z = "gpa", r = "black", g = "parinc", y = "comp",
#'                          estimator = "ri", B=250)
#'
#'   # Prescriptive analysis using DML estimator
#'   res_dml <- prescriptive(data_list = data, x = c("age", "female", "parinc"), a = "selective", z = "gpa", r = "black", g = "parinc", y = "comp",
#'                           estimator = "dml", K = 5)

#' }
prescriptive <- function(data_list, x, a, z, r, g, y,
                         weight.var = NULL,
                         estimator = "ri",
                         K = 5,
                         B = 250) {

  # -----------------------------------------
  # Checks
  # -----------------------------------------

  check_estimator <- function(estimator) {
    valid_estimators <- c("dml", "ri", "w", "dr")

    if (!estimator %in% valid_estimators) {
      stop('estimator must be one of "dml", "ri", "w", "dr')
    }
  }
  check_estimator(estimator)


  if (is.data.frame(data_list)) {
    data_list <- list(data_list)
    I <- 1
  } else if (is.list(data_list)) {
    I <- length(data_list)
  } else {
    print("Error: data must be list or dataframe")
  }


  vars_to_check <- c(a, r, g, z, y, unlist(x))
  if (any(sapply(data_list[[1]][, vars_to_check, drop = FALSE],
                 function(col) any(is.na(col) | is.nan(col))))) {
    stop("Error: Missing values (NA or NaN) detected in one or more of the variables y, r, or x.")
  }


  # -----------------------------------------
  # Setup functions
  # -----------------------------------------
  is_binary <- function(y) {
    unique_vals <- unique(y[!is.na(y)])
    length(unique_vals) == 2
  }

  trim <- function(x, min = 0.001, max = 1) {
    x[x<min] <- min
    x[x>max] <- max
    x
  }

  trimQ <- function(x, min = 0.01, max = 0.99) {
    lower <- quantile(x, min)
    upper <- quantile(x, max)
    trim(x, lower, upper)
  }

  a_form <- as.formula(paste(a, " ~ 1"))
  ar_form  <- as.formula(paste(a, " ~ r"))
  axr_form <- as.formula(paste(a, " ~ r + ", paste(c(x), collapse= "+")))
  azr_form <- as.formula(paste(a, " ~ r + ", paste(c(z), collapse= "+")))
  yaxr_form <- as.formula(paste(y, " ~ r + a + ", paste(c(x), collapse= "+")))

  # -----------------------------------------
  # Define weights function for hybrid interventions and AA-type equalization
  # -----------------------------------------

  weights_function <- function(df) {

    dfr1 <- df %>% dplyr::mutate(r = 1)
    dfr0 <- df %>% dplyr::mutate(r = 0)

    ## treatment models
    suppressWarnings({
      a_mod <- stats::glm(a_form, data = df, weights = weight, family = binomial())
      df$psA <- stats::predict(a_mod, newdata = df, type = "response")
      ar_mod <- stats::glm(ar_form, data = df, weights = weight, family = binomial())
      df$psAr0 <- stats::predict(ar_mod, newdata = dfr0, type = "response")
      df$psAr1 <- stats::predict(ar_mod, newdata = dfr1, type = "response")
      axr_mod <- stats::glm(axr_form, data = df, weights = weight, family = binomial())
      df$psAobs <- stats::predict(axr_mod, newdata = df, type = "response")
      df$psAxr0 <- stats::predict(axr_mod, newdata = dfr0, type = "response")
      df$psAxr1 <- stats::predict(axr_mod, newdata = dfr1, type = "response")
      azr_mod <- stats::glm(azr_form, data = df, weights = weight, family = binomial())
      df$psAzr0 <- stats::predict(azr_mod, newdata = dfr0, type = "response")
      df$psAzr1 <- stats::predict(azr_mod, newdata = dfr1, type = "response")
    })

    ## weights for hybrid interventions
    psA_ave <- Hmisc::wtd.mean(df$a, df$weight)
    psAxr0_ave <- Hmisc::wtd.mean(df$psAxr0, df$weight)
    psAxr1_ave <- Hmisc::wtd.mean(df$psAxr1, df$weight)
    psAzr0_ave <- Hmisc::wtd.mean(df$psAzr0, df$weight)
    psAzr1_ave <- Hmisc::wtd.mean(df$psAzr1, df$weight)
    wxr0 <- (psAxr1_ave - psA_ave)/(psAxr1_ave - psAxr0_ave)
    wxr1 <- 1 - wxr0
    wzr0 <- (psAzr1_ave - psA_ave)/(psAzr1_ave - psAzr0_ave)
    wzr1 <- 1 - wzr0

    ## AA-type equalization with and without expansion
    fA <- function(lambda, data, target){
      psAobs_new <- lambda * data$psAobs/(lambda * data$psAobs + 1 - data$psAobs)
      with(data, Hmisc::wtd.mean(psAobs_new, weight)) - target
    }

    targetA1 <- with(dplyr::filter(df, r == 0), Hmisc::wtd.mean(a, weight))
    targetA2 <- with(df, Hmisc::wtd.mean(a, weight))
    lambdaA1_r1 <- stats::uniroot(fA, interval = c(0.01, 20), data = dplyr::filter(df, r == 1),
                                  target = targetA1)$root
    lambdaA1_r0 <- 1
    lambdaA2_r1 <- stats::uniroot(fA, interval = c(0.01, 20), data = dplyr::filter(df, r == 1),
                                  target = targetA2)$root

    lambdaA2_r0 <- stats::uniroot(fA, interval = c(0.01, 20), data = dplyr::filter(df, r == 0),
                                  target = targetA2)$root

    ## AA-type equalization based on g, overall and by r
    f2A <- function(delta, data, target){
      lambda <- exp(delta[[1]] + delta[[2]] * data$g)
      psAobs_new <- lambda * data$psAobs/(lambda * data$psAobs + 1 - data$psAobs)
      zerosum <- with(data, Hmisc::wtd.mean(psAobs_new, weight)) - target
      zerocov <- with(data, Hmisc::wtd.mean(psAobs_new * g, weight) - Hmisc::wtd.mean(psAobs_new, weight) * Hmisc::wtd.mean(g, weight))
      c(zerosum, zerocov)
    }

    target_f2A <- with(df, Hmisc::wtd.mean(a, weight))
    delta_f2A <- rootSolve::multiroot(f2A, start = c(0, 0), data = df, target = target_f2A)$root
    df$lambda_f2A <- exp(delta_f2A[[1]] + delta_f2A[[2]] * df$g)
    delta_f2A_black <- rootSolve::multiroot(f2A, start = c(0, 0), data = dplyr::filter(df, r == 1),
                                            target = target_f2A)$root
    delta_f2A_white <- rootSolve::multiroot(f2A, start = c(0, 0), data = dplyr::filter(df, r == 0),
                                            target = target_f2A)$root

    df <- df %>%
      dplyr::mutate(
        lambda_f2A = exp(delta_f2A[[1]] + delta_f2A[[2]] * g),
        lambda_f2A_by_r = ifelse(r == 1,
                                 exp(delta_f2A_black[[1]] + delta_f2A_black[[2]] * g),
                                 exp(delta_f2A_white[[1]] + delta_f2A_white[[2]] * g)))
    return(list(
      df = df,
      lambdaA1_r0 = lambdaA1_r0,
      lambdaA1_r1 = lambdaA1_r1,
      lambdaA2_r0 = lambdaA2_r0,
      lambdaA2_r1 = lambdaA2_r1,
      wxr0 = wxr0,
      wxr1 = wxr1,
      wzr0 = wzr0,
      wzr1 = wzr1,
      delta_f2A = delta_f2A,
      delta_f2A_black = delta_f2A_black,
      delta_f2A_white = delta_f2A_white
    ))
  }

  # Parametric function
  parametric <- function(dfi, type = estimator) {

    boots <- matrix(NA, nrow = B, ncol = 48)

    for (b in 1:B){

      df <- dfi %>% dplyr::sample_frac(replace = TRUE)

      dfa1 <- df %>% dplyr::mutate(a = 1)
      dfa0 <- df %>% dplyr::mutate(a = 0)

      weight_output <- weights_function(df)
      df <- weight_output$df
      lambdaA1_r0 <- weight_output$lambdaA1_r0
      lambdaA1_r1 <- weight_output$lambdaA1_r1
      lambdaA2_r0 <- weight_output$lambdaA2_r0
      lambdaA2_r1 <- weight_output$lambdaA2_r1
      wxr0 <- weight_output$wxr0
      wxr1 <- weight_output$wxr1
      wzr0 <- weight_output$wzr0
      wzr1 <- weight_output$wzr1
      delta_f2A <- weight_output$delta_f2A
      delta_f2A_black <- weight_output$delta_f2A_black
      delta_f2A_white <- weight_output$delta_f2A_white

      # RI estimates of interventional disparities
      suppressWarnings({
        if (is_binary(df$y)) {
          yaxr_mod <- stats::glm(yaxr_form, data = df, weights = weight, family = binomial())

        } else {
          yaxr_mod <- stats::glm(yaxr_form, data = df, weights = weight, family = gaussian())
        } })

      df$Eyaxr <- stats::predict(yaxr_mod, newdata = df, type = "response")
      df$Eya0xr <- stats::predict(yaxr_mod, newdata = dfa0, type = "response")
      df$Eya1xr <- stats::predict(yaxr_mod, newdata = dfa1, type = "response")

      # Truncate weights at 1st and 99th percentiles
      df <- df %>%
        dplyr::mutate(
          axr_fit = psAobs,

          # Interventional propensity scores (treated as fixed)
          psAx = wxr0 * psAxr0 + wxr1 * psAxr1,
          psAz = wzr0 * psAzr0 + wzr1 * psAzr1,
          psAr0 = ifelse(r == 1, psAr0, psAobs),
          psAxr0 = ifelse(r == 1, psAxr0, psAobs),
          psAzr0 = ifelse(r == 1, psAzr0, psAobs),
          psAxr1_boost = ifelse(r == 1,
                                lambdaA1_r1 * psAobs/(lambdaA1_r1 * psAobs + 1 - psAobs),
                                lambdaA1_r0 * psAobs/(lambdaA1_r0 * psAobs + 1 - psAobs)),
          psAxr_redistr = ifelse(r == 1,
                                 lambdaA2_r1 * psAobs/(lambdaA2_r1 * psAobs + 1 - psAobs),
                                 lambdaA2_r0 * psAobs/(lambdaA2_r0 * psAobs + 1 - psAobs)),

          psAxr_f2A = lambda_f2A * psAobs/(lambda_f2A * psAobs + 1 - psAobs),
          psAxr_f2A_by_r = lambda_f2A_by_r * psAobs/(lambda_f2A_by_r * psAobs + 1 - psAobs),
          # Interventional weights (ratios of p* to p hat)
          weight_A1 = (a/axr_fit) %>% trimQ(0.01, 0.99),
          weight_Ar0 = (a * psAr0/axr_fit + (1 - a) * (1 - psAr0)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_A = (a * psA/axr_fit + (1 - a) * (1 - psA)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr0 =  (a * psAxr0/axr_fit + (1 - a) * (1 - psAxr0)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Ax = (a * psAx/axr_fit + (1 - a) * (1 - psAx)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Azr0 =  (a * psAzr0/axr_fit + (1 - a) * (1 - psAzr0)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Az =  (a * psAz/axr_fit + (1 - a) * (1 - psAz)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr1_boost = (a * psAxr1_boost/axr_fit + (1 - a) * (1 - psAxr1_boost)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr_redistr =  (a * psAxr_redistr/axr_fit + (1 - a) * (1 - psAxr_redistr)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr_f2A =  (a * psAxr_f2A/axr_fit + (1 - a) * (1 - psAxr_f2A)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr_f2A_by_r =  (a * psAxr_f2A_by_r/axr_fit + (1 - a) * (1 - psAxr_f2A_by_r)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),

          # Imputed individual-level outcomes
          y_A1 = Eya1xr,
          y_Ar0 = psAr0 * Eya1xr + (1 - psAr0) * Eya0xr,
          y_A = psA * Eya1xr + (1 - psA) * Eya0xr,
          y_Axr0 = psAxr0 * Eya1xr + (1 - psAxr0) * Eya0xr,
          y_Ax = psAx * Eya1xr + (1 - psAx) * Eya0xr,
          y_Azr0 = psAzr0 * Eya1xr + (1 - psAzr0) * Eya0xr,
          y_Az = psAz * Eya1xr + (1 - psAz) * Eya0xr,
          y_Axr1_boost = psAxr1_boost * Eya1xr + (1 - psAxr1_boost) * Eya0xr,
          y_Axr_redistr = psAxr_redistr * Eya1xr + (1 - psAxr_redistr) * Eya0xr,
          y_Axr_f2A = psAxr_f2A * Eya1xr + (1 - psAxr_f2A) * Eya0xr,
          y_Axr_f2A_by_r = psAxr_f2A_by_r * Eya1xr + (1 - psAxr_f2A_by_r) * Eya0xr,

          y_res = comp - Eyaxr,
          signal_A1 = weight_A1 * y_res + y_A1,
          signal_Ar0 = weight_Ar0 * y_res + y_Ar0,
          signal_A = weight_A * y_res + y_A,
          signal_Axr0 = weight_Axr0 * y_res + y_Axr0,
          signal_Ax = weight_Ax * y_res + y_Ax,
          signal_Azr0 = weight_Azr0 * y_res + y_Azr0,
          signal_Az = weight_Az * y_res + y_Az,
          signal_Axr1_boost = weight_Axr1_boost * y_res + y_Axr1_boost,
          signal_Axr_redistr = weight_Axr_redistr * y_res + y_Axr_redistr,
          signal_Axr_f2A = weight_Axr_f2A * y_res + y_Axr_f2A,
          signal_Axr_f2A_by_r = weight_Axr_f2A_by_r * y_res + y_Axr_f2A_by_r,

          # overwrite weights to be used in the weighting estimator
          weight_A1 = weight * (a/axr_fit) %>% trimQ(0.01, 0.99),
          weight_Ar0 = weight * (a * psAr0/axr_fit + (1 - a) * (1 - psAr0)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_A = weight * (a * psA/axr_fit + (1 - a) * (1 - psA)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr0 =  weight * (a * psAxr0/axr_fit + (1 - a) * (1 - psAxr0)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Ax = weight * (a * psAx/axr_fit + (1 - a) * (1 - psAx)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Azr0 =  weight * (a * psAzr0/axr_fit + (1 - a) * (1 - psAzr0)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Az =  weight * (a * psAz/axr_fit + (1 - a) * (1 - psAz)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr1_boost = weight * (a * psAxr1_boost/axr_fit + (1 - a) * (1 - psAxr1_boost)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr_redistr =  weight * (a * psAxr_redistr/axr_fit + (1 - a) * (1 - psAxr_redistr)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr_f2A =  weight * (a * psAxr_f2A/axr_fit + (1 - a) * (1 - psAxr_f2A)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
          weight_Axr_f2A_by_r =  weight * (a * psAxr_f2A_by_r/axr_fit + (1 - a) * (1 - psAxr_f2A_by_r)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        )

      # Weighting estimates of interventional disparities
      design_Aobs <- survey::svydesign(ids = ~ 1, data = df, weights = ~ weight)

      if (type == "w") {

        fit_svyglm <- function(data, weight_var) {
          design <- survey::svydesign(ids = ~1, data = data, weights = as.formula(paste0("~", weight_var)))
          model <- survey::svyglm(comp ~ r, design = design)
          return(model)
        }

        weight_vars <- c(
          "weight", "weight_A1", "weight_Ar0", "weight_A", "weight_Axr0",
          "weight_Ax", "weight_Azr0", "weight_Az", "weight_Axr1_boost",
          "weight_Axr_redistr", "weight_Axr_f2A", "weight_Axr_f2A_by_r")

        models <- suppressWarnings({lapply(weight_vars, function(w) fit_svyglm(df, w))})

        modexprs <- rlang::exprs(Aobs_wt, A1_wt, Ar0_wt, A_wt, Axr0_wt, Ax_wt, Azr0_wt, Az_wt,
                                 Axr1_boost_wt, Axr_redistr_wt, Axr_f2A_wt, Axr_f2A_by_r_wt)
        names <- purrr::map(modexprs, ~ stringr::str_sub(rlang::as_string(.x), end = -4))
      } else if (type == "ri") {
        fit_svyglm_ri <- function(outcome_var) {
          formula <- as.formula(paste0(outcome_var, " ~ r"))
          model <- survey::svyglm(formula, design = design_Aobs)
          return(model)
        }

        outcome_vars <- c(
          "comp", "y_A1", "y_Ar0", "y_A", "y_Axr0",
          "y_Ax", "y_Azr0", "y_Az", "y_Axr1_boost",
          "y_Axr_redistr", "y_Axr_f2A", "y_Axr_f2A_by_r")

        models <- suppressWarnings({lapply(outcome_vars, fit_svyglm_ri)})

        modexprs <- rlang::exprs(Aobs_ri, A1_ri, Ar0_ri, A_ri, Axr0_ri, Ax_ri, Azr0_ri, Az_ri,
                                 Axr1_boost_ri, Axr_redistr_ri, Axr_f2A_ri, Axr_f2A_by_r_ri)
        names <- purrr::map(modexprs, ~ stringr::str_sub(rlang::as_string(.x), end = -4))

      } else if (type == "dr") {

        fit_svyglm_dr <- function(outcome_var, design) {
          formula <- as.formula(paste0(outcome_var, " ~ r"))
          model <- survey::svyglm(formula, design = design)
          return(model)
        }

        outcome_vars_dr <- c(
          "comp", "signal_A1", "signal_Ar0", "signal_A", "signal_Axr0",
          "signal_Ax", "signal_Azr0", "signal_Az", "signal_Axr1_boost",
          "signal_Axr_redistr", "signal_Axr_f2A", "signal_Axr_f2A_by_r")

        models <- suppressWarnings({lapply(outcome_vars_dr, fit_svyglm_dr, design = design_Aobs)})

        modexprs <- rlang::exprs(Aobs_dr, A1_dr, Ar0_dr, A_dr, Axr0_dr, Ax_dr, Azr0_dr, Az_dr,
                                 Axr1_boost_dr, Axr_redistr_dr, Axr_f2A_dr, Axr_f2A_by_r_dr)
        names <- purrr::map(modexprs, ~ stringr::str_sub(rlang::as_string(.x), end = -4))
      }

      names <- purrr::map(modexprs, ~ stringr::str_sub(rlang::as_string(.x), end = -4))
      newX <- tibble::tibble(r0 = c(1, 0),
                             r1 = c(1, 1)) %>% t()

      out_df <- models %>%
        purrr::set_names(names) %>%
        tibble::enframe(name = "intervention", value = "model") %>%
        dplyr::mutate(r = replicate(nrow(.), c(0,1), simplify = FALSE),
                      est = purrr::map(model, stats::coef),
                      vcov = purrr::map(model, vcov),
                      se = purrr::map(vcov, ~ sqrt(diag(.x))),
                      pred = purrr::map(est, ~ setNames(as.double(newX %*% .x), names(est))),
                      predse = purrr::map(vcov, ~ setNames(sqrt(diag(newX %*% .x %*% t(newX))), names(est)))
        ) %>%
        dplyr::select(intervention, r, est, se, pred, predse) %>%
        tidyr::unnest(c(r, est, se, pred, predse)) %>%
        dplyr::mutate(reduced = ifelse(r == 1, est - est[intervention == "Aobs" & r == 1], NA))

      out_df <- out_df %>%
        dplyr::select(-se, -predse) %>%
        tidyr::pivot_longer(names_to = "measure", est:reduced) %>%
        dplyr::filter(!((r==0 & measure == "est")|(r==0 & measure == "reduced"))) %>%
        dplyr::arrange(estimator, intervention, measure, r) %>%
        tidyr::unite("label", c("intervention", "measure", "r"), sep = "__")

      boots[b, ] <- out_df$value
    }

    return(tibble::tibble(label = out_df$label,
                          estimate = colMeans(boots),
                          se = apply(boots, 2, sd)) %>%
             tidyr::separate(label, c("intervention", "measure", "r"), sep = "__"))
  }

  # DML function
  dml <- function(df) {

    dfa1 <- df %>% dplyr::mutate(a = 1)
    dfa0 <- df %>% dplyr::mutate(a = 0)

    weight_output <- weights_function(df)
    df <- weight_output$df
    lambdaA1_r0 <- weight_output$lambdaA1_r0
    lambdaA1_r1 <- weight_output$lambdaA1_r1
    lambdaA2_r0 <- weight_output$lambdaA2_r0
    lambdaA2_r1 <- weight_output$lambdaA2_r1
    wxr0 <- weight_output$wxr0
    wxr1 <- weight_output$wxr1
    wzr0 <- weight_output$wzr0
    wzr1 <- weight_output$wzr1
    delta_f2A <- weight_output$delta_f2A
    delta_f2A_black <- weight_output$delta_f2A_black
    delta_f2A_white <- weight_output$delta_f2A_white

    cf_fold <- caret::createFolds(df$y, K)
    main_list <- vector(mode = "list", K)

    df_xr <- model.matrix(axr_form, data = df)[, -1] %>% tibble::as_tibble()
    df_axr <- model.matrix(yaxr_form, data = df)[, -1] %>% tibble::as_tibble()
    df_a0xr <- model.matrix(yaxr_form, data = dfa0)[, -1] %>% tibble::as_tibble()
    df_a1xr <- model.matrix(yaxr_form, data = dfa1)[, -1] %>% tibble::as_tibble()

    for(k in 1:K){
      aux <- df[-cf_fold[[k]], ]
      aux_xr <- model.matrix(axr_form, data = aux)[, -1] %>% tibble::as_tibble()
      aux_axr <- model.matrix(yaxr_form, data = aux)[, -1] %>% tibble::as_tibble()

      axr_sl <- SuperLearner::SuperLearner(
        Y          = aux$a,
        X          = aux_xr,
        newX       = df_xr,
        family     = binomial(),
        obsWeights = aux$weight,
        SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"))
      df$axr_fit <- stats::predict(axr_sl, newdata = df_xr)$pred

      if (is_binary(aux$y)) {

        yaxr_sl <- SuperLearner::SuperLearner(
          Y          = aux$y,
          X          = aux_axr,
          family     = binomial(),
          obsWeights = aux$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
      } else {
        yaxr_sl <- SuperLearner::SuperLearner(
          Y          = aux$y,
          X          = aux_axr,
          family     = gaussian(),
          obsWeights = aux$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
      }
      df$Eyaxr <- stats::predict(yaxr_sl, newdata = df_axr)$pred
      df$Eya0xr <- stats::predict(yaxr_sl, newdata = df_a0xr)$pred
      df$Eya1xr <- stats::predict(yaxr_sl, newdata = df_a1xr)$pred
      main_list[[k]] <- df[cf_fold[[k]], ]
    }

    main_df <- purrr::reduce(main_list, dplyr::bind_rows) %>%
      dplyr::mutate(
        # Interventional propensity scores (treated as fixed)
        psAx = wxr0 * psAxr0 + wxr1 * psAxr1,
        psAz = wzr0 * psAzr0 + wzr1 * psAzr1,
        psAr0 = ifelse(r == 1, psAr0, psAobs),
        psAxr0 = ifelse(r == 1, psAxr0, psAobs),
        psAzr0 = ifelse(r == 1, psAzr0, psAobs),
        psAxr1_boost = ifelse(r == 1,
                              lambdaA1_r1 * psAobs/(lambdaA1_r1 * psAobs + 1 - psAobs),
                              lambdaA1_r0 * psAobs/(lambdaA1_r0 * psAobs + 1 - psAobs)),
        psAxr_redistr = ifelse(r == 1,
                               lambdaA2_r1 * psAobs/(lambdaA2_r1 * psAobs + 1 - psAobs),
                               lambdaA2_r0 * psAobs/(lambdaA2_r0 * psAobs + 1 - psAobs)),

        psAxr_f2A = lambda_f2A * psAobs/(lambda_f2A * psAobs + 1 - psAobs),
        psAxr_f2A_by_r = lambda_f2A_by_r * psAobs/(lambda_f2A_by_r * psAobs + 1 - psAobs),

        # Interventional weights (ratios of p* to p hat)
        weight_A1 = (a/axr_fit) %>% trimQ(0.01, 0.99),
        weight_Ar0 = (a * psAr0/axr_fit + (1 - a) * (1 - psAr0)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        weight_A = (a * psA/axr_fit + (1 - a) * (1 - psA)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        weight_Axr0 =  (a * psAxr0/axr_fit + (1 - a) * (1 - psAxr0)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        weight_Ax = (a * psAx/axr_fit + (1 - a) * (1 - psAx)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        weight_Azr0 =  (a * psAzr0/axr_fit + (1 - a) * (1 - psAzr0)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        weight_Az =  (a * psAz/axr_fit + (1 - a) * (1 - psAz)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        weight_Axr1_boost = (a * psAxr1_boost/axr_fit + (1 - a) * (1 - psAxr1_boost)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        weight_Axr_redistr =  (a * psAxr_redistr/axr_fit + (1 - a) * (1 - psAxr_redistr)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        weight_Axr_f2A =  (a * psAxr_f2A/axr_fit + (1 - a) * (1 - psAxr_f2A)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),
        weight_Axr_f2A_by_r =  (a * psAxr_f2A_by_r/axr_fit + (1 - a) * (1 - psAxr_f2A_by_r)/(1 - axr_fit)) %>% trimQ(0.01, 0.99),

        # Imputed individual-level outcomes
        y_A1 = Eya1xr,
        y_Ar0 = psAr0 * Eya1xr + (1 - psAr0) * Eya0xr,
        y_A = psA * Eya1xr + (1 - psA) * Eya0xr,
        y_Axr0 = psAxr0 * Eya1xr + (1 - psAxr0) * Eya0xr,
        y_Ax = psAx * Eya1xr + (1 - psAx) * Eya0xr,
        y_Azr0 = psAzr0 * Eya1xr + (1 - psAzr0) * Eya0xr,
        y_Az = psAz * Eya1xr + (1 - psAz) * Eya0xr,
        y_Axr1_boost = psAxr1_boost * Eya1xr + (1 - psAxr1_boost) * Eya0xr,
        y_Axr_redistr = psAxr_redistr * Eya1xr + (1 - psAxr_redistr) * Eya0xr,
        y_Axr_f2A = psAxr_f2A * Eya1xr + (1 - psAxr_f2A) * Eya0xr,
        y_Axr_f2A_by_r = psAxr_f2A_by_r * Eya1xr + (1 - psAxr_f2A_by_r) * Eya0xr,

        # Doubly Robust signals for interventional means
        y_res = comp - Eyaxr,
        signal_A1 = weight_A1 * y_res + y_A1,
        signal_Ar0 = weight_Ar0 * y_res + y_Ar0,
        signal_A = weight_A * y_res + y_A,
        signal_Axr0 = weight_Axr0 * y_res + y_Axr0,
        signal_Ax = weight_Ax * y_res + y_Ax,
        signal_Azr0 = weight_Azr0 * y_res + y_Azr0,
        signal_Az = weight_Az * y_res + y_Az,
        signal_Axr1_boost = weight_Axr1_boost * y_res + y_Axr1_boost,
        signal_Axr_redistr = weight_Axr_redistr * y_res + y_Axr_redistr,
        signal_Axr_f2A = weight_Axr_f2A * y_res + y_Axr_f2A,
        signal_Axr_f2A_by_r = weight_Axr_f2A_by_r * y_res + y_Axr_f2A_by_r,

        # Doubly Robust signals for changes in means
        change_Aobs = y - y,
        change_A1 = y - signal_A1,
        change_Ar0 = y - signal_Ar0,
        change_A = y - signal_A,
        change_Axr0 = y - signal_Axr0,
        change_Ax = y - signal_Ax,
        change_Azr0 = y - signal_Azr0,
        change_Az = y - signal_Az,
        change_Axr1_boost = y - signal_Axr1_boost,
        change_Axr_redistr = y - signal_Axr_redistr,
        change_Axr_f2A = y - signal_Axr_f2A,
        change_Axr_f2A_by_r = y - signal_Axr_f2A_by_r
      )

    # Estimates of interventional disparities
    design_Aobs <- survey::svydesign(ids = ~ 1, data = main_df, weights = ~ weight)
    Aobs_rm <- survey::svyglm(comp ~ black, design = design_Aobs)
    A1_rm <- survey::svyglm(signal_A1 ~ black, design = design_Aobs)
    Ar0_rm <- survey::svyglm(signal_Ar0 ~ black, design = design_Aobs)
    A_rm <- survey::svyglm(signal_A ~ black, design = design_Aobs)
    Axr0_rm <- survey::svyglm(signal_Axr0 ~ black, design = design_Aobs)
    Ax_rm <- survey::svyglm(signal_Ax ~ black, design = design_Aobs)
    Azr0_rm <- survey::svyglm(signal_Azr0 ~ black, design = design_Aobs)
    Az_rm <- survey::svyglm(signal_Az ~ black, design = design_Aobs)
    Axr1_boost_rm <- survey::svyglm(signal_Axr1_boost ~ black, design = design_Aobs)
    Axr_redistr_rm <- survey::svyglm(signal_Axr_redistr ~ black, design = design_Aobs)
    Axr_f2A_rm <- survey::svyglm(signal_Axr_f2A ~ black, design = design_Aobs)
    Axr_f2A_by_r_rm <- survey::svyglm(signal_Axr_f2A_by_r ~ black, design = design_Aobs)

    # Estimates of gap reduced
    Aobs_rd <- survey::svyglm(change_Aobs ~ black, design = design_Aobs)
    A1_rd <- survey::svyglm(change_A1 ~ black, design = design_Aobs)
    Ar0_rd <- survey::svyglm(change_Ar0 ~ black, design = design_Aobs)
    A_rd <- survey::svyglm(change_A ~ black, design = design_Aobs)
    Axr0_rd <- survey::svyglm(change_Axr0 ~ black, design = design_Aobs)
    Ax_rd <- survey::svyglm(change_Ax ~ black, design = design_Aobs)
    Azr0_rd <- survey::svyglm(change_Azr0 ~ black, design = design_Aobs)
    Az_rd <- survey::svyglm(change_Az ~ black, design = design_Aobs)
    Axr1_boost_rd <- survey::svyglm(change_Axr1_boost ~ black, design = design_Aobs)
    Axr_redistr_rd <- survey::svyglm(change_Axr_redistr ~ black, design = design_Aobs)
    Axr_f2A_rd <- survey::svyglm(change_Axr_f2A ~ black, design = design_Aobs)
    Axr_f2A_by_r_rd <- survey::svyglm(change_Axr_f2A_by_r ~ black, design = design_Aobs)

    # Output results
    modexprs <- rlang::exprs(Aobs_rm, A1_rm, Ar0_rm, A_rm, Axr0_rm, Ax_rm, Azr0_rm, Az_rm,
                             Axr1_boost_rm, Axr_redistr_rm, Axr_f2A_rm, Axr_f2A_by_r_rm)
    mods <- purrr::map(modexprs, ~ rlang::eval_tidy(.x))
    names <- purrr::map(modexprs, ~ stringr::str_sub(rlang::as_string(.x), end = -4))
    newX <- tibble::tibble(r0 = c(1, 0),
                           r1 = c(1, 1)) %>% t()

    out_df_rm <- mods %>%
      purrr::set_names(names) %>%
      tibble::enframe(name = "intervention", value = "model") %>%
      dplyr::mutate(r = replicate(nrow(.), c(0,1), simplify = FALSE),
                    est = purrr::map(model, stats::coef),
                    vcov = purrr::map(model, vcov),
                    se = purrr::map(vcov, ~ sqrt(diag(.x))),
                    pred = purrr::map(est, ~ setNames(as.double(newX %*% .x), names(est))),
                    predse = purrr::map(vcov, ~ setNames(sqrt(diag(newX %*% .x %*% t(newX))), names(est)))
      ) %>%
      dplyr::select(intervention, r, est, se, pred, predse) %>%
      tidyr::unnest(c(r, est, se, pred, predse)) %>%
      dplyr::mutate(estimand = "remaining")

    modexprs <- rlang::exprs(Aobs_rd, A1_rd, Ar0_rd, A_rd, Axr0_rd, Ax_rd, Azr0_rd, Az_rd,
                             Axr1_boost_rd, Axr_redistr_rd, Axr_f2A_rd, Axr_f2A_by_r_rd)
    mods <- purrr::map(modexprs, ~ rlang::eval_tidy(.x))
    names <- purrr::map(modexprs, ~ stringr::str_sub(rlang::as_string(.x), end = -4))

    out_df_rd <- mods %>%
      purrr::set_names(names) %>%
      tibble::enframe(name = "intervention", value = "model") %>%
      dplyr::mutate(r = replicate(nrow(.), c(0,1), simplify = FALSE),
                    est = purrr::map(model, stats::coef),
                    vcov = purrr::map(model, vcov),
                    se = purrr::map(vcov, ~ sqrt(diag(.x))),
                    pred = purrr::map(est, ~ setNames(as.double(newX %*% .x), names(est))),
                    predse = purrr::map(vcov, ~ setNames(sqrt(diag(newX %*% .x %*% t(newX))), names(est)))
      ) %>%
      dplyr::select(intervention, r, est, se, pred, predse) %>%
      tidyr::unnest(c(r, est, se, pred, predse)) %>%
      dplyr::mutate(estimand = "reduced")

    return(dplyr::bind_rows(out_df_rm, out_df_rd))
  }

  # loop over imputed datasets
  interventional_out <- vector(mode = "list", length = I)

  for (i in 1:I) {
    # Assign fixed variable names
    data <- data_list[[i]] %>%
      dplyr::mutate(
        a = !!rlang::sym(a),
        r = !!rlang::sym(r),
        g = !!rlang::sym(g),
        z = !!rlang::sym(z),
        y = !!rlang::sym(y)
      )

    # Assign weight = 1 if no sample weight indicated
    if (is.null(weight.var)) {
      data <- data %>%
        dplyr::mutate(weight = 1)
    } else {
      data <- data %>%
        dplyr::mutate(weight = !!rlang::sym(weight.var))
    }

    if (estimator == "dml") {
      interventional_out[[i]] <- dml(data)
    } else {
      interventional_out[[i]] <- parametric(data, type = estimator)
    }
  }

  if (!inherits(interventional_out, "list")) {
    interventional_out <- list(interventional_out)
  }


  if (estimator == "dml") {
    interventional_out_df <- interventional_out %>%
      purrr::imap(., ~ dplyr::mutate(.x, index = .y)) %>%
      purrr::reduce(., dplyr::bind_rows) %>%
      dplyr::group_by(intervention, r, estimand) %>%
      dplyr::summarise(coef_est = mean(est),
                       pred_est = mean(pred),
                       var_est = if (dplyr::n() > 1) var(est) else 0,
                       var_pred = if (dplyr::n() > 1) var(pred) else 0,
                       coef_se = sqrt(mean(se^2) + (1 + 1/I)*var_est),
                       pred_se = sqrt(mean(predse^2) + (1 + 1/I)*var_pred)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-var_est, -var_pred) %>%
      tidyr::pivot_longer(names_to = "measure", coef_est:pred_se) %>%
      tidyr::separate(measure, c("type", "measure")) %>%
      tidyr::pivot_wider(names_from = "measure", values_from = "value") %>%
      dplyr::filter((type=="coef" & r == 1)|(type == "pred" & estimand == "remaining")) %>%
      dplyr::mutate(
        type = dplyr::case_when(
          estimand == "reduced" ~ "reduced",
          type == "coef" ~ "gap",
          type == "pred" ~ "pred"
        ),
        est = ifelse(type == "pred", est, -est),
      ) %>%
      dplyr::select(-estimand) %>%
      dplyr::filter(intervention %in% c("Aobs", "A1",
                                        "A", "Az", "Axr_redistr",
                                        "Axr_f2A", "Axr_f2A_by_race")) %>%
      dplyr::mutate(intervention = factor(intervention,
                                          levels = c("Aobs", "A1",
                                                     "A", "Az", "Axr_redistr",
                                                     "Axr_f2A", "Axr_f2A_by_race"),
                                          labels = c("Observed Gap", "Deterministic Equalization",
                                                     "Lottery-type Marginal Equalization",
                                                     "Lottery-type Conditional Equalization",
                                                     "AA-type Equalization",
                                                     "AA-type Equalization by Class",
                                                     "AA-type Equalization by Class and Race"), ordered = TRUE)) %>%
      dplyr::arrange(factor(intervention, levels = c("Observed Gap",
                                                     "Deterministic Equalization",
                                                     "Lottery-type Marginal Equalization",
                                                     "Lottery-type Conditional Equalization",
                                                     "AA-type Equalization",
                                                     "AA-type Equalization by Class",
                                                     "AA-type Equalization by Class and Race")))

    gaps_out_df <- interventional_out_df %>%
      dplyr::filter(type == "gap") %>%
      dplyr::select(-type, -r)

    reduced_out_df <- interventional_out_df %>%
      dplyr::filter(type == "reduced") %>%
      dplyr::select(-type, -r) %>% dplyr::filter(intervention != "Observed Gap")

    preds_out_df <- interventional_out_df %>%
      dplyr::filter(type == "pred") %>%
      dplyr::select(-type)

    out <- list(gaps_out_df, reduced_out_df, preds_out_df)
    names(out) <- c("interventional_disparity", "gap_reduced", "counterfactual_means")
    return(out)

  } else {

    interventional_out_df <- interventional_out %>%
      purrr::imap(., ~ dplyr::mutate(.x, index = .y)) %>%
      purrr::reduce(., dplyr::bind_rows) %>%
      dplyr::group_by(intervention, r, measure) %>%
      dplyr::summarise(est = mean(estimate),
                       var_est = if (dplyr::n() > 1) var(estimate) else 0,
                       se = sqrt(mean(se^2) + (1 + 1/I)*var_est),
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(intervention %in% c("Aobs", "A1",
                                        "A", "Az", "Axr_redistr",
                                        "Axr_f2A", "Axr_f2A_by_race")) %>%
      dplyr::mutate(type = factor(measure, levels = c("est", "pred", "reduced"),
                                  labels = c("coef", "pred", "reduced")),
                    intervention = factor(intervention,
                                          levels = c("Aobs", "A1",
                                                     "A", "Az", "Axr_redistr",
                                                     "Axr_f2A", "Axr_f2A_by_race"),
                                          labels = c("Observed Gap", "Deterministic Equalization",
                                                     "Lottery-type Marginal Equalization",
                                                     "Lottery-type Conditional Equalization",
                                                     "AA-type Equalization",
                                                     "AA-type Equalization by Class",
                                                     "AA-type Equalization by Class and Race"), ordered = TRUE)) %>%
      dplyr::arrange(factor(intervention, levels = c("Observed Gap",
                                                     "Deterministic Equalization",
                                                     "Lottery-type Marginal Equalization",
                                                     "Lottery-type Conditional Equalization",
                                                     "AA-type Equalization",
                                                     "AA-type Equalization by Class",
                                                     "AA-type Equalization by Class and Race"))) %>%
      dplyr::select(-measure, -var_est)


    gaps_out_df <- interventional_out_df %>%
      dplyr::filter(type == "coef") %>%
      dplyr::select(-type, -r)

    reduced_out_df <- interventional_out_df %>%
      dplyr::filter(type == "reduced") %>%
      dplyr::select(-type, -r) %>% dplyr::filter(intervention != "Observed Gap")

    preds_out_df <- interventional_out_df %>%
      dplyr::filter(type == "pred") %>%
      dplyr::select(-type)

    out <- list(gaps_out_df, reduced_out_df, preds_out_df)
    names(out) <- c("interventional_disparity", "gap_reduced", "counterfactual_means")
    return(out)
  }
}
