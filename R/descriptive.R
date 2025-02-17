#' Descriptive Decomposition Analysis (Sequential or Simultaneous)
#'
#' This function implements a decomposition analysis of group differences in an outcome
#' using either a sequential decomposition or a simultaneous (Oaxaca–Blinder)
#' decomposition. The user supplies covariates via x (which may be a character vector or an ordered list of character vectors).
#'
#' When \code{type = "sequential"}, the decomposition is carried out in a sequential fashion:
#' stage 1 uses x[[1]], stage 2 uses c(x[[1]], x[[2]]), etc. For each stage the contribution
#' is computed (with bootstrap standard errors) using the chosen estimator ("ri", "w", "dr",
#' or "dml").
#'
#' When \code{type = "simultaneous"}, a standard Oaxaca–Blinder (or KOB) decomposition is performed.
#' With the default \code{estimator = "ri"}, a full decomposition is obtained—that is, the individual contributions
#' of each covariate (based on a linear regression fitted among group r) are reported, along with the total explained
#' and unexplained portions. When any other estimator is selected ("w", "dr", or "dml"), only the total explained
#' proportion is computed, as detailed decomposition is not available.
#'
#' @param data A data frame. Each data set must include the variables specified in y, r, and x.
#' @param y The outcome variable (a character string).
#' @param r The group indicator variable (binary only; a character string).
#' @param x Either an ORDERED character vector (each element is treated as one stage) or an ORDERED list of character vectors.
#'   In the sequential decomposition, the stages are built sequentially.
#' @param weight.var (Optional) A character string giving the name of the survey weight variable.
#'   If not provided, all observations receive equal weight.
#' @param estimator A character string specifying the estimator to use. Options are "ri" (regression imputation, by linear or logistic regression),
#'   "w" (weighting, by logistic regression), "dr" (doubly robust), or "dml" (double machine learning). Default is "ri".
#' @param learner For the DML estimator only: a character vector of SuperLearner library functions. Default is c("SL.mean", "SL.glmnet", "SL.ranger").
#' @param type A character string indicating the decomposition type: either "sequential" (the default) or "simultaneous".
#' @param outcome.type A character string indicating the type of outcome: either "binary" (logistic regression) or "continuous" (linear regression). Default is "binary".
#' @param B For non-DML estimators, the number of bootstrap iterations (default is 250).
#' @param K For the DML estimator, the number of cross-fitting folds (default is 5).
#'
#' @return A tibble with the estimated contributions and their standard errors. For a sequential decomposition,
#'   one row is provided for each stage (and one for the total explained gap). For a simultaneous decomposition with
#'   \code{estimator = "ri"}, additional rows for individual covariate contributions and the unexplained component are included;
#'   for other estimators ("w", "dr", or "dml"), only the total explained part is reported.
#'
#' @examples
#' \dontrun{
#'   # Sequential decomposition: each covariate is added in order (non-DML)
#'   res_seq <- descriptive(data = mydata, y = "comp", r = "black", x = c("age", "female", "parinc"),
#'                          weight.var = "weight", estimator = "ri", type = "sequential")
#'
#'   # Sequential decomposition using DML:
#'   res_seq_dml <- descriptive(data = mydata, y = "comp", r = "black", x = c("age", "female", "parinc"),
#'                              weight.var = "weight", estimator = "dml", type = "sequential")
#'
#'   # Simultaneous decomposition (Oaxaca-Blinder) with RI estimator:
#'   res_sim <- descriptive(data = mydata, y = "comp", r = "black", x = c("age", "female", "parinc"),
#'                          weight.var = "weight", estimator = "ri", type = "simultaneous")
#'
#'   # Simultaneous decomposition with a DR estimator:
#'   res_sim_dr <- descriptive(data = mydata, y = "comp", r = "black", x = c("age", "female", "parinc"),
#'                             weight.var = "weight", estimator = "dr", type = "simultaneous")
#'
#'   # Simultaneous decomposition with a DML estimator:
#'   res_sim_dml <- descriptive(data = mydata, y = "comp", r = "black", x = c("age", "female", "parinc"),
#'                              weight.var = "weight", estimator = "dml", type = "simultaneous")
#' }
#'
#' @export
descriptive <- function(data, y, r, x, weight.var = NULL,
                        estimator = "ri", learner = c("SL.mean", "SL.glm", "SL.ranger"),
                        type = "sequential", outcome.type = "binary",
                        B = 250, K = 5) {

  # Ensure data is a list
  if (!is.list(data)) {
    data <- list(data)
  }

  # Process x: if not a list, convert to list (each element a stage)
  if (!is.list(x)) {
    x <- as.list(x)
  }
  if (!all(sapply(x, is.character))) {
    stop("Argument 'x' must be either a character vector or a list of character vectors.")
  }


  # Create a vector of all variable names to check (handles x as either a vector or a list)
  vars_to_check <- c(y, r, if (is.list(x)) unlist(x) else x)
  if (any(sapply(data[, vars_to_check, drop = FALSE],
                 function(col) any(is.na(col) | is.nan(col))))) {
    stop("Error: Missing values (NA or NaN) detected in one or more of the variables y, r, or x.")
  }





  # For sequential decomposition, create sequential covariate sets.
  if (type == "sequential") {
    S <- length(x)
    cumulative_x <- lapply(seq_len(S), function(s) {
      unique(unlist(x[1:s]))
    })
  }

  # For simultaneous decomposition, use the full set of covariates.
  if (type == "simultaneous") {
    if (is.list(x)) {
      x_full <- unique(unlist(x))
    } else {
      x_full <- x
    }
  }

  # Helper: determine weights vector
  get_weights <- function(d) {
    if (!is.null(weight.var)) {
      return(d[[weight.var]])
    } else {
      return(rep(1, nrow(d)))
    }
  }


  # -----------------------------------------
  # Define bootstrap statistic functions for non-DML sequential
  # -----------------------------------------

  boot_statistic_seq <- function(data, indices, y, r, cumulative_x, weight.var, estimator, outcome.type) {
    d <- data[indices,]
    w <- get_weights(d)
    d_r1 <- d
    d_r1[[r]] <- 1  # for prediction
    mu1 <- Hmisc::wtd.mean(d[[y]], w * d[[r]])
    theta <- numeric(length(cumulative_x))
    outcome_family <- if (outcome.type == "binary") stats::binomial() else stats::gaussian()
    for (s in seq_along(cumulative_x)) {
      outcome_formula <- as.formula(paste(y, "~", paste(c(r, cumulative_x[[s]]), collapse = " + ")))
      mod_out <- stats::glm(outcome_formula, data = d, weights = w, family = outcome_family)
      ypred <- stats::predict(mod_out, newdata = d_r1, type = "response")
      theta_ri <- Hmisc::wtd.mean(ypred, w * (1 - d[[r]]))
      if (estimator == "ri") {
        theta[s] <- theta_ri
      } else {
        group_formula <- as.formula(paste(r, "~", paste(cumulative_x[[s]], collapse = " + ")))
        mod_group <- stats::glm(group_formula, data = d, weights = w, family = binomial())
        ps <- stats::predict(mod_group, newdata = d, type = "response")
        w_temp <- d[[r]] / Hmisc::wtd.mean(1 - d[[r]], w) * (1 - ps) / ps
        w_temp <- trimQ(w_temp, 0.01, 0.99)
        if (estimator == "w") {
          d$w_new <- w * w_temp
          design <- survey::svydesign(ids = ~1, data = d, weights = ~w_new)
          mod_w <- survey::svyglm(as.formula(paste(y, "~ 1")), design = design)
          theta[s] <- stats::coef(mod_w)[1]
        } else if (estimator == "dr") {
          signal <- w_temp * (d[[y]] - ypred) +
            (1 - d[[r]]) / Hmisc::wtd.mean(1 - d[[r]], w) * (ypred - theta_ri) +
            theta_ri
          design <- survey::svydesign(ids = ~1, data = d %>% dplyr::mutate(signal = signal), weights = ~w)
          mod_dr <- survey::svyglm(signal ~ 1, design = design)
          theta[s] <- stats::coef(mod_dr)[1]
        }
      }
    }
    contrib <- numeric(length(cumulative_x))
    contrib[1] <- theta[1] - mu1
    if (length(theta) > 1) {
      for (s in 2:length(theta)) {
        contrib[s] <- theta[s] - theta[s - 1]
      }
    }
    total_expl <- theta[length(theta)] - mu1
    return(c(theta, contrib, total_expl))
  }

  # -----------------------------------------
  # Sequential DML branch (non-bootstrap)
  # -----------------------------------------
  sequential_dml <- function(d, y, r, cumulative_x, weight.var, learner, outcome.type, K) {
    S <- length(cumulative_x)
    w <- get_weights(d)

    # Observed mean for group r
    mu1 <- Hmisc::wtd.mean(d[[y]], w * d[[r]])

    # Create folds for cross-fitting
    folds <- caret::createFolds(d[[y]], k = K)
    pred_list <- list()

    for (fold in folds) {
      train <- d[-fold, ]
      valid <- d[fold, ]

      # For each stage, fit outcome and group membership models via SuperLearner
      for (s in seq_along(cumulative_x)) {
        # Outcome model: fit only on train subset for observations with r == 1

        train_r <- train[train[[r]] == 1, ]

        outcome_SL <- SuperLearner::SuperLearner(
          Y = train_r[[y]],
          X = train_r[, cumulative_x[[s]], drop = FALSE],
          newX = valid[, cumulative_x[[s]], drop = FALSE],
          family = if (outcome.type == "binary") binomial() else gaussian(),
          SL.library = learner
        )
        valid[[paste0("ypredx", s)]] <- outcome_SL$SL.predict

        # Group membership model: use all train data
        ps_SL <- SuperLearner::SuperLearner(
          Y = train[[r]],
          X = train[, cumulative_x[[s]], drop = FALSE],
          newX = valid[, cumulative_x[[s]], drop = FALSE],
          family = binomial(),
          SL.library = learner
        )
        valid[[paste0("psRx", s)]] <- ps_SL$SL.predict
      }
      pred_list[[length(pred_list) + 1]] <- valid
    }

    main_df <- dplyr::bind_rows(pred_list)
    w_all <- get_weights(main_df)
    mu1_dr <- Hmisc::wtd.mean(main_df[[y]], w_all * main_df[[r]])
    theta_ri <- numeric(S)

    for (s in seq_along(cumulative_x)) {
      theta_ri[s] <- Hmisc::wtd.mean(main_df[[paste0("ypredx", s)]], w_all * (1 - main_df[[r]]))
    }

    # Compute DR signals and obtain adjusted estimates via survey regression
    theta_signal <- list()
    contrib_signal <- list()

    theta_estimates <- numeric(S)
    theta_se <- numeric(S)

    contrib_estimates <- numeric(S)
    contrib_se <- numeric(S)


    for (s in seq_along(cumulative_x)) {

      ps <- main_df[[paste0("psRx", s)]]
      ypred <- main_df[[paste0("ypredx", s)]]
      w_temp <- main_df[[r]] / Hmisc::wtd.mean(1 - main_df[[r]], w_all) * ((1 - ps) / ps)
      w_temp <- trimQ(w_temp, 0.01, 0.99)

      theta_signal[[s]] <- w_temp * (main_df[[y]] - ypred) +
        (1 - main_df[[r]]) / Hmisc::wtd.mean(1 - main_df[[r]], w_all) * (ypred - theta_ri[s]) +
        theta_ri[s]


      design_obs <- survey::svydesign(ids = ~1, data = main_df %>% dplyr::mutate(theta_sig = theta_signal[[s]]), weights = ~w_all)
      mod_theta <- survey::svyglm(theta_sig ~ 1, design = design_obs)

      theta_estimates[s] <- coef(mod_theta)[1]
      theta_se[s] <- sqrt(diag(vcov(mod_theta)))[1]

      ifelse(s == 1,
             contrib_signal[[s]] <- theta_signal[[s]] - mean(d[d[[r]]==1,y]),
             contrib_signal[[s]] <- theta_signal[[s]] - theta_signal[[s-1]]
      )

      design_obs <- survey::svydesign(ids = ~1, data = main_df %>% dplyr::mutate(contrib_sig = contrib_signal[[s]]), weights = ~w_all)
      mod_contrib <- survey::svyglm(contrib_sig ~ 1, design = design_obs)
      contrib_estimates[s] <- coef(mod_contrib)[1]
      contrib_se[s] <- sqrt(diag(vcov(mod_contrib)))[1]


    }

    total_expln_signal <- theta_signal[[length(theta_signal)]] -  mean(d[d[[r]]==1,y])
    design_obs <- survey::svydesign(ids = ~1, data = main_df %>% dplyr::mutate(total_expln_sig = total_expln_signal), weights = ~w_all)
    mod_total <- survey::svyglm(total_expln_sig ~ 1, design = design_obs)
    total_estimates <- coef(mod_total)[1]
    total_se <- sqrt(diag(vcov(mod_total)))[1]

    return(list(
      estimate = c(theta_estimates, contrib_estimates, total_estimates),
      se = c(theta_se, contrib_se, total_se)
    ))
  }


  # -----------------------------------------
  # Define bootstrap statistic functions for simultaneous estimators
  # -----------------------------------------
  ## For simultaneous decomposition with RI estimator (full Oaxaca–Blinder)
  boot_statistic_sim_ri <- function(data, indices, y, r, x_full, weight.var) {
    d <- as.data.frame(data)[indices, , drop = FALSE]
    w <- get_weights(d)
    d1 <- d[d[[r]] == 1, , drop = FALSE]
    d0 <- d[d[[r]] == 0, , drop = FALSE]
    if(nrow(d1) == 0 || nrow(d0) == 0)
      return(rep(NA, length(x_full) + 3))
    mu1 <- Hmisc::wtd.mean(d1[[y]], w[d[[r]] == 1])
    mu0 <- Hmisc::wtd.mean(d0[[y]], w[d[[r]] == 0])

    overalldiff <- mu0-mu1

    outcome_formula <- as.formula(paste(y, "~", paste(x_full, collapse = " + ")))
    mod <- survey::svyglm(outcome_formula, design  = survey::svydesign(id=~1, weights = get_weights(d1), data = d1))

    beta <- stats::coef(mod)[-1]
    mean_d1 <- apply(d1 %>% dplyr::select(x_full), 2, function(v) Hmisc::wtd.mean(v, get_weights(d1)))
    mean_d0 <- apply(d0 %>% dplyr::select(x_full), 2, function(v) Hmisc::wtd.mean(v, get_weights(d0)))


    indiv_contrib <- (mean_d0 - mean_d1) * beta # excluding intercept
    total_expl <- sum(indiv_contrib)

    unexplained <- overalldiff - total_expl

    return(c(indiv_contrib, Total_Explained = total_expl, Unexplained = unexplained, Overall_Difference = overalldiff))
  }


  ## For simultaneous decomposition with weighting estimator ("w")
  boot_statistic_sim_w <- function(data, indices, y, r, x_full, weight.var) {
    d <- data[indices, , drop = FALSE]
    w <- get_weights(d)
    d1 <- d[d[[r]] == 1, , drop = FALSE]
    if(nrow(d1) == 0)
      return(NA)
    mu1 <- Hmisc::wtd.mean(d1[[y]], w[d[[r]] == 1])
    wt_formula <- as.formula(paste(r, "~", paste(x_full, collapse = " + ")))
    mod <- stats::glm(wt_formula, data = d, weights = w, family = binomial())
    ps <- stats::predict(mod, newdata = d, type = "response")
    ps[ps < 0.01] <- 0.01
    ps[ps > 0.99] <- 0.99
    p_white <- Hmisc::wtd.mean(1 - d[[r]], w)
    w_adj <- ifelse(d[[r]] == 1,
                    trimQ(1 / p_white * ((1 - ps) / ps),0.01,0.99),
                    0)
    mu1_pred <- Hmisc::wtd.mean(d[[y]][d[[r]] == 1], w_adj[d[[r]] == 1])
    total_expl <- mu1_pred - mu1
    return(total_expl)
  }

  ## For simultaneous decomposition with DR estimator ("dr")
  boot_statistic_sim_dr <- function(data, indices, y, r, x_full, weight.var) {
    d <- data[indices, , drop = FALSE]
    w <- get_weights(d)
    d1 <- d[d[[r]] == 1, , drop = FALSE]
    d0 <- d[d[[r]] == 0, , drop = FALSE]
    if(nrow(d1) == 0 || nrow(d0) == 0)
      return(NA)
    mu1 <- Hmisc::wtd.mean(d1[[y]], w[d[[r]] == 1])

    outcome_formula <- as.formula(paste0(y,"~",r,"+",paste(x_full, collapse = " + ")))
    mod_ri <- stats::lm(outcome_formula, data = d, weights = get_weights(d))
    dmutate <- d
    dmutate[[r]] <- 1
    preds <- stats::predict(mod_ri, newdata = dmutate)
    mu1_pred_RI <- Hmisc::wtd.mean(preds[d[[r]] == 0], w[d[[r]] == 0])

    wt_formula <- as.formula(paste(r, "~", paste(x_full, collapse = " + ")))
    mod_wt <- stats::glm(wt_formula, data = d, weights = w, family = binomial())
    ps <- stats::predict(mod_wt, newdata = d, type = "response")
    ps[ps < 0.01] <- 0.01
    ps[ps > 0.99] <- 0.99
    p_white <- Hmisc::wtd.mean(1 - d[[r]], w)
    w_adj <- ifelse(d[[r]] == 1,
                    trimQ(1 / p_white * ((1 - ps) / ps),0.01,0.99),
                    0)
    score <- w_adj * (d[[y]] - preds) +
      ifelse(d[[r]] == 0, 1 / p_white * (preds - mu1_pred_RI), 0) +
      mu1_pred_RI
    mu_B_pred_DR <- Hmisc::wtd.mean(score, w)
    total_expl <- mu_B_pred_DR - mu1
    return(total_expl)
  }

  # -----------------------------------------
  # Main processing by type
  # -----------------------------------------

  if (type == "sequential") {
    if (estimator != "dml") {
        boot_obj <- boot::boot(
          data = data,
          statistic = function(data, indices) {
            boot_statistic_seq(data, indices, y, r, cumulative_x, weight.var, estimator, outcome.type)
          },
          R = B
        )
        est_vals <- boot_obj$t0  # vector of length S+1
        se_vals <- apply(boot_obj$t, 2, sd, na.rm = TRUE)

        result_df <- tibble::tibble(
          estimator  = toupper(estimator),
          type       = "sequential",
          label      = c(paste0("theta_Stage",1:length(cumulative_x)),
                         paste0("explained_Stage",1:length(cumulative_x)),
                         "total_explained"),
          estimate   = est_vals,
          se         = se_vals
        )
    } else {
      ## Sequential DML branch
        res <- sequential_dml(data, y, r, cumulative_x, weight.var, learner, outcome.type, K)
        stage_names <- if (!is.null(names(x))) names(x) else paste0("Stage_", seq_len(length(x)))
        result_df <- tibble::tibble(
          estimator  = toupper(estimator),
          type       = "sequential",
          label      = c(paste0("theta_Stage",1:length(cumulative_x)),
                         paste0("explained_Stage",1:length(cumulative_x)),
                         "total_explained"),
          estimate   = res[[1]],
          se         = res[[2]]
        )
      }
    }
   else if (type == "simultaneous") {
    if (estimator == "ri") {
      # Full Oaxaca–Blinder decomposition with detailed output via bootstrap

        boot_obj <- boot::boot(
          data = data,
          statistic = function(data, indices) {
            boot_statistic_sim_ri(data, indices, y, r, x_full, weight.var)
          },
          R = B
        )
        names_out <- c(paste0("Contrib_", setdiff(x_full, "(Intercept)")),
                       "Total_Explained", "Unexplained", "Overall_Difference")
        est_vals <- boot_obj$t0
        se_vals <- apply(boot_obj$t, 2, sd, na.rm = TRUE)
        result_df <- tibble::tibble(
          estimator  = toupper(estimator),
          type       = "simultaneous",
          label  = names_out,
          estimate   = est_vals,
          se         = se_vals
        )
    } else if (estimator == "w") {
      # Simultaneous decomposition for weighting estimator ("w") via bootstrap

        boot_obj <- boot::boot(
          data = data,
          statistic = function(data, indices) {
            boot_statistic_sim_w(data, indices, y, r, x_full, weight.var)
          },
          R = B
        )
        est_val <- boot_obj$t0
        se_val <- sd(boot_obj$t, na.rm = TRUE)
        result_df <- tibble::tibble(
          estimator  = toupper(estimator),
          type       = "simultaneous",
          label  = "Total_Explained",
          estimate   = est_val,
          se         = se_val
        )

    } else if (estimator == "dr") {
      # Simultaneous decomposition for DR estimator ("dr") via bootstrap

        boot_obj <- boot::boot(
          data = data,
          statistic = function(data, indices) {
            boot_statistic_sim_dr(data, indices, y, r, x_full, weight.var)
          },
          R = B
        )
        est_val <- boot_obj$t0
        se_val <- sd(boot_obj$t, na.rm = TRUE)
        result_df <- tibble::tibble(
          estimator  = toupper(estimator),
          type       = "simultaneous",
          label  = "Total_Explained",
          estimate   = est_val,
          se         = se_val
        )

    } else if (estimator == "dml") {
      # Simultaneous decomposition with DML estimator using cross-fitting and analytical SE.
        d <- data
        w <- get_weights(d)
        mu1 <- Hmisc::wtd.mean(d[[y]], w * d[[r]])
        folds <- caret::createFolds(d[[y]], k = K)
        phi_all <- c()
        n_total <- nrow(d)
        for (fold in folds) {
          train <- d[-fold, , drop = FALSE]
          valid <- d[fold, , drop = FALSE]
          valid_mutate <- valid
          valid_mutate[[r]] <- 1

          outcome_SL <- SuperLearner::SuperLearner(
            Y = train[[y]],
            X = train[, c(r,x_full), drop = FALSE],
            newX = valid_mutate[, c(r,x_full), drop = FALSE],
            family = if (tolower(outcome.type) == "binary") binomial() else gaussian(),
            SL.library = learner
          )

          m_hat <- outcome_SL$SL.predict

          if (sum(valid[[r]] == 0) > 0) {
            mu_pred_RI <- Hmisc::wtd.mean(m_hat[valid[[r]] == 0], get_weights(valid)[valid[[r]] == 0])
          } else {
            mu_pred_RI <- NA
          }

          ps_SL <- SuperLearner::SuperLearner(
            Y = train[[r]],
            X = train[, x_full, drop = FALSE],
            newX = valid[, x_full, drop = FALSE],
            family = binomial(),
            SL.library = learner
          )
          e_hat <- ps_SL$SL.predict
          e_hat[e_hat < 0.01] <- 0.01
          e_hat[e_hat > 0.99] <- 0.99

          p_white <- Hmisc::wtd.mean(1 - train[[r]], get_weights(train))

          W_i <- ifelse(valid[[r]] == 1, 1 / p_white * ((1 - e_hat) / e_hat), 0)
          phi <- W_i * (valid[[y]] - m_hat) +
            ifelse(valid[[r]] == 0,
                   trimQ(1 / p_white * (m_hat - mu_pred_RI),0.01,0.99),
                   0) +
            mu_pred_RI
          phi_all <- c(phi_all, phi)
        }
        mu_B_pred_DR <- mean(phi_all, na.rm = TRUE)
        total_expl <- mu_B_pred_DR - mu1
        var_phi <- var(phi_all, na.rm = TRUE)
        se <- sqrt(var_phi / n_total)
        result_df <- tibble::tibble(
          estimator  = toupper(estimator),
          type       = "simultaneous",
          label  = "Total_Explained",
          estimate   = total_expl,
          se         = se
        )

    } else {
      stop("For simultaneous decomposition, estimator must be 'ri', 'w', 'dr', or 'dml'.")
    }
  } else {
    stop("Invalid 'type'. Must be either 'sequential' or 'simultaneous'.")
  }
  return(result_df)
}
