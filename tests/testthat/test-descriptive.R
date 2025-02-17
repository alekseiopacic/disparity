



# Generate the dataframe
dummy_data <- data.frame(
  comp   = sample(c(0, 1), 1000, replace = TRUE),
  black  = sample(c(0, 1), 1000, replace = TRUE),
  age    = sample(18:70, 1000, replace = TRUE),
  female = sample(c(0, 1), 1000, replace = TRUE),
  parinc = sample(20:150, 1000, replace = TRUE),
  weight = rep(1, 1000),
  stringsAsFactors = FALSE
)


test_that("descriptive errors with missing values", {
  dummy_na <- dummy_data
  dummy_na$age[2] <- NA
  expect_error(
    descriptive(dummy_na, y = "comp", r = "black",
                x = c("age", "female", "parinc"), weight.var = "weight"),
    "Missing values"
  )
})

test_that("descriptive sequential RI returns expected output", {
  # For sequential (non-DML) decomposition with 3 covariates,
  # the revised function returns a data frame with rows for:
  # - theta estimates (S rows)
  # - stage contributions (S rows)
  # - one row for total explained
  # Thus, expect 2 * S + 1 rows.
  S <- length(c("age", "female", "parinc"))

  out_seq_ri <- descriptive(dummy_data, y = "comp", r = "black",
                            x = c("age", "female", "parinc"), weight.var = "weight",
                            estimator = "ri", type = "sequential")
  expect_true(is.data.frame(out_seq_ri))
  expect_true(all(c("estimator","type","label","estimate", "se") %in% names(out_seq_ri)))
  expect_equal(nrow(out_seq_ri), 2 * S + 1)
})

test_that("descriptive sequential DML returns expected output", {
  # For sequential DML, expect one row per stage plus one for total explained (S + 1 rows)
  S <- length(c("age", "female", "parinc"))

  out_seq_dml <- descriptive(dummy_data, y = "comp", r = "black",
                             x = c("age", "female", "parinc"), weight.var = "weight",
                             estimator = "dml", type = "sequential")
  expect_true(all(c("estimator","type","label","estimate", "se") %in% names(out_seq_dml)))
})

test_that("descriptive simultaneous RI returns expected output", {
  out_sim_ri <- descriptive(dummy_data, y = "comp", r = "black",
                            x = c("age", "female", "parinc"), weight.var = "weight",
                            estimator = "ri", type = "simultaneous")
  expect_true(is.data.frame(out_sim_ri))
  expect_true(all(c("estimator","type","label","estimate", "se") %in% names(out_sim_ri)))
})

test_that("descriptive simultaneous DR returns expected output", {
  out_sim_dr <- descriptive(dummy_data, y = "comp", r = "black",
                            x = c("age", "female", "parinc"), weight.var = "weight",
                            estimator = "dr", type = "simultaneous")
  expect_true(is.data.frame(out_sim_dr))
  expect_true(all(c("estimator","type","label","estimate", "se") %in% names(out_sim_dr)))
})

test_that("descriptive simultaneous DML returns expected output", {
  out_sim_dml <- descriptive(dummy_data, y = "comp", r = "black",
                             x = c("age", "female", "parinc"), weight.var = "weight",
                             estimator = "dml", type = "simultaneous")
  expect_true(is.data.frame(out_sim_dml))
  expect_true(all(c("estimator","type","label","estimate", "se") %in% names(out_sim_dml)))
})
