# Generate the dataframe
dummy_data <- data.frame(
  comp   = sample(c(0, 1), 1000, replace = TRUE),
  sel = sample(c(0, 1), 1000, replace = TRUE),
  black  = sample(c(0, 1), 1000, replace = TRUE),
  age    = sample(18:70, 1000, replace = TRUE),
  female = sample(c(0, 1), 1000, replace = TRUE),
  parinc = sample(20:150, 1000, replace = TRUE),
  weight = rep(1, 1000),
  stringsAsFactors = FALSE
)


test_that("prescriptive errors with missing values", {
  dummy_na <- dummy_data
  dummy_na$age[2] <- NA
  expect_error(
    prescriptive(dummy_na, x = c("age", "female", "parinc"), a = "sel", z = "age", r = "black", g = "parinc", y = "comp",
                 weight.var = "weight"),
    "Missing values"
  )
})


test_that("prescriptive ri returns expected output", {
  out_sim <-  prescriptive(dummy_data, x = c("age", "female", "parinc"), a = "sel", z = "age", r = "black", g = "parinc", y = "comp",
                           weight.var = "weight", estimator = "ri", B=20)
  expect_true(is.list(out_sim))
  expect_length(out_sim, 3)
  lapply(out_sim, function(df) expect_true(is.data.frame(df)))
  expect_true(all(c("intervention", "est", "se") %in% names(out_sim[[1]])))
  expect_true(all(c("intervention", "est", "se") %in% names(out_sim[[2]])))
  expect_true(all(c("intervention", "r", "est", "se") %in% names(out_sim[[3]])))
})

test_that("prescriptive dml returns expected output", {
  out_sim <-  prescriptive(dummy_data, x = c("age", "female", "parinc"), a = "sel", z = "age", r = "black", g = "parinc", y = "comp",
                           weight.var = "weight", estimator = "dml", K=2)
  expect_true(is.list(out_sim))
  expect_length(out_sim, 3)
  lapply(out_sim, function(df) expect_true(is.data.frame(df)))
  expect_true(all(c("intervention", "est", "se") %in% names(out_sim[[1]])))
  expect_true(all(c("intervention", "est", "se") %in% names(out_sim[[2]])))
  expect_true(all(c("intervention", "r", "est", "se") %in% names(out_sim[[3]])))
})
