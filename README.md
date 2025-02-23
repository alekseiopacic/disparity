
# disparity

<!-- badges: start -->
<!-- badges: end -->

This R package provides tools to implement the methods outlined in Opacic Wei and Zhou (2025): "Disparity Analysis: A Tale of Two Approaches". The package currently contains two functions: **descriptive** and **prescriptive**. More features are currently under development.

## Installation

You can install the development version of disparity from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("laiweisociology/disparity")
```

## Example: descriptive disparity analysis

``` r
library(disparity)

  # Sequential decomposition: each covariate is added in order (non-DML)
  res_seq <- descriptive(data = mydata, y = "comp", r = "black", x = c("age", "female", "parinc"),
                         weight.var = "weight", estimator = "ri", type = "sequential")

  # Sequential decomposition using DML:
  res_seq_dml <- descriptive(data = mydata, y = "comp", r = "black", x = c("age", "female", "parinc"),
                             weight.var = "weight", estimator = "dml", type = "sequential")

  # Simultaneous decomposition (Oaxaca-Blinder) with RI estimator:
  res_sim <- descriptive(data = mydata, y = "comp", r = "black", x = c("age", "female", "parinc"),
                         weight.var = "weight", estimator = "ri", type = "simultaneous")
```


## Example: prescriptive disparity analysis

``` r
library(disparity)

  # Prescriptive analysis using parametric (ri) estimator
  res_ri <- prescriptive(data_list = data, x = c("age", "female", "parinc"), a = "selective", z = "gpa", r = "black", g = "parinc", y = "comp",
                         estimator = "ri", B=250)

  # Prescriptive analysis using DML estimator
  res_dml <- prescriptive(data_list = data, x = c("age", "female", "parinc"), a = "selective", z = "gpa", r = "black", g = "parinc", y = "comp",
                          estimator = "dml", K = 5)
```






