lag <- dplyr::lag
select <- dplyr::select
filter <- dplyr::filter
wtd.mean <- Hmisc::wtd.mean
wtd.quantile <- Hmisc::wtd.quantile
recode <- dplyr::recode


wtd_rank <- function (x, weights = NULL, normwt = FALSE, na.rm = TRUE) {
  if (!length(weights))
    return(rank(x, na.last = if (na.rm) NA else TRUE))
  tab <- Hmisc::wtd.table(x, weights, normwt = normwt, na.rm = na.rm)
  freqs <- tab$sum.of.weights
  r <- cumsum(freqs) - 0.5 * (freqs - 1)
  approx(tab$x, r, xout = x, rule = 2)$y
}

wtd_prank <- function(x, w){
  100 * wtd_rank(x, w, normwt = TRUE)/length(x[!is.na(x)])
}

Mode <- function(x, na.rm = TRUE) {
  if(na.rm) x <- x[complete.cases(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

cumMode <- function(x) {
  map_dbl(1:length(x), ~ Mode(x[1:.x]))
}

Last <- function(x, na.rm = TRUE) {
  if(na.rm) x <- x[complete.cases(x)]
  x[length(x)]
}

# cumLast <- function(x) {
#   map_dbl(1:length(x), ~ Last(x[1:.x]))
# }

`%notin%` <- Negate(`%in%`)

`%||%` <- function(a, b) if (!is.null(a)) a else b

cummax2 <- function(x) dplyr::na_if(cummax(tidyr::replace_na(x, -Inf)), -Inf)
cummin2 <- function(x) dplyr::na_if(cummin(tidyr::replace_na(x, Inf)), Inf)


estse2star <- function(est, se){
  x <-  2 * pnorm(-abs(est/se))
  case_when(between(x, 0.05, 1) ~ "",
            between(x, 0.01, 0.05) ~ "*",
            between(x, 0.001, 0.01) ~ "**",
            between(x, 0, 0.001) ~ "***")
}

p2star <- function(x){
  case_when(between(x, 0.05, 1) ~ "",
            between(x, 0.01, 0.05) ~ "*",
            between(x, 0.001, 0.01) ~ "**",
            between(x, 0, 0.001) ~ "***")
}

combine <- function(estimate, std.error, star, format, digits){
  paste0(formatC(estimate, format = format, digits = digits), star) %>%
    paste0(., " (", formatC(std.error, format = format, digits = digits), ")")
}

regtable <- function(..., format = "f", digits = 3){
  cl <- match.call()
  cl$format <- cl$digits <-  cl[[1]] <- NULL
  eval_tidy(expr(map(list2(!!!cl), broom::tidy))) %>%
    map(~ mutate(.x,
                 star = p2star(p.value),
                 out = combine(estimate, std.error, star,
                               format = format, digits = digits))) %>%
    map(~ dplyr::select(.x, term, out)) %>%
    reduce(full_join, by = "term") %>%
    mutate_all( ~ tidyr::replace_na(.x, "")) %>%
    set_names(c("term", map_chr(ensyms(...), as_string)))
}

median2 <- function(y, x) y[[order(x)[[ceiling(length(x)/2)]]]]

rowSums2 <- function(x) ifelse(rowSums(is.na(x), na.rm = TRUE) == ncol(x), NA, rowSums(x, na.rm = TRUE))

near2 <- function(x, y) map_lgl(x, ~ any(map_lgl(y, function(j) near(.x, j))))

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

scale2 <- function(x, w) (x - Hmisc::wtd.mean(x, w))/sqrt(Hmisc::wtd.var(x, w))

predsIF <- function(mod, input){

  x <- model.matrix(mod)
  w <- mod$weights/mean(mod$weights)
  n <- nrow(x)

  term1 <- solve((t(x) %*% diag(w) %*% x)/n) # p * p matrix
  term2 <- t(diag(w * mod$residuals) %*% x) # p * n matrix

  IF <- term1 %*% term2

  input %*% IF

}
