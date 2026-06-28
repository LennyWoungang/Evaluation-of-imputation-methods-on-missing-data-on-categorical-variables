################################ Missing Data mechanism
# row-wise
create_NA <- function(data, p_row) {
  n <- nrow(data)
  p <- round(p_row * n)         
  cols <- colnames(data)
  rows_to_miss <- sample.int(n, p, replace = FALSE)
  out <- data
  for (r in rows_to_miss) {
    c <- sample(cols, 1)
    out[r, c] <- NA
  }
  out
}


### For different percentages of missingness
set.seed(6516)
M <- 500
NA5_list <- vector("list", M)
for(i in 1:M){
  NA5 <- create_NA(data = data_for_missing, p = 0.05)
  NA5_list[[i]] <- NA5
}


NA15_list <- vector("list", M)
for(i in 1:M){
  NA15 <- create_NA(data = data_for_missing, p = 0.15)
  NA15_list[[i]] <- NA15
}

NA30_list <- vector("list", M)
for(i in 1:M){
  NA30 <- create_NA(data = data_for_missing, p = 0.30)
  NA30_list[[i]] <- NA30
}

# check for missing data
count_missing_rows <- function(df) {
  sum(!complete.cases(df))
}
summary_df <- data.frame(
  five  = as.list(summary(sapply(NA5_list, count_missing_rows))),
  fifteen = as.list(summary(sapply(NA15_list, count_missing_rows))),
  thirty  = as.list(summary(sapply(NA30_list, count_missing_rows)))
)
count_total_NA <- function(df) sum(is.na(df))
total_NA_5  <- sapply(NA5_list, count_total_NA)
total_NA_15 <- sapply(NA15_list, count_total_NA)
total_NA_30 <- sapply(NA30_list, count_total_NA)


