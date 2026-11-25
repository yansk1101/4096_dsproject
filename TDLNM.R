

library(dlmtree)
library(pbs)
library(splines)

df <- read.csv("data.csv", header = TRUE, stringsAsFactors = FALSE)
df$From <- as.Date(df$From)
colnames(df) <- c("Year", "Week", "From", "temprange", "temp", "RH", "SR", "AP", "DPT", "WS", "UV", "H1", "H3", "B", "AandB", "ContainmentHealthIndex_Average", "pandemic_flag", "Population", "Log_population", "school_flag")

df$time <- seq_len(nrow(df))
df$time <- 1:nrow(df)

max_lag <- 4

# 1. Create the N x (max_lag + 1) lag matrix using base R embed()
# Note: embed() automatically drops the first 'max_lag' rows that would have NAs
my_lag_matrix <- embed(df$temprange, max_lag + 1)
colnames(my_lag_matrix) <- paste0("lag", 0:max_lag)

# 2. Align df by dropping the exact same first 'max_lag' rows
df_clean <- df[(max_lag + 1):nrow(df), ]


# Part 1: Probability of a non-zero count (Logistic TDLNM)
df_clean$y_binary <- as.numeric(df_clean$H1 > 0)

fit_part1 <- dlmtree(
    formula = y_binary ~ bs(time, df=25) + bs(ContainmentHealthIndex_Average, df=4) + pbs(Week, df=6),
    data = df_clean,
    exposure.data = my_lag_matrix,
    dlm.type = "nonlinear",
    family = "logit",
    control.diagnose = list(diagnostics = TRUE),
    control.mcmc = list(n.trees = 20, n.burn = 2000, n.iter = 5000, n.thin = 5)
)

# Part 2: Magnitude when count > 0 (Gaussian TDLNM on log scale)
pos_idx <- which(df_clean$y > 0)

fit_part2 <- dlmtree(
    formula = log(H1) ~ temprange + bs(time, df = 25) + bs(ContainmentHealthIndex_Average, df = 4) + pbs(Week, df=6),
    data = df_clean[pos_idx, ],
    exposure.data = my_lag_matrix[pos_idx, ],
    dlm.type = "nonlinear",
    family = "gaussian",
    control.diagnose = list(diagnostics = TRUE),
    control.mcmc = list(n.trees = 20, n.burn = 2000, n.iter = 5000, n.thin = 5)
)

# 1. Generate the summary object (sets the baseline exposure reference, e.g., median)
fit_sum <- summary(fit_part1, cenval = median(my_lag_matrix), mcmc=TRUE)

# 2. Print numerical summaries (cumulative effects & fixed-effect covariates)
diagnose(fit_sum)




plot(fit_part2$sigma2, type = "l", main = "MCMC Trace Plot for Sigma2", ylab = "sigma2")

# 查看固定效應參數 (gamma) 的第一欄軌跡圖
plot(fit_part2$gamma[, 1], type = "l", main = "Trace Plot for Covariate 1", ylab = "gamma")












# 1. Properly unscale fit_part2$Yhat from standardized units back to the log(y) scale
y_hat_pos_log <- fit_part2$Yhat * fit_part2$Yscale + fit_part2$Ymean

#    Verify that the unscaled mean now matches Mean log(y[pos_idx]) (~3.21)!
cat("Unscaled mean of fit_part2$Yhat:", mean(y_hat_pos_log), "\n")
cat("Actual mean of log(y[pos_idx]):", mean(log(df_clean$H1[pos_idx])), "\n")

# 2. Align the intercept of raw_log_pred to the UNSCALED y_hat_pos_log
intercept_shift <- mean(y_hat_pos_log) - mean(raw_log_pred[pos_idx])
log_count_pred  <- raw_log_pred + intercept_shift

# 3. Bound predictions to the observed range of y_hat_pos_log (not raw extreme outliers)
#    to prevent spline extrapolation spikes in zero-gaps
pred_bounds         <- quantile(y_hat_pos_log, probs = c(0.01, 0.99))
log_count_pred_safe <- pmin(pmax(log_count_pred, pred_bounds[1]), pred_bounds[2])
log_count_pred_safe[pos_idx] <- y_hat_pos_log

# 4. Back-transform with Duan's smearing factor and multiply by Part 1 P(y > 0)
smearing_factor          <- mean(exp(log(df_clean$H1[pos_idx]) - y_hat_pos_log))
expected_count_given_pos <- exp(log_count_pred_safe) * smearing_factor
combined_y_hat           <- prob_nonzero * expected_count_given_pos

# 5. Print true metrics and plot
rmse_combined <- sqrt(mean((df_clean$H1 - combined_y_hat)^2))
mae_combined  <- mean(abs(df_clean$H1 - combined_y_hat))
cor_combined  <- cor(df_clean$H1, combined_y_hat)

cat(sprintf("True Combined Hurdle -> RMSE: %.3f | MAE: %.3f | Correlation (r): %.3f\n",
            rmse_combined, mae_combined, cor_combined))

plot(df_clean$y, type = "l", col = "gray60", ylab = "Count (y)", xlab = "Time Index",
     main = "Observed vs. Combined Hurdle TDLNM (Unscaled)")
lines(combined_y_hat, col = "red", lwd = 1.5)
legend("topright", legend = c("Observed", "Predicted"), col = c("gray60", "red"), lty = 1)



