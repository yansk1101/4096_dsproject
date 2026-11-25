library(dlnm)
library(splines)
library(pbs)
library(glmmTMB)
library(DHARMa)

df <- read.csv("data.csv", header = TRUE, stringsAsFactors = FALSE)
df$From <- as.Date(df$From)
colnames(df) <- c("Year", "Week", "From", "temprange", "temp", "RH", "SR", "AP", "DPT", "WS", "UV", "H1", "H3", "B", "AandB", "ContainmentHealthIndex_Average", "pandemic_flag", "Population", "Log_population", "school_flag")
df$time <- seq_len(nrow(df))
df$time <- 1:nrow(df)


cb_temprange <- crossbasis(
    df$temprange, 
    lag = 4, 
    argvar = list(fun = "ns", df = 4), 
    arglag = list(fun = "ns", df = 3)
)


model_zinb <- glmmTMB(
    H1 ~ cb_temprange + bs(time, df = 25) + bs(ContainmentHealthIndex_Average, df = 4) + pbs(Week, df=6),
    ziformula = ~ ContainmentHealthIndex_Average, 
    family = nbinom2,
    data = df
)

# 1. Check summary, convergence, and zero-inflation intercept
print(glmmTMB:::summary.glmmTMB(model_zinb))


# 2. Check if zero-inflation is statistically needed (compare against standard NB)
model_nb <- update(model_zinb, ziformula = ~0)
print(anova(model_nb, model_zinb)) # A significant p-value (< 0.05) confirms ZINB improves fit

# 3. Check residual diagnostics (zero-inflation & dispersion) using DHARMa

sim_res <- simulateResiduals(model_zinb)
plot(sim_res)
testZeroInflation(sim_res)
testDispersion(sim_res)

# 4. Check for remaining autocorrelation in residuals (crucial for time-series)
acf(residuals(model_zinb, type = "pearson"), na.action = na.pass)

glmmTMB:::summary.glmmTMB(model_zinb)

# 4. Extract count-component parameters and run crosspred()
coef_zinb <- glmmTMB:::fixef.glmmTMB(model_zinb)$cond
vcov_zinb <- glmmTMB:::vcov.glmmTMB(model_zinb)$cond
idx <- grep("cb_temprange", names(coef_zinb))

cen_val <- median(df$temprange, na.rm = TRUE)

pred <- crosspred(
    basis = cb_temprange,
    coef = coef_zinb[idx],
    vcov = vcov_zinb[idx, idx],
    model.link = "log",
    cen = cen_val
)

# 5. Plot the clean overall cumulative IRR curve
plot(pred, "overall", 
     xlab = "Temperature Range", 
     ylab = "Cumulative IRR",
     main = "Overall Cumulative Effect (ZINB-DLNM)",
     ylim = c(0, 4))



# 1. Define exact slice values
p10 <- round(quantile(df$temprange, 0.10, na.rm = TRUE), 1)
p90 <- round(quantile(df$temprange, 0.90, na.rm = TRUE), 1)

# 2. Predict specifically at p10 and p90
pred_slices <- crosspred(
    basis = cb_temprange,
    coef = coef_zinb[idx],
    vcov = vcov_zinb[idx, idx],
    model.link = "log",
    cen = cen_val,
    at = c(p10, p90)
)

# 3. Verify what values are stored in pred_slices$predvar
print(pred_slices$predvar)

# 4. Plot using pred_slices (notice the first argument is pred_slices, not pred)
par(mfrow = c(1, 2))
plot(pred_slices, "slices", var = pred_slices$predvar[1],
     xlab = "Lag (Weeks)", ylab = "IRR",
     main = paste0("Low Temp Range (", pred_slices$predvar[1], ")"))

plot(pred_slices, "slices", var = pred_slices$predvar[2],
     xlab = "Lag (Weeks)", ylab = "IRR",
     main = paste0("High Temp Range (", pred_slices$predvar[2], ")"))
par(mfrow = c(1, 1))
