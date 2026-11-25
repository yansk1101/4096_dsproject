# ==============================================================================
# 1. LOAD REQUIRED LIBRARIES
# ==============================================================================
setwd("C:/Users/yanss/OneDrive/Codes/4096_project_new")
install.packages("dlnm")
install.packages("mgcv")
library(dlnm)     # For cross-basis matrices
library(mgcv)     # For Generalized Additive Models (GAM)

# ==============================================================================
# 2. FILE IMPORT
# ==============================================================================
# Assumes columns: date, outcome, exposure1, exposure2, instrument1, instrument2, confounder1
data <- read.csv("data.csv", header = TRUE, stringsAsFactors = FALSE)
data$From <- as.Date(data$From)

# ==============================================================================
# 3. DEFINE DLNM CROSS-BASIS MATRICES FOR EACH EXPOSURE & INSTRUMENT
# ==============================================================================
max_lag <- 30 

# --- EXPOSURE 1 ---
cb_exp1 <- crossbasis(data$Diurnal.temperature.range...C., lag = max_lag, 
                      argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_exp1 <- as.data.frame(as.matrix(cb_exp1))
colnames(df_exp1) <- paste0("exp1_col", 1:ncol(df_exp1))

# --- EXPOSURE 2 ---
cb_exp2 <- crossbasis(data$Mean.temperature...C., lag = max_lag, 
                      argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_exp2 <- as.data.frame(as.matrix(cb_exp2))
colnames(df_exp2) <- paste0("exp2_col", 1:ncol(df_exp2))

# --- EXPOSURE 3 ---
cb_exp3 <- crossbasis(data$Mean.humidity...., lag = max_lag, 
                      argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_exp3 <- as.data.frame(as.matrix(cb_exp3))
colnames(df_exp3) <- paste0("exp1_col", 1:ncol(df_exp3))

# --- EXPOSURE 4 ---
cb_exp4 <- crossbasis(data$Solar.radiation..MJ.m., lag = max_lag, 
                      argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_exp4 <- as.data.frame(as.matrix(cb_exp4))
colnames(df_exp4) <- paste0("exp1_col", 1:ncol(df_exp4))

# --- EXPOSURE5 ---
cb_exp5 <- crossbasis(data$Air.pressure..hPa., lag = max_lag, 
                      argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_exp5 <- as.data.frame(as.matrix(cb_exp5))
colnames(df_exp5) <- paste0("exp1_col", 1:ncol(df_exp5))

# --- EXPOSURE6 ---
cb_exp6 <- crossbasis(data$Dew.point.temperature...C., lag = max_lag, 
                      argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_exp6 <- as.data.frame(as.matrix(cb_exp6))
colnames(df_exp6) <- paste0("exp1_col", 1:ncol(df_exp6))

# --- EXPOSURE7 ---
cb_exp7 <- crossbasis(data$Wind.speed..km.h., lag = max_lag, 
                      argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_exp7 <- as.data.frame(as.matrix(cb_exp7))
colnames(df_exp7) <- paste0("exp1_col", 1:ncol(df_exp7))

# --- EXPOSURE8 ---
cb_exp8 <- crossbasis(data$Mean.UV, lag = max_lag, 
                      argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_exp8 <- as.data.frame(as.matrix(cb_exp8))
colnames(df_exp8) <- paste0("exp1_col", 1:ncol(df_exp8))

# --- INSTRUMENT 1 ---
cb_inst1 <- crossbasis(data$instrument1, lag = max_lag, 
                       argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_inst1 <- as.data.frame(as.matrix(cb_inst1))
colnames(df_inst1) <- paste0("inst1_col", 1:ncol(df_inst1))

# --- INSTRUMENT 2 ---
cb_inst2 <- crossbasis(data$instrument2, lag = max_lag, 
                       argvar = list(fun = "ns", df = 4), arglag = list(fun = "ns", df = 3))
df_inst2 <- as.data.frame(as.matrix(cb_inst2))
colnames(df_inst2) <- paste0("inst2_col", 1:ncol(df_inst2))

# Combine all generated basis columns back into a unified modeling dataset
model_data <- cbind(data, df_exp1, df_exp2, df_exp3, df_exp4, df_exp5, df_exp6, df_exp7, df_exp8)

# Collate formula text blocks
exp1_names  <- paste(colnames(df_exp1), collapse = " + ")
exp2_names  <- paste(colnames(df_exp2), collapse = " + ")
all_insts   <- paste(c(colnames(df_inst1), colnames(df_inst2)), collapse = " + ")

# ==============================================================================
# 4. TWO-STAGE GAM CONTROL FUNCTION WITH RESIDUAL EXTRACTION
# ==============================================================================

# --- STAGE 1: Extract Endogenous Residuals for EACH Exposure Column ---
# We loop through every column of BOTH exposure bases and regress them on instruments
res_list <- list()

# Extract residuals for Exposure 1 Basis Columns
for(col in colnames(df_exp1)) {
  f1 <- as.formula(paste(col, "~", "+ s(as.numeric(From)) + s(confounder1)"))
  m1 <- gam(f1, data = model_data)
  res_list[[paste0("res_", col)]] <- residuals(m1, type = "response")
}

# Extract residuals for Exposure 2 Basis Columns
for(col in colnames(df_exp2)) {
  f2 <- as.formula(paste(col, "~", "+ s(as.numeric(From)) + s(confounder1)"))
  m2 <- gam(f2, data = model_data)
  res_list[[paste0("res_", col)]] <- residuals(m2, type = "response")
}

# Bind all generated residual columns into the dataset
df_res <- as.data.frame(res_list)
model_data <- cbind(model_data, df_res)
res_names  <- paste(colnames(df_res), collapse = " + ")

# --- STAGE 2: Fit Final Health Outcome GAM Model ---
# The model includes: Both Exposure Bases + All Residual Terms + Smooth Confounders
gam_formula <- as.formula(paste(
  "outcome ~", exp1_names, "+", exp2_names, "+", res_names, 
  "+ s(as.numeric(From)) + s(confounder1)"
))

# Fit the model (Using quasi-poisson for overdispersed count data)
gam_iv_model <- gam(gam_formula, family = quasipoisson(), data = model_data)
summary(gam_iv_model)

# ==============================================================================
# 5. PREDICTION AND VISUALIZATION (Example: Plotting Exposure 1)
# ==============================================================================

# Isolate coefficients and vcov matrices mapping explicitly to Exposure 1
coef_exp1 <- coef(gam_iv_model)[colnames(df_exp1)]
vcov_exp1 <- vcov(gam_iv_model)[colnames(df_exp1), colnames(df_exp1)]

# Generate the DLNM causal prediction for Exposure 1
pred_exp1 <- crosspred(
  cb_exp1, 
  coef = coef_exp1, 
  vcov = vcov_exp1, 
  model = NULL, 
  at = seq(min(data$exposure1), max(data$exposure1), length.out = 50),
  cen = median(data$exposure1)
)

# Plot 3D Causal Surface for Exposure 1
plot(pred_exp1, xlab = "Exposure 1", ylab = "Lag (Days)", zlab = "Effect Size", 
     main = "Exposure 1: 3D Causal Lag Surface")

# Plot Overall Cumulative Effect for Exposure 1
plot(pred_exp1, "overall", xlab = "Exposure 1", ylab = "Cumulative Effect",
     main = "Exposure 1: Cumulative Effect Curve")
