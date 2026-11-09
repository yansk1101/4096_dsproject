# Load requried libraries
library(dlnm)
library(AER)
library(tsModel)
library(mgcv)

# File import
data <- read.csv("data2.csv", header=TRUE, stringsAsFactors=FALSE)
data$From <- as.Date(data$From)

# Determination of the value of maximum lag 
max_lag <- 4

# Cross-basis matrix construction for the exposure variables
cb_exposure <- crossbasis(
    data$'exposure_v1',
    lag=maxlag, 
    argvar=list(),
    arglag=list()
)

"cb_df <- as.data.frame(as.matrix(cb_matrix))
colnames(cb_df) <- paste0("cb_col", 1:ncol(cb_df))

model_data <- cb_data <- cbind(data, cb_data)
cb_names <- paste(colnames(cb_df), collapse=" + ")

cb_instrumental <- crossbasis(
    data$'instrumental_v1',
    lag=max_lag,
    argvar=list(),
    arglag=list()
)

cb_instrumental_df <- as.data.frame(as.matrix(cb_instrumental))
colnames(cb_instrumental_df) <- paste0("inst_col", 1:ncol(cb_instrumental_df))

model_data <- cbind(model_data, cb_instrumental_df)
instrumental_names <- paste(colnames(cb_instrumental_df), collapse=" + ")

confounders <- ""

# Final IV Formula
iv_formula <- as.formula(paste(
    "outcome ~", cb_names, "+", confounders, "|",
    instrumental_names, "+", confounders
))

iv_model <- ivreg(iv_formula, data=model_data)
summary(iv_model)"


# ==============================================================================
# REPLACED SECTION 4: TWO-STAGE INSTRUMENTAL VARIABLE REGRESSION (GAM STYLE)
# ==============================================================================
library(mgcv)

# STAGE 1: Regress the exposure cross-basis columns on the instrument basis using GAM penalization
# We predict "fitted" exposure columns stripped of endogeneity
fitted_cb_cols <- list()

for(col in colnames(cb_df)) {
  # Regress each cross-basis column against instrument bases and smooth confounders
  f1 <- as.formula(paste(col, "~", inst_names, "+ s(as.numeric(date)) + s(confounder1)"))
  stage1_mod <- gam(f1, data = model_data)
  fitted_cb_cols[[col]] <- predict(stage1_mod)
}

# Combine fitted stages into a new data frame
gam_iv_data <- model_data
for(col in colnames(cb_df)) {
  gam_iv_data[[paste0("fitted_", col)]] <- fitted_cb_cols[[col]]
}

# STAGE 2: Run the final health outcome model using the instrumented/fitted cross-basis columns
fitted_cb_names <- paste0("fitted_", colnames(cb_df), collapse = " + ")

gam_formula <- as.formula(paste(
  "outcome ~", fitted_cb_names, "+ s(as.numeric(date)) + s(confounder1)"
))

# Fit final GAM-IV model (Assuming count data, e.g., Poisson)
iv_model <- gam(gam_formula, family = poisson(), data = gam_iv_data)
summary(iv_model)

# CRITICAL JUMP FOR SECTION 5 PREDICTION:
# If using this GAM method, change Section 5's extraction naming convention to match:
coef_cb <- coef(iv_model)[paste0("fitted_cb_col", 1:ncol(cb_df))]
vcov_cb <- vcov(iv_model)[paste0("fitted_cb_col", 1:ncol(cb_df)), paste0("fitted_cb_col", 1:ncol(cb_df))]

# ==============================================================================

coef_cb <- coef(iv_model)[paste0("cb_col", 1:ncol(cb_df))]
vcov_cb <- vcov(iv_model)[paste0("cb_col", 1:ncol(cb_df)), paste0("cb_col", 1:ncol(cb_df))]
pred_iv <- crosspred(
    cb_exposure,
    coef=coef_cb,
    vcov=vcov_cb,
    model=NULL,
    at=seq(),
    cen=median(data$'exposure_v1')
)
