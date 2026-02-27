# ===== XGBoost multiclass: stratified 5-fold CV (OOF) + Pairwise ROC w/ 95% CI (SO vs O/S/H) | xgboost 3.1.3.1 =====
options(stringsAsFactors = FALSE, scipen = 99)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(xgboost)
})

# pROC (for AUC + ROC CI)
if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("缺少包 pROC。请先手动安装：install.packages('pROC')")
}
suppressPackageStartupMessages(library(pROC))

# ---- 1) paths ----
workdir <- "D:/Users/ruzuo/Desktop/菌群学习/XGBoost"
x_file  <- file.path(workdir, "V1abund_plus_species_merged_by_SampleID.csv")
outdir  <- file.path(workdir, "XGB_totalCV5plus")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- 2) read ----
dat <- readr::read_csv(x_file, show_col_types = FALSE)
stopifnot("SampleID" %in% names(dat))
sample_id <- trimws(as.character(dat$SampleID))

# ---- 3) y from SampleID ----
y <- ifelse(grepl("^SO[-_]", sample_id), "SO", sub("[-_].*$", "", sample_id))
y <- toupper(y); y[y == "SO"] <- "SO"
class_levels <- c("SO","S","O","H")
y <- factor(y, levels = class_levels)

if (any(is.na(y))) {
  bad <- unique(sample_id[is.na(y)])
  stop("以下 SampleID 无法识别分组（需形如 SO-xx / S-xx / O-xx / H-xx）：\n", paste(bad, collapse = "\n"))
}
message("Overall class counts:"); print(table(y))

# ---- 4) X matrix ----
X <- dat %>% select(-SampleID) %>% as.data.frame()
for (j in seq_along(X)) X[[j]] <- suppressWarnings(as.numeric(X[[j]]))
X <- as.matrix(X)
if (anyNA(X)) stop("X 中存在 NA（非数值/空值）。请先处理缺失或检查列。")

# ---- 5) make stratified 5 folds ----
set.seed(20260114)
K <- 5

fold_id <- integer(nrow(X))
idx_by_class <- split(seq_len(nrow(X)), y)

for (cls in names(idx_by_class)) {
  idx <- sample(idx_by_class[[cls]])
  fold_id[idx] <- rep(1:K, length.out = length(idx))
}

message("\nFold class counts:")
for (k in 1:K) {
  cat("Fold", k, ":\n")
  print(table(y[fold_id == k]))
}

# ---- helper: safe predict with best_nrounds ----
safe_predict <- function(model, dmat, best_nrounds) {
  p <- tryCatch(predict(model, dmat, iterationrange = c(1, best_nrounds)), error = function(e) NULL)
  if (is.null(p)) p <- tryCatch(predict(model, dmat, iterationrange = c(0, best_nrounds)), error = function(e) NULL)
  if (is.null(p)) p <- tryCatch(predict(model, dmat, ntreelimit = best_nrounds), error = function(e) NULL)
  if (is.null(p)) p <- predict(model, dmat)
  p
}

# ---- 6) params ----
params <- list(
  booster = "gbtree",
  objective = "multi:softprob",
  num_class = length(class_levels),
  eval_metric = "mlogloss",
  eta = 0.05,
  max_depth = 4,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  gamma = 0,
  lambda = 1,
  alpha = 0
)

label_map <- setNames(0:(length(class_levels)-1), class_levels)

# ---- 7) CV loop (OOF predictions) ----
oof_prob <- matrix(NA_real_, nrow = nrow(X), ncol = length(class_levels))
colnames(oof_prob) <- class_levels
oof_pred_class <- rep(NA_character_, nrow(X))

fold_metrics <- vector("list", K)
fold_models  <- vector("list", K)

for (k in 1:K) {
  message("\n====================")
  message("Fold ", k, "/", K)
  
  valid_idx <- which(fold_id == k)
  train_idx <- which(fold_id != k)
  
  y_train <- droplevels(y[train_idx])
  y_valid <- droplevels(y[valid_idx])
  
  message("Train class counts:"); print(table(y_train))
  message("Valid class counts:"); print(table(y_valid))
  
  # ---- 7.1) log2(x + pseudo): pseudo computed on TRAIN only ----
  pseudo <- apply(X[train_idx, , drop = FALSE], 2, function(col){
    pos <- col[is.finite(col) & col > 0]
    if (length(pos) == 0) 1e-9 else min(pos)/2
  })
  
  X_train <- log2(sweep(X[train_idx, , drop = FALSE], 2, pseudo, `+`))
  X_valid <- log2(sweep(X[valid_idx, , drop = FALSE], 2, pseudo, `+`))
  
  y_train_num <- unname(label_map[as.character(y_train)])
  y_valid_num <- unname(label_map[as.character(y_valid)])
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train_num)
  dvalid <- xgb.DMatrix(data = X_valid, label = y_valid_num)
  
  # ---- 7.2) train with early stopping ----
  fit <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 3000,
    evals = list(train = dtrain, valid = dvalid),
    early_stopping_rounds = 50,
    verbose = 1
  )
  
  best_nrounds <- as.integer(xgb.attr(fit, "best_iteration"))
  best_score   <- xgb.attr(fit, "best_score")
  if (is.na(best_nrounds) || length(best_nrounds) == 0) best_nrounds <- fit$niter
  
  message("Best nrounds = ", best_nrounds,
          if (!is.null(best_score)) paste0(" | best_score(valid mlogloss) = ", best_score) else "")
  
  # ---- 7.3) predict (OOF for this fold) ----
  pred_prob <- safe_predict(fit, dvalid, best_nrounds)
  pred_prob <- matrix(pred_prob, ncol = length(class_levels), byrow = TRUE)
  colnames(pred_prob) <- class_levels
  
  pred_class <- class_levels[max.col(pred_prob, ties.method = "first")]
  pred_class <- factor(pred_class, levels = class_levels)
  
  oof_prob[valid_idx, ] <- pred_prob
  oof_pred_class[valid_idx] <- as.character(pred_class)
  
  # ---- 7.4) fold metrics ----
  cm_k <- table(True = y_valid, Pred = pred_class)
  acc_k <- mean(pred_class == y_valid)
  
  cm_row_k <- prop.table(cm_k, 1)
  recall_k <- diag(cm_row_k)
  recall_k[!is.finite(recall_k)] <- NA_real_
  bal_acc_k <- mean(recall_k, na.rm = TRUE)
  
  fold_metrics[[k]] <- data.frame(
    Fold = k,
    n_train = length(train_idx),
    n_valid = length(valid_idx),
    best_nrounds = best_nrounds,
    best_score_valid_mlogloss = suppressWarnings(as.numeric(best_score)),
    accuracy = acc_k,
    balanced_accuracy = bal_acc_k,
    t(recall_k)
  )
  colnames(fold_metrics[[k]]) <- c(
    "Fold","n_train","n_valid","best_nrounds","best_score_valid_mlogloss",
    "accuracy","balanced_accuracy",
    paste0("recall_", class_levels)
  )
  
  readr::write_csv(as.data.frame(cm_k), file.path(outdir, paste0("fold", k, "_confusion_matrix_counts.csv")))
  
  fold_models[[k]] <- list(
    model = fit,
    pseudo = pseudo,
    best_nrounds = best_nrounds,
    best_score = best_score,
    train_idx = train_idx,
    valid_idx = valid_idx
  )
  saveRDS(fold_models[[k]], file.path(outdir, paste0("fold", k, "_model.rds")))
}

# ---- 8) CV summary ----
cv_df <- dplyr::bind_rows(fold_metrics)
print(cv_df)
readr::write_csv(cv_df, file.path(outdir, "cv5_fold_metrics.csv"))

cv_summary <- cv_df %>%
  summarise(
    folds = n(),
    accuracy_mean = mean(accuracy, na.rm = TRUE),
    accuracy_sd   = sd(accuracy, na.rm = TRUE),
    balacc_mean   = mean(balanced_accuracy, na.rm = TRUE),
    balacc_sd     = sd(balanced_accuracy, na.rm = TRUE),
    best_nrounds_mean = mean(best_nrounds, na.rm = TRUE),
    best_nrounds_median = median(best_nrounds, na.rm = TRUE)
  )
print(cv_summary)
readr::write_csv(cv_summary, file.path(outdir, "cv5_summary.csv"))

# ---- 9) OOF predictions export ----
oof_pred_class <- factor(oof_pred_class, levels = class_levels)

oof_df <- data.frame(
  SampleID = sample_id,
  Fold = fold_id,
  True = as.character(y),
  Pred = as.character(oof_pred_class),
  oof_prob,
  stringsAsFactors = FALSE
)
readr::write_csv(oof_df, file.path(outdir, "oof_predictions_by_sample.csv"))

# ---- 10) Overall OOF confusion matrix + plots ----
cm_all <- table(True = y, Pred = oof_pred_class)
print(cm_all)

acc_all <- mean(oof_pred_class == y)
recall_all <- diag(prop.table(cm_all, 1))
recall_all[!is.finite(recall_all)] <- NA_real_
bal_acc_all <- mean(recall_all, na.rm = TRUE)

cat("\nOOF Accuracy:", round(acc_all, 4), "\n")
cat("OOF Balanced Accuracy:", round(bal_acc_all, 4), "\n")
cat("OOF Recall by class:\n"); print(round(recall_all, 4))

cm_df <- as.data.frame(as.table(cm_all))
colnames(cm_df) <- c("True", "Pred", "N")
cm_df$True <- factor(cm_df$True, levels = class_levels)
cm_df$Pred <- factor(cm_df$Pred, levels = class_levels)

p_cm_counts <- ggplot(cm_df, aes(x = Pred, y = True, fill = N)) +
  geom_tile() +
  geom_text(aes(label = N), size = 4) +
  theme_classic() +
  labs(title = "Confusion Matrix (OOF, Counts)", x = "Predicted", y = "True")
ggsave(file.path(outdir, "confusion_matrix_oof_counts.pdf"), p_cm_counts, width = 6, height = 5)

cm_row <- prop.table(cm_all, 1)
cm_row_df <- as.data.frame(as.table(cm_row))
colnames(cm_row_df) <- c("True", "Pred", "P")
cm_row_df$True <- factor(cm_row_df$True, levels = class_levels)
cm_row_df$Pred <- factor(cm_row_df$Pred, levels = class_levels)

p_cm_row <- ggplot(cm_row_df, aes(x = Pred, y = True, fill = P)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", P)), size = 4) +
  theme_classic() +
  labs(title = "Confusion Matrix (OOF, Row-normalized)", x = "Predicted", y = "True")
ggsave(file.path(outdir, "confusion_matrix_oof_row_normalized.pdf"), p_cm_row, width = 6, height = 5)

# ---- 11) probability matrix mean (True x Pred prob) ----
prob_mat <- oof_df |>
  dplyr::group_by(True) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(class_levels), ~mean(.x, na.rm = TRUE)), .groups="drop")
readr::write_csv(prob_mat, file.path(outdir, "prediction_probability_matrix_oof_mean.csv"))

# ---- 12) AUC one-vs-rest on OOF ----
auc_list <- lapply(class_levels, function(cls){
  ybin <- as.integer(oof_df$True == cls)
  score <- oof_df[[cls]]
  r <- pROC::roc(response = ybin, predictor = score, quiet = TRUE)  # auto direction
  ci <- pROC::ci.auc(r)
  data.frame(
    Class = cls,
    AUC = as.numeric(pROC::auc(r)),
    CI_low = as.numeric(ci[1]),
    CI_high = as.numeric(ci[3])
  )
})
auc_df <- dplyr::bind_rows(auc_list) |>
  dplyr::mutate(AUC_macro = mean(AUC, na.rm = TRUE))
print(auc_df)
readr::write_csv(auc_df, file.path(outdir, "oof_auc_one_vs_rest.csv"))

# ---- 12C) Pairwise ROC (OOF): SO vs O / SO vs S / SO vs H + 95% CI band ----
pair_targets <- c("O", "S", "H")

ci_fill  <- "#0072B2"
ci_alpha <- 0.25
boot_n   <- 1000
ci_grid  <- seq(0, 1, by = 0.01)   # specificity grid

set.seed(20260114)

# --- helper: robust CI band builder for pROC::ci.se() (handles both matrix orientations) ---
make_ci_band_df <- function(roc_obj, ci_grid, boot_n) {
  ci_se  <- pROC::ci.se(roc_obj, specificities = ci_grid, boot.n = boot_n)
  ci_mat <- as.matrix(ci_se)
  
  # helper to parse numeric from names like "0.95", "95%", "2.5%" etc.
  parse_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9\\.]+", "", x)))
  
  # Case A: N x 3  (rows = specificities; cols = CI quantiles)
  if (nrow(ci_mat) == length(ci_grid) && ncol(ci_mat) == 3) {
    
    # specificity: prefer rownames; fallback to ci_grid
    spec_vec <- parse_num(rownames(ci_mat))
    if (all(is.na(spec_vec))) spec_vec <- ci_grid
    
    # ensure CI columns are in low-mid-high order (2.5%,50%,97.5%)
    q <- parse_num(colnames(ci_mat))
    if (all(!is.na(q))) {
      ord <- order(q)  # 2.5 < 50 < 97.5
      ci_mat <- ci_mat[, ord, drop = FALSE]
    }
    
    se_low  <- as.numeric(ci_mat[, 1])
    se_mid  <- as.numeric(ci_mat[, 2])
    se_high <- as.numeric(ci_mat[, 3])
    
    # Case B: 3 x N  (rows = CI quantiles; cols = specificities)
  } else if (nrow(ci_mat) == 3 && ncol(ci_mat) == length(ci_grid)) {
    
    # specificity: prefer colnames; fallback to ci_grid
    spec_vec <- parse_num(colnames(ci_mat))
    if (all(is.na(spec_vec))) spec_vec <- ci_grid
    
    # ensure CI rows are in low-mid-high order
    q <- parse_num(rownames(ci_mat))
    if (all(!is.na(q))) {
      ord <- order(q)
      ci_mat <- ci_mat[ord, , drop = FALSE]
    }
    
    se_low  <- as.numeric(ci_mat[1, ])
    se_mid  <- as.numeric(ci_mat[2, ])
    se_high <- as.numeric(ci_mat[3, ])
    
  } else {
    stop("ci.se() 输出矩阵维度异常：", nrow(ci_mat), " x ", ncol(ci_mat),
         "（无法识别是 3×N 还是 N×3）")
  }
  
  ci_df <- data.frame(
    spec    = as.numeric(spec_vec),
    se_low  = se_low,
    se_mid  = se_mid,
    se_high = se_high,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      FPR = 1 - spec,
      se_low  = pmax(0, se_low),
      se_high = pmin(1, se_high)
    ) %>%
    dplyr::filter(is.finite(FPR), is.finite(se_low), is.finite(se_high)) %>%
    dplyr::arrange(FPR)
  
  ci_df
}

auc_pair_list <- list()

for (cls in pair_targets) {
  
  sub <- oof_df %>%
    dplyr::filter(True %in% c("SO", cls)) %>%
    dplyr::mutate(
      ybin  = as.integer(True == "SO"),   # SO = 1, cls = 0
      score = as.numeric(.data[["SO"]])   # use P(SO)
    )
  
  if (length(unique(sub$True)) < 2) {
    warning("缺少类别：SO vs ", cls, " 中其中一类样本不存在，跳过。")
    next
  }
  
  r <- pROC::roc(response = sub$ybin, predictor = sub$score, quiet = TRUE)
  ci_auc <- pROC::ci.auc(r)
  
  auc_pair_list[[cls]] <- data.frame(
    Task = paste0("SO vs ", cls),
    n_SO = sum(sub$True == "SO"),
    n_other = sum(sub$True == cls),
    AUC = as.numeric(pROC::auc(r)),
    CI_low = as.numeric(ci_auc[1]),
    CI_high = as.numeric(ci_auc[3]),
    stringsAsFactors = FALSE
  )
  
  # ROC curve points (sorted)
  roc_df <- data.frame(
    FPR = 1 - r$specificities,
    TPR = r$sensitivities
  ) %>%
    dplyr::arrange(FPR, TPR) %>%
    dplyr::distinct(FPR, TPR, .keep_all = TRUE)
  
  # CI band
  ci_df <- make_ci_band_df(r, ci_grid = ci_grid, boot_n = boot_n)
  message("SO vs ", cls, " | CI band points = ", nrow(ci_df))
  
  p <- ggplot() +
    geom_ribbon(
      data = ci_df,
      aes(x = FPR, ymin = se_low, ymax = se_high),
      fill = ci_fill, alpha = ci_alpha,
      inherit.aes = FALSE
    ) +
    geom_step(
      data = roc_df,
      aes(x = FPR, y = TPR),
      direction = "hv",
      inherit.aes = FALSE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_classic() +
    labs(
      title = paste0("ROC (OOF): SO vs ", cls),
      subtitle = sprintf(
        "AUC=%.3f (95%% CI %.3f–%.3f), n_SO=%d, n_%s=%d | CI band: bootstrap (n=%d)",
        as.numeric(pROC::auc(r)), as.numeric(ci_auc[1]), as.numeric(ci_auc[3]),
        sum(sub$True == "SO"), cls, sum(sub$True == cls), boot_n
      ),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    )
  
  ggsave(
    filename = file.path(outdir, paste0("pairwise_roc_OOF_SO_vs_", cls, "_withCI.pdf")),
    plot = p, width = 6.5, height = 5.5
  )
  
  # export points
  readr::write_csv(ci_df,  file.path(outdir, paste0("pairwise_roc_OOF_SO_vs_", cls, "_CI_band.csv")))
  readr::write_csv(roc_df, file.path(outdir, paste0("pairwise_roc_OOF_SO_vs_", cls, "_curve_points.csv")))
}

auc_pair_df <- dplyr::bind_rows(auc_pair_list)
print(auc_pair_df)
readr::write_csv(auc_pair_df, file.path(outdir, "pairwise_auc_OOF_SO_vs_O_S_H.csv"))

# ---- 13) (optional) final model trained on ALL data using CV-chosen rounds ----
final_nrounds <- as.integer(round(cv_summary$best_nrounds_median))
message("\nTraining FINAL model on all data with nrounds = ", final_nrounds)

pseudo_all <- apply(X, 2, function(col){
  pos <- col[is.finite(col) & col > 0]
  if (length(pos) == 0) 1e-9 else min(pos)/2
})
X_all <- log2(sweep(X, 2, pseudo_all, `+`))
y_all_num <- unname(label_map[as.character(y)])
dall <- xgb.DMatrix(data = X_all, label = y_all_num)

final_fit <- xgb.train(
  params = params,
  data = dall,
  nrounds = final_nrounds,
  evals = list(train = dall),
  verbose = 1
)

final_obj <- list(
  model = final_fit,
  class_levels = class_levels,
  pseudo = pseudo_all,
  params = params,
  nrounds = final_nrounds,
  cv5_metrics = cv_df,
  fold_id = fold_id
)
saveRDS(final_obj, file.path(outdir, "xgb_model_final_trained_on_all_cv5.rds"))

message("\nDone. Outputs in: ", outdir)
# =========================
# 14) Variable importance (Top20) + SHAP importance (Top20)
# =========================

# ---- 14.0) feature names (确保不丢列名) ----
# X_all 是你 log2(pseudo) 后的矩阵；一般会保留列名
feature_names <- colnames(X_all)
if (is.null(feature_names) || any(feature_names == "")) {
  stop("X_all 的列名丢失了，无法做变量重要性。请确认输入 dat 的特征列有列名。")
}

# ---- 14.1) XGBoost built-in importance: Gain/Cover/Frequency ----
imp <- xgboost::xgb.importance(model = final_fit, feature_names = feature_names)

# 输出 Top20（按 Gain）
imp_top20 <- imp %>% dplyr::slice_max(order_by = Gain, n = 20, with_ties = FALSE)
print(imp_top20)

readr::write_csv(imp, file.path(outdir, "final_model_xgb_importance_all.csv"))
readr::write_csv(imp_top20, file.path(outdir, "final_model_xgb_importance_top20_byGain.csv"))

# 画图（Top20 by Gain）
p_imp <- ggplot(imp_top20, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    title = "XGBoost Feature Importance (Final model)",
    subtitle = "Top 20 by Gain",
    x = NULL,
    y = "Gain"
  )

ggsave(file.path(outdir, "final_model_xgb_importance_top20_byGain.pdf"),
       p_imp, width = 7, height = 6)

# ---- 14.2) SHAP: mean(|SHAP|) (robust for multiclass) ----
shap_raw <- predict(final_fit, dall, predcontrib = TRUE)

p  <- ncol(X_all)
Kc <- length(class_levels)
feature_names <- colnames(X_all)
if (is.null(feature_names) || any(feature_names == "")) {
  stop("X_all 的列名丢失了，无法做 SHAP。请确认输入特征列有列名。")
}

coerce_shap_array <- function(shap_raw, p, Kc) {
  # target: array [n, (p or p+1), K]
  if (is.array(shap_raw) && length(dim(shap_raw)) == 3) {
    d <- dim(shap_raw)
    
    # case A: n x (p or p+1) x K
    if (d[2] %in% c(p, p + 1) && d[3] == Kc) return(shap_raw)
    
    # case B: n x K x (p or p+1)  -> permute to n x (p or p+1) x K
    if (d[2] == Kc && d[3] %in% c(p, p + 1)) {
      return(aperm(shap_raw, c(1, 3, 2)))
    }
    
    stop("Unrecognized 3D SHAP dims: ", paste(d, collapse = " x "),
         " ; expected n x (p|p+1) x K OR n x K x (p|p+1).")
  }
  
  # matrix cases
  if (is.matrix(shap_raw)) {
    n <- nrow(shap_raw)
    m <- ncol(shap_raw)
    
    # common: n x ((p+1)*K)
    if (m == (p + 1) * Kc) return(array(shap_raw, dim = c(n, p + 1, Kc)))
    # sometimes: n x (p*K) (no bias)
    if (m == p * Kc) return(array(shap_raw, dim = c(n, p, Kc)))
    # binary/regression-like fallback
    if (m == (p + 1)) return(array(shap_raw, dim = c(n, p + 1, 1)))
    if (m == p) return(array(shap_raw, dim = c(n, p, 1)))
    
    stop("Unrecognized matrix SHAP shape: ncol=", m,
         " ; expected (p+1)*K=", (p+1)*Kc,
         " or p*K=", p*Kc,
         " or p+1=", p+1,
         " or p=", p, ".")
  }
  
  stop("Unsupported SHAP output type: ", paste(class(shap_raw), collapse = ", "))
}

shap_arr <- coerce_shap_array(shap_raw, p = p, Kc = Kc)
d <- dim(shap_arr)
message("[SHAP] shap_arr dims = ", paste(d, collapse = " x "), " (n x (p|p+1) x K)")

# 是否包含最后一列 BIAS
has_bias <- (d[2] == p + 1)

# 取特征贡献：n x p x K
if (has_bias) {
  shap_feat <- shap_arr[, 1:p, , drop = FALSE]
} else if (d[2] == p) {
  shap_feat <- shap_arr
} else {
  stop("After coercion, second dim is ", d[2], " but p is ", p, ".")
}

# mean(|SHAP|): across samples AND classes
mean_abs_shap <- apply(abs(shap_feat), 2, mean, na.rm = TRUE)

shap_imp <- data.frame(
  Feature = feature_names,
  MeanAbsSHAP = as.numeric(mean_abs_shap),
  stringsAsFactors = FALSE
) %>% dplyr::arrange(dplyr::desc(MeanAbsSHAP))

shap_top20 <- shap_imp %>% dplyr::slice_head(n = 20)
print(shap_top20)

readr::write_csv(shap_imp,  file.path(outdir, "final_model_shap_meanAbs_all.csv"))
readr::write_csv(shap_top20, file.path(outdir, "final_model_shap_meanAbs_top20.csv"))

if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
library(scales)

p_shap <- ggplot(shap_top20, aes(x = reorder(Feature, MeanAbsSHAP), y = MeanAbsSHAP)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.4f", MeanAbsSHAP)), hjust = -0.05, size = 3) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  expand_limits(y = max(shap_top20$MeanAbsSHAP) * 1.15) +
  labs(
    title = "SHAP Feature Importance (Final model)",
    subtitle = "Top 20 by mean(|SHAP|) across samples (and classes)",
    x = NULL,
    y = "mean(|SHAP|)"
  )

ggsave(file.path(outdir, "final_model_shap_meanAbs_top20.pdf"),
       p_shap, width = 14, height = 6)