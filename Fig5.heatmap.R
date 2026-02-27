library(ComplexHeatmap)
library(circlize)
library(grid)

## ===================== 0) 文件路径 =====================
infile <- "D:/Users/ruzuo/Desktop/菌群学习/Analysis/abund.csv"

## ===================== 1) 读入数据 =====================
dat <- read.csv(infile, check.names = FALSE)
stopifnot("SampleID" %in% colnames(dat))

sample_id <- as.character(dat$SampleID)

X <- as.matrix(dat[, setdiff(colnames(dat), "SampleID"), drop = FALSE])
rownames(X) <- sample_id
mode(X) <- "numeric"

## ===================== 2) 从 SampleID 自动分组（按前缀） =====================
prefix <- sub("-.*$", "", sample_id)

## 先判 SO，否则会被 S 吃掉
grp <- ifelse(grepl("^SO", prefix), "SO",
              ifelse(grepl("^S",  prefix), "Sarcopenia",
                     ifelse(grepl("^O", prefix), "Obesity",
                            ifelse(grepl("^H", prefix), "Healthy", NA))))

if (any(is.na(grp))) {
  warning("有样本分组没识别出来：", paste(unique(sample_id[is.na(grp)]), collapse = ", "))
}

## ✅ 组别顺序（列顺序）：SO → Sarcopenia → Obesity → Healthy
grp <- factor(grp, levels = c("SO", "Sarcopenia", "Obesity", "Healthy"))

cat("样本分组计数：\n")
print(table(grp))

## 用实际存在的组（避免空组）
group_levels <- levels(droplevels(grp))

## ===================== 3) 预处理：log10（可选） =====================
## 如果你不想 log10，改成：X_use <- X
X_use <- log10(X + 1e-12)

## ===================== 4) 计算组均值矩阵（特征/代谢物 × 组） =====================
mean_mat <- sapply(group_levels, function(g) {
  colMeans(X_use[grp == g, , drop = FALSE], na.rm = TRUE)
})
mean_mat <- as.matrix(mean_mat)
colnames(mean_mat) <- group_levels
rownames(mean_mat) <- colnames(X_use)

## ===================== 5) 选 Top30（KW + FDR） =====================
kw_p <- apply(X_use, 2, function(v) suppressWarnings(kruskal.test(v ~ grp)$p.value))
kw_fdr <- p.adjust(kw_p, method = "fdr")

N <- 30
top_met <- names(sort(kw_fdr))[1:N]

mean_top <- mean_mat[top_met, , drop = FALSE]

## 行 Z-score（每个代谢物在4组间标准化）
z_top <- t(scale(t(mean_top)))
z_top[is.na(z_top)] <- 0

## ===================== 6) 星号：每组 vs Healthy（Wilcoxon + FDR<0.05） =====================
ref_group <- "Healthy"
if (!(ref_group %in% group_levels)) stop("当前数据中没有 Healthy 组，无法进行 vs Healthy 的显著性比较。")

pair_p <- sapply(group_levels, function(g) {
  if (g == ref_group) return(rep(NA_real_, length(top_met)))
  sapply(top_met, function(m) {
    suppressWarnings(wilcox.test(X_use[grp == g, m], X_use[grp == ref_group, m])$p.value)
  })
})

pair_p <- as.matrix(pair_p)
rownames(pair_p) <- top_met
colnames(pair_p) <- group_levels

pair_fdr <- apply(pair_p, 2, function(pv) p.adjust(pv, method = "fdr"))
pair_fdr <- as.matrix(pair_fdr)
rownames(pair_fdr) <- top_met
colnames(pair_fdr) <- group_levels

sig_star <- matrix("", nrow = length(top_met), ncol = length(group_levels),
                   dimnames = list(top_met, group_levels))
sig_star[, ref_group] <- ""

other_cols <- setdiff(group_levels, ref_group)
sig_star[, other_cols] <- ifelse(pair_fdr[, other_cols] < 0.05, "*", "")

## ===================== 7) 行排序：按 SO 列从红到蓝 =====================
sort_col <- if ("SO" %in% colnames(z_top)) "SO" else group_levels[1]
ord <- order(z_top[, sort_col], decreasing = TRUE)

z_top <- z_top[ord, , drop = FALSE]
sig_star <- sig_star[rownames(z_top), , drop = FALSE]
mean_top <- mean_top[rownames(z_top), , drop = FALSE]

## ===================== 8) 顶部 signature model 柱状 =====================
model_vals <- colMeans(mean_top, na.rm = TRUE)
model_vals <- (model_vals - min(model_vals)) / (max(model_vals) - min(model_vals) + 1e-12)

grp_cols <- c(
  SO = "#EF4135",
  Sarcopenia = "#F3A332",
  Obesity = "#0055A4",
  Healthy = "#018A67"
)

ha <- HeatmapAnnotation(
  `signature model` = anno_barplot(
    model_vals,
    border = TRUE,
    gp = gpar(fill = grp_cols[group_levels]),
    height = unit(1.2, "cm")
  )
)

## ===================== 9) 画图并导出 =====================
col_fun <- colorRamp2(c(-2, 0, 2), c("#0055A4", "#FFFFFF", "#EF4135"))

ht <- Heatmap(
  z_top,
  name = "Z-score",
  col = col_fun,
  top_annotation = ha,
  cluster_rows = FALSE,     # 已手动按 SO 排序
  cluster_columns = FALSE,  # 列按组顺序固定
  rect_gp = gpar(col = "black", lwd = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 0,
  show_heatmap_legend = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (sig_star[i, j] != "") {
      grid.text(sig_star[i, j], x, y, gp = gpar(fontsize = 14))
    }
  }
)

out_pdf <- "signature_heatmap_top30_SO_to_Healthy.pdf"
pdf(out_pdf, width = 6.2, height = 9.2)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

cat("已输出：", out_pdf, "\n")
