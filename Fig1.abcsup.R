############################################
## Fig1a.Enterotype  +  Fig1b.TopGenus(伪对数)  +  Fig1c(SO比例)
## 附：导出中间结果 CSV 以用于附图/附表
############################################
while (!is.null(dev.list())) dev.off()
options(stringsAsFactors = FALSE)

# 需要的包（没有就 install.packages("包名")）
pkgs <- c("data.table","vegan","ape","cluster","dplyr","ggplot2","ellipse",
          "tidyr","readr","stringr","ggpubr","rstatix","scales","svglite")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, require, character.only = TRUE))

setwd("D:/Users/ruzuo/Desktop/菌群学习/【BGI菌群学习】/TaxonomyAnalysis/Analysis")

## ---------------- 统一颜色（按簇/ET 顺序依次使用这 4 色） ----------------
pal4 <- c("#EF4135","#F3A332","#0055A4","#018A67")  # SO/Sarc/Obe/Healthy 的同款配色
pal_group <- c(SO="#EF4135", Sarcopenia="#F3A332", Obesity="#0055A4", Healthy="#018A67")

############################################################
## 1. 读入属级相对丰度矩阵（保留 SampleID）
############################################################
g.raw <- data.table::fread("trans.genus.relative.csv", encoding = "UTF-8")

# 保存样本ID，去掉ID列
sampleID <- as.character(g.raw[[1]])
g.raw[[1]] <- NULL

# 若数据本身已是相对丰度（0-1 或 0-100），直接转矩阵
g.mat <- as.matrix(g.raw)
rownames(g.mat) <- sampleID

############################################################
## 1.5 依据“均值 ≥0.001% + 流行率 ≥10%”筛选属，并清理空样本
############################################################
# 判定是 0-1 比例还是 0-100 百分数
is_prop <- max(g.mat, na.rm = TRUE) <= 1
mean_cut <- if (is_prop) 1e-5 else 0.001   # 0.001%
prev_cut <- 0.10                            # 10%

# 计算每个属的均值与流行率（>0 视为出现）
tax_mean <- colMeans(g.mat, na.rm = TRUE)
tax_prev <- colMeans(g.mat > 0, na.rm = TRUE)

keep_cols <- (tax_mean >= mean_cut) & (tax_prev >= prev_cut)
g.use <- g.mat[, keep_cols, drop = FALSE]

# 删除筛选后变成全 0 的样本
keep_rows <- rowSums(g.use, na.rm = TRUE) > 0
g.use <- g.use[keep_rows, , drop = FALSE]
sampleID_use <- sampleID[keep_rows]

# （可选）重新按行归一化，确保行和=1
g.use <- vegan::decostand(g.use, method = "total")
rownames(g.use) <- sampleID_use

cat("筛选后保留属数：", ncol(g.use), "\n")
cat("筛选后保留样本数：", nrow(g.use), "（原始样本数：", length(sampleID), "）\n")

## —— 导出用于附表：筛选统计（均值、流行率、是否保留）
feature_stats <- data.frame(
  Genus = colnames(g.mat),
  mean  = as.numeric(tax_mean),
  prevalence = as.numeric(tax_prev),
  kept  = as.logical(colnames(g.mat) %in% colnames(g.use))
)
readr::write_csv(feature_stats, "S1_feature_filter_stats.csv")

## —— 导出用于附表：筛选后的属级矩阵（第一列 SampleID）
readr::write_csv(
  data.frame(SampleID = rownames(g.use), g.use, check.names = FALSE),
  "S2_genus_relative_filtered.csv"
)

# 基本健壮性检查
if (ncol(g.use) < 2 || nrow(g.use) < 3) {
  stop("筛选后特征数或样本数过少，无法进行后续分析。请放宽阈值或检查数据。")
}

############################################################
## 2. Bray–Curtis 距离 & PCoA
############################################################
bray  <- vegan::vegdist(g.use, method = "bray")
pcoa  <- ape::pcoa(bray)

coords <- as.data.frame(pcoa$vectors[, 1:2])
names(coords) <- c("PCoA1", "PCoA2")
rownames(coords) <- sampleID_use

############################################################
## 3. 用 Bray 距离 + PAM + silhouette 自动挑最佳 k
############################################################
k_range <- 2:min(8, nrow(g.use) - 1)   # 合理的 k 范围
sil_vec <- sapply(k_range, function(k) {
  pamk <- cluster::pam(bray, k = k, diss = TRUE)
  sil  <- cluster::silhouette(pamk$clustering, bray)
  mean(sil[, 3])  # 平均 silhouette 宽度
})

best.k <- k_range[which.max(sil_vec)]
cat("依据 Bray 距离的 silhouette，推荐 k =", best.k, "\n")

## —— 导出用于附表：不同 k 的 silhouette 宽度
readr::write_csv(data.frame(k = k_range, mean_silhouette = sil_vec),
                 "S3_silhouette_by_k.csv")

############################################################
## 4. PAM 聚类（如需固定k，可直接赋值）
############################################################
best.k <- 4
pam.res <- cluster::pam(bray, k = best.k, diss = TRUE)
coords$Cluster <- factor(pam.res$clustering)

############################################################
## 5. 给簇自动命名：取每簇平均丰度最高的属
############################################################
dominant <- sapply(unique(pam.res$clustering), function(k) {
  m <- g.use[pam.res$clustering == k, , drop = FALSE]
  names(which.max(colMeans(m)))
})

# 按簇编号顺序构造标签
levels(coords$Cluster) <- paste0(dominant, " (E", seq_along(dominant), ")")
coords$Cluster <- droplevels(coords$Cluster)  # 去掉可能的空水平

## —— 导出用于附表：每簇主导属
readr::write_csv(
  data.frame(Cluster = levels(coords$Cluster),
             DominantGenus = dominant[seq_along(levels(coords$Cluster))]),
  "S4_cluster_dominant_genus.csv"
)

############################################################
## 6. 计算 67% 置信椭圆
############################################################
ell_list <- lapply(levels(coords$Cluster), function(cl) {
  df <- coords[coords$Cluster == cl, c("PCoA1", "PCoA2")]
  e  <- ellipse::ellipse(stats::cov(df), centre = colMeans(df),
                         level = 0.67, npoints = 100)
  e <- as.data.frame(e); names(e) <- c("PCoA1","PCoA2")
  cbind(e, Cluster = cl)
})
ell <- do.call(rbind, ell_list)

## —— 导出用于附图：PCoA 坐标 + 簇标签；椭圆多边形坐标
readr::write_csv(
  data.frame(ID = rownames(coords), coords, Cluster = coords$Cluster),
  "S5_pcoa_coords_with_cluster.csv"
)
readr::write_csv(ell, "S6_pcoa_ellipse_coords.csv")

############################################################
## 7. 颜色映射（将 4 色按簇顺序分配）
############################################################
pal <- setNames(rep_len(pal4, nlevels(coords$Cluster)), levels(coords$Cluster))

############################################################
## 8. 绘图并保存（Fig1a）
############################################################
p <- ggplot(coords, aes(PCoA1, PCoA2, colour = Cluster)) +
  geom_polygon(data = ell, aes(fill = Cluster),
               colour = NA, alpha = 0.15, show.legend = FALSE) +
  geom_point(size = 2) +
  # 每簇中心（均值）用矩形表示
  geom_point(data = aggregate(cbind(PCoA1, PCoA2) ~ Cluster, coords, mean),
             shape = 15, size = 4, colour = "black") +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme_bw(base_size = 14) +
  labs(title = "Enterotype clustering (PCoA – Bray–Curtis)",
       x = sprintf("PCoA1 (%.1f%%)", pcoa$values$Rel_corr_eig[1] * 100),
       y = sprintf("PCoA2 (%.1f%%)", pcoa$values$Rel_corr_eig[2] * 100),
       colour = "Enterotype")

print(p)
ggsave("Enterotype_PCoA.png", p, width = 7, height = 6, dpi = 300)

cat("完成！图已保存为 Enterotype_PCoA.png\n")

############################################################
## 9. 导出 SampleID + Enterotype 标签表
############################################################
coords$ID <- rownames(coords)
coords$ET_code  <- stringr::str_extract(as.character(coords$Cluster), "E\\d+")
coords$ET_label <- as.character(coords$Cluster)
coords$ET_code  <- factor(coords$ET_code, levels = sort(unique(coords$ET_code)))

label.df <- data.frame(
  SampleID   = sampleID_use,
  Enterotype = as.character(coords$Cluster),
  ET_code    = coords$ET_code[match(sampleID_use, coords$ID)]
)
write.csv(label.df, "Enterotype_labels.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("已生成 Enterotype_labels.csv\n")

############################################################
## Fig1b.TopGenus  —— 伪对数轴 + 导出中间结果
############################################################
topN <- 5     # ← 想看前 5 个就改成 5

# 把丰度矩阵加上 Cluster 标签（使用筛选后 g.use 更稳）
ab.df <- as.data.frame(g.use)
ab.df$ID <- rownames(ab.df)
ab.df$Cluster <- coords$Cluster[match(ab.df$ID, coords$ID)]

# 求每簇每属均值
mean.tab <- ab.df %>%
  dplyr::select(-ID) %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), mean), .groups = "drop")

# 选每簇 Top N 属名（按均值降序）
top.list <- lapply(split(mean.tab, mean.tab$Cluster), function(df){
  df2 <- df %>% dplyr::select(-Cluster)
  vals <- unlist(df2[1, ], use.names = TRUE)
  head(names(sort(vals, decreasing = TRUE)), n = min(topN, length(vals)))
})
top.genus.all <- unique(unlist(top.list))

## —— 导出用于附表：各簇 TopN 属清单（含均值）
top_rank_df <- mean.tab |>
  tidyr::pivot_longer(-Cluster, names_to = "Genus", values_to = "MeanAbundance") |>
  dplyr::group_by(Cluster) |>
  dplyr::arrange(Cluster, dplyr::desc(MeanAbundance), .by_group = TRUE) |>
  dplyr::mutate(Rank = dplyr::row_number(),
                InTopN = Rank <= topN) |>
  dplyr::ungroup()
readr::write_csv(top_rank_df, "S7_top_genus_rank_by_cluster.csv")

# 构造成 long；转为百分比并过滤全 0 面板
plot.df <- ab.df %>%
  dplyr::select(dplyr::all_of(c(top.genus.all,"ID","Cluster"))) %>%
  tidyr::pivot_longer(cols = all_of(top.genus.all),
                      names_to = "Genus", values_to = "Abundance") %>%
  dplyr::mutate(AbundancePct = Abundance * 100) %>%
  dplyr::group_by(Genus) %>%
  dplyr::filter(sum(AbundancePct, na.rm = TRUE) > 0) %>%  # 去掉全 0 的属
  dplyr::ungroup()

## —— 导出用于附图：Fig1b 的 long 数据
readr::write_csv(plot.df, "S8_fig1b_longdata_topN_genus.csv")

# 统计检验（KW + 两两 Wilcoxon BH）
kw.res <- plot.df %>%
  dplyr::group_by(Genus) %>%
  rstatix::kruskal_test(AbundancePct ~ Cluster)
readr::write_csv(kw.res, "TopGenus_kw_test.csv")

pw.res <- plot.df %>%
  dplyr::group_by(Genus) %>%
  rstatix::pairwise_wilcox_test(AbundancePct ~ Cluster,
                                p.adjust.method = "BH",
                                exact = FALSE) %>%
  dplyr::arrange(Genus, group1, group2)
readr::write_csv(pw.res, "TopGenus_pairwise_wilcox_BH.csv")

# 基础箱线图 + 伪对数轴（sigma 决定“靠近 0 的线性区间”大小，0.1≈0.1%）
pal_cluster <- setNames(rep_len(pal4, nlevels(coords$Cluster)), levels(coords$Cluster))

p.box <- ggplot(plot.df, aes(x = Cluster, y = AbundancePct, fill = Cluster)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.35) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = pal_cluster) +
  scale_y_continuous(
    trans  = scales::pseudo_log_trans(sigma = 0.1),
    breaks = c(0, 0.01, 0.1, 1, 5, 10, 25, 50),
    labels = function(x) paste0(x, "%"),
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(title = paste0("Top ", topN, " genera of each enterotype\n",
                      paste(names(table(coords$Cluster)),
                            " n=", as.integer(table(coords$Cluster)),
                            collapse = "; ")),
       x = "Enterotype",
       y = "Relative abundance (%)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"))

# 两两比较（Wilcoxon + BH），并在各自 facet 内标注显著性
comp.list <- combn(levels(coords$Cluster), 2, simplify = FALSE)
p.box2 <- p.box +
  ggpubr::stat_compare_means(
    comparisons     = comp.list,
    method          = "wilcox.test",
    method.args     = list(exact = FALSE),
    p.adjust.method = "BH",
    label           = "p.signif",
    hide.ns         = TRUE,
    size            = 3,
    step.increase   = 0.08
  )

print(p.box2)
ggsave("Enterotype_TopGenus_Boxplot_pairwise_pseudolog.png", p.box2,
       width = 12, height = 9, dpi = 300)

cat("Fig1b 完成（伪对数轴）。\n")

############################################################
## Fig1c: Enterotype proportions by SO groups (mappingV1)
## - 将 Enterotype（E1/E2/E3...）与 mapping 合并
## - 按 SO/Sarcopenia/Obesity/Healthy 画比例堆叠柱
## - 自动做列联检验并在副标题标注 p 值（期望<5则模拟卡方）
############################################################
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(stringr); library(ggplot2); library(scales)
  library(readr); library(svglite)
})

# ---- 1) 读取 mapping，并标准化 ID / 组别标签 ----
map_file <- "mappingV1.csv"
mapping  <- data.table::fread(map_file, na.strings = c("", "NA")) %>% as.data.frame()

# 兼容大小写：把 ID 列标准化为 "ID"
idcol <- names(mapping)[tolower(names(mapping)) == "id"][1]
if (!is.na(idcol) && idcol != "ID") names(mapping)[names(mapping) == idcol] <- "ID"
mapping$ID <- trimws(as.character(mapping$ID))

if (!"group_code" %in% names(mapping)) {
  stop("mappingV1.csv 中缺少 group_code（应为 1=SO, 2=Sarcopenia, 3=Obesity, 4=Healthy）")
}

mapping <- mapping %>%
  dplyr::mutate(
    group_code  = trimws(as.character(group_code)),
    group_label = dplyr::case_when(
      group_code == "1" ~ "SO",
      group_code == "2" ~ "Sarcopenia",
      group_code == "3" ~ "Obesity",
      group_code == "4" ~ "Healthy",
      TRUE ~ NA_character_
    ),
    group_label = factor(group_label, levels = c("SO","Sarcopenia","Obesity","Healthy"))
  )

# ---- 2) 合并 Enterotype ↔ mapping ----
df_et <- coords %>%
  dplyr::select(ID, ET_code, ET_label) %>%
  dplyr::left_join(mapping %>% dplyr::select(ID, group_label), by = "ID") %>%
  dplyr::filter(!is.na(group_label), !is.na(ET_code))

## —— 导出用于附表：ID-ET-Group 对照
readr::write_csv(df_et, "S9_enterotype_with_group.csv")

# 列联表与 p 值（期望频数<5则启用模拟卡方）
tab  <- table(df_et$ET_code, df_et$group_label)
pval <- NA_real_
if (all(dim(tab) >= 2)) {
  cs <- suppressWarnings(chisq.test(tab, correct = FALSE))
  pval <- cs$p.value
  if (any(cs$expected < 5)) {
    set.seed(123)
    pval <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 9999)$p.value)
  }
}
p_to_star <- function(p) ifelse(is.na(p), "n.s.",
                                ifelse(p < 0.001, "***",
                                       ifelse(p < 0.01,  "**",
                                              ifelse(p < 0.05,  "*", "n.s."))))
subtitle_txt <- if (is.na(pval)) {
  "Group × Enterotype: insufficient sample size for p-value"
} else {
  sprintf("Group × Enterotype: chi-square test  p = %.3g  [%s]", pval, p_to_star(pval))
}

plot_df <- df_et %>%
  dplyr::count(group_label, ET_code, name = "n") %>%
  dplyr::group_by(group_label) %>%
  dplyr::mutate(pct = n / sum(n), label_txt = scales::percent(pct, accuracy = 1)) %>%
  dplyr::ungroup()

## —— 导出用于附表：列联表（宽）与组内百分比
readr::write_csv(as.data.frame(tab), "S10_enterotype_by_SO_counts.csv")
prop_tab <- prop.table(tab, margin = 2) %>% as.data.frame.matrix()
readr::write_csv(as.data.frame(prop_tab), "S11_enterotype_by_SO_colprop.csv")

# ---- 3) 颜色：将 4 色映射到 ET_code（E1/E2/E3/E4...按顺序）----
et_lvls <- levels(df_et$ET_code)
pal_et  <- setNames(rep_len(pal4, length(et_lvls)), et_lvls)

legend_map    <- df_et %>% dplyr::distinct(ET_code, ET_label) %>% dplyr::arrange(ET_code)
legend_labels <- setNames(legend_map$ET_label, legend_map$ET_code)

# ---- 4) 绘图并导出 ----
p_prop <- ggplot(plot_df, aes(x = group_label, y = pct, fill = ET_code)) +
  geom_col(width = 0.7, colour = "white") +
  geom_text(aes(label = label_txt),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = pal_et,
                    breaks = names(legend_labels),
                    labels = unname(legend_labels),
                    drop = FALSE, na.translate = FALSE) +
  labs(x = NULL, y = "Proportion", fill = "Enterotype",
       title = "Enterotype proportions by SO groups",
       subtitle = subtitle_txt) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.0, face = "bold"),
    plot.subtitle = element_text(hjust = 0, color = "gray25"),
    axis.text.x = element_text(angle = 0, vjust = 0.9, colour = pal_group[levels(mapping$group_label)])
  )

print(p_prop)
ggsave("Enterotype_proportion_by_SO.png", p_prop, width = 8.8, height = 5.2, dpi = 300)
ggsave("Enterotype_proportion_by_SO.svg", p_prop, width = 8.8, height = 5.2, device = svglite::svglite)

cat("全部完成：主图 + 伪对数 Fig1b + SO 比例图，并导出 S1–S11 多个中间结果 CSV。\n")
