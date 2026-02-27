############################################################
## Alpha 多样性 + 配对检验 + GLM(图g风格) 完整脚本
## 需求：
##  - 输入：trans.genus.relative.csv（行=样本，列=属；第一列为样本ID）
##         mappingV1.csv（至少含 ID, group_code；1=SO,2=Sarcopenia,3=Obesity,4=Healthy）
##  - 可选输入（用于图g的门水平特征）：
##         trans.phylum.relative.csv  或  genus2phylum.csv(genus,phylum)
##  - 输出：多个CSV与PNG见各 write_csv/ggsave 调用
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(broom)
  library(forcats)
})

setwd("D:\\Users\\ruzuo\\Desktop\\菌群学习\\【BGI菌群学习】\\TaxonomyAnalysis\\Analysis")

# ========= 0) 路径设置（按需修改） =========
abund_file   <- "trans.genus.relative.csv"  # 属级相对丰度
mapping_file <- "mappingV1.csv"

# ========= 1) 读入 mapping =========
mapping <- read_csv(mapping_file, guess_max = 1e6, show_col_types = FALSE)

# 统一 ID 列名
if (!"ID" %in% names(mapping)) {
  id_col <- names(mapping)[tolower(names(mapping)) == "id"][1]
  if (is.na(id_col)) stop("mappingV1.csv 中缺少 ID 列")
  mapping <- mapping %>% rename(ID = all_of(id_col))
}
mapping <- mapping %>%
  mutate(ID = as.character(ID) %>% str_trim())

if (!"group_code" %in% names(mapping)) {
  stop("mappingV1.csv 中缺少 group_code（应为 1=SO, 2=Sarcopenia, 3=Obesity, 4=Healthy）")
}
mapping <- mapping %>%
  mutate(group_code = as.character(group_code) %>% str_trim())

# ========= 2) 读入属级相对丰度 & 与肠型相同的筛选（均值与流行率）=========
g.raw <- data.table::fread(abund_file, encoding = "UTF-8")
sampleID <- as.character(g.raw[[1]])
g.raw[[1]] <- NULL
g.mat <- as.matrix(g.raw)
rownames(g.mat) <- sampleID

# 均值≥0.001% 与 流行率≥10%
is_prop <- max(g.mat, na.rm = TRUE) <= 1
mean_cut <- if (is_prop) 1e-5 else 0.001   # 0.001%
prev_cut <- 0.10                           # 10%

tax_mean <- colMeans(g.mat, na.rm = TRUE)
tax_prev <- colMeans(g.mat > 0, na.rm = TRUE)
keep_cols <- (tax_mean >= mean_cut) & (tax_prev >= prev_cut)
g.use <- g.mat[, keep_cols, drop = FALSE]

# 删除全0样本并按行归一化（总和=1）
keep_rows <- rowSums(g.use, na.rm = TRUE) > 0
g.use <- g.use[keep_rows, , drop = FALSE]
g.use <- vegan::decostand(g.use, method = "total")

cat("筛选后保留属数：", ncol(g.use), "\n")
cat("筛选后保留样本数：", nrow(g.use), "（原始样本数：", nrow(g.mat), "）\n")

# ========= 3) 计算 α 多样性（H 为 ln；Simpson 为 1−D；Pielou=H/lnS） =========
S_obs   <- rowSums(g.use > 0)                                # N0（观测到的丰富度）
H       <- vegan::diversity(g.use, index = "shannon")        # Shannon(H), ln 底
GiniSim <- vegan::diversity(g.use, index = "simpson")        # Gini–Simpson = 1 − D
InvSim  <- vegan::diversity(g.use, index = "invsimpson")     # 1/D

# Pielou（对 S<=1 的样本设 NA 以避免除零）
Pielou <- ifelse(S_obs > 1, H / log(S_obs), NA_real_)

alpha_wide <- data.frame(
  ID                  = rownames(g.use),
  N0                  = as.numeric(S_obs),
  shannon             = as.numeric(H),
  simpson             = as.numeric(GiniSim),                 # 1−D
  `Shannon diversity` = as.numeric(exp(H)),                  # Hill数 q=1
  revSimpson          = as.numeric(InvSim),                  # 1/D
  Shannonevenn        = as.numeric(exp(H) / pmax(S_obs, 1)),# 与历史文件保持一致：exp(H)/S
  Simpsonevenn        = as.numeric(InvSim / pmax(S_obs, 1)),# 与历史文件保持一致：InvSim/S
  Pielou              = as.numeric(Pielou),
  check.names         = FALSE
)
write_csv(alpha_wide, "alpha_diversity.filtered.csv")

# ========= 4) 整理为 long 并与 mapping 合并 =========
alpha_long <- alpha_wide %>%
  dplyr::rename(
    shannon_diversity = `Shannon diversity`,
    revsimpson        = revSimpson,
    shannonevenn      = Shannonevenn,
    simpsonevenn      = Simpsonevenn,
    pielou            = Pielou,
    n0                = N0
  ) %>%
  tidyr::pivot_longer(
    cols = c("shannon","simpson","pielou","n0",
             "shannon_diversity","revsimpson","shannonevenn","simpsonevenn"),
    names_to = "metric", values_to = "value"
  )

dat <- alpha_long %>%
  inner_join(mapping %>% select(ID, group_code), by = "ID") %>%
  mutate(
    group_code  = as.character(group_code),
    group_label = case_when(
      group_code == "1" ~ "SO",
      group_code == "2" ~ "Sarcopenia",
      group_code == "3" ~ "Obesity",
      group_code == "4" ~ "Healthy",
      TRUE ~ NA_character_
    ),
    group_label = factor(group_label, levels = c("SO","Sarcopenia","Obesity","Healthy")),
    group_code  = factor(group_code,  levels = c("1","2","3","4"))
  )

alpha_summary <- dat %>%
  group_by(metric, group_label, group_code) %>%
  summarise(
    n      = sum(!is.na(value)),
    mean   = mean(value, na.rm = TRUE),
    sd     = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    iqr    = IQR(value, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(metric, group_code)

write_csv(dat, "alpha_by_group_used_samples.filtered.csv")
write_csv(alpha_summary, "alpha_by_group_summary.filtered.csv")

# ========= 5) 组间两两比较（Wilcoxon，BH校正），并在图上标 p 值 =========
pal <- c("SO"="#EF4135","Sarcopenia"="#F3A332","Obesity"="#0055A4","Healthy"="#018A67")

pairwise_all <- dat %>%
  filter(!is.na(value), !is.na(group_label)) %>%
  group_by(metric) %>%
  group_modify(~{
    df <- .x
    lvls <- levels(df$group_label)
    if (is.null(lvls) || length(lvls) < 2) return(data.frame())
    prs <- t(combn(lvls, 2))
    res <- lapply(seq_len(nrow(prs)), function(i){
      g1 <- prs[i, 1]; g2 <- prs[i, 2]
      x <- df$value[df$group_label == g1]
      y <- df$value[df$group_label == g2]
      p <- if (sum(is.finite(x)) >= 3 && sum(is.finite(y)) >= 3) {
        tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
      } else NA_real_
      data.frame(g1 = g1, g2 = g2, p = p)
    }) %>% bind_rows()
    if (nrow(res)) res$padj <- p.adjust(res$p, method = "BH")
    res
  }) %>% ungroup()

write_csv(pairwise_all, "alpha_pairwise_all_wilcox_BH.csv")

# 准备注释坐标
ymax_df <- dat %>% group_by(metric, group_label) %>%
  summarise(ymax = max(value, na.rm = TRUE), .groups = "drop")
yrange_df <- dat %>% group_by(metric) %>%
  summarise(r = diff(range(value, na.rm = TRUE)), .groups = "drop")
pos_map <- data.frame(group_label = levels(dat$group_label),
                      x = seq_along(levels(dat$group_label)))

ann_pairs <- pairwise_all %>%
  filter(is.finite(p)) %>%
  left_join(pos_map, by = c("g1" = "group_label")) %>% rename(x1 = x) %>%
  left_join(pos_map, by = c("g2" = "group_label")) %>% rename(x2 = x) %>%
  left_join(ymax_df, by = c("metric", "g1" = "group_label")) %>% rename(y1 = ymax) %>%
  left_join(ymax_df, by = c("metric", "g2" = "group_label")) %>% rename(y2 = ymax) %>%
  left_join(yrange_df, by = "metric") %>%
  mutate(
    r    = ifelse(is.finite(r) & r > 0, r, 1.0),
    base = pmax(y1, y2, na.rm = TRUE) + 0.03 * r,
    span = abs(x2 - x1)
  ) %>%
  group_by(metric) %>%
  arrange(span, x1, .by_group = TRUE) %>%
  mutate(
    tier = row_number() - 1L,
    y    = base + tier * (0.06 * r),
    h    = 0.02 * r,
    xmid = (x1 + x2) / 2,
    label = sprintf("p=%.3g", p)          # 如需显示校正值改为 padj
  ) %>% ungroup()

# α多样性分组图（带p值）
p_alpha <- ggplot(dat, aes(x = group_label, y = value)) +
  geom_violin(aes(fill = group_label), trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(aes(fill = group_label), width = 0.15, outlier.shape = NA,
               alpha = 0.85, color = "white") +
  geom_jitter(aes(color = group_label), width = 0.1, alpha = 0.5, size = 0.9) +
  geom_segment(data = ann_pairs,
               aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  geom_segment(data = ann_pairs,
               aes(x = x1, xend = x1, y = y - h, yend = y),
               inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  geom_segment(data = ann_pairs,
               aes(x = x2, xend = x2, y = y - h, yend = y),
               inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  geom_text(data = ann_pairs,
            aes(x = xmid, y = y + 0.9 * h, label = label),
            inherit.aes = FALSE, size = 3.8, color = "black") +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Group", y = "Alpha diversity",
       title = "Alpha diversity by group (after genus filtering)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = pal, drop = FALSE, na.translate = FALSE) +
  scale_color_manual(values = pal, drop = FALSE, na.translate = FALSE)

print(p_alpha)
ggsave("alpha_diversity_by_group.filtered.png", p_alpha, width = 12, height = 7, dpi = 300)

# ========= 6) 图 g：GLM 回归（Shannon H, Pielou J, Gini–Simpson 1−D） =========
# 读取门水平丰度（优先）；若无则以属聚合到门；都无则用属水平直接作图
read_feature_matrix <- function(){
  if (file.exists("trans.phylum.relative.csv")) {
    m <- data.table::fread("trans.phylum.relative.csv", encoding = "UTF-8")
    level <- "Phylum"
  } else if (file.exists("genus2phylum.csv")) {
    m.genus <- data.table::fread(abund_file, encoding = "UTF-8")
    map <- readr::read_csv("genus2phylum.csv", show_col_types = FALSE)
    stopifnot(all(c("genus","phylum") %in% names(map)))
    samp <- as.character(m.genus[[1]]); m.genus[[1]] <- NULL
    mm <- as.matrix(m.genus); rownames(mm) <- samp
    gcols <- colnames(mm)
    map <- map %>% filter(genus %in% gcols)
    if (nrow(map) == 0) stop("genus2phylum.csv 中的 genus 与丰度文件列名未能匹配。")
    split_idx <- split(match(map$genus, gcols), map$phylum)
    m.mat <- sapply(split_idx, function(idx) rowSums(mm[, idx, drop = FALSE], na.rm = TRUE))
    m <- data.frame(ID = rownames(mm), m.mat, check.names = FALSE)
    level <- "Phylum"
  } else {
    m <- data.table::fread(abund_file, encoding = "UTF-8")
    level <- "Genus"
  }
  list(df = m, level = level)
}

feat <- read_feature_matrix()
feat_df <- feat$df
tax_level <- feat$level

feat_id <- as.character(feat_df[[1]])
feat_df[[1]] <- NULL
feat_mat <- as.matrix(feat_df); rownames(feat_mat) <- feat_id

# 与前述相同的均值+流行率过滤 & 行归一化
is_prop_feat <- max(feat_mat, na.rm = TRUE) <= 1
mean_cut_feat <- if (is_prop_feat) 1e-5 else 0.001
prev_cut_feat <- 0.10
tax_mean_feat <- colMeans(feat_mat, na.rm = TRUE)
tax_prev_feat <- colMeans(feat_mat > 0, na.rm = TRUE)
keep_cols_feat <- (tax_mean_feat >= mean_cut_feat) & (tax_prev_feat >= prev_cut_feat)
feat_use <- feat_mat[, keep_cols_feat, drop = FALSE]
keep_rows_feat <- rowSums(feat_use, na.rm = TRUE) > 0
feat_use <- feat_use[keep_rows_feat, , drop = FALSE]
feat_use <- vegan::decostand(feat_use, method = "total")

# 与 alpha 数据交集的样本
common_ids <- intersect(rownames(feat_use), dat$ID)
feat_use <- feat_use[common_ids, , drop = FALSE]
dat_sub <- dat %>% filter(ID %in% common_ids)

# 取丰度最高的前K个特征（避免面板拥挤）
K <- 7
top_taxa <- names(sort(colMeans(feat_use, na.rm = TRUE), decreasing = TRUE))[seq_len(min(K, ncol(feat_use)))]

# 仅使用三项：Shannon(H, ln), Pielou(J), Gini–Simpson(1−D)
metrics_for_glm <- c("shannon","pielou","simpson")
metric_labels <- c(
  shannon = "Shannon (H, ln)",
  pielou  = "Pielou's evenness (J = H/ln S)",
  simpson = "Gini–Simpson (1 − D)"
)

# 拟合 GLM： value ~ abundance
glm_rows <- list()
for (mtr in metrics_for_glm) {
  for (tx in top_taxa) {
    df <- data.frame(ID = rownames(feat_use),
                     abund = as.numeric(feat_use[, tx])) %>%
      left_join(dat_sub %>% filter(metric == mtr) %>% select(ID, value), by = "ID") %>%
      filter(is.finite(abund), is.finite(value))
    if (nrow(df) >= 8 && stats::sd(df$abund) > 0) {
      fit <- stats::glm(value ~ abund, data = df, family = gaussian())
      co <- broom::tidy(fit, conf.int = TRUE)
      co_abund <- co %>% filter(term == "abund") %>%
        mutate(metric = mtr, taxon = tx)
      glm_rows[[length(glm_rows) + 1]] <- co_abund
    }
  }
}
glm_tab <- bind_rows(glm_rows) %>%
  transmute(
    metric,
    taxon,
    beta = estimate,
    conf.low, conf.high,
    p.value
  )
write_csv(glm_tab, "Supp_Table_alpha_glm_coefficients.csv")  # 精确β与P值

# 可视化（图g样式的三联系数森林图）
glm_tab <- glm_tab %>%
  mutate(
    metric_lab = factor(metric, levels = metrics_for_glm, labels = metric_labels[metrics_for_glm]),
    sig_cat = case_when(
      p.value < 1e-4 ~ "p < 1e-4",
      p.value < 1e-3 ~ "p < 1e-3",
      p.value < 1e-2 ~ "p < 1e-2",
      TRUE ~ "NS"
    )
  ) %>%
  group_by(metric_lab) %>%
  mutate(taxon_f = forcats::fct_reorder(taxon, beta)) %>%
  ungroup()

pal_sig <- c("p < 1e-4"="#D7301F","p < 1e-3"="#FC8D59","p < 1e-2"="#3182BD","NS"="#C7B8A2")

p_g <- ggplot(glm_tab, aes(x = beta, y = taxon_f, color = sig_cat)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  facet_wrap(~ metric_lab, scales = "free_x") +
  scale_color_manual(values = pal_sig, name = "Significance") +
  labs(
    x = expression(beta*"-coefficient"),
    y = tax_level,
    title = sprintf("GLM: Shannon (H, ln), Pielou (J), and Gini–Simpson (1−D) vs %s abundance", tolower(tax_level)),
    subtitle = "Exact β and P values are provided in: Supp_Table_alpha_glm_coefficients.csv"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "right")

print(p_g)
ggsave("alpha_glm_forest_H_J_1minusD.png", p_g, width = 10, height = 6, dpi = 300)

############################################################
## 结束
############################################################

