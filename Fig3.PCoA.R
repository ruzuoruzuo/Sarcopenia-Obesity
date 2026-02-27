## =========================
## Fig3. PCoA (ONLY Bray) with confounder adjustment + p1–p4 layout
## =========================
while (!is.null(dev.list())) dev.off()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table); library(vegan); library(ape); library(dplyr)
  library(ggplot2); library(readr); library(tidyr); library(stringr)
  library(patchwork); library(multcomp); library(scales); library(grid)
  library(showtext)
})

## ---- 0) 主题 & 字体（可选）----
if (file.exists("C:/Windows/Fonts/msyh.ttc")) {
  sysfonts::font_add("CJK", regular = "C:/Windows/Fonts/msyh.ttc")
  showtext::showtext_auto(TRUE)
  theme_set(theme_bw(base_size = 12) + theme(text = element_text(family = "CJK")))
} else {
  theme_set(theme_bw(base_size = 12))
}

## ---- 1) 路径（按需修改）----
setwd("D:/Users/ruzuo/Desktop/菌群学习/【BGI菌群学习】/TaxonomyAnalysis/Analysis")
abund_file <- "trans.genus.relative.csv"   # 第1列=SampleID；不筛选，仅去全0列/行归一化
map_file   <- "mappingV1.csv"
out_png    <- "Fig3_PCoA_SO_adjBray.png"
out_pdf    <- "Fig3_PCoA_SO_adjBray.pdf"
out_dir    <- "F3_confounded"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ---- 2) 统一配色（与项目一致）----
pal_group <- c(SO="#EF4135", Sarcopenia="#F3A332", Obesity="#0055A4", Healthy="#018A67")

## ---- 3) 读 abundance：不筛选、去全0列、行归一化 ----
stopifnot("丰度文件不存在" = file.exists(abund_file))
ab_raw <- data.table::fread(abund_file, encoding = "UTF-8", check.names = FALSE)
stopifnot(ncol(ab_raw) >= 2)
sample_id <- as.character(ab_raw[[1]])
ab_raw[[1]] <- NULL
A <- as.matrix(ab_raw); rownames(A) <- sample_id
A <- A[, colSums(A, na.rm = TRUE) > 0, drop = FALSE]
A <- vegan::decostand(A, method = "total")  # 行和=1

# 附件：清洗后的丰度
readr::write_csv(
  data.frame(SampleID = rownames(A), A, check.names = FALSE),
  file.path(out_dir, "F3_input_family_relative_rowNorm.csv")
)

## ---- 4) 读 mapping & 标准化（ID / Group / 混杂）----
stopifnot("mapping文件不存在" = file.exists(map_file))
map0 <- data.table::fread(map_file, na.strings = c("", "NA")) |> as.data.frame()

# 标准化 ID 列
idcol <- names(map0)[tolower(names(map0)) == "id"][1]
if (!is.na(idcol) && idcol != "ID") names(map0)[names(map0) == idcol] <- "ID"
map0$ID <- trimws(as.character(map0$ID))

# 生成 group_label（4组）
stopifnot("mappingV1.csv 缺少 group_code（1=SO,2=Sarcopenia,3=Obesity,4=Healthy）" = "group_code" %in% names(map0))
map0 <- map0 |>
  mutate(
    group_code  = trimws(as.character(group_code)),
    Group = case_when(
      group_code == "1" ~ "SO",
      group_code == "2" ~ "Sarcopenia",
      group_code == "3" ~ "Obesity",
      group_code == "4" ~ "Healthy",
      TRUE ~ NA_character_
    ),
    Group = factor(Group, levels = c("SO","Sarcopenia","Obesity","Healthy"))
  )

#【# ---- 5) 自动挑混杂 + 构造设计矩阵（存在且有信息量才入模）----】
build_design <- function(map_df, group_var = "Group",
                         cand_num = c("age","bmi","mmse","phq9","SPPB","grip_strength"),
                         cand_fac = c("sex","hearing_impaired","teeth_loss",
                                      "constipation","yogurt","eatprob","diarrhea","laxatives")) {
  md  <- map_df
  low <- tolower(names(md))
  pick <- function(v) { idx <- which(low == tolower(v))[1]; if (length(idx)) idx else NA_integer_ }
  
  # 性别编码统一
  for (nm in c("sex","gender")) {
    idx <- pick(nm); if (is.na(idx)) next
    v <- md[[idx]]
    if (is.numeric(v)) v <- factor(ifelse(v==1,"Male", ifelse(v==2,"Female", as.character(v))))
    else               v <- factor(as.character(v))
    md[[idx]] <- v
  }
  
  keep_num <- c(); keep_fac <- c()
  for (nm in cand_num) {
    idx <- pick(nm); if (is.na(idx)) next
    v <- md[[idx]]
    if (is.numeric(v) && sum(is.finite(v)) >= 8 && sd(v, na.rm = TRUE) > 0)
      keep_num <- c(keep_num, names(md)[idx])
  }
  for (nm in cand_fac) {
    idx <- pick(nm); if (is.na(idx)) next
    v <- as.factor(md[[idx]])
    if (nlevels(droplevels(v)) > 1) keep_fac <- c(keep_fac, names(md)[idx])
  }
  
  stopifnot(group_var %in% names(md))
  grp <- droplevels(as.factor(md[[group_var]]))
  
  des <- data.frame(Sample = md$ID, Group = grp, stringsAsFactors = FALSE)
  for (nm in keep_num) des[[nm]] <- md[[nm]]
  for (nm in keep_fac) des[[nm]] <- droplevels(as.factor(md[[nm]]))
  
  cc <- complete.cases(des)
  list(design = des[cc, , drop = FALSE],
       kept_num = keep_num, kept_fac = keep_fac,
       n_dropped = sum(!cc))
}

des_list <- build_design(map0, group_var = "Group")
design   <- des_list$design
stopifnot("有效样本或分组不足" = (nrow(design) >= 6 && nlevels(design$Group) >= 2))

## ---- 6) 对齐样本、Bray 距离、PCoA（用于 p1/p2/p3 可视化）----
common <- intersect(rownames(A), design$Sample)
A1  <- A[common, , drop = FALSE]
des <- design[match(common, design$Sample), , drop = FALSE]
rownames(des) <- des$Sample
des$Sample <- NULL

bray <- vegan::vegdist(A1, method = "bray")
pcoa <- ape::pcoa(bray)

get_rel <- function(pcoa_vals) {
  if (!is.null(pcoa_vals$Relative_eig)) return(pcoa_vals$Relative_eig)
  if (!is.null(pcoa_vals$Rel_corr_eig)) return(pcoa_vals$Rel_corr_eig)
  stop("Cannot find relative eigenvalues in pcoa$values.")
}
rel_eig <- get_rel(pcoa$values)

pcoadata <- data.frame(
  Sample = rownames(pcoa$vectors),
  PC1 = pcoa$vectors[, 1],
  PC2 = pcoa$vectors[, 2],
  Group = des$Group[match(rownames(pcoa$vectors), rownames(des))]
)
pcoadata$Group <- droplevels(pcoadata$Group)

# 附件：PCoA坐标 & 解释度 & 组样本量
readr::write_csv(pcoadata, file.path(out_dir, "F3_pcoa_coords_with_group.csv"))
eig_df <- data.frame(
  Axis = paste0("PCoA", seq_along(rel_eig)),
  Relative = as.numeric(rel_eig),
  Cumulative = cumsum(as.numeric(rel_eig))
)
readr::write_csv(eig_df, file.path(out_dir, "F3_pcoa_relative_eigenvalues.csv"))
group_n <- as.data.frame(table(pcoadata$Group)); names(group_n) <- c("Group","n")
readr::write_csv(group_n, file.path(out_dir, "F3_group_sizes.csv"))

## ---- 7) PERMANOVA（仅 Bray；by="margin"：Group 的部分效应）+ PERMDISP ----
rhs <- paste(colnames(des), collapse = " + ")   # Group + 混杂
fml <- as.formula(paste("bray ~", rhs))

set.seed(1)
ad2 <- vegan::adonis2(fml, data = des, permutations = 9999, by = "margin")
adonis_tab <- as.data.frame(ad2)
readr::write_csv(adonis_tab, file.path(out_dir, "F3_PERMANOVA_adonis2_table_adjusted.csv"))
readr::write_csv(des,       file.path(out_dir, "F3_design_used_PERMANOVA.csv"))

# Group 行
rn <- rownames(adonis_tab)
term_row <- if ("Group" %in% rn) "Group" else setdiff(rn, c("Residual","Total"))[1]
df_eff <- adonis_tab[term_row, "Df"]
r2_eff <- adonis_tab[term_row, "R2"]
p_eff  <- adonis_tab[term_row, "Pr(>F)"]

# 方差齐性（PERMDISP）
bd      <- vegan::betadisper(bray, des$Group)
bd_perm <- vegan::permutest(bd, permutations = 9999)
permdisp_p <- bd_perm$tab[1,"Pr(>F)"]
# 保存 PERMDISP
bd_tab <- as.data.frame(bd_perm$tab); bd_tab$term <- rownames(bd_tab)
readr::write_csv(bd_tab, file.path(out_dir, "F3_PERMDISP_table.csv"))

## ---- 8) Tukey 字母（PC1/PC2 ~ Group；用于 p2/p3 面板）----
pad1 <- 0.08 * diff(range(pcoadata$PC1, na.rm = TRUE))
pad2 <- 0.08 * diff(range(pcoadata$PC2, na.rm = TRUE))
df1  <- aggregate(PC1 ~ Group, data = pcoadata, max); df1$y <- df1$PC1 + ifelse(is.finite(pad1), pad1, 0.1)
df2  <- aggregate(PC2 ~ Group, data = pcoadata, max); df2$y <- df2$PC2 + ifelse(is.finite(pad2), pad2, 0.1)

tuk1 <- aov(PC1 ~ Group, data = pcoadata) |>
  multcomp::glht(linfct = mcp(Group = "Tukey")) |>
  multcomp::cld()
tuk2 <- aov(PC2 ~ Group, data = pcoadata) |>
  multcomp::glht(linfct = mcp(Group = "Tukey")) |>
  multcomp::cld()

letters_df <- data.frame(
  Group = names(tuk1$mcletters$Letters),
  PC1_letter = unname(tuk1$mcletters$Letters)
) |>
  dplyr::left_join(
    data.frame(Group = names(tuk2$mcletters$Letters),
               PC2_letter = unname(tuk2$mcletters$Letters)),
    by = "Group"
  ) |>
  dplyr::left_join(df1[,c("Group","y")], by = "Group") |>
  dplyr::rename(PC1_y = y) |>
  dplyr::left_join(df2[,c("Group","y")], by = "Group") |>
  dplyr::rename(PC2_y = y)

readr::write_csv(letters_df, file.path(out_dir, "F3_Tukey_letters_PC1_PC2.csv"))

## ---- 9) p1–p4 组图（与参考布局一致）----
pal_use <- pal_group[levels(pcoadata$Group)]

# p1: PCoA散点 + 95% 椭圆 + 轴线
p1 <- ggplot(pcoadata, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 1) +
  labs(
    x = paste0("PCoA1 (", floor(rel_eig[1] * 100), "%)"),
    y = paste0("PCoA2 (", floor(rel_eig[2] * 100), "%)")
  ) +
  scale_colour_manual(values = pal_use, drop = FALSE) +
  stat_ellipse(aes(fill = Group), geom = "polygon",
               type = "t", level = 0.95, alpha = 0.25,
               color = NA, show.legend = FALSE) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 0.6, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  coord_fixed() +
  theme(
    text = element_text(size = 15),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_text(colour = 'black', size = 15),
    axis.title.y = element_text(colour = 'black', size = 15),
    axis.text     = element_text(colour = 'black', size = 12),
    legend.title  = element_blank(),
    legend.position   = 'none',
    panel.grid    = element_blank()
  ) +
  scale_fill_manual(values = pal_use, drop = FALSE)

# p2: PC1 箱线 + Tukey 字母（横置）
p2 <- ggplot(pcoadata, aes(Group, PC1)) +
  geom_boxplot(aes(fill = Group), outlier.colour = NA, alpha = 0.1) +
  scale_fill_manual(values = pal_use, drop = FALSE) +
  geom_text(data = letters_df, aes(x = Group, y = PC1_y, label = PC1_letter),
            size = 5, color = "black") +
  theme(
    panel.background = element_rect(fill='white', colour='black'),
    axis.ticks.length = grid::unit(0.4,"lines"),
    axis.ticks = element_line(color = 'black'),
    axis.line  = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(colour = 'black', size = 12),
    axis.text.x  = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  coord_flip()

# p3: PC2 箱线 + Tukey 字母
p3 <- ggplot(pcoadata, aes(Group, PC2)) +
  geom_boxplot(aes(fill = Group), outlier.colour = NA, alpha = 0.1) +
  scale_fill_manual(values = pal_use, drop = FALSE) +
  geom_text(data = letters_df, aes(x = Group, y = PC2_y, label = PC2_letter),
            size = 5, color = "black") +
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.ticks.length = grid::unit(0.4,"lines"),
    axis.ticks = element_line(color='black'),
    axis.line  = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(colour = 'black', size = 12, angle = 90,
                                vjust = 0.5, hjust = 0.5),
    axis.text.y  = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  )

# p4: PERMANOVA 注释（校正混杂后的 Group 部分效应）+ PERMDISP
xpos <- min(pcoadata$PC1, na.rm = TRUE) + 0.05 * diff(range(pcoadata$PC1, na.rm = TRUE))
ypos <- max(pcoadata$PC2, na.rm = TRUE) - 0.05 * diff(range(pcoadata$PC2, na.rm = TRUE))
txt  <- sprintf("PERMANOVA (Bray, adjusted)\ndf = %d\nR2 = %.4f\np = %.4g\nPERMDISP p = %.4g",
                as.integer(df_eff), r2_eff, p_eff, permdisp_p)

p4 <- ggplot(pcoadata, aes(PC1, PC2)) +
  annotate("text", x = xpos, y = ypos, label = txt,
           size = 3.2, family = "sans", fontface = 1, hjust = 0) +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.title  = element_blank(),
        axis.line   = element_blank(),
        axis.ticks  = element_blank(),
        axis.text   = element_blank())

# 组合：p2 + p4 | p1 + p3
p <- p2 + p4 + p1 + p3 +
  patchwork::plot_layout(heights = c(1,4), widths = c(4,1), ncol = 2, nrow = 2)
print(p)

ggsave(out_png, plot = p, dpi = 400, width = 6.2, height = 5.0)
ggsave(out_pdf, plot = p, dpi = 400, width = 6.2, height = 5.0)

