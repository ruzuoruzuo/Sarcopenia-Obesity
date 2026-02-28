# Sarcopenia-Obesity

本仓库用于复现/生成“肌少症-肥胖（Sarcopenia–Obesity）”相关的肠道菌群分析图表与机器学习结果，包含：
- Fig1：Enterotype（PAM 聚类 + PCoA）/ Top Genus（伪对数箱线图）/ Enterotype 在 4 组中的比例
- Fig2：Alpha 多样性（Shannon/Simpson/Pielou 等）+ 组间配对检验 + GLM 森林图（图 g 风格）
- Fig3：Beta 多样性（Bray-Curtis PCoA）+ 混杂因素校正的 PERMANOVA + 四联图布局
- Fig5：Top30 特征热图（组均值 Z-score）+ 各组 vs Healthy 显著性星号
- Fig6：XGBoost 四分类（SO/S/O/H）5-fold 分层交叉验证 + OOF 预测 + ROC/AUC（含 95%CI）+ Feature importance + SHAP

> 说明：当前脚本里使用了绝对路径（Windows）。使用前请先把每个脚本顶部的 `setwd()` / `workdir` / `infile` 改成你本机的路径，或改为相对路径（推荐）。

---

## 1. 环境要求

- R >= 4.1（建议 4.2/4.3）
- Windows/Mac/Linux 均可（字体部分仅在 Windows 下自动加载）
- 主要 R 包（脚本中已写明；部分脚本会自动 install，但建议手动安装更稳）：
  - 数据处理/绘图：`data.table`, `dplyr`, `tidyr`, `readr`, `stringr`, `ggplot2`, `scales`, `patchwork`
  - 多样性/距离：`vegan`, `ape`
  - 统计：`broom`, `multcomp`
  - 热图：`ComplexHeatmap`, `circlize`, `grid`
  - 机器学习：`xgboost`（脚本注明 xgboost 3.1.3.1）, `pROC`

---

## 2. 数据文件与格式

### 2.1 核心输入文件

#### `trans.genus.relative.csv`
- **用途**：Fig1 / Fig2 / Fig3
- **格式**：行=样本，列=属（Genus）相对丰度
- **第一列**：SampleID（样本编号）
- 其余列：各属的相对丰度（0~1 或比例）

#### `mappingV1.csv`
- **用途**：Fig1（Fig1c）/ Fig2 / Fig3
- **至少需要两列**：
  - `ID`：与 `trans.genus.relative.csv` 的 SampleID 对应
  - `group_code`：分组编码（脚本要求：`1=SO, 2=Sarcopenia, 3=Obesity, 4=Healthy`）
- **Fig3（可选混杂）**：若存在且有信息量，会自动纳入 PERMANOVA 校正候选变量：
  - 数值型候选：`age, bmi, mmse, phq9, SPPB, grip_strength`
  - 分类型候选：`sex, hearing_impaired, teeth_loss, constipation, yogurt, eatprob, diarrhea, laxatives`
  - 性别会尝试把 1/2 统一为 Male/Female

### 2.2 Fig2（图 g）可选输入
用于把属级丰度汇总到门（Phylum）或直接使用门水平丰度：
- `trans.phylum.relative.csv`（门水平相对丰度；第一列 SampleID）
或
- `genus2phylum.csv`（两列：`genus, phylum`）用于将 genus 聚合到 phylum

### 2.3 Fig5 输入
#### `abund.csv`
- **用途**：Fig5 热图
- **格式**：至少包含 `SampleID` 列，其余列为特征/代谢物/菌等数值矩阵
- **分组识别**：脚本用 `SampleID` 前缀自动识别组别：
  - `SO...` -> SO
  - `S...`  -> Sarcopenia
  - `O...`  -> Obesity
  - `H...`  -> Healthy  
  （注意：脚本已特意先判 SO，避免被 S 吃掉）

### 2.4 Fig6 输入
#### `V1abund_plus_species_merged_by_SampleID.csv`
- **用途**：Fig6 XGBoost
- **格式**：必须有 `SampleID` 列；其余列为数值特征（会强制转 numeric，任何 NA 会报错）
- **分组识别**：从 `SampleID` 提取类别（需要形如 `SO-xx / S-xx / O-xx / H-xx` 或用下划线分隔亦可）

---

## 3. 脚本说明与输出

### `Fig1.abcsup.R`
生成：
- **Fig1a**：Enterotype（Bray + PAM 聚类，silhouette 自动选择 k）+ PCoA 图  
  输出：`Enterotype_PCoA.png`
- **Fig1b**：Top Genus（按 Enterotype 分组箱线图，伪对数尺度）  
  输出：`Enterotype_TopGenus_Boxplot_pairwise_pseudolog.png`
- **Fig1c**：Enterotype 在 4 组（SO/Sarcopenia/Obesity/Healthy）的比例图  
  输出：`Enterotype_proportion_by_SO.png`、`Enterotype_proportion_by_SO.svg`
- 同时导出多份中间结果 CSV（前缀 `S1_` `S3_` `S6_` `S7_` `S8_` 等），便于附表/附图复核。

### `Fig2.alpha.R`
- 计算 α 多样性指标（Shannon(H, ln)、Gini–Simpson(1−D)、Pielou 等）
- 组间比较：Wilcoxon + BH 校正（输出汇总表）
- 额外：**图 g 风格 GLM 森林图**（α 指标 ~ top taxa abundance）
输出（主要）：
- `alpha_diversity_by_group.filtered.png`
- `alpha_glm_forest_H_J_1minusD.png`
- `alpha_diversity.filtered.csv`
- `alpha_by_group_used_samples.filtered.csv`
- `alpha_by_group_summary.filtered.csv`
- `alpha_pairwise_all_wilcox_BH.csv`
- `Supp_Table_alpha_glm_coefficients.csv`

### `Fig3.PCoA.R`
- Bray-Curtis 距离 + PCoA
- **PERMANOVA（adonis2）校正混杂因素**（Group 的部分效应，by="margin"）
- PERMDISP（组内离散度差异）
- p1–p4 四联图：PCoA散点 + PC1/PC2 箱线 + 统计注释
输出：
- `Fig3_PCoA_SO_adjBray.png`
- `Fig3_PCoA_SO_adjBray.pdf`
- 以及目录 `F3_confounded/` 下多份 CSV（PCoA 坐标、解释度、PERMANOVA 表、PERMDISP 表、Tukey 字母等）

### `Fig5.heatmap.R`
- 从 `abund.csv` 读入矩阵
- log10 转换（`log10(X + 1e-12)`）
- 以 **Kruskal-Wallis + FDR** 选 Top30 特征
- 计算各组均值并做行 Z-score（按 SO 从高到低排序）
- 各组 vs Healthy：Wilcoxon + FDR<0.05 标 “*”
输出：
- `signature_heatmap_top30_SO_to_Healthy.pdf`

### `Fig6.XGBOOSTSHAP.R`
- 四分类（SO/S/O/H）**分层 5 折交叉验证**
- 输出 OOF（out-of-fold）预测、混淆矩阵（计数/行归一）、AUC（one-vs-rest）、以及 SO vs (S/O/H) 的 pairwise ROC + 95%CI
- 最后用全量数据训练 final model，输出：
  - Feature importance（Gain 等）
  - SHAP mean(|SHAP|) Top20 图
输出（主要位于 `XGB_totalCV5plus/`）：
- `cv5_fold_metrics.csv`, `cv5_summary.csv`
- `oof_predictions_by_sample.csv`
- `confusion_matrix_oof_counts.pdf`
- `confusion_matrix_oof_row_normalized.pdf`
- `oof_auc_one_vs_rest.csv`
- `pairwise_auc_OOF_SO_vs_O_S_H.csv` 及各 pairwise ROC 曲线点与 CI band CSV
- `xgb_model_final_trained_on_all_cv5.rds`
- `final_model_xgb_importance_all.csv`
- `final_model_xgb_importance_top20_byGain.pdf`
- `final_model_shap_meanAbs_top20.pdf`

---

## 4. 如何运行（推荐流程）

### Step 0：放好数据
建议目录结构（仅建议，你也可以按自己习惯）：
