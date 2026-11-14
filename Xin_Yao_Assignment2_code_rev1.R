## Project 2 ###
## Yao Xin ###

# 0.prepare packages and dirs ------
library(GEOquery)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(WGCNA) 
library(RColorBrewer)

dir.create("outputs", showWarnings = FALSE)




# gds  <- getGEO("GDS2901", destdir = ".")   # downloads and caches GDS2901.soft.gz
# eset <- GDS2eSet(gds, do.log2 = FALSE)     # ExpressionSet: exprs(eset) is probes × samples
# pData(eset)                                # genotype / agent / time per sample
# data<-as.matrix(eset)


# part II - reproduce on human dataset GSE10072

## ==========1.input data, examine ============
gse  <- getGEO("GSE10072", destdir = ".")   # downloads and caches GDS2901.soft.gz
names(gse)           # e.g. "GSE10072_series_matrix.txt.gz"

eset <- gse[[1]]     # pull out the ExpressionSet
# eset                 # assayData: 22283 features, 107 samples 
# pData(eset) 
pd <- pData(eset)


# processing info: 
unique(pd$data_processing)
# "Gene expression data were processed and normalized using Bioconductor Affy package, 
# based on the Robust Multichip Average (RMA) method5 for single-channel Affymetrix chips"

experimentData(eset)              # MIAME (Minimum Information About a Microarray Experiment) object
abstract(eset)                    # GEO "Summary" text
experimentData(eset)@title
experimentData(eset)@other$overall_design
pubMedIds(eset)                   # linked paper
annotation(eset)                  # platform ID, e.g. "GPL96"

gse_soft <- getGEO("GSE10072", GSEMatrix = FALSE)
Meta(gse_soft)$title
Meta(gse_soft)$summary
Meta(gse_soft)$overall_design
Meta(gse_soft)$pubmed_id

# per-sample descriptions (tissue, smoking status, etc.)
View(pData(eset))


head(pd$title, 10)     # patient IDs

# Data inspect
colnames(pData(eset))
n_unique <- sapply(pd, function(v) length(unique(v))) 
n_unique[n_unique > 1 & n_unique < 10]    # find grouping columns
col_groups <- names(n_unique)[n_unique > 1 & n_unique < 10] 

# quick tabulate 
for (v in col_groups) {
  cat("\n**", v, "**\n")
  print(table(pd[[v]], useNA = "ifany"))
}

## =========2. sample summaries =============== 
library(gtsummary)

### ---- 2.1 Clean sample table (one row per array) ----
ch_table <- data.frame(
  patient = sub(".*_", "", pd$title),             # GTxxxxx
  sample  = pd$geo_accession,                     # GSMxxxxxx
  tissue  = pd$source_name_ch1, # Adenocarcinoma of the Lung vs normal 
  sex     = pd$`Gender:ch1`,
  smoking = pd$`Cigarette Smoking Status:ch1`,
  stage   = pd$`Stage:ch1`,
  age     = as.numeric(pd$`Age at Diagnosis:ch1`),
  row.names = pd$geo_accession,
  stringsAsFactors = FALSE
)

ch_table$tissue  <- factor(ifelse(grepl("Adenocarcinoma", ch_table$tissue), "tumor", "normal"),
                           levels = c("normal", "tumor"))
ch_table$smoking <- factor(ch_table$smoking,
                           levels = c("Never Smoked", "Former Smoker", "Current Smoker"),
                           labels = c("never", "former", "current"))
ch_table$sex     <- factor(ch_table$sex)
ch_table$stage   <- factor(ch_table$stage)
 
### ---- 2.2 Pairing: which patients contributed tumour, normal, or both ----
pair_tab <- table(ch_table$patient, ch_table$tissue)
has_t <- pair_tab[, "tumor"]  == 1
has_n <- pair_tab[, "normal"] == 1
pt_group <- ifelse(has_t & has_n, "Tumour + normal",
                   ifelse(has_t, "Tumour only", "Normal only"))
names(pt_group) <- rownames(pair_tab)

ch_table$group  <- factor(pt_group[ch_table$patient],
                          levels = c("Tumour + normal", "Tumour only", "Normal only"))
ch_table$paired <- ch_table$group == "Tumour + normal"
 
### ---- 2.3 Patient table (one row per patient) ----
pt_table <- ch_table[!duplicated(ch_table$patient),
                     c("patient", "group", "sex", "smoking", "stage", "age")]
rownames(pt_table) <- pt_table$patient

### ---- 2.4 Table 1: patients (each counted once) ----
lab <- list(sex ~ "Sex", smoking ~ "Smoking status",
            stage ~ "Stage", age ~ "Age at diagnosis (years)")

tbl_patients <- pt_table[, c("group", "sex", "smoking", "stage", "age")] |>
  tbl_summary(
    by        = group,
    type      = list(age ~ "continuous"),
    statistic = list(age ~ "{mean} ({sd}); {min}–{max}"),
    digits    = list(age ~ 1),
    label     = lab,
    missing   = "ifany"
  ) |>
  add_overall(col_label = "**All patients**  N = {N}") |>
  modify_header(all_stat_cols(stat_0 = FALSE) ~ "**{level}**  N = {n}") |>
  modify_caption("**Table 1. GSE10072 patient characteristics**") |>
  bold_labels()

tbl_patients

### ---- 2.5 Table S1: arrays by tissue (patients in 'Tumour + normal' appear in both columns) ----
tbl_samples <- ch_table[, c("tissue", "sex", "smoking", "stage", "age")] |>
  tbl_summary(
    by        = tissue,
    type      = list(age ~ "continuous"),
    statistic = list(age ~ "{mean} ({sd}); {min}–{max}"),
    digits    = list(age ~ 1),
    label     = lab,
    missing   = "ifany"
  ) |>
  modify_header(all_stat_cols() ~ "**{level} samples**  N = {n}") |>
  modify_caption("**Table S1. GSE10072 arrays by tissue (33 patients contribute to both columns)**") |>
  bold_labels()

tbl_samples

### ---- 2.6 Save: markdown (readable), rds (re-usable), csv (data) ----
foot <- c("", "n (%); Mean (SD); Min–Max")
writeLines(c(as_kable(tbl_patients, format = "pipe"), foot), "outputs/GSE10072_table1_patients.md")
writeLines(c(as_kable(tbl_samples,  format = "pipe"), foot), "outputs/GSE10072_tableS1_samples.md")
saveRDS(list(patients = tbl_patients, samples = tbl_samples), "outputs/GSE10072_tables.rds")
saveRDS(ch_table, "outputs/GSE10072_ch_table.rds")
write.csv(ch_table, "outputs/GSE10072_sample_metadata.csv", row.names = TRUE)


 
## ========= 3. QC, Expression matrices for DiffCoEx ==============

eset_gse <- gse[[1]]  
N_GENES  <- 8000   # limit    

### ---- 3.1 Expression: already RMA (log2 + quantile-normalized) per data_processing ----
x <- exprs(eset_gse)
dim(x)
identical(colnames(x), rownames(ch_table)) 

### ---- 3.2 Drop Affymetrix control probes ----
x <- x[!grepl("^AFFX", rownames(x)), ]


### ---- 3.3 Label-free gene filter: expressed (mean > cut-off), then most variable, across ALL arrays ----

#### 3.3.0 Visual check of probe means ----
m <- rowMeans(x)                                   # one mean per probe set, across all 107 arrays

png("outputs/GSE10072_probe_mean_density.png", width = 900, height = 500)
plot(density(m), main = "GSE10072 probe means (log2 RMA)", xlab = "mean expression (log2)")
abline(v = EXPR_CUT, col = "blue", lty = 2)
legend("topright", paste0("cut-off: log2 = ", EXPR_CUT, " (", sum(m > EXPR_CUT), " probes kept)"),
       col = "blue", lty = 2, bty = "n")
dev.off()

# print cutoffs 
print(sapply(c(5, 5.5, 6, 6.5, 7), function(cut) c(cutoff = cut, n_probes = sum(m > cut))))
# pick: 
EXPR_CUT <- 6          # end of the low-expression (background) shoulder in the density plot

#### 3.3.1 Apply the filter ----
expressed <- which(m > EXPR_CUT)                                   # 17,631 probe sets
x_expr    <- x[expressed, ]

# withIN-tissue variance: 
v_tum    <- rowVars(x_expr[, ch_table$tissue == "tumor"])
v_nor    <- rowVars(x_expr[, ch_table$tissue == "normal"])
v_within <- (v_tum + v_nor) / 2

keep <- order(v_within, decreasing = TRUE)[1:N_GENES]              # most variable WITHIN tissues
x_f  <- x_expr[keep, ]

stopifnot(nrow(x_f) == N_GENES, ncol(x_f) == ncol(x),
          identical(colnames(x_expr), rownames(ch_table)))         # tissue labels line up with columns

dim(x_expr); dim(x_f)



### ------ 3.4 Probe annotation (first symbol / Entrez when a probe maps to several) ----
fd  <- fData(eset_gse)[rownames(x_f), ]
ann <- data.frame(probe  = rownames(x_f),
                  symbol = sub(" ///.*", "", fd$`Gene Symbol`),
                  entrez = sub(" ///.*", "", fd$ENTREZ_GENE_ID),
                  row.names = rownames(x_f), stringsAsFactors = FALSE)
ann$symbol[ann$symbol == ""] <- NA
ann$entrez[ann$entrez == ""] <- NA

### ---- 3.5 Condition matrices: samples x genes (S1 orientation) ----
#### 3.5.1 Design "all": every array (closest to the paper: two groups compared as wholes) ----
tum_all <- rownames(ch_table)[ch_table$tissue == "tumor"]
nor_all <- rownames(ch_table)[ch_table$tissue == "normal"]

#### 3.5.2 Design "paired": 33 matched patients, rows aligned by patient ----
pairs   <- ch_table[ch_table$paired, ]
pairs   <- pairs[order(pairs$patient, pairs$tissue), ]
tum_pr  <- rownames(pairs)[pairs$tissue == "tumor"]
nor_pr  <- rownames(pairs)[pairs$tissue == "normal"]
identical(ch_table[tum_pr, "patient"], ch_table[nor_pr, "patient"])

#### 3.5.3 Design "independent": each patient contributes ONE array (mirrors the paper's design) -----
SPLIT_SEED <- 2026
set.seed(SPLIT_SEED)
pr_ids   <- unique(pairs$patient)
give_tum <- sample(pr_ids, ceiling(length(pr_ids) / 2))      # 17 paired subj give tumour, 16 give normal
ind <- ch_table[ (ch_table$group == "Tumour only") |
                   (ch_table$group == "Normal only") |
                   (ch_table$paired & ch_table$tissue == "tumor"  & ch_table$patient %in% give_tum) |
                   (ch_table$paired & ch_table$tissue == "normal" & !(ch_table$patient %in% give_tum)), ]
anyDuplicated(ind$patient)

tum_ind <- rownames(ind)[ind$tissue == "tumor"]
nor_ind <- rownames(ind)[ind$tissue == "normal"]


dc_input <- list(
  independent = list(datC1 = t(x_f[, tum_ind]), 
                     datC2 = t(x_f[, nor_ind]), 
                     seed = SPLIT_SEED),
  all    = list(datC1 = t(x_f[, tum_all]), datC2 = t(x_f[, nor_all])),
  paired = list(datC1 = t(x_f[, tum_pr]),  datC2 = t(x_f[, nor_pr]),
                patient = ch_table[tum_pr, "patient"]),
  ann    = ann,
  settings = list(N_GENES = N_GENES, EXPR_CUT = EXPR_CUT,
         filter = "mean log2 RMA > EXPR_CUT (all arrays), then top N_GENES by mean within-tissue variance")
  )


### ---- 3.6 Save + log ----
saveRDS(dc_input, "outputs/GSE10072_diffcoex_input.rds")

capture.output(
  cat("Probes after AFFX removal:", nrow(x), "\nGenes kept:", nrow(x_f), "\n\n"),
  
  cat("all   : tumour", dim(dc_input$all$datC1),    " normal", dim(dc_input$all$datC2), "\n"),
  cat("paired: tumour", dim(dc_input$paired$datC1), " normal", dim(dc_input$paired$datC2), "\n"),
  cat("paired rows aligned by patient:",
      identical(rownames(dc_input$paired$datC1), tum_pr), "\n\n"),
  cat("independent: tumour", dim(dc_input$independent$datC1),
      " normal", dim(dc_input$independent$datC2), "\n"),
  
  cat("Expression range (kept genes):", round(range(x_f), 2), "\n"),
  cat("Mean expression (kept genes): min", round(min(rowMeans(x_f)), 2),
      " median", round(median(rowMeans(x_f)), 2), "\n"),
  cat("Probes with mean >", EXPR_CUT, ":", length(expressed), "\n"),
  cat("Annotated with symbol:", sum(!is.na(ann$symbol)), "of", nrow(ann), "\n"),
  cat("Unique symbols:", length(unique(na.omit(ann$symbol))), "\n"), 
  file = "outputs/GSE10072_sec3_log.txt"
)

cols <- c(normal = "#3B5BA5", tumor = "#E8A33D")   # define colors
png("outputs/GSE10072_sec3_density.png", width = 900, height = 550, res = 120)
par(mar = c(4, 4, 3, 1))
dens <- apply(x_f, 2, density)
plot(NA, xlim = range(sapply(dens, function(d) range(d$x))),
     ylim = c(0, max(sapply(dens, function(d) max(d$y)))),
     xlab = "log2 RMA", ylab = "Density", main = "Expression distribution per array (kept genes)")
for (j in seq_along(dens))  lines(dens[[j]], col = adjustcolor(cols[as.character(ch_table$tissue[j])], alpha.f = 0.4), lwd = 1)
legend("topright", c("Normal", "Tumour"), col = cols, lwd = 2, bty = "n")
dev.off()
 
### ---- 3.7 QC figures (adapted from class code 01 and 03) ----
theme_set(theme_bw(base_size = 13))                                   # 01
PCA_N_GENES <- 2000L                                                  # 01: selected by variability, NOT by labels
cols <- c(normal = "#3B5BA5", tumor = "#E8A33D")      # one tissue palette for every figure
tissue_col <- cols[as.character(ch_table$tissue)]     # one colour per array, in column order of x

#### 3.7.1 Before/after-filter boxplots ----
# 01: "03_normalization_boxplots.pdf" (raw vs size-factor normalized).
# Microarray data are already RMA-normalized, so the two panels are all probes vs kept genes.
pdf(file.path(out_dir, "GSE10072_3.7.1_boxplots.pdf"), width = 11, height = 5)
par(mfrow = c(1, 2), mar = c(3, 4, 3, 1))                            # bottom margin 8 -> 3 (no labels)
boxplot(x, las = 2, outline = FALSE, xaxt = "n", col = tissue_col, border = "grey30",
        main = "All probes (RMA)", ylab = "log2 RMA expression")
mtext(paste(ncol(x), "arrays"), side = 1, line = 1)
boxplot(x_f, las = 2, outline = FALSE, xaxt = "n", col = tissue_col, border = "grey30",
        main = "Kept genes (filtered)", ylab = "log2 RMA expression")
mtext(paste(ncol(x_f), "arrays"), side = 1, line = 1)
legend("topright", names(cols), fill = cols, bty = "n", cex = 0.8)
dev.off()

#### 3.7.2 Sample distances ----
# 01: "05_sample_distances.pdf". Expression is already log2 RMA, so no VST step.
ann_colors <- list(tissue  = cols,
                   smoking = c(never = "grey85", former = "grey55", current = "grey20"))
pheatmap(sample_dist, annotation_col = annotation, annotation_row = annotation,
         annotation_colors = ann_colors,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Sample distances (log2 RMA)", width = 8, height = 7,
         filename = file.path(out_dir, "GSE10072_3.7.2_sample_distances.pdf"))

#### 3.7.3 PCA ----
# 01: "06_PCA.png". scale.=FALSE preserves variance differences, as in class.
variable_genes <- head(order(rowVars(x), decreasing = TRUE), PCA_N_GENES)
pca <- prcomp(t(x[variable_genes, , drop = FALSE]),
              center = TRUE, scale. = FALSE)
variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.frame(pca$x[, 1:2], ch_table, sample = rownames(ch_table))

pca_df$stage_grp <- factor(substr(as.character(pca_df$stage), 1, nchar(gsub("[AB]$", "", pca_df$stage))),
                           levels = c("I", "II", "III", "IV"))
levels(pca_df$stage_grp) <- c("I", "II", "III–IV", "III–IV")   # merge III and IV (only 3 patients with IV)
table(pca_df$stage_grp)

p <- ggplot(pca_df, aes(PC1, PC2, fill = tissue, shape = stage_grp)) +
  geom_point(size = 2.5, color = "white", stroke = 0.4) +     # thin white border
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = c(21, 24, 22)) +
  guides(shape = guide_legend(override.aes = list(fill = "grey40", color = "grey40")),
         fill  = guide_legend(override.aes = list(shape = 21))) +   # show fill colours in the legend
  labs(title = "PCA: inspect tissue and stage effects", shape = "Stage", fill = "Tissue",
       x = paste0("PC1: ", round(variance_explained[1], 1), "%"),
       y = paste0("PC2: ", round(variance_explained[2], 1), "%"))

ggsave(file.path(out_dir, "GSE10072_3.7.3_PCA.png"), p, width = 8.5, height = 5.5)
write.csv(pca_df, file.path(out_dir, "GSE10072_PCA_scores.csv"), row.names = FALSE)
# PCA clustering is descriptive. Do not remove a sample solely because it differs.   # 01

#### 3.7.4 Network QC: good samples/genes and sample clustering, per condition ----
# 03: goodSamplesGenes + "03_sample_clustering.pdf". DiffCoEx builds correlations within each
# condition, so the check is done separately for the tumour and normal arrays of the design used.
datC1 <- as.data.frame(dc_input$independent$datC1)                   # tumour, samples x genes
datC2 <- as.data.frame(dc_input$independent$datC2)                   # normal, samples x genes

gsg1 <- goodSamplesGenes(datC1, verbose = 3)
gsg2 <- goodSamplesGenes(datC2, verbose = 3)
stopifnot(gsg1$allOK, gsg2$allOK)

pdf(file.path(out_dir, "GSE10072_3.7.4_sample_clustering.pdf"), width = 10, height = 8)
par(mfrow = c(2, 1))
plot(hclust(dist(datC1), method = "average"), main = "Sample clustering for network QC: tumour",
     xlab = "", sub = "Inspect outliers", cex = 0.5,
     col.main = cols["tumor"])
plot(hclust(dist(datC2), method = "average"), main = "Sample clustering for network QC: normal",
     xlab = "", sub = "Inspect outliers", cex = 0.5,
     col.main = cols["normal"])
dev.off()







## ========= 4. DiffCoEx: modules from the difference in co-expression ===============
# Adapted from DiffCoEx S1 section 2 (Tesson et al. 2010), run on the independent design.
options(stringsAsFactors = FALSE)
disableWGCNAThreads()                           

### ---- 4.0 S1 supporting functions (from DiffCoEx_paper_Supp_file1)  ----
source("DiffCoEx_functions.R")   #  functions from the file ; 

### ---- 4.1 Input: independent design (each patient contributes one array) ----
datC1 <- dc_input$independent$datC1                     # tumour: samples x genes
datC2 <- dc_input$independent$datC2                     # normal: samples x genes
stopifnot(identical(colnames(datC1), colnames(datC2)))
dim(datC1); dim(datC2)

### ---- 4.2 Step 1: signed squared Spearman correlations within each condition ----
beta1 <- 6           # S1: user-defined soft threshold (paper used 6)
corC1 <- cor(datC1, method = "spearman")        # computed once, not twice
corC2 <- cor(datC2, method = "spearman")
AdjMatC1 <- sign(corC1) * corC1^2
AdjMatC2 <- sign(corC2) * corC2^2
diag(AdjMatC1) <- 0
diag(AdjMatC2) <- 0
rm(corC1, corC2); collectGarbage()

### ---- 4.3 Steps 2-3: adjacency difference -> topological-overlap dissimilarity ----
dissTOMC1C2 <- TOMdist((abs(AdjMatC1 - AdjMatC2) / 2)^(beta1 / 2))   # S1 line, unchanged
collectGarbage()

### ---- 4.4 Step 4: hierarchical clustering + dynamic tree cut ----
geneTreeC1C2 <- fastcluster::hclust(as.dist(dissTOMC1C2), method = "average")  # changed: was flashClust

dynamicModsHybridC1C2 <- cutreeDynamic(dendro = geneTreeC1C2, distM = dissTOMC1C2,
                                       method = "hybrid", cutHeight = .996,
                                       deepSplit = TRUE, pamRespectsDendro = FALSE,
                                       minClusterSize = 20)          # S1 settings
dynamicColorsHybridC1C2 <- labels2colors(dynamicModsHybridC1C2)

# Merge close modules (S1 cutHeight = 0.2). Changed: scale each condition first so that
# tumour-vs-normal mean shifts do not drive the merging.
colorh1C1C2 <- mergeCloseModules(rbind(scale(datC1), scale(datC2)),
                                 dynamicColorsHybridC1C2, cutHeight = .2)$colors
table(colorh1C1C2)

### ---- 4.5 Figures: dendrogram with modules; comparative correlation heatmap ----
png(file.path(out_dir, "GSE10072_4_module_dendrogram.png"), width = 1200, height = 700, res = 120)
plotDendroAndColors(geneTreeC1C2, cbind(dynamicColorsHybridC1C2, colorh1C1C2),
                    c("Dynamic tree cut", "Merged"), dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "GSE10072: gene dendrogram on TOM of co-expression difference")
dev.off()

if (any(colorh1C1C2 != "grey")) {
  hm <- plotC1C2Heatmap(colorh1C1C2, AdjMatC1, AdjMatC2, datC1, datC2,
                        file = file.path(out_dir, "GSE10072_4_C1C2_heatmap.png"))
} else {
  message("No differentially co-expressed modules found: all genes grey.")
}

### ---- 4.6 Save modules + log ----
modules <- data.frame(probe  = colnames(datC1),
                      symbol = dc_input$ann[colnames(datC1), "symbol"],
                      entrez = dc_input$ann[colnames(datC1), "entrez"],
                      module = colorh1C1C2)
write.csv(modules, file.path(out_dir, "GSE10072_4_module_assignment.csv"), row.names = FALSE)

saveRDS(list(colors = colorh1C1C2, colors_unmerged = dynamicColorsHybridC1C2,
             geneTree = geneTreeC1C2, beta1 = beta1,
             heatmap_order = if (exists("hm")) hm$module_order else NULL),
        file.path(out_dir, "GSE10072_4_diffcoex_modules.rds"))

capture.output(
  cat("Design: independent  tumour", dim(datC1), " normal", dim(datC2), "\n"),
  cat("beta1 =", beta1, " minClusterSize = 20  cutHeight = 0.996  merge cutHeight = 0.2\n\n"),
  cat("Modules before merging:", length(setdiff(unique(dynamicColorsHybridC1C2), "grey")), "\n"),
  cat("Modules after merging: ", length(setdiff(unique(colorh1C1C2), "grey")), "\n\n"),
  print(sort(table(colorh1C1C2), decreasing = TRUE)),
  cat("\nGenes in modules:", sum(colorh1C1C2 != "grey"), "of", length(colorh1C1C2), "\n"),
  if (exists("hm")) cat("Heatmap module order:", hm$module_order, "\n"),
  file = file.path(out_dir, "GSE10072_4_log.txt")
)



# ========= 5. Permutation test of differential co-expression (DiffCoEx S1 section 3) ===============
# Uses objects from section 4: datC1, datC2 (independent design), colorh1C1C2 (module colours).
SPLIT_SEED <- 2026
set.seed(SPLIT_SEED)                         
N_PERM <- 1000        # ** number of permutations (S1 used 1000)

# ---- 5.1 Permuted sample indices (S1) ----
permutations<-NULL
for (i in 1:N_PERM)   
{
  permutations<-rbind(permutations,sample(1:(nrow(datC1)+nrow(datC2)),nrow(datC1)))
}

# we scale the data in both conditions to mean 0 and variance 1.
d<-rbind(scale(datC1),scale(datC2))
# ** scaling within each condition removes tumour-vs-normal mean shifts before shuffling

# ---- 5.2 Observed dispersions and null distributions (S1 loop) ----
mods  <- setdiff(unique(colorh1C1C2), "grey")            # ** non-grey modules, in S1's order
n_mod <- length(mods)                                    # ** S1 hard-coded 8 (rat modules); use the actual count
stopifnot(n_mod >= 1)                                    # ** stop if section 4 found no modules

# compute all pairwise module to module dispersion values, and generate a null distribution from permuted scaled data
dispersionMatrix<-matrix(nrow=n_mod,ncol=n_mod)          # ** S1: length(unique(colorh1C1C2))-1, same value
dimnames(dispersionMatrix) <- list(mods, mods)           # ** label rows/cols by module colour
nullDistrib<-list()
i<-j<-0

big <- names(sort(table(colorh1C1C2[colorh1C1C2 != "grey"]), decreasing = TRUE))[1]   # ** largest module
system.time(apply(permutations[1:50, ], 1, permutationProcedureModule2Module, d, big, big, colorh1C1C2))  # ** 50 of the 1000 permutations

for (c1 in mods)  {                                    # ** S1: setdiff(unique(colorh1C1C2),"grey"), same set
  i<-i+1
  j<-0
  nullDistrib[[c1]]<-list()
  for (c2 in mods)  {
    j<-j+1
    if (j < i) next                                      # ** matrix is symmetric: skip lower triangle (halves run time)
    dispersionMatrix[i,j]<-dispersionModule2Module(c1,c2,datC1,datC2,colorh1C1C2)
    nullDistrib[[c1]][[c2]]<-apply(permutations,1,permutationProcedureModule2Module,d,c2,c1,colorh1C1C2)
  }
}

# ---- 5.3 Permutation summary (S1) ----
#We create a summary matrix indicating for each module to module
#differential coexpression the number of permuted data yielding
#an equal or higher dispersion.
permutationSummary<-matrix(nrow=n_mod,ncol=n_mod)        # ** S1: nrow=8,ncol=8
colnames(permutationSummary)<-mods
rownames(permutationSummary)<-mods
for (i in 1:n_mod) { for (j in i:n_mod) {               # ** S1: 1:8 x 1:8; upper triangle only, as computed above
  permutationSummary[i,j]<-length(which(nullDistrib[[i]][[j]] >= dispersionMatrix[i,j]))}}

# ---- 5.4 Supplementary Fig S1 (S1 plotMatrix) ----
#We plot the result (cf supplementary figure 1)
png(file.path(out_dir, "GSE10072_5_permutation_FigS1.png"), width = 700, height = 700)   # ** save to file (S1 drew to screen)
plotMatrix(permutationSummary)
dev.off()                                                # ** close the png device

# ---- 5.5 Tables: empirical p-values with BH adjustment ----
perm_tab <- data.frame(module1 = rownames(permutationSummary)[row(permutationSummary)],     # ** long format: one row per module pair
                       module2 = colnames(permutationSummary)[col(permutationSummary)],
                       dispersion = as.vector(dispersionMatrix),                             # ** observed dispersion
                       n_perm_ge  = as.vector(permutationSummary))                           # ** permutations >= observed
perm_tab <- perm_tab[!is.na(perm_tab$n_perm_ge), ]                                           # ** keep the upper triangle only
perm_tab$type   <- ifelse(perm_tab$module1 == perm_tab$module2, "within", "between")         # ** within-module vs module-to-module
perm_tab$p_perm <- (perm_tab$n_perm_ge + 1) / (N_PERM + 1)                                   # ** empirical p; +1 avoids p = 0
perm_tab$padj   <- p.adjust(perm_tab$p_perm, method = "BH")                                  # class script 03: BH across reported tests
perm_tab <- perm_tab[order(perm_tab$type != "within", perm_tab$p_perm), ]                    # ** within-module rows first, then by p
write.csv(perm_tab, file.path(out_dir, "GSE10072_5_permutation_results.csv"), row.names = FALSE)

saveRDS(list(dispersionMatrix = dispersionMatrix, permutationSummary = permutationSummary,   # ** keep everything for later sections
             nullDistrib = nullDistrib, N_PERM = N_PERM, seed = 20261005L),
        file.path(out_dir, "GSE10072_5_permutation.rds"))

capture.output(                                                                              # ** text log I can read
  cat("Permutations:", N_PERM, "  modules:", n_mod, "\n\n"),
  cat("Permutations with dispersion >= observed (upper triangle):\n"),
  print(permutationSummary),
  cat("\nWithin-module tests:\n"),
  print(perm_tab[perm_tab$type == "within", ], row.names = FALSE),
  file = file.path(out_dir, "GSE10072_5_log.txt")
)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))                  # class scripts: save sessionInfo

# --- Table 1, module enrichment (expected / observed / FDR (false discovery rate)) 


# ---- fig2b expression traces of one module’s genes in each condition
# Fig 2A, complete (heatmap + boxes + colour bar + expression bands)





