# 🧮 WSCR: Weighted Sparse Cell Regularization

This repository provides an R-based implementation of **WSCR (Weighted Sparse Cell Regularization)** — a phenotype-guided, sparsity-inducing algorithm designed to **identify TP53 mutation–associated cellular populations** by integrating single-cell RNA sequencing (scRNA-seq) and bulk transcriptomic data in **HBV-related hepatocellular carcinoma (HCC)**.

---

![WSCR Framework](https://raw.githubusercontent.com/189157/WSCR/main/WSCR_overview.png?raw=true&v=2)

---

## 🧠 Background

Tumor Protein 53 (**TP53**) mutation represents a defining hallmark of HBV-related HCC in Chinese patients, yet its **cellular heterogeneity** and **microenvironmental impact** remain poorly resolved.  
To address this, WSCR was developed to bridge **bulk phenotype labels** (e.g., TP53 mutation status) and **single-cell transcriptomic profiles**, enabling fine-grained discovery of phenotype-associated subpopulations.

---

## 📦 Requirements

Install dependencies in R:

```r
install.packages(c("Seurat", "Matrix", "ggplot2", "dplyr"))
install.packages("BiocManager")
BiocManager::install("preprocessCore")
# For adaptive sparse group lasso
# If not installed:
# install.packages("devtools")
devtools::install_github("jeffdaniel/asgl")
```

---

## 🚀 Quickstart

Below we provide **complete**, runnable code blocks for both PBIS-TP53 scoring on scRNA-seq and WSCR cell identification.  
> 💡 You can copy these sections into `.R` scripts and run them as-is after preparing your inputs.

### 1) Calculation of `PBIS-TP53` in scRNA-seq data

```r
### Calculation of PBIS.TP53 in scRNA-seq data
library(Seurat)
library(Matrix)
library(GSVA)

set.seed(123)

# -----------------------------
# 1. Load gencode.v22.annotation.txt
# -----------------------------
gencode <- read.table(
  "gencode.v22.annotation.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove entries with missing or empty gene symbols, as well as duplicated gene symbols
gencode <- gencode[!is.na(gencode$genesymbol) & gencode$genesymbol != "", ]
gencode <- gencode[!duplicated(gencode$genesymbol), ]

# Use gene symbols directly from the genesymbol column
gene_names <- gencode$genesymbol

# -----------------------------
# 2. Parameter settings
# -----------------------------
n_cells <- 1000
n_clusters <- 8
n_samples <- 10
cell_names <- paste0("Cell", seq_len(n_cells))

# -----------------------------
# 3. Randomly generate the count matrix
# -----------------------------
# To better mimic scRNA-seq data, assign different mean expression levels to different genes
gene_mu <- rgamma(length(gene_names), shape = 2, rate = 0.4)

counts_mat <- sapply(seq_len(n_cells), function(i) {
  rnbinom(length(gene_names), mu = gene_mu, size = 1)
})

counts_mat <- as.matrix(counts_mat)
rownames(counts_mat) <- gene_names
colnames(counts_mat) <- cell_names
counts_sparse <- Matrix(counts_mat, sparse = TRUE)

# -----------------------------
# 4. Construct a Seurat object
# -----------------------------
scobj <- CreateSeuratObject(
  counts = counts_sparse,
  project = "WSCR_mock"
)

# Retain only the two required metadata fields
scobj$seurat_clusters <- sample(0:(n_clusters - 1), ncol(scobj), replace = TRUE)
scobj$orig.ident <- sample(paste0("Sample", seq_len(n_samples)), ncol(scobj), replace = TRUE)

# Save the Seurat object
saveRDS(scobj, file = "scRNA-seq.rds")
```

---

### 2) WSCR (Weighted Sparse Cell Regularization)

```r
### WSCR
# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(preprocessCore)
library(Matrix)
library(asgl)

# Load single-cell dataset
sc_dataset <- readRDS("scobj.rds")

# Extract group (cluster) labels
grp.vec <- sc_dataset@meta.data$seurat_clusters
grp.vec <- as.numeric(as.character(grp.vec))

# Load bulk expression data
bulk_dataset <- read.table(file = "Bulk.txt", sep = "\t", header = TRUE, row.names = 1)
bulk_dataset <- as.matrix(bulk_dataset)

# Identify genes shared between bulk and single-cell datasets
common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)

# Compute correlation matrix between bulk and single-cell gene expression
X <- cor(bulk_dataset[common, ], sc_exprs[common, ])

# Load phenotype data
phenotype <- read.table("phenotype.txt", header = TRUE, sep = "\t")
y <- phenotype[, 2]

# Extract PBIS.TP53 from metadata
pbis_scores <- sc_dataset@meta.data$PBIS.TP53

# Compute individual weights (inverse proportional to PBIS-TP53 magnitude)
epsilon <- 1e-5
ind_weights <- 1 / (abs(pbis_scores) + epsilon)

# Compute group weights based on group-wise PBIS-TP53 distributions
uniq_grp <- sort(unique(grp.vec))
grp_weights <- numeric(length(uniq_grp))
for (i in seq_along(uniq_grp)) {
  g <- uniq_grp[i]  # True cluster ID (e.g., 0, 1, 2, ...)
  grp_scores <- pbis_scores[grp.vec == g]
  grp_weights[i] <- 1 / sqrt(sum(grp_scores^2) + epsilon)
}

# Numerical tolerance for lambda matching
tolerance <- 1e-6

# Define parameter grid
alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
           0.6, 0.7, 0.8, 0.9)
lambda_min <- c(0.002)

# Define total number of cells and sparsity threshold (20% of all cells)
total_cells <- ncol(X)
threshold <- 0.2 * total_cells

# Initialize success flag and results table
success <- FALSE
results <- data.frame(
  alpha = numeric(),
  lambda_min = numeric(),
  WSCR_pos_count = numeric(),
  WSCR_neg_count = numeric(),
  WSCR_background_count = numeric(),
  success = logical()
)

# Parameter tuning loop
for (a in alpha) {
  for (lm in lambda_min) {
    set.seed(123)

    # Fit ASGL model
    fit <- asgl(
      X, y, grp.vec,
      family = "binomial",
      ind_weights = ind_weights,
      grp_weights = grp_weights,
      alpha = a,
      standardize = FALSE,
      lambda_min = lm
    )

    # Cross-validation for optimal lambda
    source("cv_asgl.R")
    cv_asgl <- cv_asgl(
      X, y, grp.vec,
      family = "binomial",
      ind_weights = ind_weights,
      grp_weights = grp_weights,
      alpha = a,
      standardize = FALSE,
      nfolds = 10,
      lambda = fit$lambda,
      lambda_min = lm
    )

    lambda.min <- cv_asgl$lambda
    lambda_min_index <- which(abs(fit$lambda - lambda.min) < tolerance)
    coefficients <- fit$beta[, lambda_min_index]
    coefficients_df <- as.data.frame(coefficients)

    # Identify cell subsets based on coefficient signs
    WSCR_pos <- colnames(X)[which(coefficients_df$coefficients > 0)]
    WSCR_neg <- colnames(X)[which(coefficients_df$coefficients < 0)]
    WSCR_background <- colnames(X)[which(coefficients_df$coefficients == 0)]
    WSCR_total <- length(WSCR_pos) + length(WSCR_neg)

    # Record iteration results
    results <- rbind(results, data.frame(
      alpha = a,
      lambda_min = lm,
      WSCR_pos_count = length(WSCR_pos),
      WSCR_neg_count = length(WSCR_neg),
      WSCR_background_count = length(WSCR_background),
      success = (WSCR_total < threshold && WSCR_total > 0)
    ))

    # Display iteration summary
    cat("***********************************************************\n")
    cat("Alpha:", a, " Lambda_min:", lm, "\n")
    cat("WSCR_pos count:", length(WSCR_pos),
        "WSCR_neg count:", length(WSCR_neg), "\n")
    cat("***********************************************************\n")

    # Save and exit if suitable sparsity found
    if (WSCR_total < threshold && WSCR_total > 0) {
      success <- TRUE
      write.table(WSCR_pos, file = "WSCR_pos.txt", sep = "\t", quote = FALSE)
      write.table(WSCR_neg, file = "WSCR_neg.txt", sep = "\t", quote = FALSE)
      write.table(WSCR_background, file = "WSCR_background.txt", sep = "\t", quote = FALSE)
      write.table(data.frame(alpha = a, lambda_min = lm),
                  file = "final_parameters.txt",
                  sep = "\t", row.names = FALSE, quote = FALSE)
      message("Results and parameters successfully saved.")
      break
    }
  }
  if (success) break
}

# Write summary of all iterations
write.table(results, file = "iteration_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Output message if no valid configuration found
if (!success) {
  message("Could not find a configuration where 0 < (WSCR_pos + WSCR_neg) < 20% of total cells.")
}
```

---

## 📊 Output Files

| File | Description |
|------|-------------|
| `WSCR_pos.txt` | Cells with positive coefficients (TP53-mutant-like) |
| `WSCR_neg.txt` | Cells with negative coefficients (wild-type-like) |
| `WSCR_background.txt` | Cells with no association |
| `iteration_results.txt` | Summary of sparsity and convergence across parameter grid |
| `final_parameters.txt` | Final selected parameters ensuring 0 < WSCR_pos+WSCR_neg < 20% total cells |

---

## 🧩 Integration with PBIS-TP53

WSCR and PBIS-TP53 are **complementary** components:  
- **WSCR** identifies *which cells* are most strongly associated with TP53 mutation phenotype.  
- **PBIS-TP53** quantifies *how strongly each sample* reflects TP53 downstream transcriptional programs.  

You can apply WSCR-identified cell subsets (e.g., `WSCR_pos.txt`, `WSCR_neg.txt`) as input to PBIS-TP53 for pathway-based scoring and therapeutic inference.

---

## 📚 Citation

>  *Decoding TP53 Mutation at Single-Cell Resolution in Chinese HBV-Related Hepatocellular Carcinoma: From Microenvironment to Clinical Translation*  
> *[Lai et al., 2025]*

---

## 🧾 License

This project is released under the MIT License.
