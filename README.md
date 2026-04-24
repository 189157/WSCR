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

# Function definition
run_sc_weighted_gsva <- function(seurat_file,
                                 gene_length_file,
                                 group_rds,
                                 coef_file,
                                 pos_file,
                                 neg_file,
                                 output_file = "PBIS-TP53.txt") {
  # --- Step 1: Load Seurat object and extract counts ---
  scobj <- readRDS(seurat_file)
  counts <- GetAssayData(object = scobj, slot = "counts")
  # --- Step 2: TPM normalization ---
  gene_lengths <- read.table(gene_length_file, header = TRUE, sep = "\t")
  genes_in_counts <- rownames(counts)
  matched_gene_lengths <- gene_lengths[gene_lengths$genesymbol %in% genes_in_counts, ]
  counts <- counts[matched_gene_lengths$genesymbol, ]
  gene_lengths_vector <- matched_gene_lengths$length
  RPK <- counts / (gene_lengths_vector / 1000)
  TPM <- sweep(RPK, 2, colSums(RPK), FUN = "/") * 1e6
  TPM <- log2(TPM + 1)
  rownames(TPM) <- gsub("\\.", "-", rownames(TPM))
  # --- Step 3: transpose expression for downstream use ---
  X <- t(as.matrix(TPM))
  if (length(colnames(X)) != ncol(X)) {
    if (!is.null(colnames(X))) {
      colnames(X) <- colnames(X)[1:ncol(X)]
    } else {
      colnames(X) <- paste0("V", 1:ncol(X))
    }
  }
  # --- Step 4: helper functions for grouping ---
  incidenceMatrix <- function(X, group) {
    n <- nrow(X); p <- ncol(X)
    if (!is.list(group)) stop("Argument 'group' must be a list!")
    J <- length(group)
    grp.mat <- Matrix(0, nrow = J, ncol = p, sparse = TRUE)
    if (is.null(colnames(X))) colnames(X) <- paste("V", 1:p, sep = "")
    if (is.null(names(group))) names(group) <- paste("grp", 1:J, sep = "")
    if (is.numeric(group[[1]])) {
      for (i in 1:J) {
        ind <- group[[i]]
        grp.mat[i, ind] <- 1
        colnames(grp.mat)[ind] <- colnames(X)[ind]
      }
    } else {
      for (i in 1:J) {
        grp.i <- as.character(group[[i]])
        ind <- colnames(X) %in% grp.i
        grp.mat[i, ] <- 1 * ind
        colnames(grp.mat)[ind] <- colnames(X)[ind]
      }
    }
    rownames(grp.mat) <- as.character(names(group))
    if (all(grp.mat == 0)) stop("The names in X do not match the group list!")
    grp.mat
  }
  expandX <- function(X, group) {
    incidence.mat <- incidenceMatrix(X, group)
    over.mat <- Matrix(incidence.mat %*% t(incidence.mat), sparse = TRUE)
    grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat))
    X.latent <- NULL
    names <- NULL
    for (i in 1:nrow(incidence.mat)) {
      idx <- incidence.mat[i, ] == 1
      X.latent <- cbind(X.latent, X[, idx, drop = FALSE])
      names <- c(names, colnames(incidence.mat)[idx])
    }
    colnames(X.latent) <- paste("grp", grp.vec, "_", names, sep = "")
    X.latent
  }
  # --- Step 5: expand expression by groups ---
  group <- readRDS(group_rds)
  expression_data <- expandX(X, group)
  # --- Step 6: read coefficient/pos/neg files ---
  if (!file.exists(coef_file) || !file.exists(pos_file) || !file.exists(neg_file)) {
    stop("One or more input files are missing!")
  }
  weighted_expression_data <- read.table(coef_file, header = TRUE, sep = "\t")
  weighted_expression_data$Weights <- abs(weighted_expression_data$coefficients) + 1
  Pos <- read.table(pos_file, header = TRUE, sep = "\t")
  Neg <- read.table(neg_file, header = TRUE, sep = "\t")
  gene_weights <- setNames(weighted_expression_data$Weights, rownames(weighted_expression_data))
  if (!setequal(colnames(expression_data), names(gene_weights))) {
    stop("Mismatch between gene weights and expression data columns!")
  }
  gene_weights <- gene_weights[colnames(expression_data)]
  # --- Step 7: apply weights ---
  weighted_expression_matrix <- expression_data
  for (gene in colnames(expression_data)) {
    weighted_expression_matrix[, gene] <- expression_data[, gene] * gene_weights[gene]
  }
  # --- Step 8: run GSVA (ssGSEA) ---
  cc <- t(weighted_expression_matrix)
  Enrichment_score <- gsva(expr = as.matrix(cc),
                           list(Pos = Pos$x, Neg = Neg$x),
                           kcdf = "Gaussian", method = "ssgsea")
  score <- Enrichment_score["Pos", ] - Enrichment_score["Neg", ]
  Enrichment_score <- as.data.frame(t(Enrichment_score))
  Enrichment_score$PBIS.TP53 <- score
  Enrichment_score$Cell <- rownames(Enrichment_score)
  # --- Step 9: save results ---
  write.table(Enrichment_score, file = output_file, sep = "\t", quote = FALSE)
  return(Enrichment_score)
}
PBIS.TP53 <- run_sc_weighted_gsva(
  seurat_file = "scRNA-seq.rds",
  gene_length_file = "gencode.v22.annotation.txt",
  group_rds = "filtered_gene_sets.rds",
  coef_file = "coefficients_df(0.65_0.02).txt",
  pos_file = "SGL_pos(0.65_0.02).txt",
  neg_file = "SGL_neg(0.65_0.02).txt",
  output_file = "PBIS-TP53.txt"
)
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
sc_dataset <- readRDS("scRNA-seq.rds")

# Extract group labels, i.e., cluster labels
grp.vec <- sc_dataset@meta.data$seurat_clusters
grp.vec <- as.numeric(as.character(grp.vec))

# Load bulk expression data
set.seed(123)

# Load gencode.v22.annotation.txt
gencode <- read.table(
  "gencode.v22.annotation.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Retain valid gene symbols and remove duplicated gene symbols
gencode <- gencode[!is.na(gencode$genesymbol) & gencode$genesymbol != "", ]
gencode <- gencode[!duplicated(gencode$genesymbol), ]

# Extract gene symbols
gene_names <- gencode$genesymbol

# -----------------------------
# 1. Set the number of bulk samples
# -----------------------------
n_samples <- 100
sample_names <- paste0("BulkSample", seq_len(n_samples))

# -----------------------------
# 2. Randomly generate phenotype data
# The second column is coded as 0/1: 1 = mutation, 0 = wild type
# -----------------------------
phenotype_df <- data.frame(
  Sample = sample_names,
  TP53_status = sample(c(0, 1), n_samples, replace = TRUE),
  stringsAsFactors = FALSE
)

# Check the distribution of TP53 status
table(phenotype_df$TP53_status)

# -----------------------------
# 3. Generate the bulk expression matrix
# Rows = gene symbols
# Columns = samples
# -----------------------------

# Assign a baseline mean expression level to each gene
gene_mu <- rgamma(length(gene_names), shape = 2, rate = 0.3)

# Simulate bulk expression values as continuous variables
bulk_mat <- sapply(seq_len(n_samples), function(i) {
  # Baseline expression plus random noise
  rnorm(length(gene_names), mean = gene_mu, sd = gene_mu * 0.2 + 0.3)
})

bulk_mat <- as.matrix(bulk_mat)
rownames(bulk_mat) <- gene_names
colnames(bulk_mat) <- sample_names

# Set negative expression values to zero
bulk_mat[bulk_mat < 0] <- 0

# -----------------------------
# 5. Save generated files
# -----------------------------
write.table(
  bulk_mat,
  file = "Bulk.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  phenotype_df,
  file = "phenotype.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Load the generated bulk expression matrix
bulk_dataset <- read.table(file = "Bulk.txt", sep = "\t", header = TRUE, row.names = 1)
bulk_dataset <- as.matrix(bulk_dataset)

# Identify genes shared between bulk and single-cell datasets
common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))

# Extract normalized single-cell expression data
sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)

# Compute the correlation matrix between bulk samples and single cells
X <- cor(bulk_dataset[common, ], sc_exprs[common, ])

# Load phenotype data
phenotype <- read.table("phenotype.txt", header = TRUE, sep = "\t")
y <- phenotype[, 2]

# Extract PBIS-TP53 scores from metadata
pbis_scores <- sc_dataset@meta.data$PBIS.TP53

# Compute individual weights, which are inversely proportional to the magnitude of PBIS-TP53
epsilon <- 1e-5
ind_weights <- 1 / (abs(pbis_scores) + epsilon)

# Compute group weights based on group-wise PBIS-TP53 score distributions
uniq_grp <- sort(unique(grp.vec))
grp_weights <- numeric(length(uniq_grp))

for (i in seq_along(uniq_grp)) {
  g <- uniq_grp[i]  # Actual cluster ID, e.g., 0, 1, 2, ...
  grp_scores <- pbis_scores[grp.vec == g]
  grp_weights[i] <- 1 / sqrt(sum(grp_scores^2) + epsilon)
}

# Set numerical tolerance for lambda matching
tolerance <- 1e-6

# Define the parameter grid
alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
           0.6, 0.7, 0.8, 0.9)
lambda_min <- c(0.002)

# Define the total number of cells and the sparsity threshold, i.e., 20% of all cells
total_cells <- ncol(X)
threshold <- 0.2 * total_cells

# Initialize the success flag and the results table
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

    # Fit the ASGL model
    fit <- asgl(
      X, y, grp.vec,
      family = "binomial",
      ind_weights = ind_weights,
      grp_weights = grp_weights,
      alpha = a,
      standardize = FALSE,
      lambda_min = lm
    )

    # Perform cross-validation to select the optimal lambda
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

    # Extract the optimal lambda selected by cross-validation
    lambda.min <- cv_asgl$lambda

    # Match the selected lambda to the lambda sequence in the fitted model
    lambda_min_index <- which(abs(fit$lambda - lambda.min) < tolerance)

    # Extract model coefficients at the selected lambda
    coefficients <- fit$beta[, lambda_min_index]
    coefficients_df <- as.data.frame(coefficients)

    # Identify cell subsets based on the signs of model coefficients
    WSCR_pos <- colnames(X)[which(coefficients_df$coefficients > 0)]
    WSCR_neg <- colnames(X)[which(coefficients_df$coefficients < 0)]
    WSCR_background <- colnames(X)[which(coefficients_df$coefficients == 0)]

    # Calculate the total number of selected WSCR-positive and WSCR-negative cells
    WSCR_total <- length(WSCR_pos) + length(WSCR_neg)

    # Record the results of the current iteration
    results <- rbind(results, data.frame(
      alpha = a,
      lambda_min = lm,
      WSCR_pos_count = length(WSCR_pos),
      WSCR_neg_count = length(WSCR_neg),
      WSCR_background_count = length(WSCR_background),
      success = (WSCR_total < threshold && WSCR_total > 0)
    ))

    # Display the summary of the current iteration
    cat("***********************************************************\n")
    cat("Alpha:", a, " Lambda_min:", lm, "\n")
    cat("WSCR_pos count:", length(WSCR_pos),
        "WSCR_neg count:", length(WSCR_neg), "\n")
    cat("***********************************************************\n")

    # Save the results and exit the loop if a suitable sparsity level is achieved
    if (WSCR_total < threshold && WSCR_total > 0) {
      success <- TRUE

      write.table(WSCR_pos, file = "WSCR_pos.txt", sep = "\t", quote = FALSE)
      write.table(WSCR_neg, file = "WSCR_neg.txt", sep = "\t", quote = FALSE)
      write.table(WSCR_background, file = "WSCR_background.txt", sep = "\t", quote = FALSE)

      write.table(
        data.frame(alpha = a, lambda_min = lm),
        file = "final_parameters.txt",
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )

      message("Results and parameters successfully saved.")
      break
    }
  }

  if (success) break
}

# Write the summary of all iterations
write.table(results, file = "iteration_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Output a message if no valid parameter configuration is found
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

>  *Single-cell and bulk transcriptome integration maps TP53-mutant cell states in Chinese HBV-related HCC for immunotherapy stratification and Anti-gene target prioritization*  
> *[Lai et al., 2026]*

---

## 🧾 License

This project is released under the MIT License.
