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
  Enrichment_score$Score <- score
  Enrichment_score$Sample <- rownames(Enrichment_score)

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
