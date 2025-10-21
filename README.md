# 🧮 WSCR: Weighted Sparse Cell Regularization  

This repository provides an R-based implementation of **WSCR (Weighted Sparse Cell Regularization)** — a phenotype-guided, sparsity-inducing algorithm designed to **identify TP53 mutation–associated cellular populations** by integrating single-cell RNA sequencing (scRNA-seq) and bulk transcriptomic data in **HBV-related hepatocellular carcinoma (HCC)**.

---

![WSCR Framework](https://raw.githubusercontent.com/CQMU-LGC/WSCR/main/WSCR_overview.png?raw=true&v=2)

---

## 🧠 Background  

Tumor Protein 53 (**TP53**) mutation represents a defining hallmark of HBV-related HCC in Chinese patients, yet its **cellular heterogeneity** and **microenvironmental impact** remain poorly resolved.  
To address this, WSCR was developed to bridge **bulk phenotype labels** (e.g., TP53 mutation status) and **single-cell transcriptomic profiles**, enabling fine-grained discovery of phenotype-associated subpopulations.

---

## ⚙️ Algorithm Overview  

WSCR formulates phenotype-guided cell discovery as a **weighted sparse optimization** problem:

\[
\min_{\beta}\|Y - X\beta\|^2 + \lambda_1 \|\mathbf{w}_{\text{ind}}\odot \beta\|_1 + \lambda_2 \|\mathbf{w}_{\text{grp}}\odot \beta\|_2^2
\]

where:  
- **X**: correlation matrix between bulk and single-cell gene expression  
- **Y**: phenotype vector (e.g., TP53 mutation status or PBIS-TP53 score)  
- **β**: coefficients linking single-cell features to phenotype signals  
- **w_ind, w_grp**: individual- and group-level weights derived from PBIS-TP53 scores and Seurat clusters  
- **λ₁, λ₂**: sparsity and smoothness hyperparameters  

The algorithm iteratively tunes α (elasticity) and λ (regularization strength) to achieve biologically interpretable sparsity — identifying **WSCR-positive**, **WSCR-negative**, and **background** cell subsets.

---

## 📦 Requirements  

Install dependencies in R:  
```r
install.packages(c("Seurat", "Matrix", "ggplot2", "dplyr"))
install.packages("BiocManager")
BiocManager::install("preprocessCore")
# For adaptive sparse group lasso
devtools::install_github("CQMU-LGC/asgl")
