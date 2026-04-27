# =============================================================================
# Deciphering Programmed Cell Death for Gastric Cancer Prognosis: AKR1B1 as a Key Player
# GitHub repository: [your_repo_url]
# R version 4.2.1, required packages: see below
# =============================================================================

# ----------------------------- 0. Setup --------------------------------------
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

required_packages <- c(
    "survival", "survminer", "glmnet", "randomForest", "e1071", "sigFeature",
    "ggplot2", "timeROC", "GSVA", "estimate", "Seurat", "MAESTRO", "pROC", "boot"
)

for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
        if (pkg %in% c("GSVA", "estimate", "timeROC", "sigFeature", "MAESTRO"))
            BiocManager::install(pkg)
        else
            install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}

# ----------------------------- 1. Data Loading & Preprocessing ----------------
load_tcga_data <- function() {
    set.seed(42)
    n_samples <- 375
    n_genes <- 2000
    expr <- matrix(rnorm(n_samples * n_genes), nrow = n_genes, ncol = n_samples)
    rownames(expr) <- paste0("GENE", 1:n_genes)
    colnames(expr) <- paste0("TCGA-", 1:n_samples)
    clinical <- data.frame(
        sample = colnames(expr),
        OS.time = round(runif(n_samples, 30, 3650)),
        OS = rbinom(n_samples, 1, 0.4)
    )
    return(list(expr = expr, clinical = clinical))
}

# --------------------------------- 2. Cox Univariate Filtering ----------------
univariate_cox_filter <- function(expr, time, status, p_cutoff = 0.05) {
    genes <- rownames(expr)
    pvals <- numeric(length(genes))
    for (i in seq_along(genes)) {
        cox <- coxph(Surv(time, status) ~ as.numeric(expr[i, ]))
        pvals[i] <- summary(cox)$coefficients[5]
    }
    names(pvals) <- genes
    sig_genes <- names(pvals[pvals < p_cutoff])
    return(sig_genes)
}

# ---------------------------- 3. LASSO for 12 PCD Models -----------------------
build_lasso_model <- function(expr_subset, time, status) {
    x <- t(expr_subset)
    y <- Surv(time, status)
    set.seed(123)
    cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
    lambda_min <- cv_fit$lambda.min
    coefs <- as.matrix(coef(cv_fit, s = lambda_min))
    selected_genes <- rownames(coefs)[coefs[,1] != 0]
    return(list(model = cv_fit, selected_genes = selected_genes, lambda = lambda_min))
}

# ---------------------------- 4. Multivariate Cox & Risk Score -----------------
multivariate_cox_model <- function(expr_hub, time, status) {
    df <- data.frame(t(expr_hub), time = time, status = status)
    formula_str <- paste("Surv(time, status) ~", paste(rownames(expr_hub), collapse = " + "))
    cox_model <- coxph(as.formula(formula_str), data = df)
    return(cox_model)
}

compute_risk_score <- function(expr_hub, coef_vector) {
    common_genes <- intersect(rownames(expr_hub), names(coef_vector))
    expr_common <- expr_hub[common_genes, , drop = FALSE]
    beta <- coef_vector[common_genes]
    risk_scores <- colSums(expr_common * beta)
    return(risk_scores)
}

# ---------------------------- 5. Machine Learning for Key Gene Selection ------
random_forest_importance <- function(expr_subset, time, status) {
    df <- data.frame(t(expr_subset), risk = (time > median(time) & status == 0))
    df$risk <- as.factor(df$risk)
    set.seed(123)
    rf <- randomForest(risk ~ ., data = df, importance = TRUE, ntree = 500)
    importance_df <- as.data.frame(importance(rf))
    importance_df$gene <- rownames(importance_df)
    importance_df <- importance_df[order(-importance_df$MeanDecreaseGini), ]
    top5 <- head(importance_df$gene, 5)
    return(top5)
}

lasso_feature_selection <- function(expr_subset, time, status, top_k = 5) {
    response <- ifelse(time <= median(time) & status == 1, 1, 0)
    x <- t(expr_subset)
    set.seed(123)
    cv_lasso <- cv.glmnet(x, response, family = "binomial", alpha = 1, nfolds = 10)
    coefs <- as.matrix(coef(cv_lasso, s = "lambda.min"))
    selected <- rownames(coefs)[coefs[,1] != 0][-1]
    if (length(selected) > top_k) selected <- selected[1:top_k]
    return(selected)
}

svm_rfe_feature_selection <- function(expr_subset, time, status) {
    labels <- ifelse(time <= median(time) & status == 1, 1, -1)
    x <- t(expr_subset)
    set.seed(123)
    svm_model <- svm(x, labels, kernel = "linear", cost = 1, scale = FALSE)
    w <- t(svm_model$coefs) %*% svm_model$SV
    importance <- abs(w)
    names(importance) <- colnames(x)
    top5 <- names(sort(importance, decreasing = TRUE)[1:5])
    return(top5)
}

# ---------------------------- 6. External Validation --------------------------
validate_external <- function(expr_external, coef_vector, time_external, status_external, cutoff) {
    risk <- compute_risk_score(expr_external, coef_vector)
    group <- ifelse(risk > cutoff, "High", "Low")
    fit <- survfit(Surv(time_external, status_external) ~ group)
    km_plot <- ggsurvplot(fit, data = NULL, pval = TRUE, risk.table = TRUE)
    roc_data <- timeROC(T = time_external, delta = status_external,
                        marker = risk, cause = 1, times = c(365, 1095, 1825))
    return(list(km_plot = km_plot, roc = roc_data, risk_group = group))
}

# ---------------------------- 7. Immune Infiltration (ssGSEA) -------------------
compute_ssgsea <- function(expr, gene_sets) {
    ssgsea_scores <- gsva(as.matrix(expr), gene_sets, method = "ssgsea",
                          ssgsea.norm = TRUE, verbose = FALSE)
    return(ssgsea_scores)
}

cor_immune <- function(akr1b1_expr, immune_scores) {
    cor_results <- apply(immune_scores, 1, function(x) cor.test(x, akr1b1_expr, method = "spearman"))
    return(cor_results)
}

# ---------------------------- 8. Drug Sensitivity (IC50) ------------------------
cor_ic50 <- function(akr1b1_expr, ic50_values) {
    cor_test <- cor.test(akr1b1_expr, ic50_values, method = "spearman")
    return(cor_test)
}

# ---------------------------- 9. Single-cell Analysis (Seurat) -----------------
run_singlecell_analysis <- function(h5_path) {
    library(Seurat)
    data <- Read10X_h5(h5_path)
    seurat_obj <- CreateSeuratObject(counts = data, project = "GSE134520")
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, npcs = 30)
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)
    FeaturePlot(seurat_obj, features = "AKR1B1")
    return(seurat_obj)
}

# ---------------------------- 10. Statistical Tests & Figures ----------------
wilcox_test_groups <- function(expression, grouping) {
    test_res <- wilcox.test(expression ~ grouping, exact = FALSE)
    return(test_res$p.value)
}

km_analysis <- function(time, status, group) {
    fit <- survfit(Surv(time, status) ~ group)
    p_val <- surv_pvalue(fit)$pval
    return(list(fit = fit, pvalue = p_val))
}

# ---------------------------- 11. Save Results -------------------------
save_results <- function(risk_scores, cox_model, immune_scores, drug_cor, file_prefix = "results") {
    write.csv(data.frame(sample = names(risk_scores), risk_score = risk_scores),
              file = paste0(file_prefix, "_risk_scores.csv"), row.names = FALSE)
    saveRDS(cox_model, file = paste0(file_prefix, "_cox_model.rds"))
    write.csv(t(immune_scores), file = paste0(file_prefix, "_immune_scores.csv"))
    write.csv(drug_cor, file = paste0(file_prefix, "_drug_sensitivity.csv"))
}

# ---------------------------- 12. Session Info --------------------------------
sessionInfo()
