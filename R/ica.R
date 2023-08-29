run_mstd <- function(X, env_name = "ica", min_components = 10, max_components = 60, step = 2, 
                     min_mean_stability = 0.90, n_runs = 20, algorithm = "fastica_par", 
                     fun = "logcosh", max_iter = 2000, n_jobs = 1, seed = 1, plot_path = "mstd.py") {
    
    # if no reticulate venv exists, create it
    if (!reticulate::virtualenv_exists(env_name)) {
        reticulate::virtualenv_create(env_name, packages = "stabilized-ica")
    }
    
    # source the code 
    reticulate::use_virtualenv(env_name)
    reticulate::source_python("tools/ica_snf_latent_model/src/mstd.py")
    
    reticulate::py_set_seed(seed)
    n_components <- estimate_components(data.frame(scale(X, scale = FALSE)), min_components = min_components, 
                                        max_components = max_components, step = step, 
                                        min_mean_stability = min_mean_stability, n_runs = n_runs, 
                                        algorithm = algorithm, fun = fun, max_iter = max_iter, 
                                        n_jobs = n_jobs, plot_path = plot_path)
    
    return(n_components)
}

run_ica <- function(X, nc, package = "ica", method = "fast", maxit = 500, tol = 1e-6, seed = 1, ...) {
    
    set.seed(seed)
    
    # scale data
    X_centered <- scale(X, scale = FALSE)
    stopifnot(all.equal(X_centered, scale(X, scale = FALSE, center = attr(X_centered, "scaled:center"))))
    
    # calculate ICA
    if (package == "ica") {
        
        # Methods: c("fast", "imax", "jade")
        ica_res <- ica::ica(
            t(X_centered),
            nc = nc,
            maxit = maxit,
            tol = tol,
            method = method,
            ...
        )
        
        rownames(ica_res$M) <- rownames(X)
        
    } else if (package == "MineICA") {
        
        if (method == "fast") {
            method <- "fastICA"
        } else if (method == "jade") {
            method <- "JADE"
        }
        
        # Methods: c("fastICA", "JADE")
        ica_res <- MineICA::runICA(
            t(X_centered),
            nbComp = nc,
            maxit = maxit,
            tol = tol,
            method = method,
            alg.type = "parallel",
            fun = "logcosh",
            ...
        )
        
        names(ica_res)[names(ica_res) == "A"] <- "M"
        ica_res$M <- t(ica_res$M)
    }
    
    # Add additional attributes
    ica_res$X <- X
    ica_res$X_centered <- X_centered
    ica_res$parameters <- list(nc = nc, package = package, method = method, maxit = maxit, tol = tol)
    
    # name factors
    colnames(ica_res$M) <- colnames(ica_res$S) <- paste0("factor_", 1:ncol(ica_res$S))
    
    # check correct M calculation
    x_divide_s <- ica_res$X_centered %*% t(MASS::ginv(ica_res$S))
    colnames(x_divide_s) <- colnames(ica_res$M)
    stopifnot(all(bazar::almost.equal(ica_res$M, x_divide_s)))
    stopifnot(all.equal(project_ica(X, ica_res)$M, x_divide_s))
    
    # if there is skew, reorient factor
    skew <- ifelse(apply(ica_res$S, 2, moments::skewness) >= 0, 1, -1)
    for (i in 1:ncol(ica_res$S)) {
        ica_res$S[,i] <- ica_res$S[,i] * skew[i]
    }
    
    # scale factors
    ica_res$S <- scale(ica_res$S)
    
    # recalculate M after reorientation and factor scaling
    ica_res <- project_ica(X, ica_res)
    stopifnot(all.equal(ica_res$X, X))
    stopifnot(all.equal(ica_res$X_centered, X_centered))
    
    return(ica_res)
}

project_ica <- function(newdata, ica_res) {

    ica_res$X <- newdata
    ica_res$X_centered <- scale(ica_res$X, scale = FALSE, center = attr(ica_res$X_centered, "scaled:center"))
    
    # X = M * S
    # X / S = M
    # X * inv(S) = M
    ica_res$M <- ica_res$X_centered %*% t(MASS::ginv(ica_res$S))
    colnames(ica_res$M) <- colnames(ica_res$S)
    
    return(ica_res) 
}

create_ica_list <- function(ica_res, pheno_data, mart = NULL) {
    
    ica_res$sample_data <- pheno_data
    
    ica_res$feature_data <- data.frame(ensembl_id = rownames(ica_res$S))
    
    if (is.null(mart)) {
        mart <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart(dataset = "ensembl"))
    }
    
    ensembl_to_geneid <- biomaRt::getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), values = ica_res$feature_data$ensembl_id, mart = mart)
    ensembl_to_geneid <- ensembl_to_geneid[which(!duplicated(ensembl_to_geneid$ensembl_gene_id)),]
    
    ica_res$feature_data <- dplyr::left_join(ica_res$feature_data, ensembl_to_geneid, by = c("ensembl_id" = "ensembl_gene_id"))
    rownames(ica_res$feature_data) <- ica_res$feature_data$ensembl_id
    
    ica_res$feature_data$gene_id <- ica_res$feature_data$hgnc_symbol
    ica_res$feature_data$gene_id[is.na(ica_res$feature_data$gene_id)] <- ica_res$feature_data$ensembl_id[is.na(ica_res$feature_data$gene_id)]
    
    stopifnot(!any(is.na(ica_res$feature_data$ensembl_id)))
    stopifnot(!any(is.na(ica_res$feature_data$gene_id)))
    
    # checks
    stopifnot(all(colnames(ica_res$M) == colnames(ica_res$S)))
    stopifnot(all(rownames(ica_res$M) == rownames(ica_res$sample_data)))
    stopifnot(all(rownames(ica_res$S) == rownames(ica_res$feature_data)))
    stopifnot(all(rownames(ica_res$X) == rownames(ica_res$X_centered)))
    stopifnot(all(colnames(ica_res$X) == colnames(ica_res$X_centered)))
    stopifnot(all(rownames(ica_res$M) == rownames(ica_res$X)))
    stopifnot(all(rownames(ica_res$S) == colnames(ica_res$X)))
    
    return(ica_res)
}
