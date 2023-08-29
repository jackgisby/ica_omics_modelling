single_lmer <- function(data, formula_string, REML = TRUE) {
    out.model <- tryCatch(
        lmerTest::lmer(
            as.formula(formula_string),
            data = data,
            REML = REML,
            control = lme4::lmerControl(check.conv.singular = "ignore")
        ),
        warning = function(w) {
            return(lmerTest::lmer(
                as.formula(formula_string),
                data = data,
                REML = REML,
                control = lme4::lmerControl(optimizer = "Nelder_Mead", check.conv.singular = "ignore")
            ))
        }
    )
    
    
    if (class(out.model) == "lmerModLmerTest") {
        return(out.model)
    } else {
        stop("Convergence issue not caught by single_lmer")
    }
}

run_lms <- function(sample_by_feature, pheno_data, formula, method = "lmer", ret_models = FALSE) {
    
    stopifnot(rownames(sample_by_feature) == rownames(pheno_data))
    
    combined_data <- cbind(sample_by_feature, pheno_data)
    
    anovas <- data.frame()
    summaries <- data.frame()
    models <- list()
    
    for (i in 1:ncol(sample_by_feature)) {
        
        if (method == "lmer") {
            
            lmer_res <- single_lmer(combined_data, paste0(colnames(sample_by_feature)[i], formula))
            
            if (ret_models) {
                models[[i]] <- lmer_res
            }
            
            anova_res <- data.frame(anova(lmer_res))
            summary_res <- data.frame(summary(lmer_res)$coefficients)
            
        } else {
            
            lm_res <- lm(paste0(colnames(sample_by_feature)[i], formula), combined_data)
            
            if (ret_models) {
                models[[i]] <- lm_res
            }
            
            anova_res <- data.frame(car::Anova(lm_res))
            summary_res <- data.frame(summary(lm_res)$coefficients)
        }
        
        anova_res$feature <- colnames(sample_by_feature)[i]
        summary_res$feature <- colnames(sample_by_feature)[i]
        
        anova_res$term <- rownames(anova_res)
        summary_res$term <- rownames(summary_res)
        
        anovas <- rbind(anovas, anova_res)
        summaries <- rbind(summaries, summary_res)
    }
    
    anovas$p_adj <- NA
    summaries$p_adj <- NA

    for (term in unique(anovas$term)) {
        anovas$p_adj[anovas$term == term] <- p.adjust(anovas[["Pr..F."]][anovas$term == term], method = "BH")
    }
    
    for (term in unique(summaries$term)) {
        summaries$p_adj[summaries$term == term] <- p.adjust(summaries[["Pr...t.."]][summaries$term == term], method = "BH")
    }
    
    if (ret_models) {
        return(models)
    } else {
        return(list(anovas = anovas, summaries = summaries))
    }
}
