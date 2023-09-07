
factor_enrichment <- function(ica_res, method = "hyp", z_cutoff = 3, min_z_cutoff_genes = 20, adj_p_cutoff = 0.05, min_genes = 3, species = "human", universe_correction = FALSE) {
    
    enrich_res <- data.frame()
    
    for (i in 1:ncol(ica_res$S)) {
        
        if (method == "hyp") {
            
            # get entrez IDs for genes aligned with factor

            ensembl_ids <- rownames(ica_res$S)[which(abs(ica_res$S[,i]) > z_cutoff)]
            
            if (length(ensembl_ids) < min_z_cutoff_genes) {
                
                ensembl_ids <- rownames(ica_res$S)[order(abs(ica_res$S[,i]), decreasing = TRUE)][1:20]
            }
            
            de <- ica_res$feature_data$entrezgene_id[ica_res$feature_data$ensembl_id %in% ensembl_ids]

            t2g <- msigdbr::msigdbr(species = species, category = "C2") 
            t2g <- t2g[t2g$gs_subcat != "CGP" ,]
            t2g <- as.data.frame(dplyr::distinct(t2g, gs_name, entrez_gene))
            
            if (universe_correction) {
                universe <- as.character(ica_res$feature_data$entrezgene_id)
            } else {
                universe <- NULL
            }
            enrich_res_single <- clusterProfiler::enricher(de, pvalueCutoff = 1, qvalueCutoff = 1, TERM2GENE = t2g, universe = universe)
            enrich_res_single <- enrich_res_single@result
            enrich_res_single <- enrich_res_single[enrich_res_single$Count >= min_genes ,]
            
            enrich_res_single$adj_p <- p.adjust(enrich_res_single$pvalue, method = "BH")

        } else if (method == "ks") {
            
            # not implemented yet
            stopifnot(FALSE)
        }
        
        enrich_res_single <- enrich_res_single[which(enrich_res_single$adj_p < adj_p_cutoff) ,]
        
        if (nrow(enrich_res_single) >= 1) {
            enrich_res_single$factor <- i
            enrich_res_single$method <- method
            enrich_res_single$z_cutoff <- z_cutoff
            
            enrich_res <- rbind(enrich_res, enrich_res_single)
        }
    }
    
    return(enrich_res)
}
