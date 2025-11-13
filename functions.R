give_me_pretty_top_gsea <- function(top_gsea, analysis=".",context="."){
  top_gsea$`-logFDR`= -log10(top_gsea$p.adjust)
  p=ggplot(top_gsea, aes(
    x = reorder(ID, NES),
    y = NES,
    fill = NES,
    alpha=`-logFDR`,
    #size = `-logFDR`
  )) +
    geom_point(shape = 21, colour = "black",size=10) +
    scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "red",
      midpoint = 0,
      name = "Enrichment"
    ) + scale_alpha(range = c(0.3, 1),
                    name = "LogFDR") +
    
    
    
    labs(
      x = "Pathway",
      y = "NES",
      title = paste0(analysis,"\n",context)
    ) +
    my_theme +
    coord_flip()
  return(p)}
#collect genes from pathway analysis####
collect_candidate_genes= function(top_gsea){
  candidat_genes <- list()
  
  for (i in seq_along(top_gsea[[11]])) {
    path_genes <- as.character(top_gsea[i, 11])
    path_genes <- unlist(strsplit(path_genes, "/"))
    pathway_name <- as.character(top_gsea[i, 2])
    candidat_genes[[pathway_name]] <- path_genes
  }
  candidat_df <- stack(candidat_genes)
  colnames(candidat_df) <- c("gene", "pathway")
  
  return(candidat_df)
  
}
