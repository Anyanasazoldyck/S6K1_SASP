library(DESeq2)
library(dplyr)
library(gprofiler2)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(pheatmap)

library(amap)

#set -----------------------
setwd("D:/S6K1_SASP")
set.seed(1)
options(warn = 1)
theme_set(theme_classic(base_size = 10))
options(bitmapType = "cairo")

#file ---------
dir.create("results")


#loadingdata ---------
##########################################
#load the count matrix -----

counts <- read.delim("data/GSE218682_GeneCounts.txt", header=TRUE, row.names=1, check.names=FALSE)
counts = subset(counts,apply(counts, 1, mean) >= 1)
#load sample sheet -------
ss<- read.csv("data/ss.csv")

print(colnames(counts)==ss$sample_name)
#########################################





#convert names to gene ID---------------
################################################
esmbl_id<- sub("\\..*$","", rownames(counts))

# Convert Ensembl gene IDs
conv_table <- gconvert(
  query = esmbl_id,
  organism = "mmusculus",   
  target = "MGI" ,  filter_na = T
) #give me only the matched genes 

conv_table_clean <- conv_table %>%
  arrange(input, target_number) %>%   
  distinct(input, .keep_all = TRUE) %>%
  dplyr::select(input, name)                  

counts$ENSEMBL <- sub("\\..*$","", rownames(counts))

counts_prefiltered<- counts %>% left_join( conv_table_clean, by= c("ENSEMBL"="input"))#17306

colSums(is.na(counts_prefiltered)) #1241 unmapped genes 

counts_filtered<- subset(counts_prefiltered, !is.na(counts_prefiltered$name))
counts_filtered <- counts_filtered %>% distinct(name,.keep_all = T)

rownames(counts_filtered)<- counts_filtered$name #making gene names as the main 
counts_filtered <-  counts_filtered[,ss$sample_name]
colnames(counts_filtered)

print(colnames(counts_filtered)==ss$sample_name)


## Run DESeq2----------------

counts_filtered<- as.matrix(counts_filtered)
ss$sample_group<-as.factor(ss$sample_group)


dds = DESeqDataSetFromMatrix(countData=counts_filtered, colData=ss, design=~sample_group)
dds = DESeq(dds)



#diffrential analysis workflow ------------
###########################################
#get normalized expression matrix ########
# Extract the normalised counts--------------
em = as.data.frame(counts(dds, normalized=TRUE))

em = round(em,2)

em_out = cbind(row.names(em), em)
names(em_out)[1] = 'ID'
write.table(em_out, file="EM.csv",row.names=FALSE, sep="\t", quote = FALSE)

#extract dd----------------
old_wt_vs_young =as.data.frame( results(dds, c("sample_group","S6K1_WT_old","S6K1_WT_young")) )

write.table(old_wt_vs_young_sig, file="old_wt_vs_young.csv",row.names=FALSE, sep="\t", quote = FALSE)


old_KO_vs_young = as.data.frame(results(dds, c("sample_group","S6K1_KO_old","S6K1_WT_young")) )

write.table(old_KO_vs_young_sig, file="old_KO_vs_young.csv",row.names=FALSE, sep="\t", quote = FALSE)

old_KO_vs_WT = as.data.frame(results(dds, c("sample_group","S6K1_KO_old","S6K1_WT_old")) )

write.table(old_KO_vs_WT_sig, file="old_KO_vs_WT.csv",row.names=FALSE, sep="\t", quote = FALSE)


#make a master table ################
#--- Add gene column ---
#old_wt_vs_young$gene <- rownames(old_wt_vs_young)
#old_KO_vs_young$gene <- rownames(old_KO_vs_young)
#old_KO_vs_WT$gene    <- rownames(old_KO_vs_WT)

#--- Keep only key columns for merging ---
#cols <- c("gene","log2FoldChange","pvalue","padj")

#old_wt_vs_young_sub <- old_wt_vs_young[, cols]
#old_KO_vs_young_sub <- old_KO_vs_young[, cols]
#old_KO_vs_WT_sub    <- old_KO_vs_WT[, cols]

#--- Merge into a single master table ---
#master <- Reduce(function(x, y) merge(x, y, by="gene", all=TRUE),
  #               list(old_wt_vs_young_sub, old_KO_vs_young_sub, old_KO_vs_WT_sub))

#--- Rename columns for clarity ---
#colnames(master) <- c("gene",
   #                   "log2FC_WTold_vs_WTyoung", "pval_WTold_vs_WTyoung", "padj_WTold_vs_WTyoung",
    #                  "log2FC_KOold_vs_WTyoung", "pval_KOold_vs_WTyoung", "padj_KOold_vs_WTyoung",
     #                 "log2FC_KOold_vs_WTold",  "pval_KOold_vs_WTold",  "padj_KOold_vs_WTold")

#--- Add expression matrix (em) ---
# assumes em has genes as rownames
#em_df <- data.frame(gene = rownames(em), em)
#master <- merge(master, em_df, by="gene", all.x=TRUE)
#write.table(master, file="master.csv",row.names=FALSE, sep="\t", quote = FALSE)


#run PCA####################


em_scaled = na.omit(data.frame(t(scale(t(em)))))

xx = prcomp(t(em_scaled))
pca_coordinates = data.frame(xx$x)

# get % variation 
vars = apply(xx$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100

x_axis_label = paste("PC1 (" ,prop_x, "%)", sep="")
y_axis_label = paste("PC2 (" ,prop_y, "%)", sep="")


# plot
ggp = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = ss$sample_group)) + 
  geom_point() + 
  geom_label_repel(aes(label=ss$sample_name), show.legend=FALSE) + 
  labs(title = "PCA", x= x_axis_label, y= y_axis_label) 


ggsave("results/PCA_liver_samples.png", width = 6, height = 6)



#senescence #####
sen <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP") %>%
  filter(grepl("FRIDMAN_SENESCENCE_UP", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol)

#old WT vs. young WT in sen -------

genes_in_paper <- c("Smurf2","Rab13","Cryab","Igfbp5","Cdkn2a","Thbs1","Gsn","Rras","Ndn",
           "Esm1","Vim","S100a11","Irf7","Tnfaip3","Nlrg1")

hm_matrix <- em_scaled[rownames(em_scaled)%in%genes_in_paper,]

# plot
p<-pheatmap(
  hm_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)
ggsave("results/heatmap_Friedman.png", width = 10, height = 6)
png(p,"results/heatmap_Friedman.png" )
#------------------------------------path analysis ------------------------
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)

# what are the pathways unregulated in old WT vs. young wt ?######
#prepare Kegg gene set 
msi_hallmark <- msigdbr(species = "mouse", collection = "H") %>% 
  dplyr::select(gs_name,gene_symbol)

#prepare gene list
gene_list <- old_wt_vs_young$log2FoldChange
names(gene_list)<- rownames(old_wt_vs_young)
gene_list<- gene_list[!is.na(gene_list)]
head(gene_list)
gene_list<- sort(gene_list,decreasing = T)

#run gsea
gsea_results <- clusterProfiler::GSEA(geneList = gene_list,
                                      TERM2GENE    = msi_hallmark,
                                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                                      verbose      = FALSE
                                  )
#visualize 
gsea_df<- gsea_results@result
top_gsea <- gsea_df[gsea_df$p.adjust<0.05,]
p<- give_me_pretty_top_gsea(top_gsea = top_gsea, analysis = "Old WT Vs. Young WT", context = "Mouse Liver")
ggsave("results/hallmark_oldwt_vs_youngwt.png", width = 8, height = 8)



#what is the difference between KO and WT in old mice?
#prepare gene list
gene_list <- old_KO_vs_WT$log2FoldChange
names(gene_list)<- rownames(old_KO_vs_WT)
gene_list<- gene_list[!is.na(gene_list)]
head(gene_list)
gene_list<- sort(gene_list,decreasing = T)

#run gsea
gsea_results <- clusterProfiler::GSEA(geneList = gene_list,
                                      TERM2GENE    = msi_hallmark,
                                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                                      verbose      = FALSE
)
#visualize 
gsea_df<- gsea_results@result
top_gsea <- gsea_df[gsea_df$p.adjust<0.05,]
p<- give_me_pretty_top_gsea(top_gsea = top_gsea, analysis = "Old KO Vs. Old WT", context = "Mouse Liver")
ggsave("results/hallmark_oldKO_vs_oldwt.png", width = 8, height = 8)

#------------------inflammation across the three conditions ---------------------------
inflamm_pathways <- c("HALLMARK_ALLOGRAFT_REJECTION",      
"HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE",    
 "HALLMARK_TNFA_SIGNALING_VIA_NFKB" , "HALLMARK_IL6_JAK_STAT3_SIGNALING" ,
"HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_IL2_STAT5_SIGNALING",      
"HALLMARK_COMPLEMENT"       )

#get candidate genes ###############
candidat_genes <- list()

for (gs in inflamm_pathways) {

  path_genes <- as.character(top_gsea[  top_gsea$ID==gs, "core_enrichment"])
  path_genes <- unlist(strsplit(path_genes, "/"))
  pathway_name <- gs
  candidat_genes[[pathway_name]] <- path_genes
}
candidat_df <- stack(candidat_genes)
colnames(candidat_df) <- c("gene", "pathway")
candidat_df<- candidat_df%>% distinct(gene,.keep_all = TRUE)


#select the top 20 significant gene in each analysis-----
#analysis 1 : ko vs. wt
old_KO_vs_WT$gene <- rownames(old_KO_vs_WT)
top_in_old_KO_vs_old_wt <- old_KO_vs_WT[old_KO_vs_WT$gene %in% candidat_df$gene, ]
top_in_old_KO_vs_old_wt <- top_in_old_KO_vs_old_wt[top_in_old_KO_vs_old_wt$padj < 0.05, ]
top_in_old_KO_vs_old_wt <- na.omit(top_in_old_KO_vs_old_wt)

old_wt_vs_young$gene<- rownames(old_wt_vs_young)
top_in_old_vs_young_wt <- old_wt_vs_young[old_wt_vs_young$gene %in% candidat_df$gene, ]
top_in_old_vs_young_wt <- top_in_old_vs_young_wt[top_in_old_vs_young_wt$padj < 0.05, ]
top_in_old_vs_young_wt<- na.omit(top_in_old_vs_young_wt)
set1 <- rownames(top_in_old_KO_vs_old_wt)
set2 <- rownames(top_in_old_vs_young_wt)


#common genes 
# Load library
library(ggVennDiagram)

set1 <-top_in_old_KO_vs_old_wt$gene
set2<- top_in_old_vs_young_wt$gene
common <-intersect(x=set1, y=set2)

gg <- ggVennDiagram(list(
  "Old KO vs Old WT" = set1,
  "Old WT vs Young WT" = set2
))

ggsave("venn_diagram.png", gg, width = 6, height = 6, dpi = 300)


#plot a heatmap of the common genes ####
em_inflam <- em_scaled[rownames(em_scaled) %in% common,]
em_inflam<-as.matrix(em_inflam)



p<-pheatmap(
  em_inflam,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)
png("results/hm_common_inflamm_genes.png")
p
dev.off()

#plot the genes from the paper using EM #####
paper_genes <- c(
  "Mmp13","Ccl4","Lgals1","Cxcl10","Ccl2","Col5a1","Cxcl9","Napsa",
  "Angptl7","Tgfb1","Il1a","Mmp7","Mmp12","Ccl5","Gdf3","Ccl3","Nrg1",
  "Timp1","Nptx1","Ccl6","Lgals3","Il1b","Fgf13","Il6","Serpine1"
)

em_inflam <- em_scaled[rownames(em_scaled) %in% paper_genes,]
em_inflam<-as.matrix(em_inflam)



p<-pheatmap(
  em_inflam,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)
png("results/hm_common_inflamm_genes_from_paper.png", res= 300 ,width =300* 8, height = 300*8)
p
dev.off()
