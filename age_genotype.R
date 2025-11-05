library(DESeq2)
library(dplyr)
library(gprofiler2)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(msigdbr)
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
old_wt_vs_young_sig <- old_wt_vs_young[
  old_wt_vs_young$padj < 0.05 & 
  old_wt_vs_young$log2FoldChange > 1,
]

write.table(old_wt_vs_young_sig, file="old_wt_vs_young_sig.csv",row.names=FALSE, sep="\t", quote = FALSE)


old_KO_vs_young = as.data.frame(results(dds, c("sample_group","S6K1_KO_old","S6K1_WT_young")) )
old_KO_vs_young_sig <- old_KO_vs_young[
  old_KO_vs_young$padj < 0.05 & 
    old_KO_vs_young$log2FoldChange > 1,
]
write.table(old_KO_vs_young_sig, file="old_KO_vs_young_sig.csv",row.names=FALSE, sep="\t", quote = FALSE)

old_KO_vs_WT = as.data.frame(results(dds, c("sample_group","S6K1_KO_old","S6K1_WT_old")) )
old_KO_vs_WT_sig <- old_KO_vs_WT[
  old_KO_vs_WT$padj < 0.05 & 
    old_KO_vs_WT$log2FoldChange > 1,
]
write.table(old_KO_vs_WT_sig, file="old_KO_vs_WT_sig.csv",row.names=FALSE, sep="\t", quote = FALSE)


#make a master table ################
#--- Add gene column ---
old_wt_vs_young$gene <- rownames(old_wt_vs_young)
old_KO_vs_young$gene <- rownames(old_KO_vs_young)
old_KO_vs_WT$gene    <- rownames(old_KO_vs_WT)

#--- Keep only key columns for merging ---
cols <- c("gene","log2FoldChange","pvalue","padj")

old_wt_vs_young_sub <- old_wt_vs_young[, cols]
old_KO_vs_young_sub <- old_KO_vs_young[, cols]
old_KO_vs_WT_sub    <- old_KO_vs_WT[, cols]

#--- Merge into a single master table ---
master <- Reduce(function(x, y) merge(x, y, by="gene", all=TRUE),
                 list(old_wt_vs_young_sub, old_KO_vs_young_sub, old_KO_vs_WT_sub))

#--- Rename columns for clarity ---
colnames(master) <- c("gene",
                      "log2FC_WTold_vs_WTyoung", "pval_WTold_vs_WTyoung", "padj_WTold_vs_WTyoung",
                      "log2FC_KOold_vs_WTyoung", "pval_KOold_vs_WTyoung", "padj_KOold_vs_WTyoung",
                      "log2FC_KOold_vs_WTold",  "pval_KOold_vs_WTold",  "padj_KOold_vs_WTold")

#--- Add expression matrix (em) ---
# assumes em has genes as rownames
em_df <- data.frame(gene = rownames(em), em)
master <- merge(master, em_df, by="gene", all.x=TRUE)
write.table(master, file="master.csv",row.names=FALSE, sep="\t", quote = FALSE)


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
sen_WT_old_vs_WT_young <- intersect(sen, rownames(old_wt_vs_young_sig))#"Cdkn2a" "Esm1"   "Igfbp5" "Irf7"   "Ndn"    "Nrg1"  
sen_KO_old_vs_WT_young <- intersect(sen, rownames(old_KO_vs_young_sig))#"Cdkn2a" "Esm1"   "Igfbp5"
sen_KO_old_vs_WT_old <- intersect(sen, rownames(old_KO_vs_WT_sig)) #zero sig genes 

View(sen_wt_old_vs_young)
master_sen <- master[master$gene%in%sen,]
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
