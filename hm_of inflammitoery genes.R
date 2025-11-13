
#add cell_type to metadata--------------
sc_10x_flex <- SetIdent(sc_10x_flex,value = "seurat_clusters")
cell_type <- readxl::read_xlsx("results/pre_merge/flex/cell_type.xlsx",sheet = "flex")

new.cluster.ids = cell_type$cell_type
names(new.cluster.ids) = levels(sc_10x_flex)
sc_10x_flex = RenameIdents(sc_10x_flex, new.cluster.ids)
sc_10x_flex = AddMetaData(sc_10x_flex, sc_10x_flex@active.ident, col.name = "cell_types")

p=DimPlot(sc_10x_flex, label = T, label.size = 6, repel = T,pt.size = 1)
ggsave(file.path(save_file,"cell_type.png"), width = 20, height = 15)

saveRDS(sc_10x_flex,"objects/sc_10x_flex_annotated.rds")



#select the top 20 significant genes based on p adjusted value 

