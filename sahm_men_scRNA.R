library("Seurat")

### Load scRNA-Seq dataset from Sahm et al. 
#Sahm classification: ben-1, ben-2, mal, int-A, int-B
#ben = benign, mal = malignant, int = intermediate

Seurat4_filtered <- readRDS("~/Meningioma_Felix_scRNAseq/seurat4_after_filtering.rds")

View(Seurat4_filtered$MC)
table(Seurat4_filtered$MC)

### Color in UMAP by genes of interest

features <- c("ALG3", "SOX11", "HOXD13", "ABCC4", "HLA-DPB1", "HLA-DMB", "HLA-DPA1", "CD74", "NKD1", "SLC29A1", "VIT",  "PIK3R1")
FeaturePlot(Seurat4_filtered, features = features, reduction = "umap")

DimPlot(Seurat4_filtered, group.by = "MC")

features <- c("CD74", "MAGEB4", "RCN3", "KRT34")
FeaturePlot(seurat1, features = features)

#remove cells labeled as MC = int-A and int-B
seurat_sub <- subset(Seurat4_filtered, MC == "ben-1" | MC == "ben-2" | MC == "mal")
DimPlot(seurat_sub, group.by = "MC")

features <- c("ALG3", "SOX11", "HOXD13", "ABCC4", "HLA-DPB1", "HLA-DMB", "HLA-DPA1", "CD74", "NKD1", "SLC29A1", "VIT",  "PIK3R1")
FeaturePlot(seurat_sub, features = features, reduction = "umap")

### DE analysis between ben-1 and mal

seurat_sub2 <- subset(Seurat4_filtered, MC == "ben-1" | MC == "mal")
cell_selection <- SetIdent(seurat_sub2, value = "MC")
DGE_cell_selection <- FindAllMarkers(cell_selection, log2FC.threshold = 0.2, test.use = "wilcox", min.pct = 0.1, 
                                     min.diff.pct = 0.2, only.pos = TRUE, assay = "RNA", max.cells.per.ident = 1250)
write.csv(DGE_cell_selection, "~/Analysis32_scRNA_clusterB/DGE_ben-1_vs_mal_maxcells1250.csv")

### DE analysis between ben-1 and ben-2

seurat_sub2 <- subset(Seurat4_filtered, MC == "ben-1" | MC == "ben-2")
cell_selection <- SetIdent(seurat_sub2, value = "MC")
DGE_cell_selection <- FindAllMarkers(cell_selection, log2FC.threshold = 0.2, test.use = "wilcox", min.pct = 0.1, 
                                     min.diff.pct = 0.2, only.pos = TRUE, assay = "RNA")
write.csv(DGE_cell_selection, "~/Analysis32_scRNA_clusterB/DGE_ben-1_vs_ben-2.csv")

### Explore immune and non-immune cells in meningioma tumor subtypes

seurat1 <- readRDS("~/Meningioma_Felix_scRNAseq/seurat1_basic.rds")

seu_int <- Seurat::ScaleData(seurat1)
seu_int <- FindVariableFeatures(seu_int)
seu_int <- Seurat::RunPCA(seu_int, npcs = 30)
seu1_umap <- RunUMAP(seu_int, dims = 1:30)
DimPlot(seu1_umap, group.by = "MC", reduction = "umap")

seu_sub1 <- subset(seu_int, MC == "ben-1" | MC == "mal" | MC == "ben-2")
seu1_umap_1 <- RunUMAP(seu_sub1, dims = 1:30)
DimPlot(seu1_umap_1, group.by = "MC", reduction = "umap")

# add immune / non-immune labels to seu_sub1 object
nonimmune <- Cells(seurat_sub)
allcells <- Cells(seu_sub1)
seu_sub1$celltype <- ifelse(test = allcells %in% nonimmune, yes = "non-immune", no = "immune")

seu1_umap_1 <- RunUMAP(seu_sub1, dims = 1:30)
DimPlot(seu1_umap_1, group.by = "celltype", reduction = "umap", raster = FALSE)

features <- c("SOX11", "HOXD13", "HLA-DPB1", "HLA-DMB", "HLA-DPA1", "CD74", "HLA-A")
FeaturePlot(seu1_umap, features = features, reduction = "umap", raster = FALSE)


plot <- DimPlot(seu1_umap, group.by = "celltype", reduction = "umap", raster = FALSE)

plot_title_size = 20
legend_pt_size = 4
axis_text_size = 25
axis_title_size = 25
legend_text_size = 15
spacing = 1
chosen_margin = c(0.5,1,0.5,1) #top,right,bottom,left

theme_legend <- theme_blank()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size, face = "bold"),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "right",
        legend.justification = "left",
        legend.title = element_blank() )


theme_nolegend <- theme_blank()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "none",
        legend.justification = "left",
        legend.title = element_blank() )

plot + 
  theme_nolegend +
  ggtitle("Immune and non-immune cells") + 
  scale_color_manual(values = c("immune" = "gold",
                                "non-immune" = "red"
                                
  ))

setwd("~/Meningioma_Eric/Output_1298_April2023/")
ggsave("scRNA_immune_nonimmune_leg.pdf", width = 10, height = 7)
ggsave("scRNA_immune_nonimmune_noleg.pdf", width = 7, height = 7)

# only non-immune cells

seurat_sub <- subset(Seurat4_filtered, MC == "ben-1" | MC == "ben-2" | MC == "mal")
plot <- DimPlot(seurat_sub, group.by = "MC", raster = FALSE)

plot+ 
  theme_nolegend + 
  ggtitle("Non-immune cells") +
  scale_color_manual(values = c("ben-1" = "blue", 
                                "ben-2" = "red", 
                                "mal" = "forestgreen"))
ggsave("scRNA_nonimmune_leg.pdf", width = 10, height = 7)
ggsave("scRNA_nonimmune_noleg.pdf", width = 7, height = 7)

features <- c("IGKC", "IGHM", "CD37", "IGHD")
FeaturePlot(seurat_sub, features = features, reduction = "umap", raster = FALSE, order = TRUE)

### Differentiate tumor vs non-tumor (presumably immune cells) in meningioma tumors samples
### A list of genes that are highly expressed in meningioma tumors (in multiple subtypes)...
### ... compared to wildtype normal meningeal cells (Arachnoid cells) is used as marker-genes to identify tumor cells

# load markergenes 
markergenes <- as.list(markergenes)

seurat1 <- readRDS("~/HollandLabShared/Sonali/Frank/Meningioma_Felix_scRNAseq/seurat1_basic.rds")
 
# extract counts
 counts <- GetAssayData(seurat1, slot="counts", assay="RNA") %>% as.matrix(.)

# add your logic to remove meningioma specific genes here ...
 expression_threshold <- 1
subset_count <- counts[rownames(counts) %in% markergenes, ]
subset_count <- subset_count[, colSums(subset_count > expression_threshold) > 1] 

sel_cells <- colnames(subset_count)
 
# create new object with only tumor cells which have counts in markergenes. 
 counts.sub <- counts[, sel_cells]
#new_coldata = coldata[sel_cells, ] #what's this?
 new_seurat_object <- CreateSeuratObject(counts=counts.sub) 
 tumorcells <- Cells(new_seurat_object)

#generate UMAP of tumor cells
 seu_int1 <- Seurat::ScaleData(new_seurat_object)
 seu_int1 <- FindVariableFeatures(seu_int1)
 seu_int1 <- Seurat::RunPCA(seu_int1, npcs = 30)
 seutumor_umap <- RunUMAP(seu_int1, dims = 1:30)
 DimPlot(seutumor_umap, group.by = "MC", reduction = "umap", raster = FALSE)
 
 features <- c("ANXA10", "KRT34", "CCDC201", "HOXD13")
 FeaturePlot(seutumor_umap, features = features, reduction = "umap", raster = FALSE, order = TRUE)
 
 features <- c("CD19", "AIF1", "CD40", "CD14")
 FeaturePlot(seutumor_umap, features = features, reduction = "umap", raster = FALSE, order = TRUE)
 
# subset seurat1 object by sel_cells
 allcells <- Cells(seurat1)
 seurat1$celltype <- ifelse(test = allcells %in% tumorcells, yes = "tumor", no = "nontumor")
 seu_int2 <- Seurat::ScaleData(seurat1)
 seu_int2 <- FindVariableFeatures(seu_int2)
 seu_int2 <- Seurat::RunPCA(seu_int2, npcs = 30)
 seurat1_umap_2 <- RunUMAP(seu_int2, dims = 1:30)
 DimPlot(seurat1_umap_2, group.by = "celltype", reduction = "umap", raster = FALSE)
 DimPlot(seurat1_umap_2, group.by = "MC", reduction = "umap", raster = FALSE)
 features <- c("CD19", "CD86", "CD40", "CD14")
 FeaturePlot(seurat1_umap_2, features = features, reduction = "umap", raster = FALSE, order = TRUE)
 features <- c("ANXA10", "KRT34", "CCDC201", "AFP", "HOXD13")
 FeaturePlot(seurat1_umap_2, features = features, reduction = "umap", raster = FALSE, order = TRUE)
 features <- c("CD74", "HLA-DPA1")
 FeaturePlot(seurat1_umap_2, features = features, reduction = "umap", raster = FALSE, order = TRUE)
 
 
 #label seu_sub1 by tumor vs nontumor cells
 seu_sub1$celltype <- ifelse(test = allcells %in% tumorcells, yes = "tumor", no = "nontumor")
 seu1_umap_1 <- RunUMAP(seu_sub1, dims = 1:30)
 DimPlot(seu1_umap_1, group.by = "celltype", reduction = "umap", raster = FALSE)
 
 features <- c("CD19", "CD74", "CD40")
 FeaturePlot(seu1_umap_1, features = features, reduction = "umap", raster = FALSE, order = TRUE)
 
 features <- c("ANXA10", "KRT34", "CCDC201", "AFP")
 FeaturePlot(seu1_umap_1, features = features, reduction = "umap", raster = FALSE, order = TRUE)
 
 ##select ben-1, ben-2 and mal cells in Seurat1 (all cells object)
 seurat1_sub <- subset(seurat1, MC == "ben-1" | MC == "ben-2" | MC == "mal")
 seu_int3 <- Seurat::ScaleData(seurat1_sub)
 seu_int3 <- FindVariableFeatures(seu_int3)
 seu_int3 <- Seurat::RunPCA(seu_int3, npcs = 30)
 seurat1_sub_umap <- RunUMAP(seu_int3, dims = 1:30)
 DimPlot(seurat1_sub_umap, group.by = "MC")
 
 ##identify tumor cells in Seurat1
 allcells <- Cells(seurat1_sub)
 seurat1_sub$celltype <- ifelse(test = allcells %in% tumorcells, yes = "tumor", no = "nontumor")
 seu_int3 <- Seurat::ScaleData(seurat1_sub)
 seu_int3 <- FindVariableFeatures(seu_int3)
 seu_int3 <- Seurat::RunPCA(seu_int3, npcs = 30)
 seurat1_sub_umap <- RunUMAP(seu_int3, dims = 1:30)
 DimPlot(seurat1_sub_umap, group.by = "celltype", reduction = "umap", raster = FALSE)
 
 ##Color in by immune markers and tumor marker genes
 features <- c("ANXA10", "KRT34", "CCDC201", "AFP", "HOXD13")
 FeaturePlot(seurat1_sub_umap, features = features, reduction = "umap", raster = FALSE, order = TRUE)
 
 features <- c("CD19", "CD86", "CD40", "CD14")
 FeaturePlot(seurat1_sub_umap, features = features, reduction = "umap", raster = FALSE, order = TRUE)
 
 
