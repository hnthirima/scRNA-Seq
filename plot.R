
library(dplyr)
library(ggplot2)

dat = readRDS('~/HollandLabShared/Nayanga/Meningioma_Eric/Analysis33_mouse_embryonic/data.rds')
gene_module_list = colnames(dat)[1:10]

major_trajectory_color_plate = c("Neuroectoderm_and_glia"             = "#f96100",
                                 "Intermediate_neuronal_progenitors"  = "#2e0ab7",
                                 "Eye_and_other"             = "#00d450",
                                 "Ependymal_cells"           = "#b75bff",
                                 "CNS_neurons"               = "#e5c000",
                                 "Mesoderm"                  = "#bb46c5",
                                 "Definitive_erythroid"      = "#dc453e",
                                 "Epithelium"                = "#af9fb6",
                                 "Endothelium"               = "#00a34e",
                                 "Muscle_cells"              = "#ffa1f5",
                                 "Hepatocytes"               = "#185700",
                                 "White_blood_cells"         = "#7ca0ff",
                                 "Neural_crest_PNS_glia"     = "#fff167",
                                 "Adipocytes"                = "#7f3e39",
                                 "Primitive_erythroid"       = "#ffa9a1",
                                 "Neural_crest_PNS_neurons"  = "#b5ce92",
                                 "T_cells"                   = "#ff9d47",
                                 "Lung_and_airway"           = "#02b0d1",
                                 "Intestine"                 = "#ff007a",
                                 "B_cells"                   = "#01b7a6",
                                 "Olfactory_sensory_neurons" = "#e6230b",
                                 "Cardiomyocytes"            = "#643e8c",
                                 "Oligodendrocytes"          = "#916e00",
                                 "Mast_cells"                = "#005361",
                                 "Megakaryocytes"            = "#3f283d",
                                 "Testis_and_adrenal"        = "#585d3b")
for(i in 1:10){
  name = gene_module_list[i]
  print(name)
  
  df = dat[,c("major_cell_cluster", "celltype")]
  df$score = as.vector(dat[[name]])
  
  df_x = df %>% group_by(major_cell_cluster, celltype) %>% summarise(median_score = median(score))
  df_x = df_x %>%
    group_by(major_cell_cluster) %>%
    arrange(median_score, .by_group=TRUE)
  
  df$celltype = factor(df$celltype, levels = as.vector(df_x$celltype))
  
  try(ggplot(df, aes(celltype, score, fill = major_cell_cluster)) + geom_boxplot(outlier.shape = NA) + coord_flip() +
        labs(x="celltype", y="gene_module_score", title=name) +
        theme_classic(base_size = 10) +
        scale_fill_manual(values=major_trajectory_color_plate) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
        ggsave(paste0(name, ".png"),
               dpi = 300,
               height  = 15, 
               width = 12))
}

### horizontal plot, without legend
setwd("~/HollandLabShared/Nayanga/Meningioma_Eric/Analysis33_mouse_embryonic/Cluster_plots_woLegends/")
for(i in 1:10){
    name = gene_module_list[i]
    print(name)
    
    df = dat[,c("major_cell_cluster", "celltype")]
    df$score = as.vector(dat[[name]])
    
    df_x = df %>% group_by(major_cell_cluster, celltype) %>% summarise(median_score = median(score))
    df_x = df_x %>%
        group_by(major_cell_cluster) %>%
        arrange(median_score, .by_group=TRUE)
    
    df$celltype = factor(df$celltype, levels = as.vector(df_x$celltype))
    
    try(ggplot(df, aes(celltype, score, fill = major_cell_cluster)) + geom_boxplot(outlier.shape = NA) + 
            labs(x="celltype", y="gene_module_score", title=name) +
            theme_classic(base_size = 12) +
            scale_fill_manual(values=major_trajectory_color_plate) +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.text.x = element_blank(), 
                  axis.text.y = element_text(color="black", size = 12)) +
            ggsave(paste0(name, ".pdf"),
                   dpi = 300,
                   height  = 6, 
                   width = 30))
}

#color="black", angle = 90, hjust = 1, size = 17
#            theme(legend.position = "none") +

## mean of each cluster in one plot
dat2 <- dat[, c(1:7,11)]
datmean <- dat2 %>% group_by(major_cell_cluster) %>%
  summarise(across(everything(), mean),
            .groups = 'drop') %>% 
  as.data.frame()

datmeanlong <- melt(datmean)

legend_pt_size = 4
axis_text_size = 25
axis_title_size = 25
legend_text_size = 15
spacing = 1
chosen_margin = c(0.5,1,0.5,1) #top,right,bottom,left

theme_legend_em <- theme_classic()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size, face = "bold"),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "right",
        legend.justification = "left",
        legend.title = element_blank())

ggplot(datmeanlong, aes(x=variable, y=value, color = major_cell_cluster)) + geom_point(size = 2) +
  labs(x=" ", y="gene module score") +
  theme_legend_em +
  scale_color_manual(values=major_trajectory_color_plate) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color="black", size = 15 , angle = 90), axis.text.y = element_text(color="black" , size = 15), 
        axis.title = element_text(size = 15, face = "bold"))

setwd("~/HollandLabShared/Nayanga/Meningioma_Eric/Analysis33_mouse_embryonic/")
ggsave("Clusters_embryocelltype_mean_leg.pdf", width = 15, height = 7)

ggplot(datmeanlong, aes(x=variable, y=value, color = major_cell_cluster)) + geom_point(size = 2) +
  labs(x=" ", y="gene module score") +
  theme_legend_em +
  scale_color_manual(values=major_trajectory_color_plate) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color="black", size = 15 , angle = 90), axis.text.y = element_text(color="black" , size = 15), 
        axis.title = element_text(size = 15, face = "bold")) +
  theme(legend.position = "none")
ggsave("Clusters_embryocelltype_mean_noleg.pdf", width = 10, height = 7)

ggplot(datmeanlong, aes(x=variable, y=value, color = major_cell_cluster)) + geom_point(size = 2) +
  geom_text(aes(label=ifelse(value > 0.05, as.character(major_cell_cluster), '')), hjust=0, vjust=0) +
  labs(x="celltype", y="gene_module_score") +
  theme_classic(base_size = 12) +
  scale_color_manual(values=major_trajectory_color_plate) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

library(ggrepel)
ggplot(datmeanlong, aes(x=variable, y=value, color = major_cell_cluster)) + geom_point(size = 2) +
  labs(x="celltype", y="gene_module_score") +
  theme_classic(base_size = 12) +
  scale_color_manual(values=major_trajectory_color_plate) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
  geom_label_repel(data = subset(datmeanlong, value > 0.05),
                   aes(label = major_cell_cluster, size = NULL, color = NULL),
                   segment.size = 0.1, 
                   segment.color = "black", 
                   direction = "y")

## median - error in code
dat2 <- dat[, c(1:7,11)]
datmedian <- dat2 %>% group_by(major_cell_cluster) %>%
  summarise_all(median,
            .groups = 'drop') %>% 
  as.data.frame()

datmedianlong <- melt(datmedian)
ggplot(datmedianlong, aes(x=variable, y=value, color = major_cell_cluster)) + geom_point() +
  labs(x="celltype", y="gene_module_score") +
  theme_classic(base_size = 10) +
  scale_color_manual(values=major_trajectory_color_plate) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
##


##
#identifying significant cell types of each cluster

datsig <- datmean %>% filter(ClusterA < 0.05 | ClusterB < 0.05 | ClusterC < 0.05 |ClusterD < 0.05 | ClusterE < 0.05
                             | ClusterF < 0.05 | ClusterG < 0.05)

datsigA <- datmean %>% filter(ClusterA > 0.05)
datsigB <- datmean %>% filter(ClusterB > 0.05)
datsigC <- datmean %>% filter(ClusterC > 0.05)
datsigD <- datmean %>% filter(ClusterD > 0.05)
datsigE <- datmean %>% filter(ClusterE > 0.05)
datsigF <- datmean %>% filter(ClusterF > 0.05)
datsigG <- datmean %>% filter(ClusterG > 0.05)


# group dat2 by cluster and major_cell_cluster
#cluster A, by cell cluster
datA_cardio <- dat2[, c(1,8)] %>% subset(., dat2$major_cell_cluster == "Cardiomyocytes")
datA_R <- dat2[, c(1,8)] %>% subset(., !(dat2$major_cell_cluster == "Cardiomyocytes" | 
                                           dat2$major_cell_cluster == "Muscle_cells" )) %>% 
  filter(ClusterA < 0.05)
datA_t <- t.test(datA_cardio$ClusterA, datA_R$ClusterA, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)
#
datA_muscle <- dat2[, c(1,8)] %>% subset(., dat2$major_cell_cluster == "Muscle_cells")
datA_R <- dat2[, c(1,8)] %>% subset(., !(dat2$major_cell_cluster == "Cardiomyocytes" | 
                                           dat2$major_cell_cluster == "Muscle_cells" )) %>% 
  filter(ClusterA < 0.05)
datA_t <- t.test(datA_muscle$ClusterA, datA_R$ClusterA, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

#
datB_wbc <- dat2[, c(2,8)] %>% subset(., dat2$major_cell_cluster == "White_blood_cells")
datB_R <- dat2[, c(2,8)] %>% subset(., !(dat2$major_cell_cluster == "White_blood_cells")) %>% 
  filter(ClusterB < 0.05)
datB_t <- t.test(datB_wbc$ClusterB, datB_R$ClusterB, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)
#

#
datB_mast <- dat2[, c(2,8)] %>% subset(., dat2$major_cell_cluster == "Mast_cells")
datB_R <- dat2[, c(2,8)] %>% subset(., !(dat2$major_cell_cluster == "Mast_cells" | dat2$major_cell_cluster == "White_blood_cells")) %>% 
  filter(ClusterB < 0.05)
datB_t <- t.test(datB_mast$ClusterB, datB_R$ClusterB, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)
datB_t

#
datB_mast <- dat2[, c(2,8)] %>% subset(., dat2$major_cell_cluster == "Mast_cells")
datB_R <- dat2[, c(2,8)] %>% subset(., !(dat2$major_cell_cluster == "Mast_cells" )) %>% 
  filter(ClusterB < 0.05)
datB_t <- t.test(datB_mast$ClusterB, datB_R$ClusterB, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)
datB_t
#
datB_mast <- dat2[, c(2,8)] %>% subset(., dat2$major_cell_cluster == "Mesoderm")
datB_R <- dat2[, c(2,8)] %>% subset(., !(dat2$major_cell_cluster == "Mesorderm")) %>% 
  filter(ClusterB < 0.05)
datB_t <- t.test(datB_mast$ClusterB, datB_R$ClusterB, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)
datB_t

##

datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "CNS_neurons")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "CNS_neurons" |
                                           dat2$major_cell_cluster == "Intermediate_neuronal_progenitors")) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "CNS_neurons")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "CNS_neurons"
                                           )) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t
#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Oligodendrocytes")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Oligodendrocytes"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Olfactory_sensory_neurons")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Olfactory_sensory_neurons"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Neural_crest_PNS_neurons")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Neural_crest_PNS_neurons"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Lung_and_airway")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Lung_and_airway"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Eye_and_other")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Eye_and_other"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Ependymal_cells")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Ependymal_cells"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Neuroectoderm_and_glia")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Neuroectoderm_and_glia"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Lung_and_airway")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Lung_and_airway"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Epithelium")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Epithelium"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datC_intneu <- dat2[, c(3,8)] %>% subset(., dat2$major_cell_cluster == "Adipocytes")
datC_R <- dat2[, c(3,8)] %>% subset(., !(dat2$major_cell_cluster == "Adipocytes"
)) %>% 
  filter(ClusterC < 0.05)
datB_t <- t.test(datC_intneu$ClusterC, datC_R$ClusterC, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t
#
datD_de <- dat2[, c(4,8)] %>% subset(., dat2$major_cell_cluster == "Definitive_erythroid")
datD_R <- dat2[, c(4,8)] %>% subset(., !(dat2$major_cell_cluster == "Definitive_erythroid")) %>% 
  filter(ClusterD < 0.05)
datB_t <- t.test(datD_de$ClusterD, datD_R$ClusterD, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)
#

datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Epithelium")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Epithelium" |
                                           dat2$major_cell_cluster == "Intestine" |
                                           dat2$major_cell_cluster == "Lung_and_airway")) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Lung_and_airway")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Lung_and_airway" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t
#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Intestine")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Intestine" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t
#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Epithelium")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Epithelium" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Hepatocytes")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Hepatocytes" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Eye_and_other")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Eye_and_othe" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Intestine")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Intestine" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Endothelium")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Endothelium" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Neural_crest_PNS_neurons")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Neural_crest_PNS_neurons" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t #0.0001941

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Mast_cells")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Mast_cells" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t #2.2e-16

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Neural_crest_PNS_glia")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Neural_crest_PNS_glia" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Eye_and_other")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Eye_and_other" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(5,8)] %>% subset(., dat2$major_cell_cluster == "Ependymal_cells")
datC_R <- dat2[, c(5,8)] %>% subset(., !(dat2$major_cell_cluster == "Ependymal_cells" )) %>% 
  filter(ClusterE < 0.05)
datB_t <- t.test(datE_e$ClusterE, datC_R$ClusterE, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t
#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Mesoderm")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Mesoderm" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Intermediate_neuronal_progenitors")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Intermediate_neuronal_progenitors" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Neural_crest_PNS_neurons")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Neural_crest_PNS_neurons" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "B_cells")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "B_cells" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "White_blood_cells")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "White_blood_cells" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Adipocytes")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Adipocytes" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "T_cells")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "T_cells" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t
#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Oligodendrocytes")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Oligodendrocytes" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Olfactory_sensory_neurons")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Olfactory_sensory_neurons" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Neuroectoderm_and_glia")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Neuroectoderm_and_glia" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Neural_crest_PNS_glia")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Neural_crest_PNS_glia" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Muscle_cells")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Muscle_cells" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t
#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Mesoderm")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Mesoderm" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Mast_cells")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Mast_cells" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Intestine")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Intestine" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Eye_and_other")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Eye_and_other" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t
#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Epithelium")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Epithelium" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "CNS_neurons")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "CNS_neurons" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t

#
datE_e <- dat2[, c(6,8)] %>% subset(., dat2$major_cell_cluster == "Cardiomyocytes")
datC_R <- dat2[, c(6,8)] %>% subset(., !(dat2$major_cell_cluster == "Cardiomyocytes" )) %>% 
  filter(ClusterF < 0.05)
datB_t <- t.test(datE_e$ClusterF, datC_R$ClusterF, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

datB_t
# cluster G

datE_neu <- dat2[, c(7,8)] %>% subset(., dat2$major_cell_cluster == "Neural_crest_PNS_neurons")
datC_R <- dat2[, c(7,8)] %>% subset(., !(dat2$major_cell_cluster == "Neural_crest_PNS_neurons" |
                                           dat2$major_cell_cluster == "CNS_neurons" |
                                           dat2$major_cell_cluster == "Intermediate_neuronal_progenitors")) %>% 
  filter(ClusterG < 0.05)
datB_t <- t.test(datE_neu$ClusterG, datC_R$ClusterG, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

#

datE_neu <- dat2[, c(7,8)] %>% subset(., dat2$major_cell_cluster == "Neural_crest_PNS_glia")
datC_R <- dat2[, c(7,8)] %>% subset(., !(dat2$major_cell_cluster == "Neural_crest_PNS_glia")) %>% 
  filter(ClusterG < 0.05)
datB_t <- t.test(datE_neu$ClusterG, datC_R$ClusterG, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)
datB_t

#

datE_neu <- dat2[, c(7,8)] %>% subset(., dat2$major_cell_cluster == "Mesoderm")
datC_R <- dat2[, c(7,8)] %>% subset(., !(dat2$major_cell_cluster == "Mesoderm")) %>% 
  filter(ClusterG < 0.05)
datB_t <- t.test(datE_neu$ClusterG, datC_R$ClusterG, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)
datB_t

#
datE_neu <- dat2[, c(7,8)] %>% subset(., dat2$major_cell_cluster == "White_blood_cells")
datC_R <- dat2[, c(7,8)] %>% subset(., !(dat2$major_cell_cluster == "White_blood_cells")) %>% 
  filter(ClusterG < 0.05)
datB_t <- t.test(datE_neu$ClusterG, datC_R$ClusterG, alternative = "greater", 

                 
                                  paired = FALSE, conf.level = 0.95)
datB_t

#
datE_neu <- dat2[, c(7,8)] %>% subset(., dat2$major_cell_cluster == "Oligodendrocytes")
datC_R <- dat2[, c(7,8)] %>% subset(., !(dat2$major_cell_cluster == "Oligodendrocytes")) %>% 
  filter(ClusterG < 0.05)
datB_t <- t.test(datE_neu$ClusterG, datC_R$ClusterG, alternative = "greater", 
                 
                 
                 paired = FALSE, conf.level = 0.95)
datB_t

#
datE_neu <- dat2[, c(7,8)] %>% subset(., dat2$major_cell_cluster == "Neuroectoderm_and_glia")
datC_R <- dat2[, c(7,8)] %>% subset(., !(dat2$major_cell_cluster == "Neuroectoderm_and_glia")) %>% 
  filter(ClusterG < 0.05)
datB_t <- t.test(datE_neu$ClusterG, datC_R$ClusterG, alternative = "greater", 
                 
                 
                 paired = FALSE, conf.level = 0.95)
datB_t

#########
#cluster A, by cell cluster
datA_cardio <- dat2[, c(1,8)] %>% subset(., dat2$major_cell_cluster == "Cardiomyocytes")
datA_R <- dat2[, c(1,8)] %>% subset(., !(dat2$major_cell_cluster == "Cardiomyocytes" | 
                                           dat2$major_cell_cluster == "Muscle_cells" )) %>% 
  filter(ClusterA < 0.05)
datA_t <- t.test(datA_cardio$ClusterA, datA_R$ClusterA, alternative = "greater", 
                 paired = FALSE, conf.level = 0.95)

#####
dat3 <- dat[, c(1:7,11,12)]
write.csv(dat3, "~/HollandLabShared/Nayanga/Meningioma_Eric/Analysis33_mouse_embryonic/Supplmental_Table_MouseEmbryonicCells_Scores.csv")
