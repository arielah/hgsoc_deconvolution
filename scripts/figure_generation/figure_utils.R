library(scales)
library(RColorBrewer)

alpha_value <- 0.5
options(encoding = "UTF-8")

theme_set(theme_bw() +
            theme(text = element_text(size = 14),
                  strip.background = element_rect(colour = NA,
                                                  fill = "white"),
                  plot.title = element_text(hjust = 0.5)
            )
)

colors_celltypes <- c("B cells" = "#B15928",
                      "DC"="#FFFF99", 
                      "Endothelial cells"="#6A3D9A",
                      "Epithelial cells"="#CAB2D6", 
                      "Fibroblasts"="#FF7F00", 
                      "ILC"="#FDBF6F",
                      "Macrophages"="#E31A1C", 
                      "Mast cells"="#FB9A99",
                      "Monocytes"="#33A02C",
                      "NK cells"="#B2DF8A",
                      "pDC"="#1F78B4",
                      "Plasma cells"="#A6CEE3",
                      "T cells"="#808080",
                      "Myeloid cells"="#EF5A5A") #paired

colors_subtypes <- c("Mesenchymal" = "#E6AB02", "Proliferative" = "#7570B3", 
                     "Immunoreactive" = "#E7298A", "Differentiated" = "#66A61E",
                     "1" = "#A6761D", "2" = "#1B9E77", "3" = "#D95F02") #dark2

colors_bulktypes <- c("TCGA Microarray" = "#00468B", "TCGA RNA-seq" = "#925E9F",
                      "Tothill" = "#AD002A") #lancet

heatmap_scale_2d <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(20)
