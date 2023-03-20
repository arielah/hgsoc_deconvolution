library(scales)
library(RColorBrewer)

alpha_value <- 0.5
options(encoding = "UTF-8")

theme_set(theme_bw() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  text = element_text(size = 14),
                  strip.background = element_rect(colour = NA,
                                                  fill = "white"),
                  plot.title = element_text(hjust = 0.5)
            )
)

colors_celltypes <- c("B cells" = "#F8766D", "DC"="#E18A00", "Endothelial cells"="#BE9C00",
                      "Epithelial cells"="#8CAB00", "Fibroblasts"="#24B700", "ILC"="#00BE70",
                      "Macrophages"="#00C1AB", "Mast cells"="#00BBDA", "Monocytes"="#00ACFC",
                      "NK cells"="#8B93FF", "pDC"="#D575FE", "Plasma cells"="#F962DD",
                      "T cells"="#FF65AC", "Myeloid cells"="#00B8CF") #ggplot defaults

colors_subtypes <- c("Mesenchymal" = "#E6AB02", "Proliferative" = "#7570B3", "Immunoreactive" = "#E7298A", 
                     "Differentiated" = "#66A61E", "1" = "#A6761D", "2" = "#1B9E77", "3" = "#D95F02")

colors_bulktypes <- c("TCGA Microarray" = "#00468B", "TCGA RNA-seq" = "#925E9F", "Tothill" = "#AD002A") #lancet

heatmap_scale_2d <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(20)
