library(scales)
library(RColorBrewer)

alpha_value <- 0.5


theme_set(theme_bw() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  text = element_text(size = 14),
                  strip.background = element_rect(colour = "black",
                                                  fill = "white"),
                  plot.title = element_text(hjust = 0.5)
            )
)

colors_celltypes <- hue_pal()(13) #ggplot defaults

colors_subtypes <- c("Mesenchymal" = "#E6AB02", "Proliferative" = "#7570B3", "Immunoreactive" = "#E7298A", 
                     "Differentiated" = "#66A61E", "1" = "#A6761D", "2" = "#1B9E77", "3" = "#D95F02")

colors_bulktypes <- c("TCGA Microarray" = "#00468B", "TCGA RNA-seq" = "#925E9F", "Tothill" = "#AD002A") #lancet

heatmap_scale_2d <- c("#252974", "#191c4d",
                      colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(20),
                      "#660018", "#80001E")
