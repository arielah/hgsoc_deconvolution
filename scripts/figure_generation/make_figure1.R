suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path
figure_path <- params$figure_path

theme_set(theme_bw() +
              theme(text = element_text(size = 14),
                    strip.background = element_rect(colour = "black",
                                                    fill = "white"),
                    plot.title = element_text(hjust = 0.5)
              )
)

generate_barplot <- function(bulk_set, sc_set, title_tmp, taglet) {
  
  results <- fread(paste(local_data_path, "/deconvolution_output/",
                         bulk_set, "_", sc_set, "_bayesprism_results.tsv",
                         sep = ""))
  results_melt <- melt(results)
  
  g <- ggplot(results_melt, mapping = aes(x=variable, y=value, fill=cell_type, color=cell_type)) +
    geom_bar(stat = "identity") +
    labs(x="Sample", y="Proportion", fill="Cell type", colour="Cell type",
         tag = taglet, title = title_tmp) +
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank(), #remove y axis ticks
          panel.background = element_blank(),
          axis.title.x = element_text(vjust = 1)
    )
  
  g
}

pA <- generate_barplot("TCGA", "fibro", "TCGA RNA-seq", "A")
pB <- generate_barplot("microarray", "fibro", "TCGA Microarray", "B")
pC <- generate_barplot("tothill", "fibro", "Tothill", "C")
pD <- generate_barplot("TCGA", "vazquez", "TCGA RNA-seq", "D")
pE <- generate_barplot("microarray", "vazquez", "TCGA Microarray", "E")
pF <- generate_barplot("tothill", "vazquez", "Tothill", "F")

top <- pA + pB + pC + plot_layout(guides = "collect")

bottom <- pD + pE + pF + plot_layout(guides = "collect")

pdf(paste(figure_path, "figure1.pdf", sep = "/"), width = 18, height = 12, family = "sans")
top / bottom
dev.off()
