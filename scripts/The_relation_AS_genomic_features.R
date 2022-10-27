install.packages("gridExtra")
install.packages("ggeffects")
install.packages("cowplot")
library(tidyverse)
library(gridExtra)
library(ggeffects)
library(cowplot)


M82_rMATs_anno_all_gene <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno")) %>%
  mutate(size = end - str)

quantile(M82_rMATs_anno_all_gene$size, seq(0, 1, 0.025))


M82_rMATs_anno_all_gene %>%
  filter(size < 13000) %>%
  ggplot() +
  #geom_histogram(aes(size), bins = 50) +
  geom_histogram(aes(size), bins = 100) +
  scale_x_continuous(breaks = seq(0, max(M82_rMATs_anno_all_gene$size), 2000))

max(M82_rMATs_anno_all_gene$size)

AS_sig = c("AS_sig_FDR05", "AS_control")
comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
AS = c("A5S", "A3S", "RI", "SE")

AS_intersected_gene_plot_l <- vector(mode = "list", length = length(AS_sig)*length(comp)*length(AS))

for (i in seq_along(AS_sig)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(AS)) {
      path <- str_c("./data/Intersected_AS_TPM_q05_M82_anno_rMATs/", AS_sig[[i]], "_TPM_q05_", comp[[j]], "_", AS[[k]], "_intersected_M82_rMATs_gene.bed")
      path_2 <- str_c("./analysis_output/", AS_sig[[i]], "_TPM_q05_", comp[[j]], "_", AS[[k]], "_gene_size_p.jpeg")
      type_p <- str_c(AS_sig[[i]], "_", comp[[j]], "_", AS[[k]]) 
      print(path)
      
      a <- read_delim(path, delim = "\t", 
                 col_names = c("chr", "str", "end", "comp", "AS", "PI", "chr_2", "str_2", "end_2", "strand", "feature", "source", "anno")) %>%
        #summarise(count = n())
        mutate(size = (end_2 - str_2)/1000) %>%
        mutate(type = type_p) %>%
        ggplot() +
        geom_histogram(aes(size), bins = 100) +
        scale_x_continuous(breaks = seq(0, max(M82_rMATs_anno_all_gene$size)/1000, 2)) +
        facet_wrap(~type) +
        labs(x = "size (kb)", y = "counts") +
        theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
              legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
              axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))
      ggsave(path_2, a, width = 400, height = 320, units = c("mm"), dpi = 320)
      #for (l in seq_along(AS_intersected_gene_plot_l)) {
        #AS_intersected_gene_plot_l[[l]] <- a
        #ggsave(path_2, AS_intersected_gene_plot_l[[l]], width = 400, height = 320, units = c("mm"), dpi = 320)
        
      #}
    }
  }
}

AS_intersected_gene_plot_l[[1]]

ggsave("./analysis_output/test.jpeg", AS_intersected_gene_plot_l[[1]], width = 400, height = 320, units = c("mm"), dpi = 320)


test = "x"
print(test)

read_delim("./data/Intersected_AS_TPM_q05_M82_anno_rMATs/AS_sig_FDR05_TPM_q05_HS6_HS0_SE_intersected_M82_rMATs_gene.bed", delim = "\t", 
           col_names = c("chr", "str", "end", "comp", "AS", "PI", "chr_2", "str_2", "end_2", "strand", "feature", "source", "anno")) %>%
  #summarise(count = n())
  mutate(size = (end_2 - str_2)/1000) %>%
  mutate(type = str_c("HS6_HS0", "_", "RI", "_gene_size (kb)")) %>%
  ggplot() +
  geom_histogram(aes(size), bins = 100) +
  scale_x_continuous(breaks = seq(0, max(M82_rMATs_anno_all_gene$size)/1000, 2)) +
  facet_wrap(~type) +
  labs(x = "size (kb)", y = "counts") +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))

                     