install.packages("gridExtra")
install.packages("ggeffects")
install.packages("cowplot")
library(tidyverse)
library(gridExtra)
library(ggeffects)
library(cowplot)


M82_rMATs_anno_all_gene <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno")) %>%
  mutate(gene_size = end - str)

M82_rMATs_anno_all_exon_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_exon_order.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order")) %>%
  mutate(size = end - str)

quantile(M82_rMATs_anno_all_gene$size, seq(0, 1, 0.025))

M82_rMATs_anno_all_gene %>%
  filter(gene_size < 13000) %>%
  ggplot() +
  #geom_histogram(aes(size), bins = 50) +
  geom_histogram(aes(size/1000), bins = 100) +
  scale_x_continuous(breaks = seq(0, max(M82_rMATs_anno_all_gene$size)/1000, 2))

max(M82_rMATs_anno_all_gene$size)

M82_rMATs_anno_all_gene_light <- M82_rMATs_anno_all_gene %>%
  select(anno, gene_size)

M82_rMATs_anno_all_exon_order %>%
  left_join(M82_rMATs_anno_all_gene_light) %>%
  filter(gene_size >= 5000 & gene_size <= 7000) %>%
  group_by(anno) %>%
  summarise(count = n()) %>%
  group_by(count) %>%
  summarise(exon_count = n()) %>%
  mutate(sum_exon_count = sum(exon_count))

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


#combine all sig and non-sig AS with the information of intersected exons or introns
#For RI, A5S and A3S, I only imported intersected introns, and for SE, I only imported intersected exons
AS_sig = c("AS_sig_FDR05", "AS_control")
comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
AS_focus = c("A5S_intersected_M82_rMATs_intron", "A3S_intersected_M82_rMATs_intron", "RI_intersected_M82_rMATs_intron", "SE_intersected_M82_rMATs_exon")

all_AS_intersected_feature_raw <- tibble()

for (i in seq_along(AS_sig)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(AS_focus)) {
      path_AS_mark <- str_c("./data/Intersected_AS_TPM_q05_M82_anno_rMATs/", AS_sig[[i]], "_TPM_q05_", comp[[j]], "_", AS_focus[[k]], "_order.bed")
      
      a <- read_delim(path_AS_mark, delim = "\t", 
                 col_names = c("chr", "str", "end", "comp", "AS", "PI", "chr_2", "str_2", "end_2", "strand", "feature", "source", "anno", "order")) %>%
        select(feature, anno, order, comp, AS, PI) %>%
        mutate(event = AS_sig[[i]]) %>%
        mutate(event = if_else(str_detect(event, "sig")==TRUE, "sig", "non-sig")) %>%
        distinct()
      
      all_AS_intersected_feature_raw <- bind_rows(all_AS_intersected_feature_raw, a)
    
    }
  }
}

all_AS_intersected_feature_raw %>%
  group_by(event) %>%
  summarise()

all_AS_intersected_feature_raw %>%
  distinct()
 
gene_having_AS <- all_AS_intersected_feature_raw %>%
  select(feature, anno, event) %>%
  arrange(anno) %>%
  distinct()

M82_rMATs_anno_all_exon_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_exon_order.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order"))
M82_rMATs_anno_all_intron_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_intron_order.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order"))

gene_having_AS_exon_info <- M82_rMATs_anno_all_exon_order %>%
  arrange(chr, str) %>%
  left_join(gene_having_AS) %>%
  drop_na() %>%
  select(-event)

gene_having_AS_intron_info <- M82_rMATs_anno_all_intron_order %>%
  arrange(chr, str) %>%
  left_join(gene_having_AS) %>%
  drop_na() %>%
  select(-event)
  
gene_having_AS_exon_info %>%
  left_join(all_AS_intersected_feature_raw, by = c("feature", "anno", "order")) %>%
  View()
  #drop_na() %>%
  group_by(comp, PI, event) %>%
  summarise(count = n()) %>%
  View()
  
#This is to check the median size of genes containing AS sig or non sig
read_delim("./data/Intersected_AS_TPM_q05_M82_anno_rMATs/AS_control_TPM_q05_HS1_HS0_SE_intersected_M82_rMATs_exon_order.bed", delim = "\t", 
           col_names = c("chr", "str", "end", "comp", "AS", "PI", "chr_2", "str_2", "end_2", "strand", "feature", "source", "anno", "order")) %>%
  left_join(M82_rMATs_anno_all_gene_light) %>%
  #filter(order > 2) %>%
  group_by(PI) %>%
  summarise(gene_size_median = median(gene_size)) %>%
  View()
