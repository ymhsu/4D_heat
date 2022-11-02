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
  select(-event) %>%
  distinct()

gene_having_AS_intron_info <- M82_rMATs_anno_all_intron_order %>%
  arrange(chr, str) %>%
  left_join(gene_having_AS) %>%
  drop_na() %>%
  select(-event) %>%
  distinct()

#SE 
gene_having_AS_exon_info %>%
  group_by(anno) %>%
  summarise(count = n()) %>%
  group_by(count) %>%
  summarise(count_sum = sum(count)) %>%
  View()

list_for_ML_data_f <- function(data) {
  data %>%
    #remove genes with only two exons
    group_by(anno) %>%
    mutate(exon_count = n()) %>%
    filter(exon_count > 2) %>%
    select(-exon_count) %>%
    ungroup() %>%
    left_join(all_AS_intersected_feature_raw, by = c("feature", "anno", "order")) %>%
    replace_na(list(comp = "no", AS = "no", PI = "no", event = "no")) %>%
    #remove sig event in the first or second exon/intron (promoter effects)
    mutate(promoter_filter_text = str_c(order, "_", PI)) %>%
    mutate(promoter_filter = if_else(promoter_filter_text == "1_no" | promoter_filter_text == "2_no", 1, 0)) %>%
    group_by(anno) %>%
    mutate(sum_promoter_filter = sum(promoter_filter)) %>%
    filter(sum_promoter_filter == 2) %>%
    select(-promoter_filter_text, -promoter_filter, -sum_promoter_filter) %>%
    #remove genes without having PI_H or PI_L
    mutate(PI_filter = if_else(PI == "PI_H" | PI == "PI_L", 1, 0)) %>%
    group_by(anno) %>%
    mutate(sum_PI_filter = sum(PI_filter)) %>%
    filter(sum_PI_filter != 0) %>%
    select(-PI_filter, -sum_PI_filter) %>%
    #check how many genes with more exons having AS than exons lacking AS 
    #select(anno, order, AS) %>%
    #distinct() %>%
    #group_by(anno) %>%
    #mutate(count_no_AS = if_else(AS == "no", 1, 0), count_AS = if_else(AS != "no", 1, 0)) %>%
    #mutate(sum_no_AS = sum(count_no_AS), sum_AS = sum(count_AS)) %>%
    #filter(sum_no_AS < sum_AS) %>%
    #remove genes without any sig events
    mutate(sig_filter = if_else(event == "sig", 1, 0)) %>%
    mutate(sum_sig_filter = sum(sig_filter)) %>%
    filter(sum_sig_filter != 0) %>%
    select(-sig_filter, -sum_sig_filter) %>%
    split(.$anno)
}

exon_data_for_ML_list <- list_for_ML_data_f(gene_having_AS_exon_info)
intron_data_for_ML_list <- list_for_ML_data_f(gene_having_AS_intron_info)
  

make_AS_ML_data_f <- function(data, a, b){ 
#a is one of three comparisons, b is the type of PI
AS_sig_m <- data %>%
  filter(comp == a[1] & PI == b[1] & event == "sig")

AS_noAS_m <- data %>%
  filter(event == "no") 

AS_noAS_m_f <- head(AS_noAS_m[sample(1:nrow(AS_noAS_m)),], n = nrow(AS_sig_m)) 

bind_rows(AS_sig_m, AS_noAS_m_f)
}

comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
PI_type = c("PI_H", "PI_L")

exon_data_for_ML_comb_l <- vector("list", length = length(comp))
intron_data_for_ML_comb_l <- vector("list", length = length(comp))

for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    set.seed(200)
    exon_data_for_ML_comb_l[[i]][[j]] <- exon_data_for_ML_list %>%
      map(. %>% make_AS_ML_data_f(., comp[[i]], PI_type[[j]])) %>%
      bind_rows()
    
    intron_data_for_ML_comb_l[[i]][[j]] <- intron_data_for_ML_list %>%
      map(. %>% make_AS_ML_data_f(., comp[[i]], PI_type[[j]])) %>%
      bind_rows()
  }
}

for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    exon_data_l_final <- bind_rows(exon_data_for_ML_comb_l[[i]][[j]], exon_data_for_ML_comb_l[[i]][[j]]) %>%
      ungroup() %>%
      mutate(side = rep(c("five", "three"), each = nrow(exon_data_for_ML_comb_l[[i]][[j]]))) %>% 
      mutate(str_n = if_else(strand == "+" & side == "five", str - 100,
                             if_else(strand == "+" & side == "three", end,
                                     if_else(strand == "-" & side == "five", end, str - 100)))) %>%
      mutate(end_n = if_else(strand == "+" & side == "five", str,
                             if_else(strand == "+" & side == "three", end + 100,
                                     if_else(strand == "-" & side == "five", end + 100, str)))) %>%
      select(chr, str = str_n, end = end_n, strand, feature, source, anno, order, comp, AS, PI, event, seg_side = side) %>%
      arrange(chr, anno, str) %>%
      split(.$event)
    
    intro_data_l_final <- bind_rows(intron_data_for_ML_comb_l[[i]][[j]], intron_data_for_ML_comb_l[[i]][[j]]) %>%
      ungroup() %>%
      mutate(side = rep(c("five", "three"), each = nrow(intron_data_for_ML_comb_l[[i]][[j]]))) %>% 
      mutate(str_n = if_else(strand == "+" & side == "five", str - 100,
                             if_else(strand == "+" & side == "three", end,
                                     if_else(strand == "-" & side == "five", end, str - 100)))) %>%
      mutate(end_n = if_else(strand == "+" & side == "five", str,
                             if_else(strand == "+" & side == "three", end + 100,
                                     if_else(strand == "-" & side == "five", end + 100, str)))) %>%
      select(chr, str = str_n, end = end_n, strand, feature, source, anno, order, comp, AS, PI, event, seg_side = side) %>%
      arrange(chr, anno, str) %>%
      split(.$event)
    
    write_delim(exon_data_l_final$sig, str_c("./data/AS_bed_for_ML/exon_for_ML_", comp[[i]], "_", PI_type[[j]], "_", "sig.bed"), col_names = FALSE, delim = "\t")
    write_delim(exon_data_l_final$no, str_c("./data/AS_bed_for_ML/exon_for_ML_", comp[[i]], "_", PI_type[[j]], "_", "non_sig.bed"), col_names = FALSE, delim = "\t")
    
    write_delim(intro_data_l_final$sig, str_c("./data/AS_bed_for_ML/intron_for_ML_", comp[[i]], "_", PI_type[[j]], "_", "sig.bed"), col_names = FALSE, delim = "\t")
    write_delim(intro_data_l_final$no, str_c("./data/AS_bed_for_ML/intron_for_ML_", comp[[i]], "_", PI_type[[j]], "_", "non_sig.bed"), col_names = FALSE, delim = "\t")

  }
}


names(test_sig_nosig)
test_sig_nosig$sig
AS_bed_for_ML

make_AS_ML_data_f(test_exon_data_for_ML_temp_l[[2]], c("HS1_HS0"), c("PI_H"))

test_exon_data_for_ML_temp_l_small <- test_exon_data_for_ML_temp_l[1:10]

test_exon_data_for_ML_temp_l_small %>%
  map(. %>% make_AS_ML_data_f(., c("HS1_HS0"), c("PI_H"))) %>%
  bind_rows()

#This is to check the median size of genes containing AS sig or non sig
read_delim("./data/Intersected_AS_TPM_q05_M82_anno_rMATs/AS_control_TPM_q05_HS1_HS0_SE_intersected_M82_rMATs_exon_order.bed", delim = "\t", 
           col_names = c("chr", "str", "end", "comp", "AS", "PI", "chr_2", "str_2", "end_2", "strand", "feature", "source", "anno", "order")) %>%
  left_join(M82_rMATs_anno_all_gene_light) %>%
  #filter(order > 2) %>%
  group_by(PI) %>%
  summarise(gene_size_median = median(gene_size)) %>%
  View()
