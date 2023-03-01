Sys.getlocale()
Sys.setenv(LANG = "en_US.UTF-8")
Packages <- c("plyranges", "tidyverse", "doParallel", "foreach", "nullranges", "caTools", "ggpubr", "fs", "HelloRanges")
lapply(Packages, library, character.only = TRUE)

expression_level_tomato_gene <- read_csv("./data/expression_level.csv")

expression_level_average_tomato_gene <- 
expression_level_tomato_gene %>%
  dplyr::select(-...1) %>%
  filter(str_detect(ID, "Soly")) %>%
  mutate(ID = str_remove(ID, "gene:")) %>%
  mutate(level_0H = (RNA_M82_HS_0h_rep1_S4+RNA_M82_HS_0h_rep2_S5)/2,
         level_1H = (RNA_M82_HS_1h_rep1_S6+RNA_M82_HS_1h_rep2_S7)/2,
         level_6H = (RNA_M82_HS_6h_rep1_S8+RNA_M82_HS_6h_rep2_S9)/2) %>%
  dplyr::select(ID, level_0H, level_1H, level_6H) 

quantile_expression_levels <- 
function(a){  
expression_level_average_tomato_gene %>%
    dplyr::select(gene_name = ID, continuous_levels = contains(a[1])) %>%
    filter(continuous_levels > 0) %>%
    mutate(expression_discrete_level = if_else(continuous_levels <= quantile(continuous_levels, 0.25), "q1", 
                                    if_else(continuous_levels > quantile(continuous_levels, 0.25) & continuous_levels <= quantile(continuous_levels, 0.5), "q2",
                                            if_else(continuous_levels > quantile(continuous_levels, 0.75), "q4", "q3"))))
}

list_quantile_expression_levels <- 
str_c("level_", c("0H", "1H", "6H")) %>%
  map(. %>% quantile_expression_levels())

names(list_quantile_expression_levels) <- str_c("level_", c("0H", "1H", "6H"))




#import PSI combined, ctrl (no AS), DAS (sk/Inc in terms of 3 heat comp, 6 in total) as lists
Inc_AS_aft_trt <- c("inc", "sk")
comp_AS_HS0_DAS_heat_trt_label <- c("PSI_all", "ctrl", str_c("HS0_HS1_more_", Inc_AS_aft_trt),
                                    str_c("HS0_HS6_more_", Inc_AS_aft_trt), str_c("HS1_HS6_more_", Inc_AS_aft_trt))
AS_type_focused <- c("RI", "SE")
location_gr <- c("upstream_feature_gr", "feature_gr", "downstream_feature_gr")

list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression <- list()
names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression <- c()

for (i in seq_along(AS_type_focused)) {
  for (j in seq_along(comp_AS_HS0_DAS_heat_trt_label)) {
      
      if(j %in% c(1, 2)){
        
        path_input <- str_c("./data/AS_HS0_bootstrap_samples_DAS_removed/", AS_type_focused[[i]], "_", comp_AS_HS0_DAS_heat_trt_label[[j]], "_feature_gr_bootstrap_DAS_removed_1_100_combined")
        
        name_input <- str_c(AS_type_focused[[i]], "_feature_gr_", comp_AS_HS0_DAS_heat_trt_label[[j]])
        
        names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression <- append(names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression, name_input)
        
        gr_sample <- read_delim(path_input, delim = "\t", col_names = TRUE) %>%
          as_granges()
        
        list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression <- append(list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression, list(gr_sample))
      } else {
        
        path_input <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_focused[[i]], "_feature_gr_", comp_AS_HS0_DAS_heat_trt_label[[j]], "_AS_aft_trt")
        
        name_input <- str_c(AS_type_focused[[i]], "_feature_gr_", comp_AS_HS0_DAS_heat_trt_label[[j]])
        
        names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression <- append(names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression, name_input)
        
        gr_sample <- read_delim(path_input, delim = "\t", col_names = TRUE) %>%
          as_granges()
        
        list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression <- append(list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression, list(gr_sample))
      }
      
      
    
  }
}

#data < list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression[[8]]
#data2 < M82_rMATs_anno_all_intron
#data3 list_quantile_expression_levels$level_0H

group_AS_DAS_genes_by_transcription_level <- 
  function(data, data2, data3) {
    data %>%
      as_tibble() %>%
      left_join(data2) %>%
      mutate(gene_name = str_sub(gene_name, 1, 16)) %>%
      dplyr::select(-source) %>%
      distinct() %>%
      left_join(data3) %>%
      mutate(expression_discrete_level = factor(expression_discrete_level, levels = str_c("q", 1:4))) %>%
      dplyr::select(seqnames, start, end, width, strand, expression_discrete_level) %>%
      split(.$expression_discrete_level) %>%
      map(. %>% as_granges)
    
  }

group_AS_DAS_genes_by_transcription_level(list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression[[1]],
                                          M82_rMATs_anno_all_intron, 
                                          list_quantile_expression_levels$level_0H)



list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H_RI <- 
pmap(list(list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression[1:8], replicate(8, list(M82_rMATs_anno_all_intron)), replicate(8, list(list_quantile_expression_levels$level_0H))),
          group_AS_DAS_genes_by_transcription_level) %>%
  reduce(., append)

list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H_RI <- 
  pmap(list(list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression[1:8], replicate(8, list(M82_rMATs_anno_all_intron)), replicate(8, list(list_quantile_expression_levels$level_1H))),
       group_AS_DAS_genes_by_transcription_level) %>%
  reduce(., append)

list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H_RI <- 
  pmap(list(list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression[1:8], replicate(8, list(M82_rMATs_anno_all_intron)), replicate(8, list(list_quantile_expression_levels$level_6H))),
       group_AS_DAS_genes_by_transcription_level) %>%
  reduce(., append)

list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H_SE <- 
  pmap(list(list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression[9:16], replicate(8, list(M82_rMATs_anno_all_exon)), replicate(8, list(list_quantile_expression_levels$level_0H))),
       group_AS_DAS_genes_by_transcription_level) %>%
  reduce(., append)

list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H_SE <- 
  pmap(list(list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression[9:16], replicate(8, list(M82_rMATs_anno_all_exon)), replicate(8, list(list_quantile_expression_levels$level_1H))),
       group_AS_DAS_genes_by_transcription_level) %>%
  reduce(., append)

list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H_SE <- 
  pmap(list(list_gr_raw_comp_AS_HS0_DAS_heat_trt_refined_expression[9:16], replicate(8, list(M82_rMATs_anno_all_exon)), replicate(8, list(list_quantile_expression_levels$level_6H))),
       group_AS_DAS_genes_by_transcription_level) %>%
  reduce(., append)

append(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H_RI, list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H_SE)

list_M82_feature <- list(M82_rMATs_anno_all_intron, M82_rMATs_anno_all_exon)

list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H <- append(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H_RI, list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H_SE)
list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H <- append(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H_RI, list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H_SE)
list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H <- append(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H_RI, list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H_SE)

names(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H)

names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H <- c()
names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H <- c()
names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H <- c()


AS_type_focused <- c("RI")

#data > list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H

produce_list_name_AS_DAS_expression_level_group <- 
function(a, data){
name_table <- 
tibble(AS_type = rep(AS_type_focused, each = length(data)/length(AS_type_focused)),
       AS_DAS_group = rep(comp_AS_HS0_DAS_heat_trt_label, 2, each = 4),
       heat_trt = a[1],
       quantile = names(data)
       ) %>%
  mutate(names_list_final = str_c(AS_type, "_", AS_DAS_group, "_", heat_trt, "_", quantile))

name_table$names_list_final
}

names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H <- 
produce_list_name_AS_DAS_expression_level_group("0H", list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H)

names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H <- 
  produce_list_name_AS_DAS_expression_level_group("1H", list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H)

names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H <- 
  produce_list_name_AS_DAS_expression_level_group("6H", list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H)


#merge peaks of epigenetic marks with two replicates 
mark_with_more_reps <- c("Pol2", "H3K4me3")
heat_trt <- c("0H", "1H", "6H")

#make function for merging marks with two reps
merge_rep_epimark <- 
  function(a, b){
    path_peak <-
      str_c("./data/peak_epimarks/heat_trt/M82_", a[1], "_", b[1], "_rep", c(1,2), "_p0.05_peaks.bed")
    
    mark_r1 <- 
      read_delim(path_peak[[1]], delim = "\t", col_names = c("chr", "str", "end", "signal"))
    
    mark_r2 <- 
      read_delim(path_peak[[2]], delim = "\t", col_names = c("chr", "str", "end", "signal"))
    
    
    bind_rows(mark_r1, mark_r2) %>%
      arrange(chr, str) %>%
      dplyr::select(seqnames = chr, start = str, end, signal) %>%
      as_granges() %>%
      GenomicRanges::reduce()
  }


#make the output of merged epigenetic marks
for (i in seq_along(heat_trt)) {
  for (j in seq_along(mark_with_more_reps)) {
    
    merge_mark <- merge_rep_epimark(heat_trt[[i]], mark_with_more_reps[[j]])
    
    merge_mark <- merge_mark %>%
      as_tibble() %>%
      dplyr::select(-width, -strand) %>%
      mutate(signal = 0)
    
    write_delim(merge_mark, str_c("./data/peak_epimarks/heat_trt/M82_", heat_trt[[i]], "_", mark_with_more_reps[[j]], "_merged_p0.05_peaks.bed"),
                delim = "\t", col_names = FALSE)
    
  }
}

mark_heat_trt <- c("H3K9ac", "H3K18ac", "H3K27ac", "Pol2", "H3K4me3")
list_gr_mark_heat_trt <- list()
names_list_gr_mark_heat_trt <- c()


#import the files of epigenetic marks under heat treatments
for (i in seq_along(heat_trt)) {
  for (j in seq_along(mark_heat_trt)) {
    
    if(sum(str_detect(mark_with_more_reps, mark_heat_trt[[j]]))==0){
      input_gr_mark <- read_delim(str_c("./data/peak_epimarks/heat_trt/M82_", heat_trt[[i]], "_", mark_heat_trt[[j]], "_rep1_p0.05_peaks.bed"), 
                                  delim = "\t", col_names = c("seqnames", "start", "end", "signal")) %>%
        as_granges()
      
      list_gr_mark_heat_trt <- append(list_gr_mark_heat_trt, list(input_gr_mark))
      
      names_list_gr_mark_heat_trt <- append(names_list_gr_mark_heat_trt, str_c(heat_trt[[i]], "_", mark_heat_trt[[j]]))
      
    } else {
      input_gr_mark <- read_delim(str_c("./data/peak_epimarks/heat_trt/M82_", heat_trt[[i]], "_", mark_heat_trt[[j]], "_merged_p0.05_peaks.bed"), 
                                  delim = "\t", col_names = c("seqnames", "start", "end", "signal")) %>%
        as_granges()
      
      list_gr_mark_heat_trt <- append(list_gr_mark_heat_trt, list(input_gr_mark))
      
      names_list_gr_mark_heat_trt <- append(names_list_gr_mark_heat_trt, str_c(heat_trt[[i]], "_", mark_heat_trt[[j]]))
      
    }
    
  }
}

names(list_gr_mark_heat_trt) <- names_list_gr_mark_heat_trt

list_gr_mark_heat_trt %>%
  map(. %>% as_tibble) %>%
  map(. %>% summarise(median_size = median(end-start)))

#perform coenrichment analysis
comb_pair_epimark_heat_trt_raw <- combs(seq_along(mark_heat_trt), 2) 

comb_pair_epimark_heat_trt_table <- 
  tibble(mark_1 = c(comb_pair_epimark_heat_trt_raw[, 1], c(comb_pair_epimark_heat_trt_raw[, 1]+5), c(comb_pair_epimark_heat_trt_raw[, 1]+10)),
         mark_2 = c(comb_pair_epimark_heat_trt_raw[, 2], c(comb_pair_epimark_heat_trt_raw[, 2]+5), c(comb_pair_epimark_heat_trt_raw[, 2]+10)))

names_comb_pair_epimark_heat_trt_table <- tibble()

list_comb_pair_epimark_heat_trt <- vector("list", length = length(comb_pair_epimark_heat_trt_table$mark_1))

for (i in seq_along(list_comb_pair_epimark_heat_trt)) {
  
  list_comb_pair_epimark_heat_trt[[i]][[1]] <- list_gr_mark_heat_trt[[comb_pair_epimark_heat_trt_table[i, ]$mark_1]]
  list_comb_pair_epimark_heat_trt[[i]][[2]] <- list_gr_mark_heat_trt[[comb_pair_epimark_heat_trt_table[i, ]$mark_2]]
  
  name_mark_table <- 
    tibble(mark_1 = names_list_gr_mark_heat_trt[[comb_pair_epimark_heat_trt_table[i, ]$mark_1]],
           mark_2 = names_list_gr_mark_heat_trt[[comb_pair_epimark_heat_trt_table[i, ]$mark_2]])
  
  names_comb_pair_epimark_heat_trt_table <- bind_rows(names_comb_pair_epimark_heat_trt_table, name_mark_table)
}

#list_gr_raw_comp_AS_HS0_DAS_heat_trt[[1]] data

create_coenrichment_index_two_marks <- 
  function(data, a){
    mark_1_feature_intersection <- 
      join_overlap_intersect(list_comb_pair_epimark_heat_trt[[a[1]]][[1]], data) %>%
      GenomicRanges::reduce()
    
    mark_2_feature_intersection <- 
      join_overlap_intersect(list_comb_pair_epimark_heat_trt[[a[1]]][[2]], data) %>%
      GenomicRanges::reduce()
    
    sum_bp_coenrichment <- 
      join_overlap_intersect(mark_1_feature_intersection, mark_2_feature_intersection) %>%
      GenomicRanges::reduce() %>%
      as_tibble() %>%
      summarise(sum_bp_coenrichment = sum(end-start))
    
    mark1_sum_bp_intersection <- 
      as_tibble(mark_1_feature_intersection) %>%
      summarise(sum_bp_coenrichment = sum(end-start))
    
    mark2_sum_bp_intersection <- 
      as_tibble(mark_2_feature_intersection) %>%
      summarise(sum_bp_coenrichment = sum(end-start))
    
    tibble(
      mark_1 = names_comb_pair_epimark_heat_trt_table[a[1],]$mark_1,
      mark_2 = names_comb_pair_epimark_heat_trt_table[a[1],]$mark_2,
      coenrichment_index = sum_bp_coenrichment$sum_bp_coenrichment/c(mark1_sum_bp_intersection$sum_bp_coenrichment + mark2_sum_bp_intersection$sum_bp_coenrichment - sum_bp_coenrichment$sum_bp_coenrichment)
    )
  }

create_coenrichment_index_two_marks(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H[[17]], 1)
names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H


coenrichment_index_table <- tibble()
list_comb_pair_epimark_heat_trt

for (i in seq_along(list_comb_pair_epimark_heat_trt)) {
  if (i %in% c(1:10)){
  coenrichment_table_raw_0H <- map2(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H, rep(i, length(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H)), create_coenrichment_index_two_marks) %>%
    Reduce(bind_rows, .) %>%
    mutate(feature = names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_0H)
  
  coenrichment_index_table <- bind_rows(coenrichment_index_table, coenrichment_table_raw_0H)
  } else if (i %in% c(11:20)){
    coenrichment_table_raw_1H <- map2(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H, rep(i, length(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H)), create_coenrichment_index_two_marks) %>%
      Reduce(bind_rows, .) %>%
      mutate(feature = names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_1H)
    
    coenrichment_index_table <- bind_rows(coenrichment_index_table, coenrichment_table_raw_1H)  
  } else {
    coenrichment_table_raw_6H <- map2(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H, rep(i, length(list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H)), create_coenrichment_index_two_marks) %>%
      Reduce(bind_rows, .) %>%
      mutate(feature = names_list_gr_raw_comp_AS_HS0_DAS_heat_trt_expression_level_group_6H)
    
    coenrichment_index_table <- bind_rows(coenrichment_index_table, coenrichment_table_raw_6H)
  }
}

write_delim(coenrichment_index_table, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_table_combined_bootstrap_samples_transcription_level_group", delim = "\t", col_names = TRUE) 

coenrichment_index_table_expression_level_group <-
  read_delim("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_table_combined_bootstrap_samples_transcription_level_group",
             delim = "\t", col_names = TRUE)

coenrichment_index_table_expression_level_group_for_plot <- 
coenrichment_index_table_expression_level_group %>%
  mutate(AS_DAS_groups = str_sub(feature, 4, nchar(feature)-6)) %>%
  mutate(AS_heat_trt_expression_group = str_remove(feature, str_c(AS_DAS_groups, "_"))) %>%
  mutate(pair_mark = str_c(mark_1, "_", str_remove(mark_2, ".*_"))) %>%
  mutate(AS_DAS_groups = factor(AS_DAS_groups, levels = comp_AS_HS0_DAS_heat_trt_label)) %>%
  dplyr::select(pair_mark, AS_heat_trt_expression_group, coenrichment_index, AS_DAS_groups)

coenrichment_index_table_expression_level_group_for_plot_key <- 
  coenrichment_index_table_expression_level_group_for_plot %>%
  ungroup() %>%
  dplyr::select(AS_heat_trt_expression_group, pair_mark) %>%
  distinct()


produce_coenrichment_index_hist_expression_level_plots <- function(a, b){
  coenrichment_index_table_expression_level_group_for_plot %>% 
    filter(AS_heat_trt_expression_group == a[1] & pair_mark == b[1]) %>%
    ggplot(aes(AS_DAS_groups, coenrichment_index, fill=AS_DAS_groups)) +
    geom_bar(stat="identity") +
    #geom_errorbar(aes(ymin=coenrichment_index-sd_coenrichment_index, ymax=coenrichment_index+sd_coenrichment_index), width=.2) +
    scale_fill_manual(
      values = c(
        "#29f600",
        "#F6BE00",
        "#fb6a4a",
        "#9e9ac8",
        "#de2d26",
        "#756bb1",
        "#a50f15",
        "#54278f"
      )
    ) +
    facet_grid(AS_heat_trt_expression_group ~ pair_mark) +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 20), legend.text = element_text(size = 9, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 20, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
          axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle = 60, vjust = 0.5))
}

produce_coenrichment_index_hist_expression_level_plots(coenrichment_index_table_expression_level_group_for_plot$AS_heat_trt_expression_group[[4]], 
                                                       coenrichment_index_table_expression_level_group_for_plot$pair_mark[[4]])


coenrichment_index_table_expression_level_group_for_plot_list <- 
  map2(coenrichment_index_table_expression_level_group_for_plot_key$AS_heat_trt_expression_group, 
       coenrichment_index_table_expression_level_group_for_plot_key$pair_mark,
       produce_coenrichment_index_hist_expression_level_plots)

path_coenrichment_index_plots_hist_AS_DAS_expression_level <- 
  str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_analysis_heat_comb_epimark_AS_DAS_histograms_expression_level/", coenrichment_index_table_expression_level_group_for_plot_key$AS_heat_trt_expression_group, "_", coenrichment_index_table_expression_level_group_for_plot_key$pair_mark, "_coenrichment_index_AS_DAS.jpeg")

pwalk(list(path_coenrichment_index_plots_hist_AS_DAS_expression_level, coenrichment_index_table_expression_level_group_for_plot_list), ggsave,
      width = 400, height = 240, units = c("mm"), dpi = 320)
