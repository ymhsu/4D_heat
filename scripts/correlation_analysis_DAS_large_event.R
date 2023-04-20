Sys.getlocale()
Sys.setenv(LANG = "en_US.UTF-8")
Packages <- c("plyranges", "tidyverse", "doParallel", "foreach", "nullranges", "caTools", "ggpubr", "fs", "HelloRanges")
lapply(Packages, library, character.only = TRUE)

#The objective of this analysis is to focus on large DAS events for the relationship of changed coenrichment signals and PSI
#import DAS (merging H0_H1 and H0_H6 for more inc or sk) as lists
comp_AS_HS0_DAS_heat_trt_label <- c("PSI_all", "ctrl", "merged_H16_aft_H0")
AS_type_focused <- c("RI", "SE")
location_gr <- c("upstream_feature_gr", "feature_gr", "downstream_feature_gr")

list_gr_raw_comp_DAS_heat_trt <- list()
names_list_gr_raw_comp_DAS_heat_trt <- c()

for (i in seq_along(AS_type_focused)) {
  for (j in seq_along(location_gr)) {
    
    path_input <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_focused[[i]], "_", location_gr[[j]], "_merged_H16_aft_H0")
    
    name_input <- str_c(AS_type_focused[[i]], "_", location_gr[[j]], "_merged_H16_aft_H0")
    
    names_list_gr_raw_comp_DAS_heat_trt <- append(names_list_gr_raw_comp_DAS_heat_trt, name_input)
    
    gr_sample <- read_delim(path_input, delim = "\t", col_names = TRUE) %>%
      as_granges() %>%
      GenomicRanges::reduce(ignore.strand=TRUE) 
    
    list_gr_raw_comp_DAS_heat_trt <- append(list_gr_raw_comp_DAS_heat_trt, list(gr_sample))
    
    
    
  }
}

names(list_gr_raw_comp_DAS_heat_trt) <- names_list_gr_raw_comp_DAS_heat_trt

#Here we are going to separately calculate the coenrichment index for each DAS RI/SE events
#To create boxplot, each AS event should be separated for calculating co-enrichment index
#create the function to separate whole lists of DAS, AS, no-AS genes into separated events
separate_DAS_list_into_each_event <- 
  function(data){
    data %>%
      as_tibble() %>%
      mutate(AS_label = c(1:n())) %>%
      as_granges()
  }

#create the above-mentioned separated events
list_gr_raw_comp_DAS_heat_trt_separated_events <- 
  list_gr_raw_comp_DAS_heat_trt %>%
  map(. %>% separate_DAS_list_into_each_event)

names(list_gr_raw_comp_DAS_heat_trt_separated_events) <- names_list_gr_raw_comp_DAS_heat_trt


#import the files of epigenetic marks under heat treatments
heat_trt <- c("0H", "1H", "6H")
mark_heat_trt <- c("H3K9ac", "H3K18ac", "H3K27ac", "Pol2", "H3K4me3")
list_gr_mark_heat_trt <- list()
names_list_gr_mark_heat_trt <- c()

#each epigenetic mark include one replicate
for (i in seq_along(heat_trt)) {
  for (j in seq_along(mark_heat_trt)) {
    
    input_gr_mark <- read_delim(str_c("./data/peak_epimarks/heat_trt/M82_", heat_trt[[i]], "_", mark_heat_trt[[j]], "_rep1_p0.05_peaks.bed"), 
                                delim = "\t", col_names = c("seqnames", "start", "end", "signal")) %>%
      as_granges()
    
    list_gr_mark_heat_trt <- append(list_gr_mark_heat_trt, list(input_gr_mark))
    
    names_list_gr_mark_heat_trt <- append(names_list_gr_mark_heat_trt, str_c(heat_trt[[i]], "_", mark_heat_trt[[j]]))
    
  }
}







#perform co-enrichment analysis for pairs of marks and make the plots
#create the combination of two marks 
comb_pair_epimark_heat_trt_raw <- combs(seq_along(mark_heat_trt), 2) 

#since there are three heat trt for 5 epigenetic marks
#here I produce three sets of 10 combinations for 5 marks (1~5, 6~1-, 11~15)
comb_pair_epimark_heat_trt_table <- 
  tibble(mark_1 = c(comb_pair_epimark_heat_trt_raw[, 1], c(comb_pair_epimark_heat_trt_raw[, 1]+5), c(comb_pair_epimark_heat_trt_raw[, 1]+10)),
         mark_2 = c(comb_pair_epimark_heat_trt_raw[, 2], c(comb_pair_epimark_heat_trt_raw[, 2]+5), c(comb_pair_epimark_heat_trt_raw[, 2]+10)))

#create a list containing peaks of 5 epigenetic features under 3 trt (15 samples)
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

#create function for calculating coenrichment index for separated events but in a single table
#the interval position is modified according their index
#if the index is 0, the interval stays the same as their initial event interval
#if the index is above 0, the inverval will be modified based on the union of peaks of two epigenetic marks
create_coenrichment_index_two_marks_separated_events_interval_modified <- 
  function(data, a){
    #
    mark_1_feature_intersection <-
      join_overlap_intersect(list_comb_pair_epimark_heat_trt[[a[1]]][[1]], data) %>%
      as_tibble() %>%
      group_by(AS_label) %>%
      mutate(sub_AS_label = c(1:n())) %>%
      mutate(AS_label_mark_1 = str_c(AS_label, "_", sub_AS_label)) %>%
      dplyr::select(-sub_AS_label) %>%
      as_granges()
    
    #
    mark_2_feature_intersection <-
      join_overlap_intersect(list_comb_pair_epimark_heat_trt[[a[1]]][[2]], data) %>%
      as_tibble() %>%
      group_by(AS_label) %>%
      mutate(sub_AS_label = c(1:n())) %>%
      mutate(AS_label_mark_2 = str_c(AS_label, "_", sub_AS_label)) %>%
      dplyr::select(-sub_AS_label) %>%
      as_granges()
    
    
    #
    sum_bp_coenrichment <-
      join_overlap_intersect(mark_1_feature_intersection, mark_2_feature_intersection) %>%
      as_tibble() %>%
      dplyr::select(AS_label = AS_label.x, AS_label_mark_1, AS_label_mark_2, two_mark_sum_up_intersection = width) 
    
    
    
    #
    mark_1_sum_bp_intersection <-
      mark_1_feature_intersection %>%
      as_tibble() %>%
      dplyr::select(AS_label, AS_label_mark_1, mark1_sum_bp_intersection = width, start_mark_1 = start, end_mark_1 = end)
    
    
    #
    mark_2_sum_bp_intersection <-
      mark_2_feature_intersection %>%
      as_tibble() %>%
      dplyr::select(AS_label, AS_label_mark_2, mark2_sum_bp_intersection = width, start_mark_2 = start, end_mark_2 = end)
    
    
    data %>%
      as_tibble() %>%
      mutate(mark_1 = names_comb_pair_epimark_heat_trt_table[a[1],]$mark_1,
             mark_2 = names_comb_pair_epimark_heat_trt_table[a[1],]$mark_2) %>%
      left_join(sum_bp_coenrichment) %>%
      left_join(mark_1_sum_bp_intersection) %>%
      left_join(mark_2_sum_bp_intersection)%>%
      replace_na(list(start_mark_1 = 0, start_mark_2 = 0, end_mark_1 = 0, end_mark_2 = 0, AS_label_mark_1 = "no", AS_label_mark_2 = "no", two_mark_sum_up_intersection = 0, mark2_sum_bp_intersection = 0, mark1_sum_bp_intersection = 0)) %>%
      mutate(coenrichment_index = if_else(mark1_sum_bp_intersection == 0 | mark2_sum_bp_intersection == 0, 0, two_mark_sum_up_intersection/(mark1_sum_bp_intersection + mark2_sum_bp_intersection - two_mark_sum_up_intersection))) %>%
      mutate(start = if_else(start_mark_1 == 0 | start_mark_2 == 0, start, 
                             if_else(start_mark_1 <= start_mark_2, start_mark_1, start_mark_2)),
             end = if_else(end_mark_1 == 0 | end_mark_2 == 0, start, 
                           if_else(end_mark_1 >= end_mark_2, end_mark_1, end_mark_2)),
             width = end - start + 1) %>%
      dplyr::select(seqnames, start, end, width, strand, AS_label, mark_1, mark_2, AS_label_mark_1, AS_label_mark_2, coenrichment_index)
  }

list_delta_coenrichment_index_DAS <- list()
names_list_delta_coenrichment_index_DAS <- list()

for (i in seq_along(comb_pair_epimark_heat_trt_raw[,1])) {
  for (j in seq_along(list_gr_raw_comp_DAS_heat_trt_separated_events)) {
    coenrichment_index_table_0H <- 
      create_coenrichment_index_two_marks_separated_events_interval_modified(list_gr_raw_comp_DAS_heat_trt_separated_events[[j]], i) %>%
      mutate(comb_0H = str_c(mark_1, str_remove(mark_2, ".H")), coenrichment_index_0H = coenrichment_index) %>%
      dplyr::select(-mark_1, -mark_2, -coenrichment_index) %>%
      as_granges()
    
    coenrichment_index_table_1H <-
      create_coenrichment_index_two_marks_separated_events_interval_modified(list_gr_raw_comp_DAS_heat_trt_separated_events[[j]], i+10) %>%
      mutate(comb_1H = str_c(mark_1, str_remove(mark_2, ".H")), coenrichment_index_1H = coenrichment_index) %>%
      dplyr::select(-mark_1, -mark_2, -coenrichment_index) %>%
      as_granges()
    
    delta_coenrichment_index_DAS_table <- 
    join_overlap_intersect(coenrichment_index_table_0H, coenrichment_index_table_1H) %>%
      as_tibble() %>%
      filter(coenrichment_index_0H != 0 | coenrichment_index_1H != 0) %>%
      dplyr::select(seqnames, start, end, width, AS_label = AS_label.x, comb = comb_0H, coenrichment_index_0H, coenrichment_index_1H) %>%
      mutate(comb = str_remove(comb, ".H_"), delta_coenrichment_index = coenrichment_index_1H - coenrichment_index_0H) 
    
    list_delta_coenrichment_index_DAS <- append(list_delta_coenrichment_index_DAS, list(delta_coenrichment_index_DAS_table))
    name_delta_coenrichment_index_DAS_table <- str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/delta_coenrichment_index_0H_to_1H_DAS/", names(list_gr_raw_comp_DAS_heat_trt_separated_events)[[j]], str_remove(names_comb_pair_epimark_heat_trt_table$mark_1[[i]], ".H"),
                                                     str_remove(names_comb_pair_epimark_heat_trt_table$mark_2[[i]], ".H"), "_delta_coenrichment_index.txt")
    
    names_list_delta_coenrichment_index_DAS <- append(names_list_delta_coenrichment_index_DAS, list(name_delta_coenrichment_index_DAS_table))
  }
}
names(list_gr_raw_comp_DAS_heat_trt_separated_events)

pwalk(list(list_delta_coenrichment_index_DAS, names_list_delta_coenrichment_index_DAS), write_delim, delim = "\t", col_names = TRUE)

#import DAS events (splicing sites)
All_DAS_events_all_comparisons_refined_for_granges <- 
  read_delim("./data/temp/All_DAS_events_all_comparisons_refined_for_granges", delim = "\t", col_names = TRUE)

RI_SE_DAS_aft_H0 <- 
All_DAS_events_all_comparisons_refined_for_granges %>%
  filter(AS_type == "SE" | AS_type == "RI") %>%
  filter(comp != "HS1_HS6") %>%
  split(.$AS_type) %>%
  map(. %>% as_granges)

#import DAS feature (intact exon/intron)
path_input_large_size_DAS <- str_c("./data/AS_DAS_refined_grouped_by_DAS_size/", AS_type_focused, "_merged_H16_aft_H0_feature_size_group_large.bed")

RI_SE_large_DAS_aft_H0_raw <-
map(path_input_large_size_DAS, read_delim, delim = "\t", col_names = TRUE) %>%
  map(. %>% mutate(feature_event_label = c(1:n()))) %>%
  map(. %>% as_granges())

large_DAS_aft_H0_Incdf <- vector("list", length = length(AS_type_focused))
list_delta_coenrichment_index_DAS_feature <- vector("list", length = length(AS_type_focused))


for (i in seq_along(large_DAS_aft_H0_Incdf)) {
  large_DAS_aft_H0_Incdf[[i]] <-
    join_overlap_intersect(RI_SE_large_DAS_aft_H0_raw[[i]], RI_SE_DAS_aft_H0[[i]]) %>%
    as_tibble() %>%
    dplyr::select(seqnames, start, end, width, strand, Incdf) %>%
    mutate(more_retained = if_else(Incdf > 0, "yes", "no")) %>%
    as_granges()
  
  list_delta_coenrichment_index_DAS_feature[[i]] <- 
    list_delta_coenrichment_index_DAS[str_detect(names_list_delta_coenrichment_index_DAS, str_c(AS_type_focused[[i]], "_feature_gr"))] %>%
    map(. %>% as_granges())
}


list_delta_coenrichment_index_DAS_feature[[1]]

test <- 
join_overlap_intersect(list_delta_coenrichment_index_DAS_feature[[1]][[1]], large_DAS_aft_H0_Incdf[[1]]) %>%
  as_tibble() %>%
  dplyr::select(seqnames, start, end, width, strand, Incdf, delta_coenrichment_index_comb_1 = delta_coenrichment_index) %>%
  filter(Incdf > 0) %>%
  as_granges() %>%
  join_overlap_intersect(., list_delta_coenrichment_index_DAS_feature[[1]][[2]]) %>%
  as_tibble() %>%
  dplyr::select(seqnames, start, end, width, strand, Incdf, delta_coenrichment_index_comb_1, delta_coenrichment_index_comb_2 = delta_coenrichment_index) %>%
  View()

#check whether signals of 5 marks are more in DAS large events
names_list_gr_mark_heat_trt[[1]]

create_delta_epigenetic_feature_signal_plots <-
  function(a){
    DAS_event_gr <-
      large_DAS_aft_H0_Incdf[[a[2]]] %>%
      as_tibble() %>%
      mutate(DAS_label = c(1:n())) %>%
      as_granges()
    
    DAS_event_gr_epimark_0H <-
      DAS_event_gr %>%
      join_overlap_intersect(., list_gr_mark_heat_trt[[a[1]]]) %>%
      as_tibble() %>%
      mutate(peak_height_0H = signal) %>%
      dplyr::select(-signal)
    
    DAS_event_gr_epimark_1H <-
      DAS_event_gr %>%
      join_overlap_intersect(., list_gr_mark_heat_trt[[a[1] + 5]]) %>%
      as_tibble() %>%
      mutate(peak_height_1H = signal) %>%
      dplyr::select(-signal)
    
    DAS_event_gr_epimark_6H <-
      DAS_event_gr %>%
      join_overlap_intersect(., list_gr_mark_heat_trt[[a[1] + 10]]) %>%
      as_tibble() %>%
      mutate(peak_height_6H = signal) %>%
      dplyr::select(-signal)
    
    DAS_event_gr_cor_raw <-
      DAS_event_gr %>%
      as_tibble() %>%
      left_join(DAS_event_gr_epimark_0H) %>%
      left_join(DAS_event_gr_epimark_1H) %>%
      left_join(DAS_event_gr_epimark_6H) %>%
      replace_na(list(peak_height_0H = 0, peak_height_1H = 0, peak_height_6H = 0)) %>%
      mutate(delta_peak_height_0H_1H = peak_height_1H / mean(peak_height_1H) - peak_height_0H / mean(peak_height_0H),
             delta_peak_height_0H_6H = peak_height_6H / mean(peak_height_6H) - peak_height_0H / mean(peak_height_0H)) %>%
      mutate(
        AS_type = str_c("delta_", AS_type_focused[[a[2]]], "_PSI"),
        epigenetic_feature_0H_1H = str_c(
          "delta",
          str_remove(names_list_gr_mark_heat_trt[[a[1]]], ".H"),
          "_peak_height_0H_1H"
        ),
        epigenetic_feature_0H_6H = str_c(
          "delta",
          str_remove(names_list_gr_mark_heat_trt[[a[1]]], ".H"),
          "_peak_height_0H_6H"
        )
      )
    
    DAS_event_gr_cor_raw 
    
}

create_delta_epigenetic_feature_signal_plots(c(1,1))

list_table_delta_epigenetic_feature_PSI <- list()
list_plots_delta_epigenetic_features_0H_1H <- list()
list_plots_delta_epigenetic_features_0H_6H <- list()

for (i in seq_along(mark_heat_trt)) {
  for (j in seq_along(AS_type_focused)) {
    delta_epi_table <- create_delta_epigenetic_feature_signal_plots(c(i,j))
    
    delta_epi_plot_0H_1H <-
      delta_epi_table %>%
      ggplot(aes(delta_peak_height_0H_1H, Incdf)) +
      geom_point() +
      facet_grid(AS_type ~ epigenetic_feature_0H_1H) +
      theme(
        strip.text.x = element_text(
          colour = "black",
          face = "bold",
          size = 20
        ),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(
          hjust = 0.5,
          colour = "black",
          face = "bold",
          size = 18
        ),
        legend.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 20, face = "bold"),
        strip.text.y = element_text(
          colour = "black",
          face = "bold",
          size = 18
        ),
        axis.text.x = element_text(
          colour = "black",
          size = 20,
          face = "bold",
          angle = 60,
          vjust = 0.5
        )
      )
    
    delta_epi_plot_0H_6H <-
      delta_epi_table %>%
      ggplot(aes(delta_peak_height_0H_6H, Incdf)) +
      geom_point() +
      facet_grid(AS_type ~ epigenetic_feature_0H_6H) +
      theme(
        strip.text.x = element_text(
          colour = "black",
          face = "bold",
          size = 20
        ),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(
          hjust = 0.5,
          colour = "black",
          face = "bold",
          size = 18
        ),
        legend.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 20, face = "bold"),
        strip.text.y = element_text(
          colour = "black",
          face = "bold",
          size = 18
        ),
        axis.text.x = element_text(
          colour = "black",
          size = 20,
          face = "bold",
          angle = 60,
          vjust = 0.5
        )
      )
    
    list_table_delta_epigenetic_feature_PSI <- append(list_table_delta_epigenetic_feature_PSI, list(delta_epi_table))
    list_plots_delta_epigenetic_features_0H_1H <- append(list_plots_delta_epigenetic_features_0H_1H, list(delta_epi_plot_0H_1H))
    list_plots_delta_epigenetic_features_0H_6H <- append(list_plots_delta_epigenetic_features_0H_6H, list(delta_epi_plot_0H_6H))
    }
}
list_plots_delta_epigenetic_features_0H_6H[[1]]

#add the information of gene name
#import gene name
M82_rMATs_anno_all_gene <-
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "gene_name"))

M82_rMATs_anno_all_gene_gr <- 
  M82_rMATs_anno_all_gene %>%
  dplyr::select(seqnames = chr, start = str, end, feature, source, gene_name) %>%
  as_granges() 

produce_table_delta_PSI_epimark_peak_height <- 
  function(data){
  join_overlap_intersect(as_granges(data), M82_rMATs_anno_all_gene_gr) %>%
  filter(str_detect(gene_name, "MC")==FALSE) %>%
  as_tibble() %>%
  mutate(AS_type = str_remove(AS_type, "_PSI")) %>%
  mutate(AS_type = str_remove(AS_type, "delta_")) %>%
  mutate(epigenetic_feature_0H_1H = str_remove(epigenetic_feature_0H_1H, "_peak.*")) %>%
  mutate(epigenetic_feature_0H_1H = str_remove(epigenetic_feature_0H_1H, "del.*_")) %>%
  dplyr::select(seqnames, start, end, width, strand, gene_name, AS_type, epigenetic_feature = epigenetic_feature_0H_1H, Incdf, delta_peak_height_0H_1H, delta_peak_height_0H_6H,peak_height_0H, peak_height_1H, peak_height_6H) 
}  


join_overlap_intersect(as_granges(list_table_delta_epigenetic_feature_PSI[[1]]), M82_rMATs_anno_all_gene_gr) %>%
  filter(str_detect(gene_name, "MC")==FALSE) %>%
  as_tibble() %>%
  mutate(AS_type = str_remove(AS_type, "_PSI")) %>%
  mutate(AS_type = str_remove(AS_type, "delta_")) %>%
  mutate(epigenetic_feature_0H_1H = str_remove(epigenetic_feature_0H_1H, "_peak.*")) %>%
  mutate(epigenetic_feature_0H_1H = str_remove(epigenetic_feature_0H_1H, "del.*_")) %>%
  dplyr::select(seqnames, start, end, width, strand, gene_name, AS_type, epigenetic_feature = epigenetic_feature_0H_1H, Incdf, delta_peak_height_0H_1H, delta_peak_height_0H_6H,peak_height_0H, peak_height_1H, peak_height_6H) 

combined_table_delta_PSI_epimark_peak_height <- 
  map(list_table_delta_epigenetic_feature_PSI, produce_table_delta_PSI_epimark_peak_height) %>%
  bind_rows()

write_csv(combined_table_delta_PSI_epimark_peak_height, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/DAS_delta_PSI_epigenetic_feature_result/combined_table_delta_PSI_epimark_peak_height.csv", col_names = TRUE)

write_the_output_delta_PSI_epimark_table <- 
function(data){
  write_delim(data, str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/DAS_delta_PSI_epigenetic_feature_result/DAS_", str_remove(data$AS_type[[1]], "delta_"),
        "_", str_c(data$epigenetic_feature_0H_1H[[1]], "_6H"), "_table.txt"), delim = "\t", col_names = TRUE)
}

walk(list_table_delta_epigenetic_feature_PSI, write_the_output_delta_PSI_epimark_table)

write_the_output_delta_PSI_epimark_plot <- 
  function(data, data2, data3){
    ggsave(str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/DAS_delta_PSI_epigenetic_feature_result/DAS_", str_remove(data$AS_type[[1]], "delta_"),
                            "_", data$epigenetic_feature_0H_1H[[1]], ".jpeg"), data2, width = 400, height = 240, units = c("mm"), dpi = 320)
    
    ggsave(str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/DAS_delta_PSI_epigenetic_feature_result/DAS_", str_remove(data$AS_type[[1]], "delta_"),
                 "_", data$epigenetic_feature_0H_6H[[1]], ".jpeg"), data3, width = 400, height = 240, units = c("mm"), dpi = 320)
  }

pwalk(list(list_table_delta_epigenetic_feature_PSI, list_plots_delta_epigenetic_features_0H_1H, list_plots_delta_epigenetic_features_0H_6H), write_the_output_delta_PSI_epimark_plot)

list_table_delta_epigenetic_feature_PSI[[1]]$AS_type

build_model <-
  function(data){
    cor.test(data$Incdf, data$delta_signal) 
  }

list_table_delta_epigenetic_feature_PSI %>%
  map(. %>% build_model)

build_model(list_table_delta_epigenetic_feature_PSI[[1]]) 


