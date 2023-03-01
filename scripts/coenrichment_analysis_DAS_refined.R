Sys.getlocale()
Sys.setenv(LANG = "en_US.UTF-8")
Packages <- c("plyranges", "tidyverse", "doParallel", "foreach", "nullranges", "caTools", "ggpubr", "fs", "HelloRanges")
lapply(Packages, library, character.only = TRUE)

#Previously, our analysis focused on three regions (100-bp upstream/downstream junctions of exon/intron containing AS),
#Since We are more interested in the whole cassettes (intron-AS-included-exon-intron or exon-AS-included-intron-exon)
#We will modify how to define the group of upstream_gr and downstream_gr here
#Plus, we also include the DAS gene analysis, so we will extract the information of DAS gene first.
#Then the DAS info will be used to filter out PSI_all and control groups to keep these two groups containing only PSI-stable or constitutively splicing genes in response to heat

#import the file of refined DAS genes
All_DAS_events_all_comparisons_refined <-
  read_delim("./data/rMATS_out/All_DAS_events_all_comparisons_refined.tsv", delim = "\t", col_names = TRUE) 

#split Event_ID for extracting positions of AS 
All_DAS_events_all_comparisons_refined_pos_index <- 
  str_split(All_DAS_events_all_comparisons_refined$Event_ID, "_")

length(All_DAS_events_all_comparisons_refined_pos_index[[1]])

#create the function that can extract 6 positions of each DAS
extract_pos <- function(data){
  length_pos_index <- length(data)
  
  pos_index <- as.double(data[(length_pos_index-5):length_pos_index])
  
  tibble(pos_1 = pos_index[[1]],
         pos_2 = pos_index[[2]],
         pos_3 = pos_index[[3]],
         pos_4 = pos_index[[4]],
         pos_5 = pos_index[[5]],
         pos_6 = pos_index[[6]])
}

#test the function above
extract_pos(All_DAS_events_all_comparisons_refined_pos_index[[15]])

#create the table with positions of DAS
DAS_gene_pos <- 
  map(All_DAS_events_all_comparisons_refined_pos_index, extract_pos) %>%
  bind_rows()

#Since there is no strand information of DAS, I used gene annotation data to extract the strand information of each DAS gene
M82_rMATs_anno_all_gene <-
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "gene_name"))

#store DAS gene ID as a vector
DAS_geneID <- unique(All_DAS_events_all_comparisons_refined$GeneID)

#Create the function that can use the content in the DAS-gene-ID vector for producing a table
search_for_strand_DAS_gene <- function(data){
  M82_rMATs_anno_all_gene %>%
    filter(str_detect(gene_name, data)) %>%
    mutate(gene_name_f = data) %>%
    dplyr::select(GeneID = gene_name_f, strand) %>%
    distinct()
}

#create the table of DAS genes with their strand information
DAS_geneID_strand <- 
  map(DAS_geneID, search_for_strand_DAS_gene) %>%
  bind_rows()



#extract chromosome info for all DAS gene
DAS_geneID_chr <- 
  All_DAS_events_all_comparisons_refined %>%
  dplyr::select(GeneID) %>%
  mutate(chr_modified = str_remove(GeneID, "g.*")) %>%
  mutate(chr_modified = str_remove(chr_modified, "_NP.*")) %>%
  mutate(chr_modified = str_remove(chr_modified, "Solyc_")) %>%
  mutate(chr_modified = str_replace(chr_modified, "Solyc0", "chr")) %>%
  mutate(chr_modified = str_replace(chr_modified, "Solyc", "chr")) %>%
  filter(str_detect(GeneID, "Tomato")!=TRUE) %>%
  distinct()


All_DAS_events_all_comparisons_refined_raw <-   
All_DAS_events_all_comparisons_refined %>%
  left_join(DAS_geneID_strand) %>%
  left_join(DAS_geneID_chr) %>%
  bind_cols(DAS_gene_pos) %>%
  filter(str_detect(GeneID, "Tomato")!=TRUE) %>%
  mutate(str = if_else(AS_type == "A5S" & strand == "+", pos_4, 
                       if_else(AS_type == "A5S" & strand == "-", pos_1,
                               if_else(AS_type == "A3S" & strand == "+", pos_1,
                                       if_else(AS_type == "A3S" & strand == "-", pos_4, 
                                               if_else(AS_type == "SE", pos_1, pos_4)))))) %>%
  mutate(end = if_else(AS_type == "A5S" & strand == "+", pos_2, 
                       if_else(AS_type == "A5S" & strand == "-", pos_3,
                               if_else(AS_type == "A3S" & strand == "+", pos_3,
                                       if_else(AS_type == "A3S" & strand == "-", pos_2, 
                                               if_else(AS_type == "SE", pos_2, pos_5)))))) %>%
  dplyr::select(chr = chr_modified, str, end, strand, AS_type, comp, IncLevel1, IncLevel2, Incdf)

write_delim(All_DAS_events_all_comparisons_refined_raw, "./data/temp/All_DAS_events_all_comparisons_refined_raw", delim = "\t", col_names = TRUE)

IncLevel_two_trts_rep_table <- 
  tibble(
  IncLevel1_r1 = as.double(str_split(All_DAS_events_all_comparisons_refined_raw$IncLevel1, ",", simplify = TRUE)[,1]),
  IncLevel1_r2 = as.double(str_split(All_DAS_events_all_comparisons_refined_raw$IncLevel1, ",", simplify = TRUE)[,2]),
  IncLevel2_r1 = as.double(str_split(All_DAS_events_all_comparisons_refined_raw$IncLevel2, ",", simplify = TRUE)[,1]),
  IncLevel2_r2 = as.double(str_split(All_DAS_events_all_comparisons_refined_raw$IncLevel2, ",", simplify = TRUE)[,2])
)

All_DAS_events_all_comparisons_refined_for_granges <- 
All_DAS_events_all_comparisons_refined_raw %>%
  bind_cols(IncLevel_two_trts_rep_table) %>%
  mutate(IncLevel_former_tp = if_else(comp == "HS1.vs.HS6", (IncLevel1_r1 + IncLevel1_r2)/2, (IncLevel2_r1 + IncLevel2_r2)/2),
         IncLevel_latter_tp = if_else(comp == "HS1.vs.HS6", (IncLevel2_r1 + IncLevel2_r2)/2, (IncLevel1_r1 + IncLevel1_r2)/2)) %>%
  mutate(Incdf_f = IncLevel_latter_tp - IncLevel_former_tp) %>%
  mutate(more_AS = if_else(Incdf_f > 0, "yes", "no")) %>%
  mutate(comp = str_replace(comp, ".vs.", "_")) %>%
  mutate(comp = if_else(comp == "HS1_HS0", "HS0_HS1", if_else(comp == "HS6_HS0", "HS0_HS6", comp))) %>%
  select(seqnames = chr, start = str, end, strand, AS_type, comp, Incdf = Incdf_f, IncLevel_former_tp, IncLevel_latter_tp, more_AS) %>%
  arrange(seqnames, start)


write_delim(All_DAS_events_all_comparisons_refined_for_granges, "./data/temp/All_DAS_events_all_comparisons_refined_for_granges", delim = "\t", col_names = TRUE)
  
All_DAS_events_all_comparisons_refined_for_granges %>%
  as_granges()


#import bootstraped data under HS0 condition, and remove DAS genes from these data
#and stored the refined HS0 data in another subfolder
AS_type_bootstrap <- c("RI", "SE")
PSI_type_bootstrap <- c(str_c("PSI_", c("5", "5_20", "20_40", "40_60","60_80", "80_95", "95")), "PSI_all", "ctrl")
location_gr_AS_included_body <- c("feature_gr")

registerDoParallel(cores = 5)
getDoParWorkers()

 
foreach(i=seq_along(AS_type_bootstrap), .packages = c("plyranges", "tidyverse")) %:%
  foreach(j=seq_along(PSI_type_bootstrap), .packages = c("plyranges", "tidyverse")) %:%
  foreach(l=seq_along(location_gr), .packages = c("plyranges", "tidyverse")) %:%
  foreach(k=1:1000, .packages = c("plyranges", "tidyverse")) %dopar% {
    
    
    bt_sample_path <- str_c("./data/AS_HS0_bootstrap_samples/", AS_type_bootstrap[[i]], "_", PSI_type_bootstrap[[j]], "_", location_gr[[l]], "_bootstrap_", k)
    
    bt_sample_raw <- read_delim(bt_sample_path, delim = "\t", col_names = TRUE)
      
    bt_sample_output <- 
    bt_sample_raw %>%
      as_granges() %>%
      mutate(n_overlaps = count_overlaps(., All_DAS_events_all_comparisons_refined_granges)) %>%
      filter(n_overlaps == 0) %>%
      as_tibble()
    
    
    bt_sample_output_path <- str_c("./data/AS_HS0_bootstrap_samples_DAS_removed/", AS_type_bootstrap[[i]], "_", PSI_type_bootstrap[[j]], "_", location_gr[[l]], "_bootstrap_DAS_removed_", k)
    
    write_delim(bt_sample_output, bt_sample_output_path, delim = "\t", col_names = TRUE)
    
  }



#install HelloRanges
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

BiocManager::install("HelloRanges")
library(HelloRanges)

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

#Based on 1000 bootstrapped samples previously produced in coenrichment analysis for only HS0 events
#we focused on feature_gr only to combined them in one sample for following analysis
#test for one sample combined by 100 bootstrap samples for producing histograms (enrichment analysis for single mark or co-enrichment analysis for multiple marks)
AS_type_bootstrap <- c("RI", "SE")
PSI_type_simplified <- c("PSI_all", "ctrl")
location_gr_AS_included_body <- c("feature_gr")

registerDoParallel(cores = 5)
getDoParWorkers()

PSI_DAS_removed_1_1000_list <- 
foreach(i=seq_along(AS_type_bootstrap), .packages = c("plyranges", "tidyverse")) %:%
  foreach(j=seq_along(PSI_type_simplified), .packages = c("plyranges", "tidyverse")) %:%
  foreach(l=seq_along(location_gr_AS_included_body), .packages = c("plyranges", "tidyverse")) %:%
  foreach(k=1:1000, .packages = c("plyranges", "tidyverse")) %dopar% {
    
    bt_sample_output_path <- str_c("./data/AS_HS0_bootstrap_samples_DAS_removed/", AS_type_bootstrap[[i]], "_", PSI_type_simplified[[j]], "_", location_gr_AS_included_body[[l]], "_bootstrap_DAS_removed_", k)
    
    read_delim(bt_sample_output_path, delim = "\t", col_names = TRUE)
    
  }

#combined 1-1000 bootstrapped samples into one sample
PSI_DAS_removed_reduced <- 
PSI_DAS_removed_1_1000_list %>%
  map(. %>% map(. %>% map(. %>% Reduce(bind_rows, .)))) %>%
  map(. %>% map(. %>% map(. %>% distinct())))

#generate the output for combined bootstrapped samples with 1-1000 bootstraping samples
for (i in seq_along(AS_type_bootstrap)) {
  for (j in seq_along(PSI_type_simplified)) {
    for (l in seq_along(location_gr_AS_included_body)) {
      
      path_output <- str_c("./data/AS_HS0_bootstrap_samples_DAS_removed/", AS_type_bootstrap[[i]], "_", PSI_type_simplified[[j]], "_", location_gr_AS_included_body[[l]], "_bootstrap_DAS_removed_1_1000_combined")
      
      write_delim(PSI_DAS_removed_reduced[[i]][[j]][[l]], path_output, delim = "\t", col_names = TRUE)
      
    }
  }
}

#import all exons and introns of annotation data for extracting upstream/downstream AS-included intron/exon
#and extract DAS exons and introns
#import all exons
M82_rMATs_anno_all_exon <-
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_exon.bed", col_names = c("seqnames", "start", "end", "strand", "feature", "source", "gene_name")) %>%
  filter(strand == "+" | strand == "-")

M82_rMATs_anno_all_exon %>%
  summarise(median_exon = median(end-start))


#import all introns
M82_rMATs_anno_all_intron <-
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_intron.bed", col_names = c("seqnames", "start", "end", "strand", "feature", "source", "gene_name"))

M82_rMATs_anno_all_intron %>%
  summarise(median_intron = median(end-start))


#combine the annotation data of intron and exon
M82_rMATs_anno_all_intron_exon_ordered <- 
  bind_rows(M82_rMATs_anno_all_exon, M82_rMATs_anno_all_intron) %>%
  arrange(gene_name, seqnames, start) %>%
  group_by(gene_name) %>%
  mutate(genetic_feature_label = c(1:n()))

write_delim(M82_rMATs_anno_all_intron_exon_ordered, "./data/M82_annotation_data/M82_rMATs_anno_all_intron_exon_ordered.bed", delim = "\t", col_names = TRUE)

PSI_DAS_removed_reduced[[1]][[1]][[1]] %>%
  left_join(M82_rMATs_anno_all_intron_exon_ordered)
  
  

#import DAS for extracting introns or exons with AS with following analysis
All_DAS_events_all_comparisons_refined_for_granges <- 
  read_delim("./data/temp/All_DAS_events_all_comparisons_refined_for_granges", delim = "\t", col_names = TRUE) %>%
  filter(AS_type == "RI" | AS_type == "SE") %>%
  split(.$AS_type) %>%
  map(. %>% split(.$comp)) %>%
  map(. %>% map(. %>% split(.$more_AS)))


#create the function that extracts introns or exons with AS (RI or SE)
produce_feature_with_AS <-
function(data, data2){
  data %>%
    as_granges() %>%
    mutate(n_overlaps = count_overlaps(., as_granges(data2))) %>%
    filter(n_overlaps != 0) %>%
    as_tibble() %>%
    dplyr::select(seqnames, start, end, strand) %>%
    distinct()
}

#based on the DAS data, I produced exons/introns with SE/RI of DAS genes in terms of different comparisons of heat treatments and their trends (Inclusion or skipping)
heat_trt_comp <- c("HS0_HS1", "HS0_HS6", "HS1_HS6")
Inc_AS_aft_trt <- c("inc", "sk")

for (i in seq_along(AS_type_bootstrap)) {
  for (j in seq_along(heat_trt_comp)) {
    for (k in seq_along(Inc_AS_aft_trt)) {
    
      if(i == 1){
        feature_table <- produce_feature_with_AS(M82_rMATs_anno_all_intron, All_DAS_events_all_comparisons_refined_for_granges[[i]][[j]][[k]])
        
        feature_upstream_table <-
          feature_table %>%
          mutate(start = as.double(start), end = as.double(end)) %>%
          mutate(start_new = if_else(strand == "+", start - 100, end),
                 end_new = if_else(strand == "+", start, end + 100)) %>%
          dplyr::select(seqnames, start = start_new, end = end_new, strand)
        
        feature_downstream_table <-
          feature_table %>%
          mutate(start = as.double(start), end = as.double(end)) %>%
          mutate(start_new = if_else(strand == "+", end, start - 100),
                 end_new = if_else(strand == "+", end + 100, start)) %>%
          dplyr::select(seqnames, start = start_new, end = end_new, strand)
        
        path_output <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_bootstrap[[i]], "_", location_gr, "_", heat_trt_comp[[j]], "_more_", Inc_AS_aft_trt[[k]], "_AS_aft_trt")
        
        list_output <- list(feature_upstream_table, feature_table, feature_downstream_table)
        pwalk(list(list_output, path_output), write_delim, delim = "\t", col_names = TRUE)  
      } else {
        feature_table <- produce_feature_with_AS(M82_rMATs_anno_all_exon, All_DAS_events_all_comparisons_refined_for_granges[[i]][[j]][[k]])
        
        feature_upstream_table <-
          feature_table %>%
          mutate(start = as.double(start), end = as.double(end)) %>%
          mutate(start_new = if_else(strand == "+", start - 100, end),
                 end_new = if_else(strand == "+", start, end + 100)) %>%
          dplyr::select(seqnames, start = start_new, end = end_new, strand)
        
        feature_downstream_table <-
          feature_table %>%
          mutate(start = as.double(start), end = as.double(end)) %>%
          mutate(start_new = if_else(strand == "+", end, start - 100),
                 end_new = if_else(strand == "+", end + 100, start)) %>%
          dplyr::select(seqnames, start = start_new, end = end_new, strand)
        
        path_output <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_bootstrap[[i]], "_", location_gr, "_", heat_trt_comp[[j]], "_more_", Inc_AS_aft_trt[[k]], "_AS_aft_trt")
        
        list_output <- list(feature_upstream_table, feature_table, feature_downstream_table)
        pwalk(list(list_output, path_output), write_delim, delim = "\t", col_names = TRUE)
      } 
      
      
      
    }
  }
}

#import PSI combined, ctrl (no AS), DAS (sk/Inc in terms of 3 heat comp, 6 in total) as lists
comp_AS_HS0_DAS_heat_trt_label <- c("PSI_all", "ctrl", str_c("HS0_HS1_more_", Inc_AS_aft_trt),
                                    str_c("HS0_HS6_more_", Inc_AS_aft_trt), str_c("HS1_HS6_more_", Inc_AS_aft_trt))
AS_type_focused <- c("RI", "SE")
location_gr <- c("upstream_feature_gr", "feature_gr", "downstream_feature_gr")

list_gr_raw_comp_AS_HS0_DAS_heat_trt <- list()
names_list_gr_raw_comp_AS_HS0_DAS_heat_trt <- c()

for (i in seq_along(AS_type_focused)) {
  for (j in seq_along(location_gr)) {
    for (k in seq_along(comp_AS_HS0_DAS_heat_trt_label)) {
      
      if(k %in% c(1, 2)){
        
        path_input <- str_c("./data/AS_HS0_bootstrap_samples_DAS_removed/", AS_type_focused[[i]], "_", comp_AS_HS0_DAS_heat_trt_label[[k]], "_", location_gr[[j]], "_bootstrap_DAS_removed_1_100_combined")
        
        name_input <- str_c(AS_type_focused[[i]], "_", location_gr[[j]], "_", comp_AS_HS0_DAS_heat_trt_label[[k]])
        
        names_list_gr_raw_comp_AS_HS0_DAS_heat_trt <- append(names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, name_input)
        
        gr_sample <- read_delim(path_input, delim = "\t", col_names = TRUE) %>%
          as_granges()
        
        list_gr_raw_comp_AS_HS0_DAS_heat_trt <- append(list_gr_raw_comp_AS_HS0_DAS_heat_trt, list(gr_sample))
      } else {
        
        path_input <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_focused[[i]], "_", location_gr[[j]], "_", comp_AS_HS0_DAS_heat_trt_label[[k]], "_AS_aft_trt")
        
        name_input <- str_c(AS_type_focused[[i]], "_", location_gr[[j]], "_", comp_AS_HS0_DAS_heat_trt_label[[k]])
        
        names_list_gr_raw_comp_AS_HS0_DAS_heat_trt <- append(names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, name_input)
        
        gr_sample <- read_delim(path_input, delim = "\t", col_names = TRUE) %>%
          as_granges()
        
        list_gr_raw_comp_AS_HS0_DAS_heat_trt <- append(list_gr_raw_comp_AS_HS0_DAS_heat_trt, list(gr_sample))
      }
      
      
    }
  }
}

names(list_gr_raw_comp_AS_HS0_DAS_heat_trt) <- names_list_gr_raw_comp_AS_HS0_DAS_heat_trt

join_overlap_intersect(list_gr_raw_comp_AS_HS0_DAS_heat_trt[[11]], list_gr_raw_comp_AS_HS0_DAS_heat_trt[[13]])

list_gr_raw_comp_AS_HS0_DAS_heat_trt[[11]] %>%
  GenomicRanges::reduce(., ignore.strand = TRUE)

?reduce
#check the size of introns/exons for different groups
list_gr_raw_comp_AS_HS0_DAS_heat_trt %>%
  map(. %>% as_tibble()) %>%
  map(. %>% summarise(median_width = median(width),
            mean_width = mean(width))) %>%
  bind_rows() %>%
  mutate(name = names_list_gr_raw_comp_AS_HS0_DAS_heat_trt) %>%
  filter(str_detect(name, "upstream")!=TRUE) %>%
  filter(str_detect(name, "downstream")!=TRUE) %>%
  View()

#create the function for calculating single-mark enrichment index for each AS cases (PSI, ctrl and other DAS)
create_enrichment_index_single_mark <- 
function(data, data2){
  intersection_feature_mark_bp <- join_overlap_intersect(data, data2) %>%
    as_tibble() %>%
    summarise(enrichment_bp = sum(end-start))
  
  mark_sum_bp <- 
    data %>%
    as_tibble() %>%
    summarise(mark_sum_bp = sum(end-start))
  
  intersection_feature_mark_bp$enrichment_bp/mark_sum_bp$mark_sum_bp
}

#create the vector containing all indexes
single_feature_enrichment_index <- c()

for (i in seq_along(list_gr_mark_heat_trt)) {
  index_raw <- map2(replicate(length(list_gr_raw_comp_AS_HS0_DAS_heat_trt), list_gr_mark_heat_trt[[i]]), list_gr_raw_comp_AS_HS0_DAS_heat_trt, create_enrichment_index_single_mark) %>%
    unlist()
  
  single_feature_enrichment_index <- append(single_feature_enrichment_index, index_raw)
}

#make the final table with all indexes
single_feature_enrichment_index_table <- 
tibble(mark = rep(names(list_gr_mark_heat_trt), each = length(list_gr_raw_comp_AS_HS0_DAS_heat_trt)),
       feature = rep(names(list_gr_raw_comp_AS_HS0_DAS_heat_trt), length(list_gr_mark_heat_trt)), 
       enrichment_index = single_feature_enrichment_index)

write_delim(single_feature_enrichment_index_table, "./data/temp/single_feature_enrichment_index_table",
            delim = "\t", col_names = TRUE)


single_feature_enrichment_index_table <- read_delim("./data/temp/single_feature_enrichment_index_table",
                                                    delim = "\t", col_names = TRUE)
single_feature_enrichment_index_plots <- 
single_feature_enrichment_index_table %>%
  mutate(AS_DAS_groups = rep(comp_AS_HS0_DAS_heat_trt_label, nrow(single_feature_enrichment_index_table)/length(comp_AS_HS0_DAS_heat_trt_label))) %>%
  mutate(enrichment_index = if_else(AS_DAS_groups=="ctrl", enrichment_index/10, enrichment_index)) %>%
  mutate(AS_region = str_remove(feature, AS_DAS_groups)) %>%
  mutate(AS_region = str_remove(AS_region, "_gr_")) %>%
  mutate(AS_DAS_groups = if_else(AS_DAS_groups=="ctrl", "ctrl/10", AS_DAS_groups)) %>%
  mutate(AS_DAS_groups = factor(AS_DAS_groups, levels = c("PSI_all", "ctrl/10", str_c("HS0_HS1_more_", Inc_AS_aft_trt),
                                                          str_c("HS0_HS6_more_", Inc_AS_aft_trt), str_c("HS1_HS6_more_", Inc_AS_aft_trt)))) %>%
  dplyr::select(mark, AS_region, enrichment_index, AS_DAS_groups)


produce_enrichment_index_hist_plots <- function(a, b){
  single_feature_enrichment_index_plots %>% 
    filter(AS_region == a[1] & mark == b[1]) %>%
    ggplot(aes(AS_DAS_groups, enrichment_index, fill=AS_DAS_groups)) +
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
    facet_grid(AS_region ~ mark) +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 20), legend.text = element_text(size = 9, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 20, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
          axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle = 60, vjust = 0.5))
}

produce_enrichment_index_hist_plots(single_feature_enrichment_index_plots$AS_region[[1]], single_feature_enrichment_index_plots$mark[[1]])

enrichment_index_table_for_plots_key <- 
  single_feature_enrichment_index_plots %>%
  ungroup() %>%
  dplyr::select(AS_region, mark) %>%
  distinct()

enrichment_index_hist_plots_list <- 
  map2(enrichment_index_table_for_plots_key$AS_region, 
       enrichment_index_table_for_plots_key$mark,
       produce_enrichment_index_hist_plots)

enrichment_index_hist_plots_list[[1]]

path_enrichment_index_plots_hist_AS_DAS <- 
  str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/enrichment_analysis_heat_single_epimark_AS_DAS_histograms/", enrichment_index_table_for_plots_key$AS_region, "_", enrichment_index_table_for_plots_key$mark, "_enrichment_index_AS_DAS.jpeg")

pwalk(list(path_enrichment_index_plots_hist_AS_DAS, enrichment_index_hist_plots_list), ggsave,
      width = 400, height = 240, units = c("mm"), dpi = 320)


#perform co-enrichment analysis for pairs of marks and make the plots
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

length(list_gr_raw_comp_AS_HS0_DAS_heat_trt)
create_coenrichment_index_two_marks(list_gr_raw_comp_AS_HS0_DAS_heat_trt[[1]], 1)

names_list_gr_raw_comp_AS_HS0_DAS_heat_trt

coenrichment_index_table <- tibble()
list_comb_pair_epimark_heat_trt

for (i in seq_along(list_comb_pair_epimark_heat_trt)) {
  coenrichment_table_raw <- map2(list_gr_raw_comp_AS_HS0_DAS_heat_trt, rep(i, length(list_gr_raw_comp_AS_HS0_DAS_heat_trt)), create_coenrichment_index_two_marks) %>%
    Reduce(bind_rows, .) %>%
    mutate(feature = names_list_gr_raw_comp_AS_HS0_DAS_heat_trt)
  
  coenrichment_index_table <- bind_rows(coenrichment_index_table, coenrichment_table_raw)
}

write_delim(coenrichment_index_table, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_table_combined_bootstrap_samples", delim = "\t", col_names = TRUE) 

rep(comp_AS_HS0_DAS_heat_trt_label, 2)



nrow(coenrichment_index_table)/length(comp_AS_HS0_DAS_heat_trt_label)

AS_type_focused

coenrichment_index_table_for_plot <- 
coenrichment_index_table %>%
  mutate(AS_DAS_groups = rep(comp_AS_HS0_DAS_heat_trt_label, nrow(coenrichment_index_table)/length(comp_AS_HS0_DAS_heat_trt_label))) %>%
  mutate(AS_region = str_remove(feature, AS_DAS_groups)) %>%
  mutate(AS_region = str_remove(AS_region, "_gr_")) %>%
  mutate(heat_trt = if_else(str_detect(mark_1, heat_trt[[1]]), heat_trt[[1]],
                            if_else(str_detect(mark_1, heat_trt[[2]]), heat_trt[[2]], heat_trt[[3]]))) %>%
  mutate(pair_mark = str_c(mark_1, "_", str_remove(mark_2, str_c(heat_trt, "_")))) %>%
  mutate(AS_DAS_groups = factor(AS_DAS_groups, levels = comp_AS_HS0_DAS_heat_trt_label)) %>%
  dplyr::select(pair_mark, AS_region, coenrichment_index, AS_DAS_groups)


produce_coenrichment_index_hist_plots <- function(a, b){
  coenrichment_index_table_for_plot %>% 
    filter(AS_region == a[1] & pair_mark == b[1]) %>%
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
    facet_grid(AS_region ~ pair_mark) +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 20), legend.text = element_text(size = 9, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 20, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
          axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle = 60, vjust = 0.5))
}

produce_coenrichment_index_hist_plots(coenrichment_index_table_for_plot$AS_region[[1]], coenrichment_index_table_for_plot$pair_mark[[1]])

AS_region_label <- unique(coenrichment_index_table_for_plot$AS_region)
pair_mark_label <- unique(coenrichment_index_table_for_plot$pair_mark)
1440/8


coenrichment_index_table_for_plots_key <- 
  coenrichment_index_table_for_plot %>%
  ungroup() %>%
  dplyr::select(AS_region, pair_mark) %>%
  distinct()

coenrichment_index_hist_plots_list <- 
  map2(coenrichment_index_table_for_plots_key$AS_region, 
       coenrichment_index_table_for_plots_key$pair_mark,
       produce_coenrichment_index_hist_plots)

path_coenrichment_index_plots_hist_AS_DAS <- 
  str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_analysis_heat_comb_epimark_AS_DAS_histograms/", coenrichment_index_table_for_plots_key$AS_region, "_", coenrichment_index_table_for_plots_key$pair_mark, "_coenrichment_index_AS_DAS.jpeg")

pwalk(list(path_coenrichment_index_plots_hist_AS_DAS, coenrichment_index_hist_plots_list), ggsave,
      width = 400, height = 240, units = c("mm"), dpi = 320)


#perform co-enrichment analysis for three marks and make the plots
comb_three_epimark_heat_trt_raw <- combs(seq_along(mark_heat_trt), 3) 
 
comb_three_epimark_heat_trt_table <- 
  tibble(mark_1 = c(comb_three_epimark_heat_trt_raw[, 1], c(comb_three_epimark_heat_trt_raw[, 1]+5), c(comb_three_epimark_heat_trt_raw[, 1]+10)),
         mark_2 = c(comb_three_epimark_heat_trt_raw[, 2], c(comb_three_epimark_heat_trt_raw[, 2]+5), c(comb_three_epimark_heat_trt_raw[, 2]+10)),
         mark_3 = c(comb_three_epimark_heat_trt_raw[, 3], c(comb_three_epimark_heat_trt_raw[, 3]+5), c(comb_three_epimark_heat_trt_raw[, 3]+10)))

names_comb_three_epimark_heat_trt_table <- tibble()

list_comb_three_epimark_heat_trt <- vector("list", length = length(comb_three_epimark_heat_trt_table$mark_1))

for (i in seq_along(list_comb_three_epimark_heat_trt)) {
  
  list_comb_three_epimark_heat_trt[[i]][[1]] <- list_gr_mark_heat_trt[[comb_three_epimark_heat_trt_table[i, ]$mark_1]]
  list_comb_three_epimark_heat_trt[[i]][[2]] <- list_gr_mark_heat_trt[[comb_three_epimark_heat_trt_table[i, ]$mark_2]]
  list_comb_three_epimark_heat_trt[[i]][[3]] <- list_gr_mark_heat_trt[[comb_three_epimark_heat_trt_table[i, ]$mark_3]]
  
  name_mark_table <- 
    tibble(mark_1 = names_list_gr_mark_heat_trt[[comb_three_epimark_heat_trt_table[i, ]$mark_1]],
           mark_2 = names_list_gr_mark_heat_trt[[comb_three_epimark_heat_trt_table[i, ]$mark_2]],
           mark_3 = names_list_gr_mark_heat_trt[[comb_three_epimark_heat_trt_table[i, ]$mark_3]])
  
  names_comb_three_epimark_heat_trt_table <- bind_rows(names_comb_three_epimark_heat_trt_table, name_mark_table)
}


create_coenrichment_index_three_marks <- 
  function(data, a){
    mark_1_feature_intersection <- 
      join_overlap_intersect(list_comb_three_epimark_heat_trt[[a[1]]][[1]], data) %>%
      GenomicRanges::reduce()
    
    mark_2_feature_intersection <- 
      join_overlap_intersect(list_comb_three_epimark_heat_trt[[a[1]]][[2]], data) %>%
      GenomicRanges::reduce()
    
    mark_3_feature_intersection <- 
      join_overlap_intersect(list_comb_three_epimark_heat_trt[[a[1]]][[3]], data) %>%
      GenomicRanges::reduce()
    
    the_intersection_mark1_2 <-    
      join_overlap_intersect(mark_1_feature_intersection, mark_2_feature_intersection)
    
    sum_bp_coenrichment <- 
      join_overlap_intersect(the_intersection_mark1_2, mark_3_feature_intersection) %>%
      GenomicRanges::reduce() %>%
      as_tibble() %>%
      summarise(sum_bp_coenrichment = sum(end-start))
    
    mark1_sum_bp_intersection <- 
      as_tibble(mark_1_feature_intersection) %>%
      summarise(sum_bp_coenrichment = sum(end-start))
    
    mark2_sum_bp_intersection <- 
      as_tibble(mark_2_feature_intersection) %>%
      summarise(sum_bp_coenrichment = sum(end-start))
    
    mark3_sum_bp_intersection <- 
      as_tibble(mark_3_feature_intersection) %>%
      summarise(sum_bp_coenrichment = sum(end-start))
    
    tibble(
      mark_1 = names_comb_three_epimark_heat_trt_table[a[1],]$mark_1,
      mark_2 = names_comb_three_epimark_heat_trt_table[a[1],]$mark_2,
      mark_3 = names_comb_three_epimark_heat_trt_table[a[1],]$mark_3,
      coenrichment_index = sum_bp_coenrichment$sum_bp_coenrichment/c(mark1_sum_bp_intersection$sum_bp_coenrichment + mark2_sum_bp_intersection$sum_bp_coenrichment + mark2_sum_bp_intersection$sum_bp_coenrichment - 2*(sum_bp_coenrichment$sum_bp_coenrichment))
    )
  }


create_coenrichment_index_three_marks(list_gr_raw_comp_AS_HS0_DAS_heat_trt[[1]], 1)


coenrichment_index_table_3_mark <- tibble()
list_comb_pair_epimark_heat_trt

for (i in seq_along(list_comb_three_epimark_heat_trt)) {
  coenrichment_table_raw <- map2(list_gr_raw_comp_AS_HS0_DAS_heat_trt, rep(i, length(list_gr_raw_comp_AS_HS0_DAS_heat_trt)), create_coenrichment_index_three_marks) %>%
    Reduce(bind_rows, .) %>%
    mutate(feature = names_list_gr_raw_comp_AS_HS0_DAS_heat_trt)
  
  coenrichment_index_table_3_mark <- bind_rows(coenrichment_index_table_3_mark, coenrichment_table_raw)
}

write_delim(coenrichment_index_table_3_mark, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_table_3_marks_combined_bootstrap_samples", delim = "\t", col_names = TRUE) 


coenrichment_index_table_3_mark_for_plot <- 
  coenrichment_index_table_3_mark %>%
  mutate(AS_DAS_groups = rep(comp_AS_HS0_DAS_heat_trt_label, nrow(coenrichment_index_table_3_mark)/length(comp_AS_HS0_DAS_heat_trt_label))) %>%
  mutate(AS_region = str_remove(feature, AS_DAS_groups)) %>%
  mutate(AS_region = str_remove(AS_region, "_gr_")) %>%
  mutate(heat_trt = if_else(str_detect(mark_1, heat_trt[[1]]), heat_trt[[1]],
                            if_else(str_detect(mark_1, heat_trt[[2]]), heat_trt[[2]], heat_trt[[3]]))) %>%
  mutate(pair_mark = str_c(mark_1, "_", str_remove(mark_2, str_c(heat_trt, "_")), "_", str_remove(mark_3, str_c(heat_trt, "_")))) %>%
  mutate(AS_DAS_groups = factor(AS_DAS_groups, levels = comp_AS_HS0_DAS_heat_trt_label)) %>%
  dplyr::select(pair_mark, AS_region, coenrichment_index, AS_DAS_groups)


produce_coenrichment_index_3_marks_hist_plots <- function(a, b){
  coenrichment_index_table_3_mark_for_plot %>% 
    filter(AS_region == a[1] & pair_mark == b[1]) %>%
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
    facet_grid(AS_region ~ pair_mark) +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 20), legend.text = element_text(size = 9, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 20, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
          axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle = 60, vjust = 0.5))
}

produce_coenrichment_index_3_marks_hist_plots(coenrichment_index_table_3_mark_for_plot$AS_region[[1]], coenrichment_index_table_3_mark_for_plot$pair_mark[[1]])

AS_region_label <- unique(coenrichment_index_table_for_plot$AS_region)
pair_mark_label <- unique(coenrichment_index_table_for_plot$pair_mark)
1440/8


coenrichment_index_table_3_mark_for_plots_key <- 
  coenrichment_index_table_3_mark_for_plot %>%
  ungroup() %>%
  dplyr::select(AS_region, pair_mark) %>%
  distinct()

coenrichment_index_3_mark_hist_plots_list <- 
  map2(coenrichment_index_table_3_mark_for_plots_key$AS_region, 
       coenrichment_index_table_3_mark_for_plots_key$pair_mark,
       produce_coenrichment_index_3_marks_hist_plots)

path_coenrichment_index_3_mark_plots_hist_AS_DAS <- 
  str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_analysis_heat_three_epimark_AS_DAS_histograms/", coenrichment_index_table_3_mark_for_plots_key$AS_region, "_", coenrichment_index_table_3_mark_for_plots_key$pair_mark, "_coenrichment_index_3_mark_AS_DAS.jpeg")

pwalk(list(path_coenrichment_index_3_mark_plots_hist_AS_DAS, coenrichment_index_3_mark_hist_plots_list), ggsave,
      width = 400, height = 240, units = c("mm"), dpi = 320)
