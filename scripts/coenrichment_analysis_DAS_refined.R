Sys.getlocale()
Sys.setenv(LANG = "en_US.UTF-8")
Packages <- c("plyranges", "tidyverse", "doParallel", "foreach", "nullranges", "caTools", "ggpubr", "fs", "HelloRanges")
lapply(Packages, library, character.only = TRUE)


#Previously, our analysis focused on three regions (100-bp upstream/downstream junctions of exon/intron containing AS),
#Since We are more interested in the whole cassettes (intron-AS-included-exon-intron or exon-AS-included-intron-exon)
#We will modify how to define the group of upstream_gr and downstream_gr here
#Plus, we also include the DAS gene analysis, so we will extract the information of DAS gene first.
#Then the DAS info will be used to filter out PSI_all and control groups to keep these two groups containing only PSI-stable or constitutively spliced genes in response to heat treatment
#Since this task is focusing on DAS, I only use 5 epigenetic marks that we have under heat treatment for the analysis

#import the file of refined DAS genes
All_DAS_events_all_comparisons_refined <-
  read_delim("./data/rMATS_out/All_DAS_events_all_comparisons_refined.tsv", delim = "\t", col_names = TRUE) 

All_DAS_events_all_comparisons_refined %>%
  View()

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

DAS_gene <- 
  tibble(
    geneID = DAS_geneID
    )

write_delim(DAS_gene, "./analysis/DAS_list_gene_ID.txt", delim = "\t", col_names = FALSE)

#Create the function that can use the content in the DAS-gene-ID vector for producing a table
search_for_strand_DAS_gene <- function(a){
  M82_rMATs_anno_all_gene %>%
    filter(str_detect(gene_name, a[1])) %>%
    mutate(gene_name_f = a[1]) %>%
    dplyr::select(GeneID = gene_name_f, strand) %>%
    distinct()
}

#test the function above
search_for_strand_DAS_gene(DAS_geneID[[1]])

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

#define str and end of DAS 
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

#extract inclusion levels of DAS
IncLevel_two_trts_rep_table <- 
  tibble(
  IncLevel1_r1 = as.double(str_split(All_DAS_events_all_comparisons_refined_raw$IncLevel1, ",", simplify = TRUE)[,1]),
  IncLevel1_r2 = as.double(str_split(All_DAS_events_all_comparisons_refined_raw$IncLevel1, ",", simplify = TRUE)[,2]),
  IncLevel2_r1 = as.double(str_split(All_DAS_events_all_comparisons_refined_raw$IncLevel2, ",", simplify = TRUE)[,1]),
  IncLevel2_r2 = as.double(str_split(All_DAS_events_all_comparisons_refined_raw$IncLevel2, ",", simplify = TRUE)[,2])
)

#produce the table for granges
All_DAS_events_all_comparisons_refined_for_granges <- 
All_DAS_events_all_comparisons_refined_raw %>%
  bind_cols(IncLevel_two_trts_rep_table) %>%
  mutate(IncLevel_former_tp = if_else(comp == "HS1.vs.HS6", (IncLevel1_r1 + IncLevel1_r2)/2, (IncLevel2_r1 + IncLevel2_r2)/2),
         IncLevel_latter_tp = if_else(comp == "HS1.vs.HS6", (IncLevel2_r1 + IncLevel2_r2)/2, (IncLevel1_r1 + IncLevel1_r2)/2)) %>%
  mutate(Incdf_f = IncLevel_latter_tp - IncLevel_former_tp) %>%
  mutate(more_AS = if_else(Incdf_f > 0, "yes", "no")) %>%
  mutate(comp = str_replace(comp, ".vs.", "_")) %>%
  mutate(comp = if_else(comp == "HS1_HS0", "HS0_HS1", if_else(comp == "HS6_HS0", "HS0_HS6", comp))) %>%
  dplyr::select(seqnames = chr, start = str, end, strand, AS_type, comp, Incdf = Incdf_f, IncLevel_former_tp, IncLevel_latter_tp, more_AS) %>%
  arrange(seqnames, start)


write_delim(All_DAS_events_all_comparisons_refined_for_granges, "./data/temp/All_DAS_events_all_comparisons_refined_for_granges", delim = "\t", col_names = TRUE)
  
All_DAS_events_all_comparisons_refined_for_granges %>%
  as_granges()

#-----delete this part-----#
#import bootstrap samples under HS0 condition, and remove DAS genes from these samples
#and stored the refined HS0 bootstrap samples in another sub-folder
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
#AS_DAS_refined_bootstrap_samples
#-----delete this part-----#

#Based on 1000 bootstrap samples previously produced in co-enrichment analysis for only HS0 events
#we focused only on feature_gr to combine them in one sample from which I extract their upstream and downstream genomic features
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

#combined 1-1000 bootstrap samples into one sample
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

#calculate the medium of exon sizes
M82_rMATs_anno_all_exon %>%
  summarise(median_exon = median(end-start))


#import all introns
M82_rMATs_anno_all_intron <-
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_intron.bed", col_names = c("seqnames", "start", "end", "strand", "feature", "source", "gene_name"))

#calculate the medium of intron sizes
M82_rMATs_anno_all_intron %>%
  summarise(median_intron = median(end-start))


#combine the annotation data of intron and exon
M82_rMATs_anno_all_intron_exon_ordered <- 
  bind_rows(M82_rMATs_anno_all_exon, M82_rMATs_anno_all_intron) %>%
  arrange(gene_name, seqnames, start) %>%
  group_by(gene_name) %>%
  mutate(genetic_feature_label = if_else(strand == "+", c(1:n()), c(n():1))) %>%
  ungroup()

M82_rMATs_anno_all_intron_exon_ordered %>%
  View()

write_delim(M82_rMATs_anno_all_intron_exon_ordered, "./data/M82_annotation_data/M82_rMATs_anno_all_intron_exon_ordered.bed", delim = "\t", col_names = TRUE)

#create the function that can create genomic features upstream or downstream AS-included introns/exons
#this function below can be used for extracting flanking genomic features of PSI-stable, control and PSI-changing exons/introns
produce_flanking_genomic_feature <- function(data) {
  genetic_feature_raw <- 
    data %>%
    left_join(M82_rMATs_anno_all_intron_exon_ordered) %>%
    dplyr::select(strand, gene_name, genetic_feature_label)

  upstream_genetic_feature_raw <-
    genetic_feature_raw %>%
    mutate(genetic_feature_label = if_else(strand=="+", genetic_feature_label - 1, genetic_feature_label + 1)) %>%
    mutate(location_gr = "upstream")
  
  downstream_genetic_feature_raw <-   
    genetic_feature_raw %>%
    mutate(genetic_feature_label = if_else(strand=="+", genetic_feature_label + 1, genetic_feature_label - 1)) %>%
    mutate(location_gr = "downstream")


#create a function that can create upstream/downstream genomic feature flanking AS-included exon/intron at the same time
produce_flanking_genomic_feature_inside <- function(data) {

  flanking_genomic_feature_raw <- 
    M82_rMATs_anno_all_intron_exon_ordered %>%
    left_join(data) %>%
    drop_na() %>%
    dplyr::select(seqnames, start, end, strand, feature) %>%
    distinct() %>%  
    as_granges() %>%
    GenomicRanges::reduce(ignore.strand=TRUE) 
  
  seqlevels(flanking_genomic_feature_raw) <- str_c("chr", c(1:12))
  
  
  sort(flanking_genomic_feature_raw) %>%
    as_tibble()

}

map(list(upstream_genetic_feature_raw, downstream_genetic_feature_raw), produce_flanking_genomic_feature_inside)

}

#create upstream/downstream genomic features for PSI_all and controls (no AS)
for (i in seq_along(AS_type_bootstrap)) {
  for (j in seq_along(PSI_type_simplified)) {
    for (l in seq_along(location_gr_AS_included_body)) {
      
      location_gr_flanking_regions <- c("upstream_feature_gr", "downstream_feature_gr")
      
      path_output <- str_c("./data/AS_HS0_bootstrap_samples_DAS_removed/", AS_type_bootstrap[[i]], "_", PSI_type_simplified[[j]], "_", location_gr_flanking_regions, "_bootstrap_DAS_removed_1_1000_combined")
      
      flanking_genomic_feature_output <- produce_flanking_genomic_feature(PSI_DAS_removed_reduced[[i]][[j]][[l]])
      
      pwalk(list(flanking_genomic_feature_output, path_output), write_delim, delim = "\t", col_names = TRUE)
      
    }
  }
}

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
#then I extract their flanking regions and store them  
heat_trt_comp <- c("HS0_HS1", "HS0_HS6", "HS1_HS6")
Inc_AS_aft_trt <- c("inc", "sk")

for (i in seq_along(AS_type_bootstrap)) {
  for (j in seq_along(heat_trt_comp)) {
    for (k in seq_along(Inc_AS_aft_trt)) {
      location_gr <- c("upstream_feature_gr", "downstream_feature_gr", "feature_gr")
      if(i == 1){
        feature_table <- produce_feature_with_AS(M82_rMATs_anno_all_intron, All_DAS_events_all_comparisons_refined_for_granges[[i]][[j]][[k]])
        
        list_flanking_feature_table <- produce_flanking_genomic_feature(feature_table)
        
        path_output <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_bootstrap[[i]], "_", location_gr, "_", heat_trt_comp[[j]], "_more_", Inc_AS_aft_trt[[k]], "_AS_aft_trt")
        
        list_output <- append(list_flanking_feature_table, list(feature_table))
          
        pwalk(list(list_output, path_output), write_delim, delim = "\t", col_names = TRUE)  
      } else {
        feature_table <- produce_feature_with_AS(M82_rMATs_anno_all_exon, All_DAS_events_all_comparisons_refined_for_granges[[i]][[j]][[k]])
        
        list_flanking_feature_table <- produce_flanking_genomic_feature(feature_table)
        
        path_output <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_bootstrap[[i]], "_", location_gr, "_", heat_trt_comp[[j]], "_more_", Inc_AS_aft_trt[[k]], "_AS_aft_trt")
        
        list_output <- append(list_flanking_feature_table, list(feature_table))
        
        pwalk(list(list_output, path_output), write_delim, delim = "\t", col_names = TRUE)
      } 
      
      
      
    }
  }
}

#Even though I separated DAS as different groups previously for their outputs in terms of the trend of changing PSI and heat treatment
#here I combined H0/H1 and H0/H6 together for the simplicity of co-enrichment index calculation and their visualization
location_gr <- c("upstream_feature_gr", "downstream_feature_gr", "feature_gr")
comp_DAS_heat_trt_label <- c(str_c("HS0_HS6_more_", Inc_AS_aft_trt), str_c("HS1_HS6_more_", Inc_AS_aft_trt))

DAS_combined_table_PSI_trend_heat_trt <- tibble()

#import separate files into a merged table
for (i in seq_along(AS_type_bootstrap)) {
  for (j in seq_along(location_gr)) {
    for (k in seq_along(comp_DAS_heat_trt_label)) {
    
      path_input <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_bootstrap[[i]], "_", location_gr[[j]], "_", comp_DAS_heat_trt_label[[k]], "_AS_aft_trt")
      
      table_input <- read_delim(path_input, delim = "\t", col_names = TRUE) %>%
        mutate(AS_type = AS_type_bootstrap[[i]], location = location_gr[[j]])
      
      DAS_combined_table_PSI_trend_heat_trt <- bind_rows(DAS_combined_table_PSI_trend_heat_trt, table_input)
    }
  }
}

#produce the output
DAS_combined_table_PSI_trend_heat_trt %>%
  filter(AS_type == AS_type_bootstrap[[1]] & location == location_gr[[1]])

for (i in seq_along(AS_type_bootstrap)) {
  for (j in seq_along(location_gr)) {
      path_output <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_bootstrap[[i]], "_", location_gr[[j]], "_merged_H16_aft_H0")
      
      table_output <- 
        DAS_combined_table_PSI_trend_heat_trt %>%
        filter(AS_type == AS_type_bootstrap[[i]] & location == location_gr[[j]]) %>%
        dplyr::select(-AS_type, -location)
      
      write_delim(table_output, path_output, delim = "\t", col_names = TRUE)
    
  }
}

#install HelloRanges(this package can allow us to use r function to simulate functions in bedtools)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

BiocManager::install("HelloRanges")
library(HelloRanges)

#since we only have two marks with two replicates
#we focus on replicate 1 for all analysis
#but we still merge two replicates of two epigenetic marks for the future analysis
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
#this is the old one that includes merged epigenetic marks
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

#the new one to include only replicate 1 of epigenetic marks
for (i in seq_along(heat_trt)) {
  for (j in seq_along(mark_heat_trt)) {
    
      input_gr_mark <- read_delim(str_c("./data/peak_epimarks/heat_trt/M82_", heat_trt[[i]], "_", mark_heat_trt[[j]], "_rep1_p0.05_peaks.bed"), 
                                  delim = "\t", col_names = c("seqnames", "start", "end", "signal")) %>%
        as_granges()
      
      list_gr_mark_heat_trt <- append(list_gr_mark_heat_trt, list(input_gr_mark))
      
      names_list_gr_mark_heat_trt <- append(names_list_gr_mark_heat_trt, str_c(heat_trt[[i]], "_", mark_heat_trt[[j]]))
    
  }
}


names(list_gr_mark_heat_trt) <- names_list_gr_mark_heat_trt


#import PSI combined, ctrl (no AS), DAS (merging H0_H1 and H0_H6 for more inc or sk) as lists
comp_AS_HS0_DAS_heat_trt_label <- c("PSI_all", "ctrl", "merged_H16_aft_H0")
AS_type_focused <- c("RI", "SE")
location_gr <- c("upstream_feature_gr", "feature_gr", "downstream_feature_gr")

list_gr_raw_comp_AS_HS0_DAS_heat_trt <- list()
names_list_gr_raw_comp_AS_HS0_DAS_heat_trt <- c()

for (i in seq_along(AS_type_focused)) {
  for (j in seq_along(location_gr)) {
    for (k in seq_along(comp_AS_HS0_DAS_heat_trt_label)) {
      
      if(k %in% c(1, 2)){
        
        path_input <- str_c("./data/AS_HS0_bootstrap_samples_DAS_removed/", AS_type_focused[[i]], "_", comp_AS_HS0_DAS_heat_trt_label[[k]], "_", location_gr[[j]], "_bootstrap_DAS_removed_1_1000_combined")
        
        name_input <- str_c(AS_type_focused[[i]], "_", location_gr[[j]], "_", comp_AS_HS0_DAS_heat_trt_label[[k]])
        
        names_list_gr_raw_comp_AS_HS0_DAS_heat_trt <- append(names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, name_input)
        
        gr_sample <- read_delim(path_input, delim = "\t", col_names = TRUE) %>%
          as_granges() %>%
          GenomicRanges::reduce(ignore.strand=TRUE) 
        
        list_gr_raw_comp_AS_HS0_DAS_heat_trt <- append(list_gr_raw_comp_AS_HS0_DAS_heat_trt, list(gr_sample))
      } else {
        
        path_input <- str_c("./data/DAS_refined_bootstrap_and_combined/DAS_", AS_type_focused[[i]], "_", location_gr[[j]], "_", comp_AS_HS0_DAS_heat_trt_label[[k]])
        
        name_input <- str_c(AS_type_focused[[i]], "_", location_gr[[j]], "_", comp_AS_HS0_DAS_heat_trt_label[[k]])
        
        names_list_gr_raw_comp_AS_HS0_DAS_heat_trt <- append(names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, name_input)
        
        gr_sample <- read_delim(path_input, delim = "\t", col_names = TRUE) %>%
          as_granges() %>%
          GenomicRanges::reduce(ignore.strand=TRUE) 
        
        list_gr_raw_comp_AS_HS0_DAS_heat_trt <- append(list_gr_raw_comp_AS_HS0_DAS_heat_trt, list(gr_sample))
      }
      
      
    }
  }
}

names(list_gr_raw_comp_AS_HS0_DAS_heat_trt) <- names_list_gr_raw_comp_AS_HS0_DAS_heat_trt

list_gr_raw_comp_AS_HS0_DAS_heat_trt[[5]]

#make the gene list of DAS, AS, and no-AS gene
str_match_all()

M82_rMATs_anno_all_gene_modified <- 
M82_rMATs_anno_all_gene %>%
  filter(source != "REGARN") %>%
  dplyr::select(seqnames = chr, start = str, end, gene_name) %>%
  mutate(gene_name_m = if_else(str_detect(gene_name, "Solyc")==TRUE, str_sub(gene_name, 1, 16), str_sub(gene_name, 1, 100))) %>%
  dplyr::select(seqnames, start, end, gene_name = gene_name_m) %>%
  as_granges()

extract_gene_name_AS_DAS <- function(data){
  find_overlaps(M82_rMATs_anno_all_gene_modified, data) %>%
    as_tibble() %>%
    dplyr::select(gene_name) %>%
    distinct()
}

#create gene list of changed and stable PSI
gene_list_RI_stable_PSI_all <- 
  extract_gene_name_AS_DAS(list_gr_raw_comp_AS_HS0_DAS_heat_trt[[4]])

gene_list_RI_changed_PSI <- 
  extract_gene_name_AS_DAS(list_gr_raw_comp_AS_HS0_DAS_heat_trt[[6]])


gene_list_SE_stable_PSI_all <- 
  extract_gene_name_AS_DAS(list_gr_raw_comp_AS_HS0_DAS_heat_trt[[13]])

gene_list_SE_changed_PSI <- 
  extract_gene_name_AS_DAS(list_gr_raw_comp_AS_HS0_DAS_heat_trt[[15]])


write_delim(gene_list_RI_stable_PSI_all, "./analysis/gene_list_RI_stable_PSI_all.txt", delim = "\t", col_names = FALSE)  
write_delim(gene_list_RI_changed_PSI, "./analysis/gene_list_RI_changed_PSI.txt", delim = "\t", col_names = FALSE)  
write_delim(gene_list_SE_stable_PSI_all, "./analysis/gene_list_SE_stable_PSI_all.txt", delim = "\t", col_names = FALSE)  
write_delim(gene_list_SE_changed_PSI, "./analysis/gene_list_SE_changed_PSI.txt", delim = "\t", col_names = FALSE)  


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

#since we would like to check if size of features could influence profiles of epigenetic features
#check quantiles of size for these AS events
AS_event_size_check <- 
function(data, data2) {
data %>%
  as_tibble() %>%
  mutate(name = data2) %>%
  mutate(size_group = if_else(width <= quantile(width, 0.33), "small",
                              if_else(width > quantile(width, 0.66), "large", "medium"))) %>%
  group_by(name, size_group) %>%
  summarise(median_size_group = median(width))
}


map2(list_gr_raw_comp_AS_HS0_DAS_heat_trt, names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, AS_event_size_check) %>%
  bind_rows() %>%
  View()

AS_event_size_check_quantile <- 
function(data, data2) {
  data_t <- data %>%
    as_tibble()
  
  tibble(name = data2, 
         small = quantile(data_t$width, 0.33),
         large = quantile(data_t$width, 0.66))
}

map2(list_gr_raw_comp_AS_HS0_DAS_heat_trt, names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, AS_event_size_check_quantile) %>%
  bind_rows() %>%
  View()

transform_AS_DAS_granges_into_tibble <- 
  function(data, data2) {
    data_t <- data %>%
      as_tibble() %>%
      mutate(name = data2)
  }

#combine the list of feature_gr of RI/SE as a list
combined_table_AS_ctrl_DAS_feature_body <-
  map2(list_gr_raw_comp_AS_HS0_DAS_heat_trt, names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, transform_AS_DAS_granges_into_tibble) %>%
  bind_rows() %>%
  filter(str_detect(name, "RI_feature_gr") | str_detect(name, "SE_feature_gr"))

#produce list of AS/DAS groups based on size selection using two cutting points making three quantiles of DAS
produce_list_AS_DAS_size_selection <- 
function(a){
DAS_feature_body_tibble <- 
  combined_table_AS_ctrl_DAS_feature_body %>%
  filter(str_detect(name, a[1]) & str_detect(name, "merge")) 

combined_table_AS_ctrl_DAS_feature_body %>%
  filter(str_detect(name, a[1])) %>%
  mutate(size_first_point = quantile(DAS_feature_body_tibble$width, 0.33),
         size_second_point = quantile(DAS_feature_body_tibble$width, 0.67)) %>%
  mutate(size_group = if_else(width <= size_first_point, "small",
                              if_else(width > size_second_point, "large", "medium"))) %>%
  mutate(AS_DAS_group = str_remove(name, str_c(str_sub(name, 1, 2), "_feature_gr_"))) %>%
  group_by(AS_DAS_group, size_group) %>%
  ungroup() %>%
  split(.$AS_DAS_group) %>%
  map(. %>% split(.$size_group))
}

list_RI_AS_DAS_size_selection <- produce_list_AS_DAS_size_selection("RI")
list_SE_AS_DAS_size_selection <- produce_list_AS_DAS_size_selection("SE")

list_AS_DAS_size_selection <- 
  map(AS_type_focused, produce_list_AS_DAS_size_selection)

for (i in seq_along(list_AS_DAS_size_selection)) {
  for (j in seq_along(list_AS_DAS_size_selection[[1]])) {
    for (k in seq_along(list_AS_DAS_size_selection[[1]][[1]])) {
      
      path_output <- str_c("./data/AS_DAS_refined_grouped_by_DAS_size/", AS_type_focused[[i]], "_", names(list_AS_DAS_size_selection[[i]])[[j]], "_feature_size_group_", 
            names(list_AS_DAS_size_selection[[i]][[j]])[[k]], ".bed")
      print(path_output)
      #write_delim(list_AS_DAS_size_selection[[i]][[j]][[k]], path_output, delim = "\t", col_names = FALSE)
      
    }
  }
}

#based on the similar idea using deeptool to check profiles of epigenetic marks of DAS, AS (stable PSI), and no-AS genes
#we focused on genomic features in large genes, and make profiles according to the junction of splicing sites instead of the boundaries of each genomic feature
#import all AS (SE/RI) found in HS0

path_input_AS_HS0 <- str_c("./data/AS_events_HS0/AS_HS0_PSI_segment/AS_HS0_", AS_type_focused, ".bed")

AS_HS0_raw <- 
  map(path_input_AS_HS0, read_delim, delim = "\t", col_names = c("seqnames", "start", "end", "trt", "AS_type", "PSI_type")) %>%
  map(. %>% as_granges)

#create a table with all DAS (from line 315)
All_DAS_events_all_comparisons_refined_for_granges_combined <- 
All_DAS_events_all_comparisons_refined_for_granges %>%
  map(. %>% map(. %>% bind_rows)) %>%
  map(. %>% bind_rows) %>%
  bind_rows() %>%
  as_granges()

#remove DAS from all AS found in HS0
AS_HS0_PSI_all_DAS_removed <- 
  AS_HS0_raw %>%
  map(. %>% mutate(n_overlaps = count_overlaps(., All_DAS_events_all_comparisons_refined_for_granges_combined))) %>%
  map(. %>% filter(n_overlaps == 0))

#extract only large-size RI or SE containing stable-PSI AS 
list_AS_PSI_all_DAS_removed_feature <- list(list_AS_DAS_size_selection[[1]]$PSI_all$large, list_AS_DAS_size_selection[[2]]$PSI_all$large)

#create the function that can extract splicing sites of AS events from large genes (specifically, large exons or introns)
extract_AS_splicing_sites_large_genomic_feature <- 
  function(data, data2){
    data %>%
      mutate(n_overlaps = count_overlaps(., as_granges(data2))) %>%
      filter(n_overlaps != 0)
  }

#produce splicing sites of stable-PSI AS events from large genomic features
AS_splicing_sites_HS0_PSI_all_DAS_removed_granges <- 
  map2(AS_HS0_PSI_all_DAS_removed, list_AS_PSI_all_DAS_removed_feature, extract_AS_splicing_sites_large_genomic_feature)


#extract only large-size RI or SE containing changing-PSI AS 
list_DAS_changed_PSI_feature <- list(list_AS_DAS_size_selection[[1]]$merged_H16_aft_H0$large, list_AS_DAS_size_selection[[2]]$merged_H16_aft_H0$large)

#produce splicing sites of changed-PSI AS events from large genomic features
DAS_splicing_sites_merged_H16_aft_H0_granges <-
  map2(list(All_DAS_events_all_comparisons_refined_for_granges_combined, All_DAS_events_all_comparisons_refined_for_granges_combined), list_DAS_changed_PSI_feature, extract_AS_splicing_sites_large_genomic_feature) %>%
  map(. %>% filter(comp != "HS1_HS6"))


#combine AS and DAS with the information of splicing sites together
list_AS_DAS_splicing_site_large_events <- 
  append(AS_splicing_sites_HS0_PSI_all_DAS_removed_granges, DAS_splicing_sites_merged_H16_aft_H0_granges) %>%
  map(. %>% as_tibble())

#produce the output of AS splicing sites
AS_splicing_site_path_output <- str_c("./data/AS_DAS_refined_grouped_by_DAS_size/", rep(AS_type_focused, 2), "_splicing_sites", rep(c("_PSI_all_feature_size_group_large", "_merged_H16_aft_H0_feature_size_group_large"), each = 2), ".bed")



pwalk(list(list_AS_DAS_splicing_site_large_events, AS_splicing_site_path_output), write_delim,
      delim = "\t", col_names = FALSE)


#AS_DAS_refined_grouped_by_DAS_size

#produce bootstrap samples for the following analysis (sample size = 0.5 * the size of initial data)

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

#list_gr_raw_comp_AS_HS0_DAS_heat_trt[[1]] data

#the idea of coenrichment index is based on the intersection of each epigenetic mark and AS event (followingly designated as AS-mark intersection)
#we used the intersection between two AS-mark intersections divided by the union of them
create_coenrichment_index_two_marks <- 
function(data, a){
mark_1_feature_intersection <- 
join_overlap_intersect(list_comb_pair_epimark_heat_trt[[a[1]]][[1]], data) %>%
  GenomicRanges::reduce()

nrow_mark_1_feature_intersection <- nrow(as_tibble(mark_1_feature_intersection))

mark_2_feature_intersection <- 
join_overlap_intersect(list_comb_pair_epimark_heat_trt[[a[1]]][[2]], data) %>%
  GenomicRanges::reduce()

nrow_mark_2_feature_intersection <- nrow(as_tibble(mark_2_feature_intersection))

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
  coenrichment_index = if_else(nrow_mark_1_feature_intersection==0 | nrow_mark_2_feature_intersection==0, 0, sum_bp_coenrichment$sum_bp_coenrichment/c(mark1_sum_bp_intersection$sum_bp_coenrichment + mark2_sum_bp_intersection$sum_bp_coenrichment - sum_bp_coenrichment$sum_bp_coenrichment))
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

#To create boxplot, each AS event should be separated for calculating co-enrichment index
#create the function to separate whole lists of DAS, AS, no-AS genes into separated events
separate_AS_DAS_ctrl_list_into_each_event <- 
  function(data){
    data %>%
      as_tibble() %>%
      mutate(AS_label = c(1:n())) %>%
      as_granges()
  }

#create the above-mentioned separated events
list_gr_raw_comp_AS_HS0_DAS_heat_trt_separated_events <- 
  list_gr_raw_comp_AS_HS0_DAS_heat_trt %>%
  map(. %>% separate_AS_DAS_ctrl_list_into_each_event)

#create function for calculating coenrichment index for separated events but in a single table
create_coenrichment_index_two_marks_separated_events <- 
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
      dplyr::select(AS_label, AS_label_mark_1, mark1_sum_bp_intersection = width)
    
    
    #
    mark_2_sum_bp_intersection <-
      mark_2_feature_intersection %>%
      as_tibble() %>%
      dplyr::select(AS_label, AS_label_mark_2, mark2_sum_bp_intersection = width)
    
    
    data %>%
      as_tibble() %>%
      mutate(mark_1 = names_comb_pair_epimark_heat_trt_table[a[1],]$mark_1,
             mark_2 = names_comb_pair_epimark_heat_trt_table[a[1],]$mark_2) %>%
      left_join(sum_bp_coenrichment_test) %>%
      left_join(mark_1_sum_bp_intersection_test) %>%
      left_join(mark_2_sum_bp_intersection_test)%>%
      replace_na(list(AS_label_mark_1 = "no", AS_label_mark_2 = "no", two_mark_sum_up_intersection = 0, mark2_sum_bp_intersection = 0, mark1_sum_bp_intersection = 0)) %>%
      mutate(coenrichment_index = if_else(mark1_sum_bp_intersection == 0 | mark2_sum_bp_intersection == 0, 0, two_mark_sum_up_intersection/(mark1_sum_bp_intersection + mark2_sum_bp_intersection - two_mark_sum_up_intersection)))
  }

for (i in seq_along(list_comb_pair_epimark_heat_trt)) {
  for (j in seq_along(list_gr_raw_comp_AS_HS0_DAS_heat_trt_separated_events)) {
    coenrichment_table_raw <- create_coenrichment_index_two_marks_separated_events(list_gr_raw_comp_AS_HS0_DAS_heat_trt_separated_events[[j]], i) %>%
      mutate(feature = names_list_gr_raw_comp_AS_HS0_DAS_heat_trt[[j]])
    
    path_output <- str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/AS_DAS_separated_events_refined_coenrichment_result_2/", names_list_gr_raw_comp_AS_HS0_DAS_heat_trt[[j]], "_mark_comb_", i, ".txt")
    
    write_delim(coenrichment_table_raw, path_output, delim = "\t", col_names = TRUE)
  }
}

#import single table of co-enrichment index for separated events
coenrichment_index_separated_events_path_input <- dir_ls("./analysis/AS_RI_SE_epimark_coenrichment_analysis/AS_DAS_separated_events_refined_coenrichment_result_2/", glob = "*.txt")

length(coenrichment_index_separated_events_path_input)

list_coenrichment_index_separated_events_path_input <- list()

str_remove(coenrichment_index_separated_events_path_input, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/AS_DAS_separated_events_refined_coenrichment_result_2/")

grep("RIupstream_feature_gr", str_remove(coenrichment_index_separated_events_path_input, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/AS_DAS_separated_events_refined_coenrichment_result_2/")
)

test <- 
read_delim(coenrichment_index_separated_events_path_input[[1]], delim = "\t", col_names = TRUE)



test %>%
  mutate(comb = str_c(mark_1, str_remove(mark_2, ".H"))) %>%
  ggplot(aes(x=comb, y=coenrichment_index)) + 
  geom_boxplot()

for (i in seq_along(coenrichment_index_separated_events_path_input)) {
  
  read_delim()
  
}

comp_AS_HS0_DAS_heat_trt_label <- c("PSI_all", "ctrl", "merged_H16_aft_H0")
AS_type_focused <- c("RI", "SE")
location_gr <- c("upstream_feature_gr", "feature_gr", "downstream_feature_gr")


coenrichment_index_violin_plot_list <- list()

for (i in seq_along(AS_type_focused)) {
  for (j in seq_along(location_gr)) {
    for (k in seq_along(1:30)) {
      
      path_input <- str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/AS_DAS_separated_events_refined_coenrichment_result_2/", AS_type_focused[[i]], "_", location_gr[[j]], "_", comp_AS_HS0_DAS_heat_trt_label, "_mark_comb_", k,".txt")
      
      combined_table <-
        map(path_input, read_delim, delim = "\t", col_names = TRUE) %>%
        Reduce(bind_rows, .)
      
      final_table <- 
      combined_table %>%
        mutate(comb = str_c(mark_1, str_remove(mark_2, ".H"))) %>%
        mutate(AS_region = str_remove(str_c(AS_type_focused[[i]], "_", location_gr[[j]]), "_gr")) %>%
        mutate(AS_DAS_groups = str_remove(feature, str_c(AS_region, "_gr_"))) %>%
        mutate(AS_DAS_groups = if_else(AS_DAS_groups=="PSI_all", "AS_HS0_stable_PSI",
                                       if_else(AS_DAS_groups=="merged_H16_aft_H0", "DAS_changed_PSI_aft_HS0", "ctrl"))) %>%
        mutate(AS_DAS_groups = factor(AS_DAS_groups, levels = c("AS_HS0_stable_PSI", "ctrl", "DAS_changed_PSI_aft_HS0"))) 
      
      final_plot <- 
        final_table %>%
        ggplot(aes(x=AS_DAS_groups, y=coenrichment_index, fill = AS_DAS_groups)) + 
        #geom_boxplot() +
        geom_violin() +
        facet_grid(AS_region ~ comb) +
        scale_fill_manual(
          values = c(
            "#29f600",
            "#F6BE00",
            "#fb6a4a"
          )
        ) +
        theme(strip.text.x = element_text(colour = "black", face = "bold", size = 20), legend.text = element_text(size = 9, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
              legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
              axis.text.y = element_text(color = "black", size = 20, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
              axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle = 60, vjust = 0.5))
      
      path_coenrichment_index_plots_boxplot_AS_DAS <- 
        str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_analysis_heat_comb_epimark_AS_combined_DAS_boxplots/", AS_type_focused[[i]], "_", location_gr[[j]], "_", final_table$comb[[1]], "_coenrichment_index_AS_DAS_boxplot.jpeg")
      
      
      ggsave(path_coenrichment_index_plots_boxplot_AS_DAS, final_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
      
    }
  }
}

?ggsave



#perform bootstrap (50% for 1000 times, for creating error bars)
registerDoParallel(cores = 5)
getDoParWorkers()

names_list_gr_raw_comp_AS_HS0_DAS_heat_trt[[1]]

foreach(i=seq_along(list_comb_pair_epimark_heat_trt), .packages = c("plyranges", "tidyverse")) %:%
  foreach(k=1:1000, .packages = c("plyranges", "tidyverse")) %dopar% {
    
    
    list_bt_sample <- map(list_gr_raw_comp_AS_HS0_DAS_heat_trt, produce_bootstrap_sample)
    
    bt_sample_path <- str_c("./data/AS_DAS_refined_bootstrap_samples/", names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, "_bootstrap_sample_", k)
    
    pwalk(list(list_bt_sample, bt_sample_path), write_delim, delim = "\t", col_names = TRUE)
    
    #coenrichment_table_raw <- map2(list_bt_sample, rep(i, length(list_bt_sample)), create_coenrichment_index_two_marks) %>%
      #Reduce(bind_rows, .) %>%
      #mutate(feature = names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, bt = k)
    
    #write_delim(coenrichment_table_raw, str_c("./analysis/AS_DAS_refined_bootstrap_result/coenrichment_index_bt_",k), delim = "\t", col_names = TRUE)
    
  }


foreach(i=seq_along(list_comb_pair_epimark_heat_trt), .packages = c("plyranges", "tidyverse")) %:%
  foreach(k=1:1000, .packages = c("plyranges", "tidyverse")) %dopar% {
    
    bt_sample_input_path <- str_c("./data/AS_DAS_refined_bootstrap_samples/", names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, "_bootstrap_sample_", k)
    
    list_bt_sample <- map(bt_sample_input_path, read_delim, delim = "\t", col_names = TRUE) %>%
      map(. %>% as_granges)
    
    coenrichment_table_raw <- 
    map2(list_bt_sample, rep(i, length(list_bt_sample)), create_coenrichment_index_two_marks) %>%
      Reduce(bind_rows, .) %>%
      mutate(feature = names_list_gr_raw_comp_AS_HS0_DAS_heat_trt, bt = k)
    
    path_output <- str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/AS_DAS_refined_coenrichment_result/coenrichment_index_pair_mark_", i, "_bt_", k, ".txt")
    
    write_delim(coenrichment_table_raw, path_output, delim = "\t", col_names = TRUE)
    
  }


#import coenrichment index result

coenrichment_index_result_list <- 
foreach(i=seq_along(list_comb_pair_epimark_heat_trt), .packages = c("plyranges", "tidyverse")) %:%
  foreach(k=1:1000, .packages = c("plyranges", "tidyverse")) %dopar% {
    
    
    path_input <- str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/AS_DAS_refined_coenrichment_result/coenrichment_index_pair_mark_", i, "_bt_", k, ".txt")
    
    read_delim(path_input, delim = "\t", col_names = TRUE)
    
  }

#with error bar (using bt samples)
coenrichment_index_result_for_plot <- 
coenrichment_index_result_list %>%
  map(. %>% bind_rows) %>%
  bind_rows() %>%
  group_by(feature, mark_1, mark_2) %>%
  summarise(mean_coenrichment_index = mean(coenrichment_index), sd_coenrichment_index = sd(coenrichment_index), count = n())

write_delim(coenrichment_index_result_for_plot, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_table_AS_combined_DAS_1000_bootstrap_samples.txt", delim = "\t", col_names = TRUE)

coenrichment_index_result_for_plot_f <- 
coenrichment_index_result_for_plot %>%
  ungroup() %>%
  mutate(AS_type = if_else(str_detect(feature, "RI"), "RI", "SE"),
         region = if_else(str_detect(feature, "upstream_feature"), "upstream_feature", 
                          if_else(str_detect(feature, "downstream_feature"), "downstream_feature", "feature"))) %>%
  mutate(AS_DAS_groups = str_remove(feature, str_c(AS_type, "_", region, "_gr_"))) %>%
  mutate(AS_region = str_c(AS_type, "_", region)) %>%
  mutate(pair_mark = str_c(mark_1, str_sub(mark_2, 3, nchar(mark_2)))) %>%
  mutate(AS_DAS_groups = if_else(AS_DAS_groups=="PSI_all", "AS_HS0_stable_PSI",
                                 if_else(AS_DAS_groups=="merged_H16_aft_H0", "DAS_changed_PSI_aft_HS0", "ctrl"))) %>%
  mutate(AS_DAS_groups = factor(AS_DAS_groups, levels = c("AS_HS0_stable_PSI", "ctrl", "DAS_changed_PSI_aft_HS0"))) %>%
  dplyr::select(pair_mark, AS_region, coenrichment_index = mean_coenrichment_index, sd_coenrichment_index, AS_DAS_groups)

#no error bar (all events in one sample)
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

#produce function for creating error bars
produce_bt_coenrichment_index_hist_plots <- function(a, b){
  coenrichment_index_result_for_plot_f %>% 
    filter(AS_region == a[1] & pair_mark ==b[1]) %>%
    ggplot(aes(AS_DAS_groups, coenrichment_index, fill=AS_DAS_groups)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=coenrichment_index-sd_coenrichment_index, ymax=coenrichment_index+sd_coenrichment_index), width=.2) +
    scale_fill_manual(
      values = c(
        "#29f600",
        "#F6BE00",
        "#fb6a4a"
      )
    ) +
    facet_grid(AS_region ~ pair_mark) +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 20), legend.text = element_text(size = 9, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 20, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
          axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle = 60, vjust = 0.5))
}

coenrichment_index_table_for_bt_plots_key <- 
  coenrichment_index_result_for_plot_f %>%
  ungroup() %>%
  dplyr::select(AS_region, pair_mark) %>%
  distinct()

bt_error_bars_coenrichment_index_hist_plots_list <- 
  map2(coenrichment_index_table_for_bt_plots_key$AS_region, 
       coenrichment_index_table_for_bt_plots_key$pair_mark,
       produce_bt_coenrichment_index_hist_plots)

produce_bt_coenrichment_index_hist_plots(coenrichment_index_table_for_bt_plots_key$AS_region[[1]],
                                         coenrichment_index_table_for_bt_plots_key$pair_mark[[1]])

path_bt_coenrichment_index_plots_hist_AS_combined_DAS <- 
  str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_analysis_heat_comb_epimark_AS_combined_DAS_histograms/", coenrichment_index_table_for_bt_plots_key$AS_region, "_", coenrichment_index_table_for_bt_plots_key$pair_mark, "_coenrichment_index_AS_combined_DAS.jpeg")

pwalk(list(path_bt_coenrichment_index_plots_hist_AS_combined_DAS, bt_error_bars_coenrichment_index_hist_plots_list), ggsave,
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
