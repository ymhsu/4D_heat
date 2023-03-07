Sys.getlocale()
Sys.setenv(LANG = "en_US.UTF-8")
Packages <- c("plyranges", "tidyverse", "doParallel", "foreach", "nullranges", "caTools", "ggpubr", "fs")
lapply(Packages, library, character.only = TRUE)

#The initial idea for performing bootstrap analysis is based on block bootstrap from the nature communication paper
#However, we found that this way may lead to many comparisons between introns/exons with AS with irrelevant regions
#So in the end, we developed another way for comparisons based on introns/exons with and without AS

##This part is based on block bootstrap, so here I just keep these codes for a record (until line 235)
#produce all grange files from bed files for performing bootstraping test
produce_grange_data <- function(a){
  raw_data <- read_delim(a[1], delim = "\t", col_names = c("seqnames", "start", "end", "strand", "range_id")) %>%
    select(-range_id)
  
  raw_data <- raw_data %>%
    distinct()
  
  data_gr <- makeGRangesFromDataFrame(raw_data, starts.in.df.are.0based = TRUE)
  
  #specify chromosome lengths of tomato genome
  seqlengths(data_gr) <- c(92851076, 55085968, 66835077, 66250971, 68507990, 52260739, 68464519, 66735172, 68734460, 66285255, 59556480, 70392938)
  
  data_gr %>%
    sort() %>%
    mutate(id = seq_along(.))
}

path_AS_features_files <- list.files("./data/AS_bed_for_coenrichment_HS0/", pattern = "*.bed",
                               full.names = TRUE)


#produce grange data as a list 
grange_list <- lapply(path_AS_features_files, produce_grange_data)


#categorize all data into 4 lists for following analysis depending on their AS type and location
feature <- c("H3K36me3", "H3K4me1", "H3K4me3", "H3K9ac", "H4K12ac", "H4K20ac", "Pol2")
AS_type <- c("RI", "SE")
location <- c("upstream", "downstream")

coenrichment_AS_data <- vector("list", length = length(AS_type))

for (i in seq_along(AS_type)) {
  for (j in seq_along(location)) {
    
    coenrichment_AS_data[[i]][[j]] <- grange_list[which(str_detect(path_AS_features_files, AS_type[[i]]) & str_detect(path_AS_features_files, location[[j]]))]

  }
}


registerDoParallel(cores = 7)
getDoParWorkers()

bootstrap_list_multithreads <-
  foreach(i=c(1:2), .packages = c("nullranges")) %:%
  foreach(j=c(1:2), .packages = c("nullranges")) %:%
  foreach(k=1:7, .packages = c("nullranges")) %dopar% {
    bootRanges(coenrichment_AS_data[[i]][[j]][[k]], blockLength=500000, R=10000)
  }

data_temp <- bootstrap_list_multithreads[[1]][[1]][[7]] %>%
  as_tibble() 

AS_epi_feature_bootstrap_data <- list()

for (i in seq_along(location)) {
  for (j in seq_along(AS_type)) {
    for (k in seq_along(feature)) {
      
      data_temp <- bootstrap_list_multithreads[[i]][[j]][[k]] %>%
        as_tibble()
      
      write_delim(data_temp, str_c("./data/bootstrap_data/bootstrap_10000_temp_", location[[i]], "_", AS_type[[j]], "_", feature[[k]], ".bed"), delim = "\t", col_names = TRUE) 
    }
  }
}

#produce the table for creating p-value for each pair of epigenetic feature
produce_p_value_table_boot_obs <- 
function(a, b, c){
location_index <- which(str_detect(location, a[1]))
AS_type_index <- which(str_detect(AS_type, b[1]))
feature_index <- seq_along(feature)[-c[1]]

output_table <- tibble()

for (i in feature_index) {
  obs_count_table <- coenrichment_AS_data[[location_index]][[AS_type_index]][[c[1]]] %>% join_overlap_inner(coenrichment_AS_data[[location_index]][[AS_type_index]][[i]]) %>%
    summarize(n_overlaps = n())
  
  null <- coenrichment_AS_data[[location_index]][[AS_type_index]][[c[1]]] %>% join_overlap_inner(bootstrap_list_multithreads[[location_index]][[AS_type_index]][[i]]) %>%
    group_by(iter) %>%
    summarize(n_overlaps = n()) %>%
    as_tibble() %>%
    complete(iter=as.factor(c(1:10000)), fill=list(n_overlaps = 0)) %>%
    mutate(count_boot_larger_obs = if_else(n_overlaps > obs_count_table$n_overlaps, 1, 0)) 
  
  mean_bootstrap_overlaps <- mean(null$n_overlaps)
  sum_bootstrap_larger_obs <- sum(test_table$count_boot_larger_obs)
  
  temp_output_table <- tibble(
    location = location[[location_index]],
    AS_type = AS_type[[AS_type_index]],
    query = feature[[c[1]]],
    subject = feature[[i]],
    obs = obs_count_table$n_overlaps,
    mean_bootstrap = mean_bootstrap_overlaps,
    n_bt_lr_obs = sum_bootstrap_larger_obs
  )
  
  output_table <- bind_rows(output_table, temp_output_table)
}

output_table

}

seq_along(feature)[-2]

#produce the final table
output_final <- tibble()

for (i in seq_along(location)) {
  for (j in seq_along(AS_type)) {
    for (k in seq_along(feature)) {
      
      table_temp <- produce_p_value_table_boot_obs(location[[i]], AS_type[[j]], k)
      
      output_final <- bind_rows(output_final, table_temp)
    }
  }
}

write_delim(output_final, "./analysis/bootstrap_result/bootstrap_result_10000_times", delim = "\t", col_names = TRUE)

#perform the same analysis (block bootstrap) in the euchromatic regions of chr1#

grange_list_chr1 <- vector("list", length = length(grange_list))

for (i in seq_along(grange_list_chr1)) {
  grange_list_chr1[[i]] <- grange_list[[i]] %>%
    filter(seqnames == "chr1") %>%
    as_tibble() %>%
    filter(end <= 4000000 | start >= 70000000) %>%
    as_granges() %>%
    sort() %>%
    mutate(id = seq_along(.))
  
  seqlevels(grange_list_chr1[[i]], pruning.mode="coarse") <- c("chr1")
  seqlengths(grange_list_chr1[[i]]) <- c(92851076)
}

grange_list_chr1[[1]]


coenrichment_AS_data

coenrichment_AS_data_chr1 <- vector("list", length = length(AS_type))

for (i in seq_along(AS_type)) {
  for (j in seq_along(location)) {
    
    coenrichment_AS_data_chr1[[i]][[j]] <- grange_list_chr1[which(str_detect(path_AS_features_files, AS_type[[i]]) & str_detect(path_AS_features_files, location[[j]]))]
  }
}

registerDoParallel(cores = 7)
getDoParWorkers()

bootstrap_list_chr1_multithreads <-
  foreach(i=c(1:2), .packages = c("nullranges")) %:%
  foreach(j=c(1:2), .packages = c("nullranges")) %:%
  foreach(k=1:7, .packages = c("nullranges")) %dopar% {
    bootRanges(coenrichment_AS_data_chr1[[i]][[j]][[k]], blockLength=500000, R=100)
  }

#produce the table for creating p-value for each pair of epigenetic feature
produce_p_value_table_boot_obs_chr1 <- 
  function(a, b, c){
    location_index <- which(str_detect(location, a[1]))
    AS_type_index <- which(str_detect(AS_type, b[1]))
    feature_index <- seq_along(feature)[-c[1]]
    
    output_table <- tibble()
    
    for (i in feature_index) {
      obs_count_table <- coenrichment_AS_data_chr1[[location_index]][[AS_type_index]][[c[1]]] %>% join_overlap_inner(coenrichment_AS_data_chr1[[location_index]][[AS_type_index]][[i]]) %>%
        summarize(n_overlaps = n())
      
      null <- coenrichment_AS_data_chr1[[location_index]][[AS_type_index]][[c[1]]] %>% join_overlap_inner(bootstrap_list_chr1_multithreads[[location_index]][[AS_type_index]][[i]]) %>%
        group_by(iter) %>%
        summarize(n_overlaps = n()) %>%
        as_tibble() %>%
        complete(iter=as.factor(c(1:10000)), fill=list(n_overlaps = 0)) %>%
        mutate(count_boot_larger_obs = if_else(n_overlaps > obs_count_table$n_overlaps, 1, 0)) 
      
      mean_bootstrap_overlaps <- mean(null$n_overlaps)
      sum_bootstrap_larger_obs <- sum(null$count_boot_larger_obs)
      
      temp_output_table <- tibble(
        location = location[[location_index]],
        AS_type = AS_type[[AS_type_index]],
        query = feature[[c[1]]],
        subject = feature[[i]],
        obs = obs_count_table$n_overlaps,
        mean_bootstrap = mean_bootstrap_overlaps,
        n_bt_lr_obs = sum_bootstrap_larger_obs
      )
      
      output_table <- bind_rows(output_table, temp_output_table)
    }
    
    output_table
    
  }

seq_along(feature)[-2]

#produce the final table
output_final_chr1 <- tibble()

for (i in seq_along(location)) {
  for (j in seq_along(AS_type)) {
    for (k in seq_along(feature)) {
      
      table_temp <- produce_p_value_table_boot_obs_chr1(location[[i]], AS_type[[j]], k)
      
      output_final_chr1 <- bind_rows(output_final_chr1, table_temp)
    }
  }
}

write_delim(output_final_chr1, "./analysis/bootstrap_result/bootstrap_result_100_times_chr1_euchromatic_regions", delim = "\t", col_names = TRUE)



##Here I started to use a different way to perform bootstarp analysis
#produce a raw data set for bootstrap with following criteria
#genes with at least three introns or exons are considered in the analysis
#the first stage is to define all AS events in from the output of rMATs analysis
#import all AS segments (four AS types), and remove any of feature overlapping AS events
#take the code from "AS_HS0_PSI.R"
#The modification of the header for all output files of rMATs
col_names_uni <- c("ID_1", "GeneID", "geneSymbol", "chr", "strand", "pos_1", "pos_2", "pos_3", "pos_4", "pos_5", 
                   "pos_6", "ID_2", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncFormLen",
                   "SkipFormLen", "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference")

#create the function of producing HS0 AS events with PSI
produce_AS_HS0_PSI <- function(a){
  #import raw data of all .MATS.JC.txt files
  #variable "a" is different types of AS
  data_raw <- read_delim(str_c("./data/AS_events_HS0/", a[1], ".MATS.JC.txt"), col_names = col_names_uni, skip = 1) %>%
    select(chr, strand, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, IncLevel1) %>%
    filter(str_detect(IncLevel1, "NA")==FALSE) %>%
    mutate(AS_type = a[1])
  
  #split IncLevel1 (there are two replicates which will be used for calculating average PSI)
  IncLevel1_m <- str_split(data_raw$IncLevel1, ",", simplify = TRUE)
  
  #add PSI data from two replicates to compute "PSI" which is the average PSI of r1 and r2 
  data_raw %>%
    mutate(
      IncLevel_r1 = as.double(IncLevel1_m[,1]),
      IncLevel_r2 = as.double(IncLevel1_m[,2])
    ) %>%
    mutate(PSI = (IncLevel_r1 + IncLevel_r2)/2) %>%
    select(-IncLevel_r1, -IncLevel_r2)
}

#assign the variable for 4 types of AS
AS_type <- c("RI", "SE", "A5SS", "A3SS")

#create AS data with the average of PSI
AS_PSI_raw <- 
  AS_type %>%
  map(., produce_AS_HS0_PSI) %>%
  bind_rows() %>%
  mutate(AS_type = if_else(AS_type == "A5SS" | AS_type == "A3SS", str_remove(AS_type, "S"), AS_type))

#produce the list of 4 types of AS with 7 classes of PSI (nested list with two layers)
#In brief, PSI_5: PSI <= 0.05, PSI_5_20:  0.05 < PSI <= 0.2, etc.
#PSI with some given values will be removed from the list of AS
#SE with PSI as 1, RI with PSI as 0, A5S/A3S with PSI as 1 or 0, these cases are recognized as constitutive splicing 
AS_HSO_PSI_all <- 
  AS_PSI_raw %>%
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
  mutate(PI = if_else(PSI <= 0.05, "PSI_5",
                      if_else(PSI > 0.05 & PSI <= 0.2, "PSI_5_20",
                              if_else(PSI > 0.2 & PSI <= 0.4, "PSI_20_40", 
                                      if_else(PSI > 0.4 & PSI <= 0.6, "PSI_40_60",
                                              if_else(PSI > 0.6 & PSI <= 0.8, "PSI_60_80", 
                                                      if_else(PSI > 0.8 & PSI <= 0.95, "PSI_80_95", "PSI_95"))))))) %>%
  mutate(trt = "HS0") %>%
  select(chr, str, end, strand, trt, AS_type, PSI, PI) %>%
  mutate(selection_index = if_else(AS_type == "RI" & PSI == 0, 1,
                                   if_else(AS_type == "SE" & PSI == 1, 1,
                                           if_else(AS_type == "A5S" & PSI == 0, 1, 
                                                   if_else(AS_type == "A5S" & PSI == 1, 1,
                                                           if_else(AS_type == "A3S" & PSI == 0, 1,
                                                                   if_else(AS_type == "A3S" & PSI == 1, 1, 0))))))) %>%
  filter(selection_index == 0 & chr != "chrctg00000403") %>%
  select(seqnames = chr, start = str, end = end, strand, AS_type, PSI, PI) %>%
  as_granges()




#produce the function that can not only create control set universal for exon and intron from the genes without any AS events
#but RI/SE occurred features and their flanking regions depending on the need
produce_list_of_features_in_gr <- 
  function(a, b, c, data){
    #only keep introns or exons at the third one or following ones
    feature_raw <- 
      read_delim(str_c("./data/M82_annotation_data/M82_rMATs_anno_all_", a[1], ".bed"), col_names = c("seqnames", "start", "end", "strand", "feature", "source", "anno")) %>%
      group_by(anno) %>%
      mutate(count_feature = n()) %>%
      filter(count_feature >= 3) %>%
      mutate(order = if_else(strand == "+", c(1:n()), c(n():1))) %>%
      filter(order != 1 & order != 2) %>%
      ungroup() %>%
      select(seqnames, start, end, feature, strand) %>%
      distinct() %>%
      select(seqnames, start, end, strand) %>%
      as_granges()
    
    #extract the control features
    feature_ctrl_gr <- 
    feature_raw %>% 
      mutate(n_overlaps = count_overlaps(., data)) %>%
      filter(n_overlaps == 0) %>%
      mutate(type = "ctrl") %>%
      select(-n_overlaps)
    
    #extract the table with the interested AS type
    
    PSI_filter <- c[1]
    
    if (str_detect(PSI_filter, "PSI")==TRUE) {
      AS_event_gr <- 
        AS_HSO_PSI_all %>%
        filter(AS_type == b[1] & PI == c[1])
    } else {
      AS_event_gr <- 
        AS_HSO_PSI_all %>%
        filter(AS_type == b[1])
    }
    
    #extract the target features
    feature_AS_gr <- 
      feature_raw %>% 
      mutate(n_overlaps = count_overlaps(., AS_event_gr)) %>%
      filter(n_overlaps != 0) %>%
      mutate(type = b[1]) %>%
      select(-n_overlaps)
    
    #utilize the control/target feature to extract two sides of 100-bp flanking regions
    #upstream
    upstream_feature_gr <- 
      bind_rows(as_tibble(feature_ctrl_gr), as_tibble(feature_AS_gr)) %>%
      mutate(start = as.double(start), end = as.double(end)) %>%
      mutate(start_new = if_else(strand == "+", start - 100, end),
             end_new = if_else(strand == "+", start, end + 100)) %>%
      split(.$type) %>%
      map(. %>% select(seqnames, start = start_new, end = end_new, strand, type)) %>%
      map(. %>% as_granges())
      
    
    #downstream
    downstream_feature_gr <- 
      bind_rows(as_tibble(feature_ctrl_gr), as_tibble(feature_AS_gr)) %>%
      mutate(start = as.double(start), end = as.double(end)) %>%
      mutate(start_new = if_else(strand == "+", end, start - 100),
             end_new = if_else(strand == "+", end + 100, start)) %>%
      split(.$type) %>%
      map(. %>% select(seqnames, start = start_new, end = end_new, strand, type)) %>%
      map(. %>% as_granges())
      
    #make lists containing granges from two flanking regions at the junction of features with or without AS
    #and features as well
    
    #the the for AS and ctrl
    if (str_detect(unique(as_tibble(upstream_feature_gr[[1]])$type), "ctrl")==TRUE) {
      feature_ctrl_gr_l <- list(upstream_feature_gr[[1]], feature_ctrl_gr, downstream_feature_gr[[1]])
      feature_target_gr_l <- list(upstream_feature_gr[[2]], feature_AS_gr, downstream_feature_gr[[2]])
    } else {
      feature_ctrl_gr_l <- list(upstream_feature_gr[[2]], feature_ctrl_gr, downstream_feature_gr[[2]])
      feature_target_gr_l <- list(upstream_feature_gr[[1]], feature_AS_gr, downstream_feature_gr[[1]])
    }
    
    names(feature_target_gr_l) <- c("upstream_feature_gr", "feature_gr", "downstream_feature_gr")
    
    #name both lists
    names(feature_target_gr_l) <- c("upstream_feature_gr", "feature_gr", "downstream_feature_gr")
    
    names(feature_ctrl_gr_l) <- c("upstream_feature_gr", "feature_gr", "downstream_feature_gr")
    
    #make the final list including AS and ctrl
    feature_gr_final_l <- list(feature_target_gr_l, feature_ctrl_gr_l)
    
    names(feature_gr_final_l) <- c(str_c(b[1], "_feature"), str_c(b[1], "_ctrl"))
    
    feature_gr_final_l
  }

#make the list of grange files based on RI or SE (all PSI)
gr_RI_all <- produce_list_of_features_in_gr("intron", "RI", "all", AS_HSO_PSI_all)

gr_SE_all <- produce_list_of_features_in_gr("exon", "SE", "all", AS_HSO_PSI_all)

#make the list of grange files
PSI_type <- str_c("PSI_", c("5", "5_20", "20_40", "40_60","60_80", "80_95", "95"))

gr_RI_separate_PSI <- 
  pmap(list(rep("intron", length(PSI_type)), rep("RI", length(PSI_type)), PSI_type, replicate(length(PSI_type), as_granges(AS_HSO_PSI_all))), produce_list_of_features_in_gr)

gr_SE_separate_PSI <- 
  pmap(list(rep("exon", length(PSI_type)), rep("SE", length(PSI_type)), PSI_type, replicate(length(PSI_type), as_granges(AS_HSO_PSI_all))), produce_list_of_features_in_gr)

gr_RI_separate_PSI[[1]]

produce_coenrichment_index(gr_RI_separate_PSI[[7]], 2, 1, c(comb_pair_epimark[1,]))

#import the peak of 7 epigenetic features
epi_feature <- c("H3K36me3", "H3K4me1", "H3K4me3", "H3K9ac", "H4K12ac", "H4K20ac", "Pol2")

path_epi_feature <- str_c("./data/peak_epimarks/no_trt/peak_pileup_norm_file_simplied_name/", epi_feature, "_pileup_norm.bed")

epi_feature_gr <- 
lapply(path_epi_feature, read_delim, delim = "\t", col_names = c("seqnames", "start", "end", "signal")) %>%
  map(. %>% as_granges())

names(epi_feature_gr) <- epi_feature

#create the function for the "co-enrichment index"
#the idea of this index is to use the intersection between peaks overlapping features, 
#and divide it by the union

gr_SE_all$SE_feature
dplyr::n()
gr_RI_all$RI_feature$feature_gr %>% 
  join_overlap_inner(epi_feature_gr$H3K36me3) %>%
  summarize(n_overlaps = n()) 

produce_coenrichment_index <-
  function(data, a, b, c){
    
    data_target <- data[[a]][[b]]
    
    produce_epimark_intersection_f <-  
      function(b, data){
        
        data %>%
          mutate(label_gr = c(1:plyranges::n())) %>%
          mutate(n_overlaps = count_overlaps(., epi_feature_gr[[b]])) %>%
          filter(n_overlaps != 0)
      }
    
    feature_epimark_list <- lapply(c(c), produce_epimark_intersection_f, data_target)

    feature_epimark1_2_intersection <-
      feature_epimark_list[[1]] %>%
      mutate(n_overlaps = count_overlaps(., feature_epimark_list[[2]])) %>%
      filter(n_overlaps != 0) %>%
      as_tibble()
    
    feature_epimark2_1_intersection <-
      feature_epimark_list[[2]] %>%
      mutate(n_overlaps = count_overlaps(., feature_epimark_list[[1]])) %>%
      filter(n_overlaps != 0) %>%
      as_tibble() 


    if (nrow(feature_epimark1_2_intersection) > nrow(feature_epimark2_1_intersection)) {
      coenrichment_index <-
        nrow(feature_epimark1_2_intersection) / c(nrow(as_tibble(feature_epimark_list[[2]])) - nrow(feature_epimark2_1_intersection) + nrow(as_tibble(feature_epimark_list[[1]])))
    } else {
      coenrichment_index <-
        nrow(feature_epimark2_1_intersection) / c(nrow(as_tibble(feature_epimark_list[[2]])) - nrow(feature_epimark1_2_intersection) + nrow(as_tibble(feature_epimark_list[[1]])))
    }
    
    tibble(
      target = names(data)[[a]],
      location = names(data[[a]])[[b]],
      epimark_1 = names(epi_feature_gr)[c[1]],
      epimark_2 = names(epi_feature_gr)[c[2]],
      coenrichment_index = coenrichment_index
    )
    
  }

feature
names(gr_RI_all)[[1]]

system.time(
nrow(as_tibble(gr_RI_all[[1]]$upstream_feature_gr))
)

nrow(gr_RI_all[[1]]$upstream_feature_gr)

system.time(
length(gr_RI_all[[1]]$upstream_feature_gr$type)
)

system.time(
produce_coenrichment_index_v2(gr_RI_all, 1, 2, c(2, 7))
)

system.time(
pmap(list(replicate(3, gr_RI_all, simplify = FALSE), rep(1, 3), rep(2, 3), replicate(3, c(2, 7), simplify = FALSE)), produce_coenrichment_index)
)

comb_pair_epimark

#produce coenrichment index tables for SE and RI
#all PSI combined
comb_pair_epimark <- combs(seq_along(epi_feature_gr), 2) #seven epigenetic features

coenrichment_all_comb_RI <- tibble()
coenrichment_all_comb_SE <- tibble()


for (i in seq_along(gr_RI_all)) {
  for (j in seq_along(gr_RI_all[[1]])) {
    for (k in seq_along(comb_pair_epimark[,1])) {
      
      single_result_RI <- produce_coenrichment_index(gr_RI_all, i, j, c(comb_pair_epimark[k,]))
      
      single_result_SE <- produce_coenrichment_index(gr_SE_all, i, j, c(comb_pair_epimark[k,]))
      
      
      coenrichment_all_comb_RI <- bind_rows(coenrichment_all_comb_RI, single_result_RI)
      
      coenrichment_all_comb_SE <- bind_rows(coenrichment_all_comb_SE, single_result_SE)
    }
  }
}




#compare roughly the average of coenrichment index between feature with AS and without AS
coenrichment_all_comb_RI %>%
  mutate(x_label = rep(c(1:21), 6)) %>%
  ggplot() +
  geom_point(aes(x_label, coenrichment_index, color = target)) +
  facet_wrap(~location)

coenrichment_all_comb_SE %>%
  mutate(x_label = rep(c(1:21), 6)) %>%
  ggplot() +
  geom_point(aes(x_label, coenrichment_index, color = target)) +
  facet_wrap(~location)

coenrichment_all_comb_AS_all_PSI <- 
bind_rows(coenrichment_all_comb_RI, coenrichment_all_comb_SE)

coenrichment_all_comb_AS_all_PSI %>%
  View()

write_delim(coenrichment_all_comb_AS_all_PSI, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_all_comb_AS_all_PSI.txt", delim = "\t", col_names = TRUE)


#PSI separated
coenrichment_all_comb_RI_PSI_separated <- tibble()
coenrichment_all_comb_SE_PSI_separated <- tibble()


for (i in seq_along(PSI_type)) {
  for (j in seq_along(gr_RI_separate_PSI[[1]])) {
    for (k in seq_along(gr_RI_separate_PSI[[1]][[1]])) {
      for (l in seq_along(comb_pair_epimark[,1])) {
        print(str_c(PSI_type[[i]], "_", str_remove(names(gr_RI_separate_PSI[[1]]), "RI_")[[j]], "_", location_gr[[k]], "_", epi_feature[comb_pair_epimark[l,][[1]]], "_", epi_feature[comb_pair_epimark[l,][[2]]]))
        
        single_result_RI_PSI_separated <-
          produce_coenrichment_index(gr_RI_separate_PSI[[i]], j, k, c(comb_pair_epimark[l, ])) %>%
          mutate(PSI = PSI_type[[i]])
        
        single_result_SE_PSI_separated <-
          produce_coenrichment_index(gr_SE_separate_PSI[[i]], j, k, c(comb_pair_epimark[l, ]))%>%
          mutate(PSI = PSI_type[[i]])
        
        coenrichment_all_comb_RI_PSI_separated <-
          bind_rows(coenrichment_all_comb_RI_PSI_separated, single_result_RI_PSI_separated)
        
        coenrichment_all_comb_SE_PSI_separated <-
          bind_rows(coenrichment_all_comb_SE_PSI_separated, single_result_SE_PSI_separated)
      }
    }
  }
}


#compare df PSI with ctrl and PSI combined
ctrl_coenrichment_index <- 
coenrichment_all_comb_AS_all_PSI %>%
  filter(str_detect(target, "ctrl")) %>%
  mutate(target = str_remove(target, "_ctrl")) %>%
  select(target, location, epimark_1, epimark_2, ctrl = coenrichment_index)

coenrichment_index_all_separate_PSI_ctrl <- 
coenrichment_all_comb_AS_all_PSI %>% 
  mutate(PSI = "all") %>%
  bind_rows(coenrichment_all_comb_RI_PSI_separated, coenrichment_all_comb_SE_PSI_separated) %>%
  filter(str_detect(target, "ctrl")!=TRUE) %>%
  spread(PSI, coenrichment_index) %>%
  mutate(target = str_remove(target, "_feature")) %>%
  left_join(ctrl_coenrichment_index) %>%
  select(target, location, epimark_1, epimark_2, all, ctrl, PSI_5, PSI_5_20, PSI_20_40, PSI_40_60, PSI_60_80, PSI_80_95, PSI_95) 

write_delim(coenrichment_index_all_separate_PSI_ctrl, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_all_separate_PSI_ctrl.txt.", delim = "\t", col_names = TRUE)

tibble(
  target = test$target,
  epimark_1 = test$epimark_1,
  epimark_2 = test$epimark_2,
  k = test$PSI,
  v = test$coenrichment_index
  ) %>%
  arrange(epimark_1, epimark_2) %>%
  View()


#create the function for bootstraping
gr_bootstrap_f <- 
  function(data){
    #transfer granges data into tibble
    all_data <- 
      data %>%
      mutate(label = c(1:plyranges::n())) %>%
      as_tibble()
    
    #randomly sample 10% of data
    all_data_granges <-
      all_data[sample(c(1:nrow(all_data)), round(nrow(all_data)/10)),] %>%
      as_granges() %>%
      sort()
    
    all_data_granges %>%
      as_tibble()
  }

##here I am going to create the loop for creating bootstrap samples for combined PSI, separate PSI, and control for RI and SE
#create the function for extracting features with AS from separate PSI lists

extract_separate_PSI_list <-
  function(data){
  data_label <- which(str_detect(names(data), "ctrl", negate = TRUE))
  
  #trick to remove a level of list
  #https://stackoverflow.com/questions/17164715/how-to-remove-a-level-of-lists-from-a-list-of-lists
  data_list <- unlist(data[data_label], recursive = FALSE)
  
  names(data_list) <- str_remove(names(data_list), ".*feature\\.")
  
  data_list
}

#merge all PSI, separate PSI and ctrl into one list with the order (PSI_type, PSI_combined, ctrl)
RI_list_for_bootstrap <- map(gr_RI_separate_PSI, extract_separate_PSI_list)

RI_list_for_bootstrap <- 
  append(append(RI_list_for_bootstrap, list(gr_RI_all$RI_feature)), list(gr_RI_all$RI_ctrl))

SE_list_for_bootstrap <- map(gr_SE_separate_PSI, extract_separate_PSI_list)

SE_list_for_bootstrap <- 
  append(append(SE_list_for_bootstrap, list(gr_SE_all$SE_feature)), list(gr_SE_all$SE_ctrl))

AS_bootstrap_rawdata_list <- list(RI_list_for_bootstrap, SE_list_for_bootstrap)

#create bootstraped data
AS_type_bootstrap <- c("RI", "SE")
PSI_type_bootstrap <- c(PSI_type, "PSI_all", "ctrl")

set.seed(200)

#the below is the choice for using different types of cores
#single core
for (i in seq_along(AS_type_bootstrap)) {
  for (j in seq_along(PSI_type_bootstrap)) {
    for (k in c(1:10)) {
     
      bt_sample <- map(AS_bootstrap_rawdata_list[[i]][[j]], gr_bootstrap_f)
      
      bt_sample_path <- str_c("./data/AS_HS0_bootstrap_samples/", AS_type_bootstrap[[i]], "_", PSI_type_bootstrap[[j]], "_", names(test), "_bootstrap_", k)
      
      pwalk(list(bt_sample, bt_sample_path), write_delim, delim = "\t", col_names = TRUE)
       
    }
  }
}

#multiple cores
set.seed(200)

registerDoParallel(cores = 4)
getDoParWorkers()

foreach(i=seq_along(AS_type_bootstrap), .packages = c("plyranges", "tidyverse")) %:%
  foreach(j=seq_along(PSI_type_bootstrap), .packages = c("plyranges", "tidyverse")) %:%
  foreach(k=1:1000, .packages = c("plyranges", "tidyverse")) %dopar% {
    bt_sample <- map(AS_bootstrap_rawdata_list[[i]][[j]], gr_bootstrap_f)
    
    bt_sample_path <- str_c("./data/AS_HS0_bootstrap_samples/", AS_type_bootstrap[[i]], "_", PSI_type_bootstrap[[j]], "_", names(bt_sample), "_bootstrap_", k)
    
    pwalk(list(bt_sample, bt_sample_path), write_delim, delim = "\t", col_names = TRUE)
    
  }

#read all bootstrap files for producing histograms with error bars and pvalue (comparing features with and without AS) 
AS_HS0_bootstrap_bootstrap_list <- fs::dir_ls("./data/AS_HS0_bootstrap_samples/", glob = "*_bootstrap_*")

location <- c("upstream_feature", "feature", "downstream_feature")



registerDoParallel(cores = 10)
getDoParWorkers()

foreach(i = seq_along(AS_type_bootstrap),
        .packages = c("plyranges", "tidyverse")) %:%
  foreach(j = seq_along(PSI_type_bootstrap),
          .packages = c("plyranges", "tidyverse")) %:%
  foreach(k = seq_along(location),
          .packages = c("plyranges", "tidyverse")) %:%
  foreach(l = seq_along(comb_pair_epimark[, 1]),
          .packages = c("plyranges", "tidyverse")) %:%
  foreach(m = 1:1000,
          .packages = c("plyranges", "tidyverse")) %dopar% {
            path_bootstrap <-
              str_c(
                "./data/AS_HS0_bootstrap_samples/",
                AS_type_bootstrap[[i]],
                "_",
                PSI_type_bootstrap[[j]],
                "_",
                location,
                "_gr_bootstrap_",
                m
              )
            
            #read bootstrap sample
            bootstrap_sample <-
              map(path_bootstrap,
                  read_delim,
                  delim = "\t",
                  col_names = TRUE)
            
            #transmform tibble into granges
            bootstrap_sample <-
              bootstrap_sample %>%
              map(. %>% as_granges)
            
            #add another level of list for running the function of producing coenrichment index
            bootstrap_sample <-
              list(bootstrap_sample)
            
            #add the name of the first level
            names(bootstrap_sample) <-
              str_c(AS_type_bootstrap[[i]], "_", PSI_type_bootstrap[[j]])
            
            #add the name of the second level
            names(bootstrap_sample[[1]]) <- location
            
            #produce the output file of coenrichment index
            output <-
              produce_coenrichment_index(bootstrap_sample, 1, k, c(comb_pair_epimark[l,])) %>%
              mutate(bootstrap_label = m)
            
            path_output <- 
              str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/bootstrap_coenrichment_result/", AS_type_bootstrap[[i]], "_", PSI_type_bootstrap[[j]], "_", location[[k]], "_epimark_", c(comb_pair_epimark[l,])[1], "_", c(comb_pair_epimark[l,])[2], "_bt_", m, "_coenrichment_index.txt")
            
            write_delim(output, path_output, delim = "\t", col_names = TRUE)
            
          }


registerDoParallel(cores = 5)
getDoParWorkers()

seq_along(AS_type_bootstrap)
seq_along(PSI_type_bootstrap)
seq_along(location)
coenrichment_index_1000_bootstrap_list <- 
foreach(i = 1,
        .packages = c("plyranges", "tidyverse")) %:%
  foreach(j = 1,
          .packages = c("plyranges", "tidyverse")) %:%
  foreach(k = 1,
          .packages = c("plyranges", "tidyverse")) %:%
  foreach(l = seq_along(comb_pair_epimark[, 1]),
          .packages = c("plyranges", "tidyverse")) %:%
  foreach(m = 1:5,
          .packages = c("plyranges", "tidyverse")) %dopar% {
            
            path_output <- 
              str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/bootstrap_coenrichment_result/", AS_type_bootstrap[[i]], "_", PSI_type_bootstrap[[j]], "_", location[[k]], "_epimark_", c(comb_pair_epimark[l,])[1], "_", c(comb_pair_epimark[l,])[2], "_bt_", m, "_coenrichment_index.txt")
            
            read_delim(path_output, delim = "\t", col_names = TRUE)
          }


#import the table with all coenrichment index including 1000 bootstraping
coenrichment_index_1000_bootstrap_total <-
  read_delim("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_1000_bootstrap_total.txt", delim = "\t", col_names = c("target", "location", "epimark_1", "epimark_2", "coenrichment_index", "bootstrap_label"))

coenrichment_index_1000_bootstrap_total %>%
  select(target) %>%
  distinct()

#import the table with coenrichment index of separate PSI
coenrichment_index_all_separate_PSI_ctrl <- 
read_delim("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_all_separate_PSI_ctrl.txt", delim = "\t", col_names = TRUE)

coenrichment_index_all_separate_PSI_ctrl_for_pvalue <- 
coenrichment_index_all_separate_PSI_ctrl %>%
  gather(key = "PSI", value = "coenrichment_index_all", c(5:13)) %>%
  mutate(target = str_c(target, "_", PSI), location = str_remove(location, "_gr")) %>%
  select(-PSI) %>%
  mutate(ctrl = str_c(str_sub(target, 1, 2), "_ctrl")) %>%
  select(target, location, epimark_1, epimark_2, ctrl, coenrichment_index_all)

coenrichment_index_1000_bootstrap_total_ctrl <- 
  coenrichment_index_1000_bootstrap_total %>%
  filter(str_detect(target, "ctrl")) %>%
  select(location, epimark_1, epimark_2, ctrl = target, coenrichment_index, bootstrap_label) 

coenrichment_index_all_separate_PSI_ctrl_for_larger_count <- 
coenrichment_index_all_separate_PSI_ctrl_for_pvalue %>%
  left_join(coenrichment_index_1000_bootstrap_total_ctrl) %>%
  filter(coenrichment_index > coenrichment_index_all) %>%
  group_by(target, location, epimark_1, epimark_2) %>%
  summarise(count = n())

coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table <- 
coenrichment_index_all_separate_PSI_ctrl_for_pvalue %>%
  left_join(coenrichment_index_all_separate_PSI_ctrl_for_larger_count) %>%
  replace_na(list(count = 0)) %>%
  mutate(pvalue = count/1000) 

write_delim(coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table.txt", delim = "\t", col_names = TRUE)

coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table <- 
read_delim("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table.txt", delim = "\t", col_names = TRUE)

coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table_simplified <- 
coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table %>%
  select(-ctrl, -coenrichment_index_all, -count) %>%
  spread(key = location, value = pvalue) %>%
  select(AS_PSI_group = target, epimark_1, epimark_2, upstream_feature, feature, downstream_feature) %>%
  filter(str_detect(AS_PSI_group, "ctrl")==FALSE)

write_delim(coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table_simplified, "./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table_simplified.txt", delim = "\t", col_names = TRUE)

coenrichment_index_all_separate_PSI_ctrl_pvalue_plot_comb_epimark_raw <- 
coenrichment_index_all_separate_PSI_ctrl_pvalue_final_table_simplified %>%
  gather(key = peak_location, value = pvalue, c(4:6)) %>%
  mutate(peak_location = factor(peak_location, levels = c("upstream_feature", "feature", "downstream_feature"))) %>%
  mutate(comb_epimark = str_c(epimark_1, "_", epimark_2)) %>%
  mutate(AS_group = str_sub(AS_PSI_group, 1, 2), PSI_group = str_remove(AS_PSI_group, str_c(AS_group, "_"))) %>%
  mutate(PSI_group = factor(PSI_group, levels = c("all", str_c(
    "PSI_",
    c("5", "5_20", "20_40", "40_60", "60_80", "80_95", "95")
  )))) 

making_plot_comb_epimark_pvalue <- 
function(a, b){
  coenrichment_index_all_separate_PSI_ctrl_pvalue_plot_comb_epimark_raw %>%
  filter(AS_group == a[1], comb_epimark == b[1]) %>%
  ggplot(aes(peak_location, pvalue, color = PSI_group)) +
  geom_point(size = 2) +
  geom_line(aes(group = PSI_group), linewidth = 1) +
  scale_color_manual(
    values = c(
      "#29f600",
      "#22a7f0",
      "#63bff0",
      "#a7d5ed",
      "#ac00f6",
      "#e1a692",
      "#de6e56",
      "#e14b31"
    )
  ) +
  facet_grid(AS_group ~ comb_epimark) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 9, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        legend.title = element_blank(), 
        axis.text.y = element_text(color = "black", size = 12, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 60, vjust = 0.5))
}

making_plot_comb_epimark_pvalue("RI", "H3K4me1_H4K20ac")

rep("RI", 5)

RI_coenrichment_rejected_comb_epimark <-
  c("H3K4me1_H3K4me3", "H3K4me1_H4K12ac", "H3K4me1_H4K20ac", "H3K4me3_H4K12ac", "H3K4me3_H4K20ac")

SE_coenrichment_rejected_comb_epimark <-
  c("H3K4me1_H3K4me3", "H3K4me1_H4K12ac", "H3K4me1_H4K20ac", "H3K9ac_Pol2", "H4K12ac_Pol2")

RI_plot_coenrichment_trend_separated_PSI <- 
  map2(rep("RI", 5), RI_coenrichment_rejected_comb_epimark, making_plot_comb_epimark_pvalue)

SE_plot_coenrichment_trend_separated_PSI <- 
  map2(rep("SE", 5), SE_coenrichment_rejected_comb_epimark, making_plot_comb_epimark_pvalue)

path_RI_plot_coenrichment_trend_separated_PSI <- str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_analysis_rejected_comb_epimark_plots/RI_plot_coenrichment_trend_separated_PSI_", c(1:5), ".jpeg")

path_SE_plot_coenrichment_trend_separated_PSI <- str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_analysis_rejected_comb_epimark_plots/SE_plot_coenrichment_trend_separated_PSI_", c(1:5), ".jpeg")


ggsave(path_RI_plot_coenrichment_trend_separated_PSI[[1]], RI_plot_coenrichment_trend_separated_PSI[[1]], width = 400, height = 240, units = c("mm"), dpi = 320)

pwalk(list(path_RI_plot_coenrichment_trend_separated_PSI, RI_plot_coenrichment_trend_separated_PSI), ggsave, width = 400, height = 240, units = c("mm"), dpi = 320)
pwalk(list(path_SE_plot_coenrichment_trend_separated_PSI, SE_plot_coenrichment_trend_separated_PSI), ggsave, width = 400, height = 240, units = c("mm"), dpi = 320)


#produce histograms of coenrichment of separated PSI with error bars
coenrichment_index_all_separate_PSI_ctrl <- 
read_delim("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_all_separate_PSI_ctrl.txt", delim = "\t", col_names = TRUE) %>%
  gather(key = PSI_group, value = coenrichment_index, c(5:13)) %>%
  mutate(PSI_group = if_else(PSI_group == "all", "PSI_all", PSI_group)) %>%
  mutate(target = str_c(target, "_", PSI_group), location = str_remove(location, "_gr"))

coenrichment_index_1000_bootstrap_total_hist_plots_raw <- 
coenrichment_index_1000_bootstrap_total %>%
  group_by(target, location, epimark_1, epimark_2) %>%
  summarise(sd_coenrichment_index = sd(coenrichment_index), bt_coenrichment_index_mean = mean(coenrichment_index)) %>%
  left_join(coenrichment_index_all_separate_PSI_ctrl) %>%
  mutate(AS_group = str_sub(target, 1, 2)) %>%
  mutate(PSI_group = factor(PSI_group, levels = c("PSI_all", "ctrl", str_c(
    "PSI_",
    c("5", "5_20", "20_40", "40_60", "60_80", "80_95", "95")
  )))) %>%
  mutate(AS_region = str_c(AS_group, "_", location), pair_mark = str_c(epimark_1, "_", epimark_2))

coenrichment_index_1000_bootstrap_total_hist_plots_raw %>%
  View()

coenrichment_index_1000_bootstrap_total_hist_plots_key <- 
coenrichment_index_1000_bootstrap_total_hist_plots_raw %>%
  ungroup() %>%
  select(AS_region, pair_mark) %>%
  distinct()

produce_coenrichment_index_hist_plots <- function(a, b){
coenrichment_index_1000_bootstrap_total_hist_plots_raw %>% 
  filter(AS_region == a[1] & pair_mark == b[1]) %>%
  ggplot(aes(PSI_group, bt_coenrichment_index_mean, fill=PSI_group)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=bt_coenrichment_index_mean-sd_coenrichment_index, ymax=bt_coenrichment_index_mean+sd_coenrichment_index), width=.2) +
  scale_fill_manual(
    values = c(
      "#29f600",
      "#F6BE00",
      "#22a7f0",
      "#63bff0",
      "#a7d5ed",
      "#ac00f6",
      "#e1a692",
      "#de6e56",
      "#e14b31"
    )
  ) +
  facet_grid(AS_region ~ pair_mark) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 20), legend.text = element_text(size = 9, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 20, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
        axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle = 60, vjust = 0.5))
}

coenrichment_index_1000_bootstrap_total_hist_plots_raw %>% 
  filter(AS_region == "SE_feature" & pair_mark == "H3K4me1_H3K4me3") %>%
  View()

coenrichment_index_1000_bootstrap_total_hist_plots_key %>%
  View()

produce_coenrichment_index_hist_plots(coenrichment_index_1000_bootstrap_total_hist_plots_key$AS_region[[91]],
                                      coenrichment_index_1000_bootstrap_total_hist_plots_key$pair_mark[[91]])

coenrichment_index_hist_plots_list <- 
map2(coenrichment_index_1000_bootstrap_total_hist_plots_key$AS_region, 
     coenrichment_index_1000_bootstrap_total_hist_plots_key$pair_mark,
     produce_coenrichment_index_hist_plots)



path_coenrichment_index_plots_hist_separated_PSI <- 
str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_analysis_comb_epimark_histgram_separated_PSI/", coenrichment_index_1000_bootstrap_total_hist_plots_key$AS_region, "_", coenrichment_index_1000_bootstrap_total_hist_plots_key$pair_mark, "_coenrichment_index_separated_PSI.jpeg")

pwalk(list(path_coenrichment_index_plots_hist_separated_PSI, coenrichment_index_hist_plots_list), ggsave,
      width = 400, height = 240, units = c("mm"), dpi = 320)

coenrichment_all_comb_RI
length(gr_ctrl_final_l)
AS_pvalue_raw_tibble <- tibble()

AS_coenrichment_all_comb <- list(coenrichment_all_comb_RI, coenrichment_all_comb_SE)

names(gr_ctrl_final_l)

AS_type <- c("RI", "SE")
comb_pair_epimark[1,][[1]]

epi_feature[comb_pair_epimark[l,][[2]]]
location_gr <- c("upstream_feature_gr", "feature_gr", "downstream_feature_gr")


for (i in seq_along(gr_ctrl_final_l)) {
  for (j in seq_along(gr_ctrl_final_l[[1]])) {
    for (k in seq_along(c("upstream_feature_gr", "feature_gr", "downstream_feature_gr"))) {
      for (l in seq_along(comb_pair_epimark[,1])) {
      
        print(str_c("The_iteration_", i, "_", AS_type[[j]], "_", location_gr[[k]], "_", epi_feature[comb_pair_epimark[l,][[1]]], "_", epi_feature[comb_pair_epimark[l,][[2]]]))
        
        AS_single_ctrl_result <-
          produce_coenrichment_index(gr_ctrl_final_l[[i]], j, k, c(comb_pair_epimark[l,]))
        
        #print(AS_single_ctrl_result)
        
        AS_feature_table <-
          AS_coenrichment_all_comb[[j]] %>%
          filter(
            location == AS_single_ctrl_result$location[[1]] &
              epimark_1 == AS_single_ctrl_result$epimark_1[[1]] &
              epimark_2 == AS_single_ctrl_result$epimark_2[[1]]
          ) %>%
          filter(str_detect(target, "feature") == TRUE)
        
        #print(AS_feature_table)
        
        AS_single_ctrl_result <-
          AS_single_ctrl_result %>%
          mutate(target_coenrichment_index = AS_feature_table$coenrichment_index,
                 iter = i) %>%
          filter(coenrichment_index > target_coenrichment_index)
        
        #print(AS_single_ctrl_result)
        
        AS_pvalue_raw_tibble <-
          bind_rows(AS_pvalue_raw_tibble, AS_single_ctrl_result)
        
        print(AS_pvalue_raw_tibble)
      
      }
    }
  }
}

write_delim(AS_pvalue_raw_tibble, "./analysis/AS_pvalue_raw_tibble_100_bootstrap", delim = "\t", col_names = TRUE)


AS_pvalue_raw_simplified <- 
AS_pvalue_raw_tibble %>%
  group_by(target, location, epimark_1, epimark_2) %>%
  summarise(count = dplyr::n()) %>%
  mutate(count = as.double(count))

p_value_result_two_epimark_RI_SE <- 
bind_rows(AS_coenrichment_all_comb) %>%
  filter(str_detect(target, "_feature")==FALSE) %>%
  left_join(AS_pvalue_raw_simplified) %>%
  replace_na(list(count = 0)) %>%
  mutate(pvalue = count/100)

p_value_result_two_epimark_RI_SE_f <- 
p_value_result_two_epimark_RI_SE %>%
  select(-coenrichment_index, -count) %>%
  spread(location, pvalue) %>%
  mutate(target = str_remove(target, "_ctrl"))
   
write_delim(p_value_result_two_epimark_RI_SE_f, "./analysis/p_value_result_two_epimark_RI_SE_100_bootstrap.txt", delim = "\t", col_names = TRUE)
              

#make plots for coenrichment index for combined samples and separate PSI
library(tidyverse)

coenrichment_index_all_separate_PSI_ctrl <- 
  read_delim("./analysis/AS_RI_SE_epimark_coenrichment_analysis/coenrichment_index_all_separate_PSI_ctrl.txt", delim = "\t", col_names = TRUE)

p_value_result_two_epimark_RI_SE_f <- 
  read_delim("./analysis/p_value_result_two_epimark_RI_SE_100_bootstrap.txt", delim = "\t", col_names = TRUE)

p_value_result_two_epimark_RI_SE_sig <- 
  p_value_result_two_epimark_RI_SE_f %>%
  filter(upstream_feature_gr < 0.05)

coenrichment_index_data_plot_raw <-
  coenrichment_index_all_separate_PSI_ctrl %>%
  gather(key = PSI_group, value = coenrichment_index, c("all", "ctrl", str_c(
    "PSI_", c("5", "5_20", "20_40", "40_60", "60_80", "80_95", "95")
  ))) %>%
  group_by(target, PSI_group, epimark_1, epimark_2) %>%
  summarise(mean_coenrich_index = mean(coenrichment_index)) %>%
  left_join(p_value_result_two_epimark_RI_SE_sig) %>%
  drop_na() %>%
  arrange(target,-mean_coenrich_index) %>%
  mutate(mark_pair = str_c(epimark_1, "/", epimark_2)) %>%
  split(.$target) 

theme_ym <- function(data){
  data +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 9, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 12, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
          axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 60, vjust = 0.5)) 
}

coenrichment_index_plot_PSI_combined <- function(data){
  data %>%
    ggplot() +
    geom_bar(aes(reorder(mark_pair, -mean_coenrich_index), mean_coenrich_index, fill=mean_coenrich_index),
             col="red", alpha = .2, stat="identity") +
    scale_fill_gradient("Percentage", low = "blue", high = "red") +
    #scale_y_continuous(labels = percent_format()) +
    labs(title="Histogram for coenrichment index") +
    labs(x="Pairs of marks", y="Coenrichment index") +
    facet_wrap(~target)
}

coenrichment_index_RI_PSI_combined <-
  subset(coenrichment_index_data_plot_raw$RI, PSI_group == "all") %>%
  coenrichment_index_plot_PSI_combined() %>%
  theme_ym() + theme(axis.text.x = element_blank(), legend.position = "bottom")

coenrichment_index_SE_PSI_combined <-
  subset(coenrichment_index_data_plot_raw$SE, PSI_group == "all") %>%
  coenrichment_index_plot_PSI_combined() %>%
  theme_ym() + theme(axis.text.x = element_blank(), legend.position = "bottom")

#make the order of pairs of epimarks for generating the same order of epimarks in plots with separate PSI groups
pair_mark_fct_RI <- 
reorder(subset(coenrichment_index_data_plot_raw$RI, PSI_group == "all")$mark_pair, -subset(coenrichment_index_data_plot_raw$RI, PSI_group == "all")$mean_coenrich_index)

pair_mark_fct_SE <- 
  reorder(subset(coenrichment_index_data_plot_raw$SE, PSI_group == "all")$mark_pair, -subset(coenrichment_index_data_plot_raw$SE, PSI_group == "all")$mean_coenrich_index)

coenrichment_index_plot_PSI_separated <- function(data, data2) {
  data %>%
    filter(PSI_group != "all") %>%
    mutate(mark_pair = factor(mark_pair, levels = data2)) %>%
    mutate(PSI_group = factor(PSI_group, levels = c("ctrl", str_c(
      "PSI_",
      c("5", "5_20", "20_40", "40_60", "60_80", "80_95", "95")
    )))) %>%
    #arrange(desc(mean_coenrich_index)) %>%
    ggplot(aes(x = mark_pair, y = mean_coenrich_index, fill = PSI_group)) +
    #ggplot(aes(reorder(mark_pair, -mean_coenrich_index), mean_coenrich_index, fill=PSI_group))+
    geom_bar(stat = "identity", position = position_identity()) +
    scale_fill_manual(
      values = c(
        "#F6BE00",
        "#22a7f0",
        "#63bff0",
        "#a7d5ed",
        "#e2e2e2",
        "#e1a692",
        "#de6e56",
        "#e14b31"
      )
    )
}

coenrichment_index_data_plot_raw$RI %>%
  group_by(mark_pair) %>%
  summarise()

coenrichment_index_RI_PSI_separated <-
  coenrichment_index_data_plot_raw$RI %>%
  coenrichment_index_plot_PSI_separated(., data2 = pair_mark_fct_RI) %>%
  theme_ym() + theme(legend.position = "bottom")

coenrichment_index_SE_PSI_separated <-
  coenrichment_index_data_plot_raw$SE %>%
  coenrichment_index_plot_PSI_separated(., data2 = pair_mark_fct_SE) %>%
  theme_ym() + theme(legend.position = "bottom")
  
library(ggpubr)
coenrichment_index_RI_SE_PSI_all_separated_plot <-
  ggarrange(
    coenrichment_index_RI_PSI_combined,
    coenrichment_index_SE_PSI_combined,
    coenrichment_index_RI_PSI_separated,
    coenrichment_index_SE_PSI_separated,
    ncol = 2,
    nrow = 2
  )

ggsave(str_c("./analysis/AS_RI_SE_epimark_coenrichment_analysis/", "coenrichment_index_RI_SE_PSI_all_separated_plot", ".jpeg"), coenrichment_index_RI_SE_PSI_all_separated_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
