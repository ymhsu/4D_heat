install.packages("pacman")
library(pacman)
Packages <- c("tidyverse", "doParallel", "foreach", "gridExtra", "ggeffects", "cowplot", "caTools", "grid", "ggplotify", "Boruta")
p_load(Packages, character.only = TRUE)


#import all genes
M82_rMATs_anno_all_gene <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno")) %>%
  mutate(gene_size = end - str)

#import all exons with their order in genes
M82_rMATs_anno_all_exon_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_exon_order.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order")) %>%
  mutate(size = end - str)

#check the quantile of gene sizes for all genes (the median gene size is about 3 kb)
quantile(M82_rMATs_anno_all_gene$gene_size, seq(0, 1, 0.025))

#quick check the distribution of gene size of all genes
M82_rMATs_anno_all_gene %>%
  filter(gene_size < 13000) %>%
  ggplot() +
  #geom_histogram(aes(size), bins = 50) +
  geom_histogram(aes(gene_size/1000), bins = 100) +
  scale_x_continuous(breaks = seq(0, max(M82_rMATs_anno_all_gene$gene_size)/1000, 2))

#check the maximum gene size of all genes
max(M82_rMATs_anno_all_gene$gene_size)

#make table of gene list with only their names and gene size
M82_rMATs_anno_all_gene_light <- M82_rMATs_anno_all_gene %>%
  select(anno, gene_size)

#get to know how many exons acquired by genes with size from 5 to 7 kb
#since most genes with AS sig events have size at this range
#This part is to make sure whether we can use similar strategy for ML as done by nat com paper (They remove the first and second exons)
M82_rMATs_anno_all_exon_order %>%
  left_join(M82_rMATs_anno_all_gene_light) %>%
  filter(gene_size >= 5000 & gene_size <= 7000) %>%
  group_by(anno) %>%
  summarise(count = n()) %>%
  group_by(count) %>%
  summarise(exon_count = n()) %>%
  mutate(sum_exon_count = sum(exon_count))

#create the histogram plots for sig (FDR < 0.05) and control (FDR > 0.05) sets
AS_sig = c("AS_sig_FDR05", "AS_control")
comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
AS = c("A5S", "A3S", "RI", "SE")

AS_intersected_gene_plot_l <- vector(mode = "list", length = length(AS_sig)*length(comp)*length(AS))

for (i in seq_along(AS_sig)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(AS)) {
      path <- str_c("./data/Intersected_AS_TPM_q05_M82_anno_rMATs/", AS_sig[[i]], "_TPM_q05_", comp[[j]], "_", AS[[k]], "_intersected_M82_rMATs_gene.bed")
      path_2 <- str_c("./analysis_output/AS_intersected_gene_size_plots/", AS_sig[[i]], "_TPM_q05_", comp[[j]], "_", AS[[k]], "_gene_size_p.jpeg")
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


#Combine all sig and non-sig AS with the information of intersected exons or introns
#Import all intron and exon intersected with different types of AS,
#Then combine all of these files into one file to intersect with the file of intron and exon with order information
AS_sig = c("AS_sig_FDR05", "AS_control")
comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
AS_type = c("A5S", "A3S", "SE", "RI")
feature = c("exon", "intron")

454436
all_AS_intersected_feature_raw <- tibble()

for (i in seq_along(AS_sig)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(AS_type)) {
      for (l in seq_along(feature)) {
        path_AS_mark <- str_c("./data/Intersected_AS_TPM_q05_M82_anno_rMATs/", AS_sig[[i]], "_TPM_q05_", comp[[j]], "_", AS_type[[k]], "_intersected_M82_rMATs_", feature[[l]], "_order.bed")
        
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
}

all_AS_intersected_feature_raw

all_AS_intersected_feature_raw %>%
  filter(anno == "Solyc03g019850.3.1") %>%
  View()

#combine all genes intersected with AS (four types)  
gene_having_AS <- all_AS_intersected_feature_raw %>%
  select(feature, anno, event) %>%
  arrange(anno) %>%
  distinct()

all_AS_intersected_feature_raw %>%
  select(anno) %>%
  distinct() %>%
  #group_by(feature) %>%
  summarise(count = n())

#Import all exon and intron with the order into
M82_rMATs_anno_all_exon_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_exon_order.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order"))
M82_rMATs_anno_all_intron_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_intron_order.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order"))

#combine all exon and intron as a list
gene_having_AS_info_l <- list(M82_rMATs_anno_all_exon_order, M82_rMATs_anno_all_intron_order)
names(gene_having_AS_info_l) <- c("exon", "intron")

#In this step, we keep all information of genes if they are intersected with any kind of AS (sig, non-sig)
#And we also produce a geneset without overlapping AS, this is for producing another kind of control (gene without AS) 
#gene with AS
gene_having_AS_info_l_f <- gene_having_AS_info_l %>%
  map(. %>% arrange(chr, str)) %>%
  map(. %>% left_join(gene_having_AS)) %>%
  map(. %>% drop_na()) %>%
  map(. %>% select(-event)) %>%
  map(. %>% distinct())

gene_having_AS_light <- gene_having_AS %>%
  select(anno) %>%
  distinct() %>%
  mutate(having_AS = "yes")

#gene without AS (external ctrl)
#since this is for choosing external ctrl, the first/second exon and intron are removed

gene_lacking_AS_info_l <- bind_rows(gene_having_AS_info_l) %>%
  left_join(gene_having_AS_light) %>%
  replace_na(list(having_AS = "no")) %>%
  filter(having_AS == "no") %>%
  select(-having_AS) %>%
  left_join(M82_rMATs_anno_all_gene_light) %>%
  mutate(gene_size_dis_500bp = gene_size %/% 500) %>%
  group_by(gene_size_dis_500bp) %>%
  mutate(nrow = n()) %>%
  mutate(nrow_remaining = nrow %% 5) %>%
  #making five groups with similar amount of genes with similar sizes
  mutate(sub_group = c(rep(c(1:5), each = n() %/% 5), sample(c(1:5), n() %% 5))) %>%
  mutate(sub_group_sh = sample(sub_group)) %>%
  mutate(sub_group_sh = str_c("sh", sub_group_sh)) %>%
  select(-nrow, -nrow_remaining, -sub_group) %>%
  filter(order != 1 & order != 2) %>%
  ungroup()

gene_lacking_AS_info_l %>%
  group_by(anno) %>%
  summarise(count = n())

write_delim(gene_lacking_AS_info_l, "./data/important_features_comp/external_ctrl.bed", delim = "\t", col_names = FALSE)

#have a quick look how many exons do genes with SE have in general 
gene_having_AS_info_l_f$exon %>%
  group_by(anno) %>%
  summarise(count = n()) %>%
  group_by(count) %>%
  summarise(count_sum = sum(count)) %>%
  View()

gene_lackng_AS_info_l_f$exon

#create the function for deciding which gene should be kept for the following ML analysis
#1. remove promoter effects (the first/second exon/intron are removed)
#2. only focusing on PI_H and PI_L (modified: focus on every categories)
#3. only focusing on genes with significant difference bw two comp (modified: extract non-sig and sig events)
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
    #mutate(PI_filter = if_else(PI == "PI_H" | PI == "PI_L", 1, 0)) %>%
    #group_by(anno) %>%
    #mutate(sum_PI_filter = sum(PI_filter)) %>%
    #filter(sum_PI_filter != 0) %>%
    #select(-PI_filter, -sum_PI_filter) %>%
    #check how many genes with more exons having AS than exons lacking AS 
    #select(anno, order, AS) %>%
    #distinct() %>%
    #group_by(anno) %>%
    #mutate(count_no_AS = if_else(AS == "no", 1, 0), count_AS = if_else(AS != "no", 1, 0)) %>%
    #mutate(sum_no_AS = sum(count_no_AS), sum_AS = sum(count_AS)) %>%
    #filter(sum_no_AS < sum_AS) %>%
    #remove genes without any sig events
    #mutate(sig_filter = if_else(event == "sig", 1, 0)) %>%
    #mutate(sum_sig_filter = sum(sig_filter)) %>%
    #filter(sum_sig_filter != 0) %>%
    #select(-sig_filter, -sum_sig_filter) %>%
    split(.$anno)
}

#create the list genes with sig events for the random selection of event without AS in the same gene
exon_data_for_ML_list <- list_for_ML_data_f(gene_having_AS_info_l_f$exon)
intron_data_for_ML_list <- list_for_ML_data_f(gene_having_AS_info_l_f$intron)

exon_intron_data_for_ML_list <- gene_having_AS_info_l_f %>%
  map(. %>% list_for_ML_data_f())

exon_intron_data_for_ML_list$exon$Solyc03g019850.3.1 %>%
  View()

length(exon_intron_data_for_ML_list$intron)
bind_rows(exon_intron_data_for_ML_list$exon) %>%
  group_by(PI) %>%
  summarise(count = n())

exon_intron_data_for_ML_list$exon$Solyc_chr1_NPC00337.MC1
gene_having_AS %>%
  filter(anno == "Solyc_chr1_NPC00337.MC1")
exon_intron_data_for_ML_list$intron$Solyc_chr1_NPC00337.MC1

exon_intron_data_for_ML_list$exon$Solyc_chr1_NPC00337.MC1


#create the function for random choosing the same amount of exon/intron as that of AS events occurred in that gene 
#(internal control)
make_AS_ML_data_f <- function(data, a){ 
  data_a <- data %>%
    ungroup()
  AS_sig_m <- data_a %>%
    filter(comp == a[1] & PI == a[2] &  AS == a[3] & event == "sig") 
  
  AS_non_sig_m <- data_a %>%
    filter(comp == a[1] & PI == a[2] &  AS == a[3] & event == "non-sig") 
  
  if(nrow(AS_sig_m)==0){
    AS_noAS_sig_m <- tibble()
  } else {
    AS_noAS_sig_m <- data_a %>%
      filter(event == "no") %>%
      filter(!order %in% c(c(unique(AS_sig_m$order)-1), unique(AS_sig_m$order), unique(AS_sig_m$order)+1))
  }
  
  if(nrow(AS_non_sig_m)==0){
    AS_noAS_non_sig_m <- tibble()
  } else {
    AS_noAS_non_sig_m <- data_a %>%
      filter(event == "no") %>%
      filter(!order %in% c(c(unique(AS_non_sig_m$order)-1), unique(AS_non_sig_m$order), unique(AS_non_sig_m$order)+1))
  }
  
  AS_noAS_sig_m_f <- head(AS_noAS_sig_m[sample(1:nrow(AS_noAS_sig_m)),], n = nrow(AS_sig_m)) 
  
  AS_sig_pair_f <- bind_rows(AS_sig_m, AS_noAS_sig_m_f) %>%
    mutate(pair = "sig")
  
  AS_noAS_non_sig_m_f <- head(AS_noAS_non_sig_m[sample(1:nrow(AS_noAS_non_sig_m)),], n = nrow(AS_non_sig_m)) 
  
  AS_non_sig_pair_f <- bind_rows(AS_non_sig_m, AS_noAS_non_sig_m_f) %>%
    mutate(pair = "non-sig")
  
  bind_rows(AS_sig_pair_f, AS_non_sig_pair_f)
}


#make list for parallelism 
comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
PI_type = c("PI_L", "PI_ML", "PI_M", "PI_MH", "PI_H")  
AS_type = c("A5S", "A3S", "SE", "RI")

comb_comp_PI_AS_list_raw <- list()

for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    for (k in seq_along(AS_type)) {
      comb_raw <- c(comp[[i]], PI_type[[j]], AS_type[[k]])
      comb_comp_PI_AS_list_raw <- append(comb_comp_PI_AS_list_raw, list(comb_raw))
    }
  }
}

#I will use 10 cores, so here I made the length of list as 6 (60/10 = 6)
comb_comp_PI_AS_list <- vector("list", length = 5)

for (i in seq_along(comb_comp_PI_AS_list)) {
  print(seq(12*(i-1)+1, 12*(i-1)+12, 1))
  comb_comp_PI_AS_list[[i]] <- comb_comp_PI_AS_list_raw[seq(12*(i-1)+1, 12*(i-1)+12, 1)]
}

#start to use 10 cores
registerDoParallel(cores = 12)
getDoParWorkers()

#create all control set from all genes with AS event whether these events are sig or non-sig
#internal control
exon_intron_data_for_ML_list_v2 <- 
  foreach(i=c(1:5), .packages = c("tidyverse")) %:%
  foreach(j=1:12, .packages = c("tidyverse")) %dopar% {
    exon_intron_data_for_ML_list %>%
      map(. %>% map(. %>% make_AS_ML_data_f(., comb_comp_PI_AS_list[[i]][[j]]))) %>%
      map(. %>% bind_rows())
    
  }

#modify the content in the columns of "comp", "AS" and "PI" for randomly selected events (they are noted as "no" previously)
exon_intron_data_for_ML_list_v3 <- 
  foreach(i=c(1:5), .packages = c("tidyverse")) %:%
  foreach(j=1:12, .packages = c("tidyverse")) %dopar% {
    exon_intron_data_for_ML_list_v2[[i]][[j]] %>%
      map(. %>% mutate(comp = if_else(comp == "no", comb_comp_PI_AS_list[[i]][[j]][[1]], comp),
                       AS = if_else(AS == "no", comb_comp_PI_AS_list[[i]][[j]][[3]], AS),
                       PI = if_else(PI == "no", comb_comp_PI_AS_list[[i]][[j]][[2]], PI))) 
    
  }

#reorganize the list into another form of list based on "comp", "PI", and "AS"
#plus, exon or intron were combined
#Since what I wanted to focus on are AS sig/non-sig and no-event segment
#whether these events are located in exon or intron is not an importnat issue
#one AS event can be located in an exon of a gene but an intron of the isoform of the same genes
exon_intron_data_for_ML_list_v4 <- exon_intron_data_for_ML_list_v3 %>%
  map(. %>% map(. %>% bind_rows())) %>%
  map(. %>% bind_rows()) %>%
  bind_rows() %>%
  split(.$comp) %>%
  map(. %>% split(.$PI)) %>%
  map(. %>% map(. %>% split(.$AS)))

#produce two sides (100 bp extension) of each exon or intron with or without AS for the following analysis (target and internal ctrl)
#(the intersection between these events and histone marks)
for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    for (k in seq_along(AS_type)) {
      #produce 100-bp extension on two sides of "target" and "internal ctrl"
      exon_intron_data_for_ML_f_combined <- bind_rows(exon_intron_data_for_ML_list_v4[[i]][[j]][[k]], exon_intron_data_for_ML_list_v4[[i]][[j]][[k]]) %>%
        ungroup() %>%
        mutate(side = rep(c("five", "three"), each = nrow(exon_intron_data_for_ML_list_v4[[i]][[j]][[k]]))) %>% 
        mutate(str_n = if_else(strand == "+" & side == "five", str - 100,
                               if_else(strand == "+" & side == "three", end,
                                       if_else(strand == "-" & side == "five", end, str - 100)))) %>%
        mutate(end_n = if_else(strand == "+" & side == "five", str,
                               if_else(strand == "+" & side == "three", end + 100,
                                       if_else(strand == "-" & side == "five", end + 100, str)))) %>%
        select(chr, str = str_n, end = end_n, strand, feature, source, anno, order, comp, AS, PI, event, pair, seg_side = side) %>%
        arrange(chr, str)
      
      #extract internal "ctrl"
      exon_intron_data_for_ML_f_ctrl <- exon_intron_data_for_ML_f_combined %>%
        filter(event == "no")
      
      #extract "target"
      exon_intron_data_for_ML_f_target <- exon_intron_data_for_ML_f_combined %>%
        filter(event != "no")
      
      #make the labels for the name of output file
      label_comp <- unique(exon_intron_data_for_ML_list_v4[[i]][[j]][[k]]$comp)
      label_PI <- unique(exon_intron_data_for_ML_list_v4[[i]][[j]][[k]]$PI)
      label_AS <- unique(exon_intron_data_for_ML_list_v4[[i]][[j]][[k]]$AS)
      
      #write output files
      write_delim(exon_intron_data_for_ML_f_target, str_c("./data/AS_bed_for_ML/segment_for_ML_", label_comp, "_", label_PI, "_", label_AS, "_", "target.bed"), col_names = FALSE, delim = "\t")
      write_delim(exon_intron_data_for_ML_f_ctrl, str_c("./data/AS_bed_for_ML/segment_for_ML_", label_comp, "_", label_PI, "_", label_AS, "_", "ctrl.bed"), col_names = FALSE, delim = "\t")
      
    }
  }
}

#produce two sides of 100-bp extension (external ctrl)
exon_intron_data_for_ML_f_ex_ctrl <- bind_rows(gene_lacking_AS_info_l, gene_lacking_AS_info_l) %>%
  mutate(side = rep(c("five", "three"), each = nrow(gene_lacking_AS_info_l))) %>% 
  mutate(str_n = if_else(strand == "+" & side == "five", str - 100,
                         if_else(strand == "+" & side == "three", end,
                                 if_else(strand == "-" & side == "five", end, str - 100)))) %>%
  mutate(end_n = if_else(strand == "+" & side == "five", str,
                         if_else(strand == "+" & side == "three", end + 100,
                                 if_else(strand == "-" & side == "five", end + 100, str)))) %>%
  mutate(comp = "no", AS = "no", PI = sub_group_sh, event = "no", pair = "external_ctrl") %>%
  select(chr, str = str_n, end = end_n, strand, feature, source, anno, order, comp, AS, PI, event, pair, seg_side = side) %>%
  arrange(chr, str) %>%
  split(.$PI)


pwalk(list(exon_intron_data_for_ML_f_ex_ctrl, str_c("./data/AS_bed_for_ML/segment_for_ML_", names(exon_intron_data_for_ML_f_ex_ctrl), "_ex_ctrl.bed")),
      write_delim, col_names = FALSE, delim = "\t" )

##after using bedtools intersect for target/ctrl with mark signals from all bases or peaks
##we can produce data for running ML
#stage 1 (based on each individual PSI)
#1. produce the information of histone modification marks including read count and name (read counts are used for calculating RPKM)
histone_mark_list <- read_delim("./data/Chip_seq/no_treatment/bam/mark_list", delim = "\t", col_names = c("name_raw"))

histone_mark_read_count <- c()

for (i in seq_along(histone_mark_list$name_raw)) {
  read_count_file <- read_delim(str_c("./data/Chip_seq/no_treatment/bam/", histone_mark_list$name_raw[[i]], "_read_count"), delim = "\t", col_names = c("read_count"))
  histone_mark_read_count <- append(histone_mark_read_count, read_count_file$read_count)
}

#this data is one of necessary file for running the function "histone_RPKM_f"
histone_mark_list_m <- histone_mark_list %>%
  mutate(read_count = histone_mark_read_count) %>%
  mutate(name_light = str_remove_all(name_raw, "0h_")) %>%
  mutate(name_light = str_remove_all(name_light, "_.*"))

#2.import all segments of exon/intron as the anchor for the following tables
#necessary for producing datasets for running ML
comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
PI_type = c("PI_H", "PI_MH", "PI_M", "PI_ML", "PI_L")
AS_type = c("A5S", "A3S", "SE", "RI")
target = c("target", "ctrl")


AS_exon_intron_sig_no_event_l <- list()
names_AS_exon_intron_sig_no_event_l <- c()

for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    for (k in seq_along(AS_type)) {
      for (l in seq_along(target)) {
        
        AS_table <- read_delim(str_c("./data/AS_bed_for_ML/segment_for_ML_", comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]], ".bed"), delim = "\t",
                                                                   col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order", "comp", "AS", "PI", "event", "pair", "seg_side")) %>%
          #the reason to have order_t is to separate five and three end of the same gene for taking the information of histone modification marks
          mutate(order_t = c(1:n()))
        names_AS_exon_intron_sig_no_event_l <- append(names_AS_exon_intron_sig_no_event_l, str_c("segment_for_ML_", comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]]))
        AS_exon_intron_sig_no_event_l <- append(AS_exon_intron_sig_no_event_l, list(AS_table))
      }
    }
  }
}

length(AS_exon_intron_sig_no_event_l)
AS_exon_intron_sig_no_event_l[[1]]
exon_intron_data_for_ML_f_ex_ctrl[[1]]
names(AS_exon_intron_sig_no_event_l) <- names_AS_exon_intron_sig_no_event_l
names_AS_exon_intron_sig_no_event_l
str_c("./data/AS_bed_for_ML/", histone_mark_list_m$name_raw[[1]], "_", names_AS_exon_intron_sig_no_event_l[[1]], "_12chr.bed")

#3. calculate the signal of marks in events containing (target) or not containing AS (ctrl)
#since we use RPKM to normalize the signal of 19 features for each exon/intron events from all samples
#RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 ) (from: https://www.metagenomics.wiki/pdf/qc/RPKM)
#geneLength is replaced by the exon or intron length

#create the function for the list of length of each intron/exon
feature_size_f <- function(data){
  data %>%
    group_by(anno, feature, order) %>%
    mutate(str_n = min(end), end_n = max(str)) %>%
    mutate(feature_size = end_n - str_n) %>%
    ungroup() %>%
    select(anno, feature, order, seg_side, order_t, feature_size) %>%
    distinct() %>%
    mutate(order_t = c(1:n()))
}

#make the list of feature size for all target or ctrl segments
AS_exon_intron_sig_no_event_size_l <- 
  AS_exon_intron_sig_no_event_l %>% 
  map(., feature_size_f)

#create the function for calculating RPKM for each mark from each sample
histone_RPKM_f <- function(a, b, data, data2) {
  #a is the order in histone modification mark
  #b is the middle name of file for the intersection between feature segments and marks
  #data is names_AS_exon_intron_sig_no_event_l
  #data2 is AS_exon_intron_sig_no_event_size_l
  histone_RPKM <- read_delim(
    str_c(
      "./data/AS_bed_for_ML/all_segments/",
      histone_mark_list_m$name_raw[a[1]],
      b[1],
      str_remove(data, "segment_"),
      "_12chr.bed.gz"
    ),
    delim = "\t",
    col_names = c(
      "chr",
      "str",
      "end",
      "strand",
      "feature",
      "source",
      "anno",
      "order",
      "comp",
      "AS",
      "PI",
      "event",
      "pair",
      "seg_side",
      "AS_signal"
    )
  ) %>%
    group_by(anno, feature, order, comp, AS, PI, event, seg_side) %>%
    summarise(mean_signal = sum(AS_signal) / n()) %>%
    mutate(read_count = histone_mark_list_m$read_count[a[1]]) %>%
    mutate(RPM = mean_signal / read_count * 10 ^ 6) %>%
    left_join(data2) %>%
    mutate(RPKM = RPM / feature_size) %>%
    ungroup() %>%
    select(order_t, RPKM) %>%
    distinct()
  
  names(histone_RPKM) <-
    c("order_t", str_c(histone_mark_list_m$name_light[a[1]], "_RPKM"))
  histone_RPKM
}

#a) based on all detected signal of histone modification marks

for (i in seq_along(AS_exon_intron_sig_no_event_l)) {
  #generate the list with 19 files of feature size
  feature_size_list <- AS_exon_intron_sig_no_event_size_l[[i]][rep(1:nrow(AS_exon_intron_sig_no_event_size_l[[i]]),19),] %>%
    mutate(list_group = rep(c(1:19), each = nrow(AS_exon_intron_sig_no_event_size_l[[i]]))) %>%
    split(.$list_group)
  
  #generate the RPKM file for each "AS_exon_intron_sig_no_event_l" based on signals of all bases
  RPKM_list <-
    pmap(list(seq(1:19), rep("_", 19), 
              rep(names_AS_exon_intron_sig_no_event_l[[i]], 19), feature_size_list), 
         histone_RPKM_f) %>%
    reduce(., full_join)
  
  #produce the output file ready for ML
  ready_ML <- AS_exon_intron_sig_no_event_l[[i]] %>%
    left_join(RPKM_list) %>%
    select(-order_t)
  
  write_delim(ready_ML, str_c("./data/AS_bed_for_ML/all_segments/table_ready_", str_remove(names_AS_exon_intron_sig_no_event_l[[i]], "segment_")),
              delim = "\t", col_names = TRUE)
}


#b) based on all detected peaks of histone modification marks

for (i in seq_along(AS_exon_intron_sig_no_event_l)) {
  #generate the list with 19 files of feature size
  feature_size_list <- AS_exon_intron_sig_no_event_size_l[[i]][rep(1:nrow(AS_exon_intron_sig_no_event_size_l[[i]]),19),] %>%
    mutate(list_group = rep(c(1:19), each = nrow(AS_exon_intron_sig_no_event_size_l[[i]]))) %>%
    split(.$list_group)
  
  #generate the RPKM file for each "AS_exon_intron_sig_no_event_l" based on peaks
  RPKM_list <-
    pmap(list(seq(1:19), rep("_p0.05_peaks_pileup_", 19), 
              rep(names_AS_exon_intron_sig_no_event_l[[i]], 19), feature_size_list), 
         histone_RPKM_f) %>%
    reduce(., full_join)
  
  #produce the output file ready for ML
  ready_ML <- AS_exon_intron_sig_no_event_l[[i]] %>%
    left_join(RPKM_list) %>%
    select(-order_t)
  
  write_delim(ready_ML, str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_", str_remove(names_AS_exon_intron_sig_no_event_l[[i]], "segment_")),
              delim = "\t", col_names = TRUE)
}

##Run ML
#In the previous part, the output data from method (a) and (b) still contain some information
#which have to be removed from the table if one wants to perform ML
#import RPKM raw data
comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
PI_type = c("PI_H", "PI_MH", "PI_M", "PI_ML", "PI_L")
AS_type = c("A5S", "A3S", "SE", "RI")
target = c("target", "ctrl")

path_RPKM_AS_all_seg <- str_c("./data/AS_bed_for_ML/all_segments/table_ready_for_ML_", 
                              comp[[1]], "_", PI_type[[1]], "_", AS_type[[1]], "_", target[[1]])
#all signals of marks
table_ready_for_ML_raw <- list()
name_table_ready_for_ML_raw <- c()

for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    for (k in seq_along(AS_type)) {
      for (l in seq_along(target)) {
        #print(str_c(comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", event[[l]]))
        path_RPKM_AS_all_seg <- str_c("./data/AS_bed_for_ML/all_segments/table_ready_for_ML_", 
                                      comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]])
        raw_table <- read_delim(path_RPKM_AS_all_seg, delim = "\t") %>%
          mutate(type = as.factor(str_c(comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]])))
        
        raw_table <- raw_table[c(9:11, (ncol(raw_table)-19):ncol(raw_table))] 
        table_ready_for_ML_raw <- append(table_ready_for_ML_raw, list(raw_table))
        name_table_ready_for_ML_raw <- append(name_table_ready_for_ML_raw, str_c(comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]]))
      }
    }
  }
}

names(table_ready_for_ML_raw) <- name_table_ready_for_ML_raw

table_ready_for_ML_raw_target <- table_ready_for_ML_raw[which(str_detect(name_table_ready_for_ML_raw, "target"))]

#only peaks of marks
table_ready_for_ML_raw_peak <- list()
name_table_ready_for_ML_raw_peak <- c()

for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    for (k in seq_along(AS_type)) {
      for (l in seq_along(target)) {
        print(str_c(comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]]))
        path_RPKM_AS_all_seg <- str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_for_ML_", 
                                      comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]])
        raw_table <- read_delim(path_RPKM_AS_all_seg, delim = "\t") %>%
          mutate(type = as.factor(str_c(comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]])))
        
        raw_table <- raw_table[c(9:11, (ncol(raw_table)-19):ncol(raw_table))] 
        table_ready_for_ML_raw_peak <- append(table_ready_for_ML_raw_peak, list(raw_table))
        name_table_ready_for_ML_raw_peak <- append(name_table_ready_for_ML_raw_peak, str_c(comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]]))
      }
    }
  }
}

names(table_ready_for_ML_raw_peak) <- name_table_ready_for_ML_raw_peak

table_ready_for_ML_raw_peak_target <- table_ready_for_ML_raw_peak[which(str_detect(name_table_ready_for_ML_raw_peak, "target"))]
table_ready_for_ML_raw_peak_ctrl <- table_ready_for_ML_raw_peak[which(str_detect(name_table_ready_for_ML_raw_peak, "ctrl"))]


##sort the list of target table by the order of variables followed: AS, comp, PI
#all signals of marks
table_ready_for_ML_raw_target_sorted <- table_ready_for_ML_raw_target %>%
  map(. %>% as_tibble()) %>%
  bind_rows() %>%
  split(.$AS) %>%
  map(. %>% split(.$comp)) %>%
  map(. %>% map(. %>% split(.$PI))) %>%
  map(. %>% map(. %>% map(. %>% select(-AS, -comp, -PI))))

#peaks of marks
#target
table_ready_for_ML_raw_peak_target_sorted <- table_ready_for_ML_raw_peak_target %>%
  bind_rows() %>%
  replace_na(list(Pol2_RPKM = 0, H3K14ac_RPKM = 0, H3K4ac_RPKM = 0, K18ac_RPKM = 0, K27ac_RPKM = 0,
                  K27me3_RPKM = 0, K4me1_RPKM = 0, K4me3_RPKM = 0, K9ac_RPKM = 0, H3K36me3_RPKM = 0,
                  H3K79ac_RPKM = 0, H4K12ac_RPKM = 0, H4K16ac_RPKM = 0, H4K20ac_RPKM = 0, H4K5ac_RPKM = 0,
                  H4K8ac_RPKM = 0, K36ac_RPKM = 0, K4me2_RPKM = 0, K56ac_RPKM = 0)) %>%
  split(.$AS) %>%
  map(. %>% split(.$comp)) %>%
  map(. %>% map(. %>% split(.$PI))) %>%
  map(. %>% map(. %>% map(. %>% select(-AS, -comp, -PI))))

#ctrl  
table_ready_for_ML_raw_peak_ctrl_sorted <- table_ready_for_ML_raw_peak_ctrl %>%
  bind_rows() %>%
  replace_na(list(Pol2_RPKM = 0, H3K14ac_RPKM = 0, H3K4ac_RPKM = 0, K18ac_RPKM = 0, K27ac_RPKM = 0,
                  K27me3_RPKM = 0, K4me1_RPKM = 0, K4me3_RPKM = 0, K9ac_RPKM = 0, H3K36me3_RPKM = 0,
                  H3K79ac_RPKM = 0, H4K12ac_RPKM = 0, H4K16ac_RPKM = 0, H4K20ac_RPKM = 0, H4K5ac_RPKM = 0,
                  H4K8ac_RPKM = 0, K36ac_RPKM = 0, K4me2_RPKM = 0, K56ac_RPKM = 0)) %>%
  split(.$AS) %>%
  map(. %>% split(.$comp)) %>%
  map(. %>% map(. %>% split(.$PI))) %>%
  map(. %>% map(. %>% map(. %>% select(-AS, -comp, -PI))))
  
#the order of names in the sorted list
names(table_ready_for_ML_raw_sorted)
names(table_ready_for_ML_raw_sorted[[1]])
names(table_ready_for_ML_raw_sorted[[1]][[1]])

table_ready_for_ML_raw_target_sorted[[1]][[1]][[1]] %>%
  #View()
  group_by(type) %>%
  summarise()


length(table_ready_for_ML_raw_sorted[[1]])

length(table_ready_for_ML_raw_sorted[[1]][[1]])

##make pairwise table for running Boruta 
comb_pair_AS_seg_type <- combs(c(1:5), 2) #five PSI type

#all signals
table_ready_for_ML_stage_1_raw_target <- list()
name_pairwise_table_for_ML_stage_1_raw_target <- list()

for (i in seq_along(AS_type)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(comb_pair_AS_seg_type[,1])) {
      a <- table_ready_for_ML_raw_target_sorted[[i]][[j]][comb_pair_AS_seg_type[k,]] %>%
        bind_rows()
      
      type_check <- a %>%
        select(type) %>%
        unique()
      print(type_check)
      a <- a %>%
        as.data.frame()
      table_ready_for_ML_stage_1_raw_target <- append(table_ready_for_ML_stage_1_raw_target, list(a))
      name_pairwise_table_for_ML_stage_1_raw_target <- append(name_pairwise_table_for_ML_stage_1_raw_target, list(type_check))
    }
  }
}


for (i in seq_along(table_ready_for_ML_stage_1_raw_target)) {
  print(nrow(table_ready_for_ML_stage_1_raw_target[[i]]))
}

table_ready_for_ML_stage_1_target_boruta <- vector("list", length = length(table_ready_for_ML_stage_1_raw_target))

for (i in seq_along(table_ready_for_ML_stage_1_target_boruta)) {
  print(str_c("start boruta for the ", i, "th dataset"))
  table_ready_for_ML_stage_1_target_boruta[[i]] <- Boruta(type~., table_ready_for_ML_stage_1_raw_target[[i]], doTrace = 2)
}


#peaks
#target
table_ready_for_ML_stage_1_raw_peak_target <- list()
name_pairwise_table_for_ML_stage_1_raw_peak_target <- list()
comb_pair_AS_seg_type[1,]
for (i in seq_along(AS_type)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(comb_pair_AS_seg_type[,1])) {
      a <- table_ready_for_ML_raw_peak_target_sorted[[i]][[j]][comb_pair_AS_seg_type[k,]] %>%
        bind_rows()
      
      type_check <- a %>%
        select(type) %>%
        unique()
      print(type_check)
      a <- a %>%
        as.data.frame()
      table_ready_for_ML_stage_1_raw_peak_target <- append(table_ready_for_ML_stage_1_raw_peak_target, list(a))
      name_pairwise_table_for_ML_stage_1_raw_peak_target <- append(name_pairwise_table_for_ML_stage_1_raw_peak_target, list(type_check))
    }
  }
}
unique(table_ready_for_ML_stage_1_raw_peak_target[[1]]$type)

table_ready_for_ML_stage_1_peak_target_boruta <- vector("list", length = length(table_ready_for_ML_stage_1_raw_peak_target))

for (i in seq_along(table_ready_for_ML_stage_1_peak_target_boruta)) {
  print(str_c("start boruta for the ", i, "th dataset"))
  table_ready_for_ML_stage_1_peak_target_boruta[[i]] <- Boruta(type~., table_ready_for_ML_stage_1_raw_peak_target[[i]], doTrace = 2)
}

#ctrl
table_ready_for_ML_stage_1_raw_peak_ctrl <- list()
name_pairwise_table_for_ML_stage_1_raw_peak_ctrl <- list()

for (i in seq_along(AS_type)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(comb_pair_AS_seg_type[,1])) {
      a <- table_ready_for_ML_raw_peak_ctrl_sorted[[i]][[j]][comb_pair_AS_seg_type[k,]] %>%
        bind_rows()
      
      type_check <- a %>%
        select(type) %>%
        unique()
      print(type_check)
      a <- a %>%
        as.data.frame()
      table_ready_for_ML_stage_1_raw_peak_ctrl <- append(table_ready_for_ML_stage_1_raw_peak_ctrl, list(a))
      name_pairwise_table_for_ML_stage_1_raw_peak_ctrl <- append(name_pairwise_table_for_ML_stage_1_raw_peak_ctrl, list(type_check))
    }
  }
}

unique(table_ready_for_ML_stage_1_raw_peak_ctrl[[1]]$type)

table_ready_for_ML_stage_1_peak_ctrl_boruta <- vector("list", length = length(table_ready_for_ML_stage_1_raw_peak_ctrl))

for (i in c(1:110)) {
  print(str_c("start boruta for the ", i+10, "th dataset"))
  table_ready_for_ML_stage_1_peak_ctrl_boruta[[i+10]] <- Boruta(type~., table_ready_for_ML_stage_1_raw_peak_ctrl[[i+10]], doTrace = 2)
}

##prepare the list of names for producing plots by ggplot

name_ML_stage_1_target_for_plot_f <- function(data){
  test_name_list <- str_split(data$type, "_")
  str_c(test_name_list[[1]][[1]], "_", test_name_list[[1]][[2]], 
        "_", test_name_list[[1]][[3]], "_", test_name_list[[1]][[4]], 
        "_vs_", test_name_list[[2]][[3]], "_", test_name_list[[2]][[4]],
        "_", test_name_list[[1]][[5]], "_", test_name_list[[1]][[6]])
}

#all signals
name_ML_stage_1_target_for_plot <- name_pairwise_table_for_ML_stage_1_raw_target %>%
  map(. %>% name_ML_stage_1_target_for_plot_f())

#peaks
name_ML_stage_1_target_peaks_for_plot <- name_pairwise_table_for_ML_stage_1_raw_peak_target %>%
  map(. %>% name_ML_stage_1_target_for_plot_f())

name_ML_stage_1_ctrl_peaks_for_plot <- name_pairwise_table_for_ML_stage_1_raw_peak_ctrl %>%
  map(. %>% name_ML_stage_1_target_for_plot_f())

##reference modify Boruta's result
##use this function to modify the output of boruta and produce plots
#https://stackoverflow.com/questions/73415232/how-to-use-ggplot2-to-plot-box-plots-from-borutas-results-in-r
process_the_Boruta_data <- function(x, whichShadow=c(TRUE,TRUE,TRUE),
                                    colCode=c('green','yellow','red','blue'),
                                    col=NULL) {
  if(is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.')
  
  #Removal of -Infs and conversion to a list
  lz <- lapply(1:ncol(x$ImpHistory),
               function(i) x$ImpHistory[is.finite(x$ImpHistory[,i]),i])
  colnames(x$ImpHistory) -> names(lz)
  
  #Selection of shadow meta-attributes
  numShadow <- sum(whichShadow)
  lz[c(rep(TRUE,length(x$finalDecision)),whichShadow)] -> lz
  
  generateCol<-function(x,colCode,col,numShadow){
    #Checking arguments
    if(is.null(col) & length(colCode)!=4)
      stop('colCode should have 4 elements.')
    #Generating col
    if(is.null(col)){
      rep(colCode[4],length(x$finalDecision)+numShadow)->cc
      cc[c(x$finalDecision=='Confirmed',rep(FALSE,numShadow))]<-colCode[1]
      cc[c(x$finalDecision=='Tentative',rep(FALSE,numShadow))]<-colCode[2]
      cc[c(x$finalDecision=='Rejected',rep(FALSE,numShadow))]<-colCode[3]
      col=cc
    }
    return(col)
  }
  
  #Generating color vector
  col <- generateCol(x, colCode, col, numShadow)
  
  #Ordering boxes due to attribute median importance
  ii<-order(sapply(lz,stats::median))
  lz[ii] -> lz
  col <- col[ii]
  lz_df <- do.call(rbind.data.frame, lz)
  df <- as.data.frame(t(lz_df))
  names(df) <- names(lz)
  rownames(df) <- NULL
  return(df)
}
#the function of raw ggplot
plot_boruta <- function(data){
  data %>%
    pivot_longer(everything()) %>%
    mutate(name = str_remove(name, "_RPKM")) %>%
    mutate(mark = if_else(str_detect(name, "shadow")==TRUE, "shadowmarks", "marks")) %>%
    ggplot(aes(x = fct_reorder(name, value, median), y = value)) +
    geom_boxplot(aes(fill = mark))
}

#ggplot customized theme function
theme_ym <- function(data){
  data +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 36),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 28, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
          axis.text.x = element_text(colour = "black", size = 32, face = "bold", angle = 45, vjust = 0.5))
}

#all signals
for (i in seq_along(name_ML_stage_1_target_for_plot)) {
  output_table <- process_the_Boruta_data(table_ready_for_ML_stage_1_target_boruta[[i]])
  write_delim(output_table, str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_1/", name_ML_stage_1_target_for_plot[[i]]), delim = "\t", col_names = TRUE)
  
  g <- process_the_Boruta_data(table_ready_for_ML_stage_1_target_boruta[[i]]) %>%
    plot_boruta() +
    ggtitle(name_ML_stage_1_target_for_plot[[i]])
  
  output_plot <- theme_ym(g)
  
  ggsave(str_c("./analysis/AS_all_seg_ML_stage_1/", name_ML_stage_1_target_for_plot[[i]], ".jpeg"), output_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
  
}


#peaks
#target
temp_ml_s1_target <- list()

for (i in seq_along(name_ML_stage_1_target_peaks_for_plot)) {
  #output_table <- process_the_Boruta_data(table_ready_for_ML_stage_1_peak_target_boruta[[i]])
  #write_delim(output_table, str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_1_peaks/", name_ML_stage_1_target_peaks_for_plot[[i]], "_peaks"), delim = "\t", col_names = TRUE)
  
  #g <- process_the_Boruta_data(table_ready_for_ML_stage_1_peak_target_boruta[[i]]) %>%
    #plot_boruta() +
    #ggtitle(str_c(name_ML_stage_1_target_peaks_for_plot[[i]], "_peaks")) 
  
  raw_table <- read_delim(str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_1_peaks/", name_ML_stage_1_target_peaks_for_plot[[i]], "_peaks"), delim = "\t", col_names = TRUE)
  
  temp_ml_s1_target <- append(temp_ml_s1_target, list(raw_table))
  g <- raw_table %>%
  plot_boruta() +
    ggtitle(str_c(name_ML_stage_1_target_peaks_for_plot[[i]], "_peaks"))
  
  output_plot <- theme_ym(g)
  
  
  ggsave(str_c("./analysis/AS_all_seg_ML_stage_1_peaks/", name_ML_stage_1_target_peaks_for_plot[[i]], "_peaks.jpeg"), output_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
  
}

bind_rows(temp_ml_s1_target) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(q = mean(value)) %>%
  arrange(q) %>%
  View()

#ctrl
for (i in c(39)) {
  print(str_c(i, "th"))
  print(nrow(process_the_Boruta_data(table_ready_for_ML_stage_1_peak_ctrl_boruta[[i]])))
  output_table <- process_the_Boruta_data(table_ready_for_ML_stage_1_peak_ctrl_boruta[[i]]) 
  #print(warnings())
  #write_delim(output_table, str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_1_peaks/", name_ML_stage_1_ctrl_peaks_for_plot[[i]], "_peaks"), delim = "\t", col_names = TRUE)
}  

##stage 2: To use randomized PSI to perform feature selection
#write the function suitable for randomize PSI (eg. data = table_ready_for_ML_raw_peak_target_sorted$A3S$HS6_HS0$PS_H)
regroup_ramdom_PSI_f <- function(data){
  
  str_loc_remove <- str_locate_all(unique(data$type), "_") 
  
  str_loc_remove_f <- unlist(str_loc_remove)[c(2,4)]
  
  data %>%
    mutate(group = c(rep(c(1:5), each = nrow(.) %/% 5), rep(1, nrow(.) %% 5))) %>%
    mutate(group = sample(group)) %>%
    mutate(type_2 = str_sub(type, str_loc_remove_f[[1]] + 1, str_loc_remove_f[[2]])) %>%
    mutate(type = str_remove(type, type_2)) %>%
    mutate(type = str_c(type, "_", group)) %>%
    select(-type_2, -group)
}

str_locate_all(unique(table_ready_for_ML_raw_peak_target_sorted$A3S$HS6_HS0$PI_H$type), "_")
regroup_ramdom_PSI_f(table_ready_for_ML_raw_peak_target_sorted$A3S$HS6_HS0$PI_H)

#sort table_ready_for_ML_raw_peak_target_sorted into random PSI
#target
table_ready_for_ML_raw_peak_target_sorted_ran_PSI <- table_ready_for_ML_raw_peak_target_sorted %>%
  map(. %>% map(. %>% map(. %>% regroup_ramdom_PSI_f))) %>%
  map(. %>% map(. %>% bind_rows)) %>%
  map(. %>% map(. %>% split(.$type))) 

#make pairwise table for running Boruta (stage 2: based on randomized PSI)
comb_pair_AS_seg_type <- combs(c(1:5), 2) #five PSI type

#peaks
table_ready_for_ML_stage_2_raw_target <- list()
name_pairwise_table_for_ML_stage_2_raw_target <- list()

for (i in seq_along(AS_type)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(comb_pair_AS_seg_type[,1])) {
      a <- table_ready_for_ML_raw_peak_target_sorted_ran_PSI[[i]][[j]][comb_pair_AS_seg_type[k,]] %>%
        bind_rows() %>%
        mutate(type = as.factor(type))
      
      type_check <- a %>%
        select(type) %>%
        unique()
      print(type_check)
      a <- a %>%
        as.data.frame()
      table_ready_for_ML_stage_2_raw_target <- append(table_ready_for_ML_stage_2_raw_target, list(a))
      name_pairwise_table_for_ML_stage_2_raw_target <- append(name_pairwise_table_for_ML_stage_2_raw_target, list(type_check))
    }
  }
}
name_pairwise_table_for_ML_stage_2_raw_target[[63]]
for (i in seq_along(table_ready_for_ML_stage_2_raw_target)) {
  print(nrow(table_ready_for_ML_stage_2_raw_target[[i]]))
}

table_ready_for_ML_stage_2_target_boruta <- vector("list", length = length(table_ready_for_ML_stage_2_raw_target))

for (i in seq_along(table_ready_for_ML_stage_2_target_boruta)) {
  print(str_c("start boruta for the ", i, "th dataset"))
  table_ready_for_ML_stage_2_target_boruta[[i]] <- Boruta(type~., table_ready_for_ML_stage_2_raw_target[[i]], doTrace = 2)
}
table_ready_for_ML_stage_2_raw_target[[1]]

str_split(name_pairwise_table_for_ML_stage_2_raw_target[[1]]$type, "_")

setdiff(name_pairwise_table_for_ML_stage_2_raw_target[[1]]$type[[1]],
        name_pairwise_table_for_ML_stage_2_raw_target[[1]]$type[[2]])

for (i in seq_along(name_pairwise_table_for_ML_stage_2_raw_target)) {
  #output_table <- process_the_Boruta_data(table_ready_for_ML_stage_2_target_boruta[[i]])
  #write_delim(output_table, str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_2_peaks/random_AS_PSI_peaks_", i), delim = "\t", col_names = TRUE)
  #print(str_c(i, "th boruta result"))
  #print(sum(attStats(table_ready_for_ML_stage_2_target_boruta[[i]])[,6]=="Rejected"))
  
  title <- str_split(name_pairwise_table_for_ML_stage_2_raw_target[[i]]$type, "_")
  print(str_c(title[[1]][[1]], "_", title[[1]][[2]], "_", title[[1]][[3]], "_sh", title[[1]][[5]], "_vs_sh", title[[2]][[5]]))
  
  title_f <- str_c(title[[1]][[1]], "_", title[[1]][[2]], "_", title[[1]][[3]], "_PI_sh", title[[1]][[5]], "_vs_PI_sh", title[[2]][[5]])
  
  g <- read_delim(str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_2_peaks/random_AS_PSI_peaks_", i), delim = "\t", col_names = TRUE) %>%
    plot_boruta() +
    ggtitle(str_c(title_f, "_peaks")) 

    #g <- process_the_Boruta_data(table_ready_for_ML_stage_1_peak_target_boruta[[i]]) %>%
    #plot_boruta() +
    #ggtitle(str_c(name_ML_stage_1_target_peaks_for_plot[[i]], "_peaks")) 
  
  output_plot <- theme_ym(g)
  
  ggsave(str_c("./analysis/AS_all_seg_ML_stage_2_peaks/", title_f, "_peaks.jpeg"), output_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
  
}


print(sum(attStats(table_ready_for_ML_stage_2_target_boruta[[1]])[,6]=="Rejected"))
attStats(table_ready_for_ML_stage_2_target_boruta[[100]])

##stage 3: non-sig vs sig
#produce data for running ML
table_ready_for_ML_raw_peak_s3 <- list()
name_table_ready_for_ML_raw_peak_s3 <- c()

for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    for (k in seq_along(AS_type)) {
      for (l in seq_along(target)) {
        print(str_c(comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]]))
        path_RPKM_AS_all_seg <- str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_for_ML_", 
                                      comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]])
        raw_table <- read_delim(path_RPKM_AS_all_seg, delim = "\t") %>%
          replace_na(list(Pol2_RPKM = 0, H3K14ac_RPKM = 0, H3K4ac_RPKM = 0, K18ac_RPKM = 0, K27ac_RPKM = 0,
                          K27me3_RPKM = 0, K4me1_RPKM = 0, K4me3_RPKM = 0, K9ac_RPKM = 0, H3K36me3_RPKM = 0,
                          H3K79ac_RPKM = 0, H4K12ac_RPKM = 0, H4K16ac_RPKM = 0, H4K20ac_RPKM = 0, H4K5ac_RPKM = 0,
                          H4K8ac_RPKM = 0, K36ac_RPKM = 0, K4me2_RPKM = 0, K56ac_RPKM = 0))
        
        raw_table <- raw_table[c(13, (ncol(raw_table)-18):ncol(raw_table))] %>%
          as.data.frame()
        table_ready_for_ML_raw_peak_s3 <- append(table_ready_for_ML_raw_peak_s3, list(raw_table))
        name_table_ready_for_ML_raw_peak_s3 <- append(name_table_ready_for_ML_raw_peak_s3, str_c(comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", target[[l]]))
      }
    }
  }
}

table_ready_for_ML_raw_peak_s3 <- table_ready_for_ML_raw_peak_s3 %>%
  map(. %>% mutate(pair = as.factor(pair)))

names(table_ready_for_ML_raw_peak_s3) <- name_table_ready_for_ML_raw_peak_s3

#running ML
table_ready_for_ML_raw_peak_s3_boruta <- vector("list", length = length(table_ready_for_ML_raw_peak_s3))

for (i in seq_along(table_ready_for_ML_raw_peak_s3_boruta)) {
  print(str_c("start boruta for the ", i, "th dataset"))
  table_ready_for_ML_raw_peak_s3_boruta[[i]] <- Boruta(pair~., table_ready_for_ML_raw_peak_s3[[i]], doTrace = 2)
}
AS_all_seg_ML_s3_peaks
ML_results_table_stage_2_peaks

#produce the output files
for (i in seq_along(table_ready_for_ML_raw_peak_s3_boruta)) {
  output_table <- process_the_Boruta_data(table_ready_for_ML_raw_peak_s3_boruta[[i]])
  write_delim(output_table, str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_3_peaks/", name_table_ready_for_ML_raw_peak_s3[[i]], "_sig_vs_non_sig"), delim = "\t", col_names = TRUE)
  
  g <- process_the_Boruta_data(table_ready_for_ML_raw_peak_s3_boruta[[i]]) %>%
    plot_boruta() +
  ggtitle(str_c(name_table_ready_for_ML_raw_peak_s3[[i]], "_peaks_sig_vs_non_sig")) 
  
  output_plot <- theme_ym(g)
  
  ggsave(str_c("./analysis/AS_all_seg_ML_s3_peaks/", name_table_ready_for_ML_raw_peak_s3[[i]], "_peaks_sig_vs_non_sig.jpeg"), output_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
  
}

table_ready_for_ML_raw_peak_s3_boruta %>%
  map(. %>% process_the_Boruta_data()) %>%
  bind_rows() %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(mean_value = mean(value)) %>%
  arrange(-mean_value) %>%
  View()

warnings()


#stage 4 PL_M compared with the extremity of PI_H and PI_L
#import all segments of exon/intron as the anchor for the following tables
#necessary for producing datasets for running ML
#in this stage, the strategy to produce files is to extract the extremity of PI_L (0) & PI_H (1) from all PI_H & PI_L

comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
PI_type_s4_ex = c("PI_H", "PI_L")
AS_type = c("A5S", "A3S", "SE", "RI")
target_s4 = c("target")


AS_exon_intron_sig_no_event_l_s4 <- list()
names_AS_exon_intron_sig_no_event_l_s4 <- c()

for (i in seq_along(comp)) {
  for (j in seq_along(PI_type_s4_ex)) {
    for (k in seq_along(AS_type)) {
      for (l in seq_along(target_s4)) {
        
        AS_table <- read_delim(str_c("./data/AS_bed_for_ML/segment_for_ML_", comp[[i]], "_", PI_type_s4_ex[[j]], "_", AS_type[[k]], "_", target_s4[[l]], ".bed"), delim = "\t",
                               col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order", "comp", "AS", "PI", "event", "pair", "seg_side")) %>%
          mutate(order_t = c(1:n()))
        names_AS_exon_intron_sig_no_event_l_s4 <- append(names_AS_exon_intron_sig_no_event_l_s4, str_c("segment_for_ML_", comp[[i]], "_", PI_type_s4_ex[[j]], "_", AS_type[[k]], "_", target_s4[[l]]))
        AS_exon_intron_sig_no_event_l_s4 <- append(AS_exon_intron_sig_no_event_l_s4, list(AS_table))
      }
    }
  }
}

names(AS_exon_intron_sig_no_event_l_s4) <- names_AS_exon_intron_sig_no_event_l_s4

names_AS_exon_intron_sig_no_event_l_s4

#import all AS with PI extremity
#this extremity data was produced by "all_AS_events_TPM.R"
all_AS_events_bed_TPM_extremity_PI_sorted <- read_delim("./data/all_AS_events_bed_TPM_extremity_PI", delim = "\t") %>%
  select(chr, comp, AS = AS_type, PI, str_event, end_event) %>%
  arrange(chr, str_event) %>%
  group_by(comp, AS, PI, chr) %>%
  mutate(label_event = c(1:n())) 
  


#sort the integrated file (all_AS_events_bed_TPM_extremity_PI) with the extremity of PI_H/PI_L based on comp/PI_type_s4_ex/AS_type
all_AS_events_bed_TPM_extremity_PI_sorted_l = list()
names_all_AS_events_bed_TPM_extremity_PI_sorted_l = list()


for (i in seq_along(comp)) {
  for (j in seq_along(PI_type_s4_ex)) {
    for (k in seq_along(AS_type)) {
      var <- c(comp[[i]], PI_type_s4_ex[[j]], AS_type[[k]])
      a <- all_AS_events_bed_TPM_extremity_PI_sorted %>%
        filter(comp == var[1] & PI == var[2] & AS == var[3])
      
      b <- str_c("PI_ex_", comp[[i]], "_", PI_type_s4_ex[[j]], "_", AS_type[[k]])
      all_AS_events_bed_TPM_extremity_PI_sorted_l <- append(all_AS_events_bed_TPM_extremity_PI_sorted_l, list(a))
      names_all_AS_events_bed_TPM_extremity_PI_sorted_l <- append(names_all_AS_events_bed_TPM_extremity_PI_sorted_l, b)
    }
  }
}

names(all_AS_events_bed_TPM_extremity_PI_sorted_l) <- names_all_AS_events_bed_TPM_extremity_PI_sorted_l


#produce the function to extract the extremity of PI_H/PI_L from the files with all segments
#since these segments are the 100-bp extension of events (intron or exon)
#this function is to keep events containing the extreme PI_H/PI_L
filter_AS_ex_f <- function(data, data2){
 
 #make the list of each extreme PI_H/PI_L event
 a <- data2 %>% 
    split(.$label_event)
 
 #make the empty tibble for followinly merging
 b <- tibble()
 
 #run loops to keep and merge all events which have extreme PI_H/PI_L
 for (i in seq_along(c(1:length(a)))) {
   c <- data %>%
     arrange(anno, chr, str) %>%
     group_by(anno, order) %>%
     mutate(str_segment = min(end), end_segment = max(str)) %>%
     left_join(a[[i]]) %>%
     filter(str_event >= str_segment & end_event <= end_segment)
   b <- bind_rows(b, c)
 }
 b
}

test <- head(all_AS_events_bed_TPM_extremity_PI_sorted_l[[1]], 10)
filter_AS_ex_f(AS_exon_intron_sig_no_event_l_s4[[1]], test)

#make the files with all extreme events
AS_exon_intron_sig_no_event_l_s4_v2 <- list()

for (i in seq_along(AS_exon_intron_sig_no_event_l_s4)) {
  a <- filter_AS_ex_f(AS_exon_intron_sig_no_event_l_s4[[i]], all_AS_events_bed_TPM_extremity_PI_sorted_l[[i]])
  write_delim(a, str_c("./data/temp/temp", i), delim = "\t", col_names = TRUE)
  AS_exon_intron_sig_no_event_l_s4_v2 <- append(AS_exon_intron_sig_no_event_l_s4_v2, list(a))
  
}
AS_exon_intron_sig_no_event_l_s4_v2[[1]]
#import PI extremity data (modified PI_H & PI_L) 

AS_exon_intron_sig_no_event_l_s4_v2 <- list()

for (i in seq_along(AS_exon_intron_sig_no_event_l_s4)) {
  a <- read_delim(str_c("./data/temp/temp", i), delim = "\t", col_names = TRUE)
  
  AS_exon_intron_sig_no_event_l_s4_v2 <- append(AS_exon_intron_sig_no_event_l_s4_v2, list(a))
  
}

#produce the size of each segments from modified PI_H & PI_L
AS_exon_intron_sig_no_event_size_l_s4_v2 <- 
  AS_exon_intron_sig_no_event_l_s4_v2 %>% 
  map(., feature_size_f)

for (i in seq_along(AS_exon_intron_sig_no_event_l_s4_v2)) {
  #generate the list with 19 files of feature size
  feature_size_list <- AS_exon_intron_sig_no_event_size_l_s4_v2[[i]][rep(1:nrow(AS_exon_intron_sig_no_event_size_l_s4_v2[[i]]),19),] %>%
    mutate(list_group = rep(c(1:19), each = nrow(AS_exon_intron_sig_no_event_size_l_s4_v2[[i]]))) %>%
    split(.$list_group)
  
  #generate the RPKM file for each "AS_exon_intron_sig_no_event_l_s4_v2" based on peaks
  RPKM_list <-
    pmap(list(seq(1:19), rep("_p0.05_peaks_pileup_", 19), 
              rep(names_AS_exon_intron_sig_no_event_l_s4[[i]], 19), feature_size_list), 
         histone_RPKM_f) %>%
    reduce(., full_join)
  
  #produce the output file ready for ML
  ready_ML <- AS_exon_intron_sig_no_event_l_s4_v2[[i]] %>%
    left_join(RPKM_list) %>%
    select(-order_t)
  
  write_delim(ready_ML, str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_", str_remove(names_AS_exon_intron_sig_no_event_l_s4[[i]], "segment_"), "_ext"),
              delim = "\t", col_names = TRUE)
}

#create the table of PI_H_ext and PI_L_ext separately compared with PI_M
table_ready_for_ML_raw_peak_s4 <- list()

for (i in seq_along(AS_exon_intron_sig_no_event_l_s4_v2)) {
  raw_table <- read_delim(str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_", str_remove(names_AS_exon_intron_sig_no_event_l_s4[[i]], "segment_"), "_ext"), delim = "\t", col_names = TRUE) %>%
    replace_na(list(Pol2_RPKM = 0, H3K14ac_RPKM = 0, H3K4ac_RPKM = 0, K18ac_RPKM = 0, K27ac_RPKM = 0,
                    K27me3_RPKM = 0, K4me1_RPKM = 0, K4me3_RPKM = 0, K9ac_RPKM = 0, H3K36me3_RPKM = 0,
                    H3K79ac_RPKM = 0, H4K12ac_RPKM = 0, H4K16ac_RPKM = 0, H4K20ac_RPKM = 0, H4K5ac_RPKM = 0,
                    H4K8ac_RPKM = 0, K36ac_RPKM = 0, K4me2_RPKM = 0, K56ac_RPKM = 0)) %>%
    mutate(type = as.factor(names_AS_exon_intron_sig_no_event_l_s4[[i]]))
  
  raw_table <- raw_table[c((ncol(raw_table)-19):ncol(raw_table))]
  
  if(str_detect(names_AS_exon_intron_sig_no_event_l_s4[[i]], "PI_H")==TRUE){
    path_PI_M <- str_replace(names_AS_exon_intron_sig_no_event_l_s4[[i]], "PI_H", "PI_M")
  } else {
    path_PI_M <- str_replace(names_AS_exon_intron_sig_no_event_l_s4[[i]], "PI_L", "PI_M")
  }
  
  raw_table_m <- read_delim(str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_", str_remove(path_PI_M, "segment_")), delim = "\t", col_names = TRUE) %>%
    replace_na(list(Pol2_RPKM = 0, H3K14ac_RPKM = 0, H3K4ac_RPKM = 0, K18ac_RPKM = 0, K27ac_RPKM = 0,
                    K27me3_RPKM = 0, K4me1_RPKM = 0, K4me3_RPKM = 0, K9ac_RPKM = 0, H3K36me3_RPKM = 0,
                    H3K79ac_RPKM = 0, H4K12ac_RPKM = 0, H4K16ac_RPKM = 0, H4K20ac_RPKM = 0, H4K5ac_RPKM = 0,
                    H4K8ac_RPKM = 0, K36ac_RPKM = 0, K4me2_RPKM = 0, K56ac_RPKM = 0)) %>%
    mutate(type = as.factor(path_PI_M))
  
  raw_table_m <- raw_table_m[c((ncol(raw_table_m)-19):ncol(raw_table_m))]
  
  raw_table_f <- bind_rows(raw_table, raw_table_m) %>%
    as.data.frame()
  
  table_ready_for_ML_raw_peak_s4 <- append(table_ready_for_ML_raw_peak_s4, list(raw_table_f))
}

#run ML for comparing M with PI_H_ext, PI_L_ext
table_ready_for_ML_raw_peak_s4_boruta <- vector("list", length = length(table_ready_for_ML_raw_peak_s4))

for (i in seq_along(table_ready_for_ML_raw_peak_s4_boruta)) {
  print(str_c("start boruta for the ", i, "th dataset"))
  table_ready_for_ML_raw_peak_s4_boruta[[i]] <- Boruta(type~., table_ready_for_ML_raw_peak_s4[[i]], doTrace = 2)
}

#make output file and figures
for (i in seq_along(table_ready_for_ML_raw_peak_s4_boruta)) {
  output_table <- process_the_Boruta_data(table_ready_for_ML_raw_peak_s4_boruta[[i]])
  
  names_output <- str_c(str_remove(names_AS_exon_intron_sig_no_event_l_s4[[i]], "segment_for_ML_"), "_vs_PI_M")
  
  write_delim(output_table, str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_4_peaks/", names_output), delim = "\t", col_names = TRUE)
  
  g <- process_the_Boruta_data(table_ready_for_ML_raw_peak_s4_boruta[[i]]) %>%
    plot_boruta() +
    ggtitle(names_output) 
  
  output_plot <- theme_ym(g)
  
  ggsave(str_c("./analysis/AS_all_seg_ML_s4_peaks/", names_output, "_peaks.jpeg"), output_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
  
}

#stage 4 v2: to identify whether there is a difference between PI_H_ext and PI_L_ext
#create the table of the combination of PI_H_ext and PI_L_ext (compare PI_H_ext with PI_L_ext)
table_ready_for_ML_raw_peak_s4_v2 <- list()

for (i in c(1:4, 9:12, 17:20)) {
  raw_table_H <- read_delim(str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_", str_remove(names_AS_exon_intron_sig_no_event_l_s4[[i]], "segment_"), "_ext"), delim = "\t", col_names = TRUE) %>%
    replace_na(list(Pol2_RPKM = 0, H3K14ac_RPKM = 0, H3K4ac_RPKM = 0, K18ac_RPKM = 0, K27ac_RPKM = 0,
                    K27me3_RPKM = 0, K4me1_RPKM = 0, K4me3_RPKM = 0, K9ac_RPKM = 0, H3K36me3_RPKM = 0,
                    H3K79ac_RPKM = 0, H4K12ac_RPKM = 0, H4K16ac_RPKM = 0, H4K20ac_RPKM = 0, H4K5ac_RPKM = 0,
                    H4K8ac_RPKM = 0, K36ac_RPKM = 0, K4me2_RPKM = 0, K56ac_RPKM = 0)) %>%
    mutate(type = as.factor(names_AS_exon_intron_sig_no_event_l_s4[[i]]))
  
  raw_table_L <- read_delim(str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_", str_remove(names_AS_exon_intron_sig_no_event_l_s4[[i+4]], "segment_"), "_ext"), delim = "\t", col_names = TRUE) %>%
    replace_na(list(Pol2_RPKM = 0, H3K14ac_RPKM = 0, H3K4ac_RPKM = 0, K18ac_RPKM = 0, K27ac_RPKM = 0,
                    K27me3_RPKM = 0, K4me1_RPKM = 0, K4me3_RPKM = 0, K9ac_RPKM = 0, H3K36me3_RPKM = 0,
                    H3K79ac_RPKM = 0, H4K12ac_RPKM = 0, H4K16ac_RPKM = 0, H4K20ac_RPKM = 0, H4K5ac_RPKM = 0,
                    H4K8ac_RPKM = 0, K36ac_RPKM = 0, K4me2_RPKM = 0, K56ac_RPKM = 0)) %>%
    mutate(type = as.factor(names_AS_exon_intron_sig_no_event_l_s4[[i+4]]))

  raw_table_f <- bind_rows(raw_table_H, raw_table_L) %>%
    as.data.frame()
  
  raw_table_f <- raw_table_f[c((ncol(raw_table_f)-19):ncol(raw_table_f))]
  
  table_ready_for_ML_raw_peak_s4_v2 <- append(table_ready_for_ML_raw_peak_s4_v2, list(raw_table_f))
}

#stage 4 version 2 (run ML)

table_ready_for_ML_raw_peak_s4_boruta_v2 <- vector("list", length = length(table_ready_for_ML_raw_peak_s4_v2))

for (i in seq_along(table_ready_for_ML_raw_peak_s4_boruta_v2)) {
  print(str_c("start boruta for the ", i, "th dataset"))
  table_ready_for_ML_raw_peak_s4_boruta_v2[[i]] <- Boruta(type~., table_ready_for_ML_raw_peak_s4_v2[[i]], doTrace = 2)
}

#produce output
for (i in seq_along(table_ready_for_ML_raw_peak_s4_boruta_v2)) {
  output_table <- process_the_Boruta_data(table_ready_for_ML_raw_peak_s4_boruta_v2[[i]])
  
  names_output <- str_c(str_remove(unique(table_ready_for_ML_raw_peak_s4_v2[[i]]$type)[1], "segment_for_ML_"), "_vs_PI_L_ext")
  
  write_delim(output_table, str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_4_peaks/", names_output), delim = "\t", col_names = TRUE)
  
  g <- process_the_Boruta_data(table_ready_for_ML_raw_peak_s4_boruta_v2[[i]]) %>%
    plot_boruta() +
    ggtitle(names_output) 
  
  output_plot <- theme_ym(g)
  
  ggsave(str_c("./analysis/AS_all_seg_ML_s4_peaks/", names_output, "_peaks.jpeg"), output_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
  
}

#data check for missing in our data
#if there are missing data, one can check here (https://www.datacamp.com/tutorial/feature-selection-R-boruta)
for (i in c(1:4, 9:12, 17:20)) {
  test <- AS_exon_intron_sig_no_event_l[[i]] %>%
    nrow()
  test_na <- AS_exon_intron_sig_no_event_l[[i]] %>%
    drop_na() %>%
    nrow()
  print(str_c(test, ", ", test_na))
}

AS_exon_intron_sig_no_event_l[[1]]


names_AS_exon_intron_sig_no_event_l

?append

make_AS_ML_data_f(test_exon_data_for_ML_temp_l[[2]], c("HS1_HS0"), c("PI_H"))

test_exon_data_for_ML_temp_l_small <- test_exon_data_for_ML_temp_l[1:10]

test_exon_data_for_ML_temp_l_small %>%
  map(. %>% make_AS_ML_data_f(., c("HS1_HS0"), c("PI_H"))) %>%
  bind_rows()

#stage 5 external ctrl
#the logic here is to simulate previous ML using shuffled PSI
#So I made 5 subgroups of external ctrl in the part of #gene without AS
#Then, I can make ten pairwise comparisons

#stage 5 v1, external ctrl vs external ctrl
#produce ML data of external ctrl
names_exon_intron_data_for_ML_f_ex_ctrl <- str_c("segment_for_ML_sh", seq(1:5), "_ex_ctrl")

#add the label of three and five prime side of each external ctrl segment
exon_intron_data_for_ML_f_ex_ctrl_l <- exon_intron_data_for_ML_f_ex_ctrl %>%
  map(. %>% mutate(order_t = c(1:n()))) 

exon_intron_data_for_ML_f_ex_ctrl_size_l <- 
  exon_intron_data_for_ML_f_ex_ctrl_l %>% 
  map(., feature_size_f)

#import and calculate the intersection between external ctrl and each histone modification marks
#and merge them into one table for each external ctrl subgroup sample
for (i in seq_along(exon_intron_data_for_ML_f_ex_ctrl_l)) {
  #generate the list with 19 files of feature size
  feature_size_list <- exon_intron_data_for_ML_f_ex_ctrl_size_l[[i]][rep(1:nrow(exon_intron_data_for_ML_f_ex_ctrl_size_l[[i]]),19),] %>%
    mutate(list_group = rep(c(1:19), each = nrow(exon_intron_data_for_ML_f_ex_ctrl_size_l[[i]]))) %>%
    split(.$list_group)
  
  #generate the RPKM file for each "exon_intron_data_for_ML_f_ex_ctrl" based on peaks
  RPKM_list <-
    pmap(list(seq(1:19), rep("_p0.05_peaks_pileup_", 19), 
              rep(names_exon_intron_data_for_ML_f_ex_ctrl[[i]], 19), feature_size_list), 
         histone_RPKM_f) %>%
    reduce(., full_join)
  
  #produce the output file ready for ML
  ready_ML <- exon_intron_data_for_ML_f_ex_ctrl_l[[i]] %>%
    left_join(RPKM_list) %>%
    select(-order_t)
  
  write_delim(ready_ML, str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_", str_remove(names_exon_intron_data_for_ML_f_ex_ctrl[[i]], "segment_")),
              delim = "\t", col_names = TRUE)
}

#import ML ready data and merge them into the form for pairwise comparison
table_ready_for_ML_stage_5_ex_ctrl_raw <- list()
name_table_ready_for_ML_stage_5_ex_ctrl_raw <- c()

for (i in seq_along(names_exon_intron_data_for_ML_f_ex_ctrl)) {
  path_RPKM_AS_all_seg <- str_c("./data/AS_bed_for_ML/all_segments/table_ready_p005_peaks_pileup_", str_remove(names_exon_intron_data_for_ML_f_ex_ctrl[[i]], "segment_"))
  
  a <- read_delim(path_RPKM_AS_all_seg, delim = "\t") %>%
    replace_na(list(Pol2_RPKM = 0, H3K14ac_RPKM = 0, H3K4ac_RPKM = 0, K18ac_RPKM = 0, K27ac_RPKM = 0,
                    K27me3_RPKM = 0, K4me1_RPKM = 0, K4me3_RPKM = 0, K9ac_RPKM = 0, H3K36me3_RPKM = 0,
                    H3K79ac_RPKM = 0, H4K12ac_RPKM = 0, H4K16ac_RPKM = 0, H4K20ac_RPKM = 0, H4K5ac_RPKM = 0,
                    H4K8ac_RPKM = 0, K36ac_RPKM = 0, K4me2_RPKM = 0, K56ac_RPKM = 0)) %>%
    mutate(type = as.factor(str_c(pair, "_", PI)))
  
  
  raw_table <- a[c(ncol(a)-19):ncol(a)]
  
  table_ready_for_ML_stage_5_ex_ctrl_raw <- append(table_ready_for_ML_stage_5_ex_ctrl_raw, list(raw_table))
  name_table_ready_for_ML_stage_5_ex_ctrl_raw <- append(name_table_ready_for_ML_stage_5_ex_ctrl_raw, raw_table$type[[1]])
}

comb_pair_AS_seg_type <- combs(c(1:5), 2) #five sh external type

table_ready_for_ML_stage_5_ex_ctrl <- list()
name_table_ready_for_ML_stage_5_ex_ctrl <- list()

for (i in seq_along(comb_pair_AS_seg_type[,1])) {
  
  a <- table_ready_for_ML_stage_5_ex_ctrl_raw[comb_pair_AS_seg_type[i,]] %>%
    bind_rows()
  
  type_check <- a %>%
    select(type) %>%
    unique()
  
  print(type_check)
  
  a <- a %>%
    as.data.frame()
  
  table_ready_for_ML_stage_5_ex_ctrl <- append(table_ready_for_ML_stage_5_ex_ctrl, list(a))
  name_table_ready_for_ML_stage_5_ex_ctrl <- append(name_table_ready_for_ML_stage_5_ex_ctrl, list(type_check))
}

#run boruta
table_ready_for_ML_stage_5_ex_ctrl_boruta <- vector("list", length = length(table_ready_for_ML_stage_5_ex_ctrl))

for (i in seq_along(table_ready_for_ML_stage_5_ex_ctrl_boruta)) {
  print(str_c("start boruta for the ", i, "th dataset"))
  table_ready_for_ML_stage_5_ex_ctrl_boruta[[i]] <- Boruta(type~., table_ready_for_ML_stage_5_ex_ctrl[[i]], doTrace = 2)
}

#since the result showed the these 5 external ctrl subgroups are not the same
#I combined them together as the whole external ctrl, and compared it with df PSI target
table_ready_for_ML_stage_5_ex_ctrl_combined <- bind_rows(table_ready_for_ML_stage_5_ex_ctrl_raw) %>%
  mutate(type = as.factor("external_ctrl"))


table_ready_for_ML_stage_5_raw_target <- list()
name_pairwise_table_for_ML_stage_5_raw_target <- c()

for (i in seq_along(AS_type)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(PI_type)) {
      #combine df PSI and external ctrl
      a <- bind_rows(table_ready_for_ML_raw_peak_target_sorted[[i]][[j]][[k]], 
                     table_ready_for_ML_stage_5_ex_ctrl_combined)
      
      a <- a %>%
        as.data.frame()
      
      #produce the name for the following output and ggplot title
      name_a <- str_c(table_ready_for_ML_raw_peak_target_sorted[[i]][[j]][[k]]$type[[1]], "_vs_", "external_ctrl")
      
      #append data with df and external ctrl into a list
      table_ready_for_ML_stage_5_raw_target <- append(table_ready_for_ML_stage_5_raw_target, list(a))
      
      #append names into a vector
      name_pairwise_table_for_ML_stage_5_raw_target <- append(name_pairwise_table_for_ML_stage_5_raw_target, name_a)
    }
  }
}

#run boruta for df PSI and external ctrl
table_ready_for_ML_stage_5_PSI_vs_ex_ctrl_boruta <- vector("list", length = length(table_ready_for_ML_stage_5_raw_target))

for (i in seq_along(table_ready_for_ML_stage_5_PSI_vs_ex_ctrl_boruta)) {
  print(str_c("start boruta for the ", i, "th dataset", "for ", name_pairwise_table_for_ML_stage_5_raw_target[[i]]))
  table_ready_for_ML_stage_5_PSI_vs_ex_ctrl_boruta[[i]] <- Boruta(type~., table_ready_for_ML_stage_5_raw_target[[i]], doTrace = 2)
}

#write the boruta output file and plots
for (i in seq_along(name_pairwise_table_for_ML_stage_5_raw_target)) {
  output_table <- process_the_Boruta_data(table_ready_for_ML_stage_5_PSI_vs_ex_ctrl_boruta[[i]])
  write_delim(output_table, str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_5_peaks/", name_pairwise_table_for_ML_stage_5_raw_target[[i]]), delim = "\t", col_names = TRUE)
  
  g <- process_the_Boruta_data(table_ready_for_ML_stage_5_PSI_vs_ex_ctrl_boruta[[i]]) %>%
    plot_boruta() +
    ggtitle(name_pairwise_table_for_ML_stage_5_raw_target[[i]])
  
  output_plot <- theme_ym(g)
  
  ggsave(str_c("./analysis/AS_all_seg_ML_s5_peaks/", name_pairwise_table_for_ML_stage_5_raw_target[[i]], ".jpeg"), output_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
  
}

table_ready_for_ML_stage_5_PSI_vs_ex_ctrl_boruta %>%
  map(. %>% process_the_Boruta_data()) %>%
  bind_rows() %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(q = mean(value)) %>%
  arrange(q) %>%
  View()

#stage 6 target vs internal ctrl
table_ready_for_ML_stage_6_raw_peak_target <- list()
name_pairwise_table_for_ML_stage_6_raw_peak_target <- list()
table_ready_for_ML_raw_peak_target_sorted[[1]][[1]][[2]]$type[[1]]
table_ready_for_ML_raw_peak_ctrl_sorted[[1]][[1]][[2]]$type[[1]]
PI_type

for (i in seq_along(AS_type)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(PI_type)) {
      a <- table_ready_for_ML_raw_peak_target_sorted[[i]][[j]][[k]]
      b <- table_ready_for_ML_raw_peak_ctrl_sorted[[i]][[j]][[k]]
      
      c <- bind_rows(a, b)
      
      type_check <- c %>%
        select(type) %>%
        unique()
      print(type_check)
      c <- c %>%
        as.data.frame()
      table_ready_for_ML_stage_6_raw_peak_target <- append(table_ready_for_ML_stage_6_raw_peak_target, list(c))
      name_pairwise_table_for_ML_stage_6_raw_peak_target <- append(name_pairwise_table_for_ML_stage_6_raw_peak_target, list(type_check))
    }
  }
}





unique(table_ready_for_ML_stage_1_raw_peak_ctrl[[1]]$type)

table_ready_for_ML_stage_6_raw_peak_boruta <- vector("list", length = length(table_ready_for_ML_stage_6_raw_peak_target))

for (i in seq_along(table_ready_for_ML_stage_6_raw_peak_target)) {
  print(str_c("start boruta for the ", i, "th dataset"))
  table_ready_for_ML_stage_6_raw_peak_boruta[[i]] <- Boruta(type~., table_ready_for_ML_stage_6_raw_peak_target[[i]], doTrace = 2)
}

#write the boruta output file and plots
for (i in seq_along(name_pairwise_table_for_ML_stage_5_raw_target)) {
  output_table <- process_the_Boruta_data(table_ready_for_ML_stage_6_raw_peak_boruta[[i]])
  write_delim(output_table, str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_6_peaks/", name_pairwise_table_for_ML_stage_6_raw_peak_target[[i]]$type[[1]], "_vs_internal_ctrl"), delim = "\t", col_names = TRUE)
  
  g <- process_the_Boruta_data(table_ready_for_ML_stage_6_raw_peak_boruta[[i]]) %>%
    plot_boruta() +
    ggtitle(str_c(name_pairwise_table_for_ML_stage_6_raw_peak_target[[i]]$type[[1]], "_vs_internal_ctrl"))
  
  output_plot <- theme_ym(g)
  
  ggsave(str_c("./analysis/AS_all_seg_ML_s6_peaks/", name_pairwise_table_for_ML_stage_6_raw_peak_target[[i]]$type[[1]], "_vs_internal_ctrl.jpeg"), output_plot, width = 400, height = 240, units = c("mm"), dpi = 320)
  
}

table_ready_for_ML_stage_6_raw_peak_boruta %>%
  map(. %>% process_the_Boruta_data) %>%
  bind_rows() %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(mean = mean(value)) %>%
  arrange(-mean)


#This is to check the median size of genes containing AS sig or non sig
read_delim("./data/Intersected_AS_TPM_q05_M82_anno_rMATs/AS_control_TPM_q05_HS1_HS0_SE_intersected_M82_rMATs_exon_order.bed", delim = "\t", 
           col_names = c("chr", "str", "end", "comp", "AS", "PI", "chr_2", "str_2", "end_2", "strand", "feature", "source", "anno", "order")) %>%
  left_join(M82_rMATs_anno_all_gene_light) %>%
  #filter(order > 2) %>%
  group_by(PI) %>%
  summarise(gene_size_median = median(gene_size)) %>%
  View()
