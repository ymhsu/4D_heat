install.packages("pacman")
library(pacman)
Packages <- c("tidyverse", "doParallel", "foreach", "gridExtra", "ggeffects", "cowplot", "caTools", "grid", "ggplotify")
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

#For RI, A5S and A3S, I combined them into intron, and for SE, it is belong to the file of intersected exons
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

all_AS_intersected_feature_raw %>%
  filter(anno == "Solyc03g019850.3.1") %>%
  group_by(feature, event) %>%
  summarise(count = n())

all_AS_intersected_feature_raw %>%
  filter(anno == "Solyc03g019850.3.1") %>%
  View()

#combine all genes intersected with AS (four types)  
gene_having_AS <- all_AS_intersected_feature_raw %>%
  select(feature, anno, event) %>%
  arrange(anno) %>%
  distinct()

#Import all exon and intron with the order into
M82_rMATs_anno_all_exon_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_exon_order.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order"))
M82_rMATs_anno_all_intron_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_intron_order.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order"))

#combine all exon and intron as a list
gene_having_AS_info_l <- list(M82_rMATs_anno_all_exon_order, M82_rMATs_anno_all_intron_order)
names(gene_having_AS_info_l) <- c("exon", "intron")

#In this step, we keep all information of genes if they are intersected with any kind of AS (sig, non-sig)
gene_having_AS_info_l_f <- gene_having_AS_info_l %>%
  map(. %>% arrange(chr, str)) %>%
  map(. %>% left_join(gene_having_AS)) %>%
  map(. %>% drop_na()) %>%
  map(. %>% select(-event)) %>%
  map(. %>% distinct())

#have a quick look how many exons do genes with SE have in general 
gene_having_AS_info_l_f$exon %>%
  group_by(anno) %>%
  summarise(count = n()) %>%
  group_by(count) %>%
  summarise(count_sum = sum(count)) %>%
  View()

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


#create the function for random choosing the same amount exon/intron as that of AS events occurred in that gene
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
registerDoParallel(cores = 10)
getDoParWorkers()

#create all control set from all genes with AS event whehter these events are sig or non-sig
exon_intron_data_for_ML_list_v2 <- 
  foreach(i=c(1:5), .packages = c("tidyverse")) %:%
  foreach(j=1:12, .packages = c("tidyverse")) %dopar% {
    exon_intron_data_for_ML_list %>%
      map(. %>% map(. %>% make_AS_ML_data_f(., comb_comp_PI_AS_list[[i]][[j]]))) %>%
      map(. %>% bind_rows())
    
  }

#modify the content in the columns of "comp", "AS" and "PI" for randomed selected events (they are noted as "no" previously)
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
#whether these events are located in exon or intron are not an importnat issue
#one AS event can be located in an exon of a gene but an intron of the isoform of the same genes
exon_intron_data_for_ML_list_v4 <- exon_intron_data_for_ML_list_v3 %>%
  map(. %>% map(. %>% bind_rows())) %>%
  map(. %>% bind_rows()) %>%
  bind_rows() %>%
  split(.$comp) %>%
  map(. %>% split(.$PI)) %>%
  map(. %>% map(. %>% split(.$AS)))

#produce two sides (100 bp extension) of each exon or intron with or without AS for the following analysis 
#(the intersection between these events and histone marks)
for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    for (k in seq_along(AS_type)) {
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
      
      exon_intron_data_for_ML_f_ctrl <- exon_intron_data_for_ML_f_combined %>%
        filter(event == "no")
      
      exon_intron_data_for_ML_f_target <- exon_intron_data_for_ML_f_combined %>%
        filter(event != "no")
      
      label_comp <- unique(exon_intron_data_for_ML_list_v4[[i]][[j]][[k]]$comp)
      label_PI <- unique(exon_intron_data_for_ML_list_v4[[i]][[j]][[k]]$PI)
      label_AS <- unique(exon_intron_data_for_ML_list_v4[[i]][[j]][[k]]$AS)
      
      write_delim(exon_intron_data_for_ML_f_target, str_c("./data/AS_bed_for_ML/segment_for_ML_", label_comp, "_", label_PI, "_", label_AS, "_", "target.bed"), col_names = FALSE, delim = "\t")
      write_delim(exon_intron_data_for_ML_f_ctrl, str_c("./data/AS_bed_for_ML/segment_for_ML_", label_comp, "_", label_PI, "_", label_AS, "_", "ctrl.bed"), col_names = FALSE, delim = "\t")
      
    }
  }
}


#import bedgraph of exon/intron with sig AS or no event to calculate mark signals
histone_mark_list <- read_delim("./data/Chip_seq/no_treatment/bam/mark_list", delim = "\t", col_names = c("name_raw"))

histone_mark_read_count <- c()

for (i in seq_along(histone_mark_list$name_raw)) {
  read_count_file <- read_delim(str_c("./data/Chip_seq/no_treatment/bam/", histone_mark_list$name_raw[[i]], "_read_count"), delim = "\t", col_names = c("read_count"))
  histone_mark_read_count <- append(histone_mark_read_count, read_count_file$read_count)
}

histone_mark_list_m <- histone_mark_list %>%
  mutate(read_count = histone_mark_read_count) %>%
  mutate(name_light = str_remove_all(name_raw, "0h_")) %>%
  mutate(name_light = str_remove_all(name_light, "_.*"))

#segment_for_ML_HS6_HS0_PI_MH_A5S_ctrl
comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
PI_type = c("PI_H", "PI_MH", "PI_M", "PI_ML", "PI_L")
AS_type = c("A5S", "A3S", "SE", "RI")
event = c("target", "ctrl")


AS_exon_intron_sig_no_event_l <- list()
names_AS_exon_intron_sig_no_event_l <- c()

for (i in seq_along(comp)) {
  for (j in seq_along(PI_type)) {
    for (k in seq_along(AS_type)) {
      for (l in seq_along(event)) {
        
        AS_table <- read_delim(str_c("./data/AS_bed_for_ML/segment_for_ML_", comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", event[[l]], ".bed"), delim = "\t",
                                                                   col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order", "comp", "AS", "PI", "event", "pair", "seg_side"))
        names_AS_exon_intron_sig_no_event_l <- append(names_AS_exon_intron_sig_no_event_l, str_c("segment_for_ML_", comp[[i]], "_", PI_type[[j]], "_", AS_type[[k]], "_", event[[l]]))
        AS_exon_intron_sig_no_event_l <- append(AS_exon_intron_sig_no_event_l, list(AS_table))
      }
    }
  }
}

length(AS_exon_intron_sig_no_event_l)

names(AS_exon_intron_sig_no_event_l)
str_c("./data/AS_bed_for_ML/", histone_mark_list_m$name_raw[[1]], "_", names_AS_exon_intron_sig_no_event_l[[1]], "_12chr.bed")

read_delim
<<<<<<< HEAD
AS_exon_intron_sig_no_event_l[[1]]
=======

#produce data for ML
#make a list for operating 10 cores
#I will use 10 cores, so here I made the length of list as 12 (120/10 = 12)

AS_exon_intron_sig_no_event_l_10cores <- vector("list", length = length(AS_exon_intron_sig_no_event_l)/10)
names_AS_exon_intron_sig_no_event_l_10cores <- vector("list", length = length(AS_exon_intron_sig_no_event_l)/10)


for (i in seq_along(AS_exon_intron_sig_no_event_l_10cores)) {
  print(seq(10*(i-1)+1, 10*(i-1)+10, 1))
  #comb_comp_PI_AS_list[[i]] <- comb_comp_PI_AS_list_raw[seq(12*(i-1)+1, 12*(i-1)+12, 1)]
  AS_exon_intron_sig_no_event_l_10cores[[i]] <- AS_exon_intron_sig_no_event_l[seq(10*(i-1)+1, 10*(i-1)+10, 1)]
  names_AS_exon_intron_sig_no_event_l_10cores[[i]] <- names_AS_exon_intron_sig_no_event_l[seq(10*(i-1)+1, 10*(i-1)+10, 1)]
}

names_AS_exon_intron_sig_no_event_l_10cores[[1]][[1]]

AS_all_seg_his_signal_raw <- vector("list", length = length(AS_exon_intron_sig_no_event_l))
length(AS_all_seg_his_signal_raw)

>>>>>>> bb73b9b34b4fa9a023e91d4f83d0bb1ffc9a7318
for (i in seq_along(histone_mark_list_m$name_raw)) {
  for (j in seq_along(AS_exon_intron_sig_no_event_l)) {
    feature_size <- AS_exon_intron_sig_no_event_l[[j]] %>%
      group_by(anno, feature, order) %>%
      mutate(str_n = min(end), end_n = max(str)) %>%
      mutate(feature_size = end_n - str_n) %>%
      ungroup() %>%
      select(anno, feature, order, feature_size) %>%
      distinct()
    
    test <-
      read_delim(
        str_c(
          "./data/AS_bed_for_ML/all_segments/",
          histone_mark_list_m$name_raw[[i]],
          "_",
          str_remove(names_AS_exon_intron_sig_no_event_l[[j]], "segment_"),
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
      mutate(read_count = histone_mark_list_m$read_count[[i]]) %>%
      mutate(RPM = mean_signal / read_count * 10 ^ 6) %>%
      left_join(feature_size) %>%
      mutate(RPKM = RPM / feature_size) %>%
      select(-mean_signal, -read_count, -RPM, -feature_size)
    
    names(test) <-
      c(
        "anno",
        "feature",
        "order",
        "comp",
        "AS",
        "PI",
        "event",
        "seg_side",
        str_c(histone_mark_list_m$name_light[[i]], "_RPKM")
      )
    if(i!=19){
    AS_exon_intron_sig_no_event_l[[j]] <- AS_exon_intron_sig_no_event_l[[j]] %>%
      left_join(test)
    } else {
      AS_exon_intron_sig_no_event_l[[j]] <- AS_exon_intron_sig_no_event_l[[j]] %>%
        left_join(test)
      write_delim(AS_exon_intron_sig_no_event_l[[j]], str_c("./data/AS_bed_for_ML/all_segments/table_ready_", str_remove(names_AS_exon_intron_sig_no_event_l[[j]], "segment_")),
                  delim = "\t", col_names = TRUE)
    }
  }
}




for (i in seq_along(AS_exon_intron_sig_no_event_l)) {
  if(i==1){
    print("test")
  } else {
    print("testv2")
  }
}


  foreach(i=c(1:19), .packages = c("tidyverse")) %:%
  foreach(j=c(1:12), .packages = c("tidyverse")) %:%
  foreach(k=1:10, .packages = c("tidyverse")) %dopar% {
    
    feature_size <- AS_exon_intron_sig_no_event_l_10cores[[j]][[k]] %>%
      group_by(anno, feature, order) %>%
      mutate(str_n = min(end), end_n = max(str)) %>%
      mutate(feature_size = end_n - str_n) %>%
      ungroup() %>%
      select(anno, feature, order, feature_size) %>%
      distinct()
    
    
  }

AS_all_seg_his_signal_raw_v2 <- vector("list", length = length(AS_exon_intron_sig_no_event_l))

for (i in seq_along(histone_mark_list_m$name_raw)) {
  for (j in seq_along(AS_exon_intron_sig_no_event_l_10cores)) {
    for (k in seq_along(AS_exon_intron_sig_no_event_l_10cores[[1]])) {
      print(c(10*(j-1)+k))
      AS_all_seg_his_signal_raw_v2[[c(10*(j-1)+k)]] <- append(AS_all_seg_his_signal_raw_v2[[c(10*(j-1)+k)]], list(AS_all_seg_his_signal_raw[[i]][[j]][[k]]))
      AS_all_seg_his_signal_raw_v2[[c(10*(j-1)+k)]] %>%
        reduce(full_join)
    }
  }
}


#write the data to the output (RPKM of each feature are ready)
path_AS_exon_intron_sig_no_event_l <- str_c("./data/AS_bed_for_ML/ML_data_each_AS_comp_PI_sig/", names_AS_exon_intron_sig_no_event_l, "_data_ready_raw")

pwalk(list(
  AS_exon_intron_sig_no_event_l,
  path_AS_exon_intron_sig_no_event_l,
  delim = "\t",
  col_names = FALSE
),
write_delim)

#import RPKM raw data

for (i in seq_along(path_AS_exon_intron_sig_no_event_l)) {
  AS_exon_intron_sig_no_event_l[[i]] <- read_delim(path_AS_exon_intron_sig_no_event_l[[i]], 
                                                   delim = "\t", col_names = c(colnames(AS_exon_intron_sig_no_event_l[[1]])[1:13], str_c(histone_mark_list_m$name_light, "_RPKM")))
}

colnames(AS_exon_intron_sig_no_event_l[[1]])[1:13]
ncol(read_delim(path_AS_exon_intron_sig_no_event_l[[1]], delim = "\t", col_names = c(colnames(AS_exon_intron_sig_no_event_l[[1]]), str_c(histone_mark_list_m$name_light, "_RPKM"))))

AS_exon_intron_sig_no_event_l_v2 <- AS_exon_intron_sig_no_event_l %>%
  map(. %>% mutate(type = as.factor(str_c(feature, "_", comp, "_", PI, "_", event)))) %>%
  map(. %>% select(-chr, -str, -end, -strand, -feature, -source, -anno, -order, -comp, -AS, -PI, -event, -seg_side))

for (i in seq_along(AS_exon_intron_sig_no_event_l_v2)) {
  AS_exon_intron_sig_no_event_l_v2[[i]] <- AS_exon_intron_sig_no_event_l_v2[[i]] %>%
    mutate(type = as.factor(names_AS_exon_intron_sig_no_event_l[[i]]))
}

pair_1_2 <- vector("list", length = length(names_AS_exon_intron_sig_no_event_l)/2)
for (i in seq_along(pair_1_2)) {
  pair_1_2[[i]] <- as.data.frame(bind_rows(AS_exon_intron_sig_no_event_l_v2[[(i-1)*2+1]], AS_exon_intron_sig_no_event_l_v2[[(i-1)*2+2]]))

}

pair_2_4 <- vector("list", length = length(names_AS_exon_intron_sig_no_event_l)/4)

for (i in seq_along(pair_2_4)) {
    pair_2_4[[i]] <- as.data.frame(bind_rows(AS_exon_intron_sig_no_event_l_v2[[(i-1)*4+2]], AS_exon_intron_sig_no_event_l_v2[[(i-1)*4 + 4]]))
  
}

pair_1_2_3_4 <- vector("list", length = length(names_AS_exon_intron_sig_no_event_l)/4)
for (i in seq_along(pair_1_2_3_4)) {
  pair_1_2_3_4[[i]] <- as.data.frame(bind_rows(AS_exon_intron_sig_no_event_l_v2[[(i-1)*4+1]], AS_exon_intron_sig_no_event_l_v2[[(i-1)*4+2]],
                                               AS_exon_intron_sig_no_event_l_v2[[(i-1)*4+3]], AS_exon_intron_sig_no_event_l_v2[[(i-1)*4+4]]))
  
}


ML_l_f[[1]]
ML_l_f <- list()
ML_l_f_reshuffled <- list()
ML_l_f <- append(ML_l_f, pair_1_2)
ML_l_f <- append(ML_l_f, pair_2_4)
ML_l_f <- append(ML_l_f, pair_1_2_3_4)

path_ML_l_f_temp_l <- str_c("./data/ML_running_temp/ML_l_f_temp_", seq(1:24))

pwalk(list(
  ML_l_f,
  path_ML_l_f_temp_l,
  delim = "\t",
  col_names = TRUE
),
write_delim)

for (i in seq_along(c(1:24))) {
  ML_l_f[[i]] <- read_delim(path_ML_l_f_temp_l[[i]], delim = "\t", col_names = TRUE) %>%
    mutate(type = as.factor(type))
  ML_l_f_reshuffled[[i]] <- ML_l_f[[i]] %>%
    mutate(type = sample(type))
  ML_l_f[[i]] <- as.data.frame(ML_l_f[[i]])
  ML_l_f_reshuffled[[i]] <- as.data.frame(ML_l_f_reshuffled[[i]])
}

#test the first data for ML
ML_l_result[[1]] <- Boruta(type~.,data=ML_l_f[[1]],doTrace = 2)
?plot.Boruta
ML_l_result[[1]]

ML_l_resh_result <- Boruta(type~.,data=ML_l_f_reshuffled[[1]],doTrace = 2)
str_
ML_l_f_1_combined <- bind_rows(as_tibble(ML_l_f[[1]]), as_tibble(ML_l_f[[2]])) %>%
  mutate(type = str_replace(type, "PI_H", "major")) %>%
  mutate(type = str_replace(type, "PI_L", "major")) %>%
  mutate(type = as.factor(type)) %>%
  as.data.frame()
  
ML_l_f_1_combined_res <- Boruta(type~.,data=ML_l_f_1_combined,doTrace = 2)

plot(ML_l_resh_result)

as.grob(densityplot(~mpg|cyl, data=mtcars))
plot.Boruta(ML_l_result[[1]])

test_plot <- as.ggplot(~plot(ML_l_result[[1]], las = 2))

test_plot + 
  theme(plot.background = element_rect(color = 1,
                                       size = 1),
        plot.margin = margin(t = 1,  # Top margin
                             r = 1,  # Right margin
                             b = 3,  # Bottom margin
                             l = 2,  # Left margin
                             unit = "cm")) 


as_tibble(attStats(ML_l_result[[3]])) %>%
  mutate(mark = row.names(attStats(ML_l_result[[3]]))) %>%
  arrange(medianImp)

ML_l_f[[3]]
attStats(ML_l_result[[3]])

split <- sample.split(ML_l_f[[3]],SplitRatio = 0.75)
train <- subset(ML_l_f[[3]],split==T)
test <- subset(ML_l_f[[3]],split==F)

rfmodel <- randomForest(type ~., data = train)
pred_full <- predict(rfmodel, test)
confusionMatrix(table(pred_full, test$type))

rfmodel_without_K4me1 <- randomForest(type ~ H3K4ac_RPKM, data = train)
rfmodel_without_K4me1

K4me1_RPKM

Pol2_RPKM + H3K14ac_RPKM + H3K4ac_RPKM + K18ac_RPKM + K27ac_RPKM + K27me3_RPKM + K4me3_RPKM +
  K9ac_RPKM + H3K36me3_RPKM + H3K79ac_RPKM + H4K12ac_RPKM + H4K16ac_RPKM + H4K20ac_RPKM + H4K5ac_RPKM + H4K8ac_RPKM +
  K36ac_RPKM + K4me2_RPKM + K56ac_RPKM

pred_full_without_K4me1 <- predict(rfmodel_without_K4me1, test)
confusionMatrix(table(pred_full_without_K4me1, test$type))

colnames(train)

ML_l_result <- vector("list", length = length(c(1:12)))

for (i in seq_along(c(1:12))) {
  ML_l_result[[i]] <- Boruta(type~.,data=ML_l_f[[i]],doTrace = 2)
}

test_plot <- as.ggplot(function()plot(ML_l_result[[1]], las=2))
par(mar = c(12, 5, 4, 2)+ 0.1)
plot(ML_l_result[[1]], las=2)

plot(ML_l_result[[7]], las=2)
attStats(ML_l_result[[2]]) %>%
  

as_tibble(ML_l_f[[7]]) %>%
  #group_by(type) %>%
  #summarise(count = n())
  gather(key = "feature", value = "signal", c(1:19)) %>%
  group_by(type, feature) %>%
  #mutate(average = mean(signal)) %>%
  summarise(average = mean(signal)) %>%
  #View()
  ggplot() + geom_point(aes(x=type, y=average, color = type)) +
  #geom_jitter(aes(x=type, y=log(signal)), position = position_jitter(width = 0.4)) +
  #geom_line(aes(x=type,y=log(average),color=type)) +
  facet_wrap(~feature)

as_tibble(ML_l_f[[2]]) %>%
  #group_by(type) %>%
  #summarise(count = n())
  gather(key = "feature", value = "signal", c(1:19)) %>%  
  #filter(feature == "K4me1_RPKM") %>%
  ggplot() +
  #geom_violin(aes(type, log(signal))) +
  geom_histogram(aes(log(signal), fill = type), position = "identity") +
  facet_wrap(~feature, scales="free_x")

?geom_violin

names_AS_exon_intron_sig_no_event_l

as_tibble(attStats(ML_l_result[[1]])) %>%
  mutate(feature = row.names(attStats(ML_l_result[[1]]))) %>%
  arrange(-meanImp)

ML_l_f_test <- as_tibble(ML_l_f[[1]]) %>%
  mutate(type_n = type) %>%
  mutate(type = runif(nrow(ML_l_f[[1]]))) %>%
  as.data.frame()
  

print(ML_l_result[[24]])

plotImpHistory(ML_l_result[[24]])

length(ML_l_f)
as_tibble(ML_l_f[[24]]) %>%
  select(type) %>%
  distinct()
names_AS_exon_intron_sig_no_event_l

test_ML_ex_in_PI_H <- as.data.frame(bind_rows(AS_exon_intron_sig_no_event_l_v2[[1]], AS_exon_intron_sig_no_event_l_v2[[13]]))

ML_l_f_test_f <- Boruta(type_n~.,data=ML_l_f_test,doTrace = 2)
plot(ML_l_f_test_f)

#data check for missing in our data
#if there are missing data, one can check here (https://www.datacamp.com/tutorial/feature-selection-R-boruta)
for (i in seq_along(AS_exon_intron_sig_no_event_l)) {
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

#This is to check the median size of genes containing AS sig or non sig
read_delim("./data/Intersected_AS_TPM_q05_M82_anno_rMATs/AS_control_TPM_q05_HS1_HS0_SE_intersected_M82_rMATs_exon_order.bed", delim = "\t", 
           col_names = c("chr", "str", "end", "comp", "AS", "PI", "chr_2", "str_2", "end_2", "strand", "feature", "source", "anno", "order")) %>%
  left_join(M82_rMATs_anno_all_gene_light) %>%
  #filter(order > 2) %>%
  group_by(PI) %>%
  summarise(gene_size_median = median(gene_size)) %>%
  View()
