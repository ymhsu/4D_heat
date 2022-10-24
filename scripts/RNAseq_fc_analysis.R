install.packages("pacman")
library(pacman)
Packages <- c("tidyverse", "fc", "vroom", "gt", "plotly", "modelr")
p_load(Packages, character.only = TRUE)


#depending on all types of AS to perform the analysis similar to the workflow in the file of "The intersection_RI_marks.R"
#Analysis will be based on DAS and AS
#DAS (Differentially alternative splicing) was to check the trend of marks overlapping these DAS
#AS was to generate control set (meaning they don't have AS) for producing plots based on deeptools

all_AS_events_temp <- read_delim("./data/rMATS_out/All_events_all_comparisons.tsv", delim = "\t")
all_AS_events_temp %>%
  View()

all_AS_events <- read_delim("./data/rMATS_out/All_events_all_comparisons.tsv", delim = "\t") %>%
  mutate(ID_modified = str_remove(ID, GeneID))

#make the header of this file using "All_DAS_events_all_comparisons.tsv"
#the reason to extract this header is that the header of "all_DAS_events" cannot be parsed correctly
#cat All_DAS_events_all_comparisons.tsv | awk -F'\t' 'BEGIN { OFS="," } { $1=$1 } 1' | awk -F' ' 'BEGIN { OFS="," } { $1=$1 } 1' | head -n 1 > header_All_DAS_events_all_comparisons
#https://www.thegeekstuff.com/2010/01/8-powerful-awk-built-in-variables-fs-ofs-rs-ors-nr-nf-filename-fnr/

header_All_DAS_events_all_comparisons <- read_csv("./data/rMATS_out/header_All_DAS_events_all_comparisons")

all_DAS_events <- read_delim("./data/rMATS_out/All_DAS_events_all_comparisons.tsv", skip = 1, delim = "\t", col_names = str_c("X", c(1:24)))

colnames(all_DAS_events) = colnames(header_All_DAS_events_all_comparisons)

#in the file "all_AS_events", ID contains GeneID and the coordinates of other events (eg. RI: exonStart, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE)
#the six positions of other events are listed in https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md
#the info that we want to get is from upstreamEE to downstreamES
#the file "All_DAS_events_all_comparisions.tsv) is based on FDR < 0.01
all_DAS_events_005 <- all_AS_events %>%
  filter(FDR < 0.05)

#check the data format
all_DAS_events_005 %>%
  View()

head(all_DAS_events_005$ID_modified)

#split ID_modified to get 6 position of the necessary information for alternative splicing
DAS_005_coordinate_all <- str_split(all_DAS_events_005$ID_modified, "_", simplify = TRUE)
AS_coordinate_all <- str_split(all_AS_events$ID_modified, "_", simplify = TRUE)

head(DAS_005_coordinate_all)

#check how many types of AS in this data
all_DAS_events_005 %>%
  group_by(AS_type) %>%
  summarise()

#generate the raw table for the following intersection files between marks and AS
#also generate the raw table for all events for producing the control set for comparing the profiles of different marks

#for DAS 
all_DAS_events_005_bed_raw_v1 <- tibble(
    chr = all_DAS_events_005$chr,
    GeneID = all_DAS_events_005$GeneID,
    strand = all_DAS_events_005$strand,
    pos_1 = as.double(DAS_005_coordinate_all[,2]),
    pos_2 = as.double(DAS_005_coordinate_all[,3]),
    pos_3 = as.double(DAS_005_coordinate_all[,4]),
    pos_4 = as.double(DAS_005_coordinate_all[,5]),
    pos_5 = as.double(DAS_005_coordinate_all[,6]),
    pos_6 = as.double(DAS_005_coordinate_all[,7]),
    AS_type = all_DAS_events_005$AS_type,
    comp = all_DAS_events_005$comp,
    FDR = all_DAS_events_005$FDR,
    IncLevel1 = all_DAS_events_005$IncLevel1,
    IncLevel2 = all_DAS_events_005$IncLevel2,
    Incdf = all_DAS_events_005$IncLevelDifference
  ) 

#for all AS
all_AS_events %>%
  View()


all_AS_events_bed_raw_v1 <- tibble(
  chr = all_AS_events$chr,
  GeneID = all_AS_events$GeneID,
  strand = all_AS_events$strand,
  pos_1 = as.double(AS_coordinate_all[,2]),
  pos_2 = as.double(AS_coordinate_all[,3]),
  pos_3 = as.double(AS_coordinate_all[,4]),
  pos_4 = as.double(AS_coordinate_all[,5]),
  pos_5 = as.double(AS_coordinate_all[,6]),
  pos_6 = as.double(AS_coordinate_all[,7]),
  AS_type = all_AS_events$AS_type,
  comp = all_AS_events$comp,
  FDR = all_AS_events$FDR,
  IncLevel1 = all_AS_events$IncLevel1,
  IncLevel2 = all_AS_events$IncLevel2,
  Incdf = all_AS_events$IncLevelDifference
) 


#based on the characteristics of different AS, different positions from these 6 ones will be different
#more information about how to choose correct position can be check from the website (https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md#output)
#Briefly, RI/SE are simple, the interval between upstreamEE (pos_4) and downstreamES (pos_5) are used, so I didn't include segments outside the event 
#For A5S/A3S, I used shortEE and flankingES, but the orientation was concerned (+/-)
#And for producing bed files, while end position was used for the starting position, I shifted 1 base to the left
#If I used starting position for the end position, I shifted 1 base to the right

all_DAS_events_005_bed_raw <- all_DAS_events_005_bed_raw_v1 %>%
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
  select(chr, str, end, AS_type, comp, FDR, IncLevel1, IncLevel2, Incdf, GeneID, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6) %>%
  filter(Incdf != 0) %>%
  mutate(Inclusion = if_else(Incdf > 0, "Inc", "Skip")) %>%
  arrange(FDR) %>%
  group_by(comp, AS_type, Inclusion) %>%
  mutate(event = c(1:n())) %>%
  ungroup() %>%
  arrange(chr, str) 

all_AS_events_bed_raw <- all_AS_events_bed_raw_v1 %>%
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
  select(chr, str, end, AS_type, comp, FDR, IncLevel1, IncLevel2, Incdf, GeneID, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6) %>%
  filter(Incdf != 0) %>%
  mutate(Inclusion = if_else(Incdf > 0, "Inc", "Skip")) %>%
  arrange(FDR) %>%
  group_by(comp, AS_type, Inclusion) %>%
  mutate(event = c(1:n())) %>%
  ungroup() %>%
  arrange(chr, str) 

write_delim(all_AS_events_bed_raw, "./analysis/all_AS_events_bed_raw", delim = "\t", col_names = FALSE)



#Calculate the difference of inclusion levels between treatments using IncLevel1 and IncLevel2
all_DAS005_IncLevel1_l <- str_split(all_DAS_events_005_bed_raw$IncLevel1, ",", simplify = TRUE)
all_DAS005_IncLevel2_l <- str_split(all_DAS_events_005_bed_raw$IncLevel2, ",", simplify = TRUE)
 
all_DAS_events_005_bed_full <- all_DAS_events_005_bed_raw %>%
   mutate(Inclvl1_r1 = as.double(all_DAS005_IncLevel1_l[,1]), Inclvl1_r2 = as.double(all_DAS005_IncLevel1_l[,2]),
          Inclvl2_r1 = as.double(all_DAS005_IncLevel2_l[,1]), Inclvl2_r2 = as.double(all_DAS005_IncLevel2_l[,2])) %>%
   replace_na(list(Inclvl1_r1 = 10^5, Inclvl1_r2 = 10^5, Inclvl2_r1 = 10^5, Inclvl2_r2 = 10^5)) %>%
   mutate(group = c(1:n())) %>%
   group_by(group) %>%
   mutate(Inclvl1_m = if_else(Inclvl1_r1 == 10^5 || Inclvl1_r2 == 10^5, Inclvl1_r1 + Inclvl1_r2 - 10^5, (Inclvl1_r1 + Inclvl1_r2)/2)) %>%
   mutate(Inclvl2_m = if_else(Inclvl2_r1 == 10^5 || Inclvl2_r2 == 10^5, Inclvl2_r1 + Inclvl2_r2 - 10^5, (Inclvl2_r1 + Inclvl2_r2)/2)) %>%
   filter(Inclvl1_m != Inclvl2_m) %>%
   mutate(Inc_max = if_else(Inclvl1_m > Inclvl2_m, Inclvl1_m, Inclvl2_m)) %>%
   ungroup() %>%
   filter(Inc_max >= 0.01) 

#making list for the following intersection (based on comparison and AS types)
all_DAS_events_005_bed_l <- all_DAS_events_005_bed_full %>%
  select(chr, str, end, FDR, Incdf, Inclusion, event, comp, AS_type) %>%
  split(.$comp) %>%
  map(. %>% split(.$AS_type))

comp_name_full = c("H1_H0", "H1_H6", "H6_H0")

for (i in seq_along(comp_name_full)) {
  a <- str_c(names(all_DAS_events_005_bed_l[[i]]), "_", comp_name_full[[i]], "_FDR005")
  
  b <-
    str_c(
      "data/AS_sig_FDR_bed/",
      a,
      ".bed"
    )
  
  pwalk(list(
    all_DAS_events_005_bed_l[[i]],
    b,
    delim = "\t",
    col_names = FALSE
  ),
  write_delim)
}

#I also used all events to intersect with epigenetic signals without using any list
all_DAS_events_005_bed_simplified <- all_DAS_events_005_bed_l %>%
  map(. %>% bind_rows(.)) %>%
  bind_rows() %>%
  select(chr, str, end, comp, AS_type) %>%
  distinct() %>%
  mutate(significance = "yes") %>%
  spread(key = comp, value = significance) %>%
  replace_na(list(HS1.vs.HS0 = "no", HS1.vs.HS6 = "no", HS6.vs.HS0 = "no")) %>%
  mutate(event = c(1:n()))

write_delim(all_DAS_events_005_bed_simplified, "./data/AS_sig_FDR_bed/all_DAS_events_005_bed_simplified", delim = "\t", col_names = FALSE)

#run bash scripts to first get the intersection between H1 vs H0 and H6 vs H0
#then get the intersection between the peaks of each mark and H610 intersected file
#make violin plot for the intensity of all marks intersected with H610 intervals 
#depending on the name of AS files to upload them together
allAS_FDR005_marks_intersected_files <- fs::dir_ls("./data/Intersected_marks_AS/", glob = "*FDR005_all.bed")

#import Inc bed files together (IS means intersected) 
#produce the function to create files which can be used to plot
allAS_mark_IS_combined_rep_f <- function(data){
  a <- vroom(data, col_types = cols(), id = "name", delim = "\t", col_names = c("chr", "str", "end", "intensity", "Inclusion", "label_H1_H0", "label_H6_H0", "AS_type"))
  
  b <- str_split(a$name, "_", simplify = TRUE)
  
  a %>%
    mutate(trt = b[,6], mark = b[,7], rep = b[,8]) %>%
    select(-name) %>%
    split(.$rep)
}

#produce df reps of Inc/Skip 
allAS_FDR005_mark_IS_combined_rep <- allAS_mark_IS_combined_rep_f(allAS_FDR005_marks_intersected_files)
allAS_FDR005_mark_IS_combined_rep$rep1

#produce Inclusion/Skip isoforms for rep1
allAS_FDR005_mark_IS_combined_r1_raw <- allAS_FDR005_mark_IS_combined_rep$rep1 %>%
  mutate(trt = str_c("trt_H", str_sub(trt, 1, 1))) %>%
  spread(key = trt, value = intensity) %>%
  #here we need to think about whether to keep na
  #since not all DAS overlap the peaks of marks from two comp
  replace_na(list(trt_H0 = 1, trt_H1 = 1, trt_H6 = 1)) %>%
  #drop_na() %>%
  gather(key = "trt", value = "intensity", c("trt_H0", "trt_H1", "trt_H6"))

#produce intervals with all AS in the comparisons (H1 vs.H0 and H6 vs.H0)
AS_coordinate_all <- str_split(all_AS_events$ID_modified, "_", simplify = TRUE)

all_AS_events_H60_H10 <- tibble(
  chr = all_AS_events$chr,
  strand = all_AS_events$strand,
  pos_1 = as.double(AS_coordinate_all[,2]),
  pos_2 = as.double(AS_coordinate_all[,3]),
  pos_3 = as.double(AS_coordinate_all[,4]),
  pos_4 = as.double(AS_coordinate_all[,5]),
  pos_5 = as.double(AS_coordinate_all[,6]),
  pos_6 = as.double(AS_coordinate_all[,7]),
  AS_type = all_AS_events$AS_type,
  comp = all_AS_events$comp
) %>%
  mutate(str = if_else(AS_type == "A5S" & strand == "+", pos_4, 
                       if_else(AS_type == "A5S" & strand == "-", pos_6,
                               if_else(AS_type == "A3S" & strand == "+", pos_6,
                                       if_else(AS_type == "A3S" & strand == "-", pos_4, pos_4))))) %>%
  mutate(end = if_else(AS_type == "A5S" & strand == "+", pos_5, 
                       if_else(AS_type == "A5S" & strand == "-", pos_3,
                               if_else(AS_type == "A3S" & strand == "+", pos_3,
                                       if_else(AS_type == "A3S" & strand == "-", pos_5, pos_5))))) %>%
  select(chr, str, end, AS_type, comp) %>%
  mutate(comp = str_remove_all(comp, "S")) %>%
  mutate(comp = str_replace(comp, ".vs.", "_")) %>%
  filter(comp == "H1_H0" | comp == "H6_H0") %>%
  arrange(chr, str)
  
write_delim(all_AS_events_H60_H10, "./data/shuffled_intervals_for_control/all_AS_events_H60_H10.bed", delim = "\t", col_names = FALSE)

#make plots

allAS_FDR005_mark_IS_combined_r1_raw %>%
  View()
  filter(intensity < 0)

allAS_FDR005_mark_IS_r1_plot <- allAS_FDR005_mark_IS_combined_r1_raw %>%
  split(.$AS_type) %>%
  map(. %>% split(.$Inclusion)) %>%
  map(. %>% map(. %>% violin_plot_f()))

allAS_FDR005_mark_IS_r1_plot$RI$Inc

trt_order <- c("A3S", "A5S", "RI", "SE")
Inclusion_order <- c("Inc", "Skip")


for (i in seq_along(trt_order)) {
  for (j in seq_along(Inclusion_order)) {
    ggsave(str_c("./analysis/", trt_order[[i]], "_", Inclusion_order[[j]], "_H610_r1.jpeg"), allAS_FDR005_mark_IS_r1_plot[[i]][[j]], width = 300, height = 180, units = c("mm"), dpi = 600)
  }
}

#take the 10% of Incdf from either H1 vs H0 or H6 vs H0 to check their plots and perform linear regression
all_DAS_events_005_bed_full_Incdf_q90 <- all_DAS_events_005_bed_full %>%
  group_by(AS_type, Inclusion) %>%
  filter(abs(Incdf) >= quantile(abs(Incdf), 0.9)) %>%
  ungroup() %>%
  select(AS_type, Inclusion, comp, Incdf, event) %>%
  mutate(comp = str_remove_all(comp, "S")) %>%
  mutate(comp = str_replace(comp, ".vs.", "_"))

#make plots
allAS_FDR005_mark_IS_r1_Incdf_q90_plot <- allAS_FDR005_mark_IS_combined_r1_raw %>%
  gather(key = comp, value = event, c("label_H1_H0", "label_H6_H0")) %>%
  mutate(comp = str_remove(comp, "label_")) %>%
  left_join(all_DAS_events_005_bed_full_Incdf_q90) %>%
  drop_na() %>%
  split(.$AS_type) %>%
  map(. %>% split(.$Inclusion)) %>%
  map(. %>% map(. %>% violin_plot_f()))


for (i in seq_along(trt_order)) {
  for (j in seq_along(Inclusion_order)) {
    ggsave(str_c("./analysis/", trt_order[[i]], "_", Inclusion_order[[j]], "_Incdf_q90_H610_r1.jpeg"), allAS_FDR005_mark_IS_r1_Incdf_q90_plot[[i]][[j]], width = 300, height = 180, units = c("mm"), dpi = 600)
  }
}

#based on the intersection between all AS events (mixed in one file) and 21 epigenetic mark files
#I calculated the signal of each mark in each event
all_DAS_005_marks_files <- fs::dir_ls("./data/Intersected_marks_AS/", regexp = "all_DAS_005_M82")

all_DAS_005_marks_files_c <- vroom(all_DAS_005_marks_files, col_types = cols(), id = "name", delim = "\t", col_names = c(colnames(all_DAS_events_005_bed_simplified), "pileup", "intersected_bp"))

file_name_split <- str_split(all_DAS_005_marks_files_c$name, "_", simplify = TRUE)

all_DAS_005_mark_mean_signal <- all_DAS_005_marks_files_c %>%
  mutate(
    trt = file_name_split[, 7],
    mark = file_name_split[, 8],
    rep = str_sub(file_name_split[, 9], 1, 4)
  ) %>%
  mutate(mark_trt_rep = str_c(mark, "_", trt, "_", rep)) %>%
  #select(-name, -trt, -mark, -rep) %>%
  mutate(total_signal = intersected_bp * pileup) %>%
  group_by(event, mark_trt_rep, trt, mark, rep, str, end) %>%
  summarise(sum_total_signal = sum(total_signal)) %>%
  mutate(mean_signal = sum_total_signal / (end - str)) %>%
  ungroup() %>%
  select(-str, -end, -sum_total_signal) %>%
  arrange(event, mark, trt) %>%
  select(-trt, -mark, -rep) %>%
  spread(mark_trt_rep, mean_signal) %>%
  replace_na(
    list(
      H3K18ac_0h_rep1 = 0,
      H3K18ac_1h_rep1 = 0,
      H3K18ac_6h_rep1 = 0,
      H3K27ac_0h_rep1 = 0,
      H3K27ac_1h_rep1 = 0,
      H3K27ac_6h_rep1 = 0,
      H3K4me3_0h_rep1 = 0,
      H3K4me3_0h_rep2 = 0,
      H3K4me3_1h_rep1 = 0,
      H3K4me3_1h_rep2 = 0,
      H3K4me3_6h_rep1 = 0,
      H3K4me3_6h_rep2 = 0,
      H3K9ac_0h_rep1 = 0,
      H3K9ac_1h_rep1 = 0,
      H3K9ac_6h_rep1 = 0,
      Pol2_0h_rep1 = 0,
      Pol2_0h_rep2 = 0,
      Pol2_1h_rep1 = 0,
      Pol2_1h_rep2 = 0,
      Pol2_6h_rep1 = 0,
      Pol2_6h_rep2 = 0
    )
  ) %>%
  select(
    event,
    H3K18ac_0h_rep1,
    H3K18ac_1h_rep1,
    H3K18ac_6h_rep1,
    H3K27ac_0h_rep1,
    H3K27ac_1h_rep1,
    H3K27ac_6h_rep1,
    H3K4me3_0h_rep1,
    H3K4me3_0h_rep2,
    H3K4me3_1h_rep1,
    H3K4me3_1h_rep2,
    H3K4me3_6h_rep1,
    H3K4me3_6h_rep2,
    H3K9ac_0h_rep1,
    H3K9ac_1h_rep1,
    H3K9ac_6h_rep1,
    Pol2_0h_rep1,
    Pol2_0h_rep2,
    Pol2_1h_rep1,
    Pol2_1h_rep2,
    Pol2_6h_rep1,
    Pol2_6h_rep2
  )



all_DAS_events_005_bed_mark_signal <- all_DAS_events_005_bed_simplified %>%
  left_join(all_DAS_005_mark_mean_signal) %>%
  select(-event)

all_DAS_events_005_bed_mark_final <- all_DAS_events_005_bed_raw %>%
  select(-event) %>%
  left_join(all_DAS_events_005_bed_mark_signal, by = c("chr", "str", "end", "AS_type")) %>%
  replace_na(
    list(
      H3K18ac_0h_rep1 = 0,
      H3K18ac_1h_rep1 = 0,
      H3K18ac_6h_rep1 = 0,
      H3K27ac_0h_rep1 = 0,
      H3K27ac_1h_rep1 = 0,
      H3K27ac_6h_rep1 = 0,
      H3K4me3_0h_rep1 = 0,
      H3K4me3_0h_rep2 = 0,
      H3K4me3_1h_rep1 = 0,
      H3K4me3_1h_rep2 = 0,
      H3K4me3_6h_rep1 = 0,
      H3K4me3_6h_rep2 = 0,
      H3K9ac_0h_rep1 = 0,
      H3K9ac_1h_rep1 = 0,
      H3K9ac_6h_rep1 = 0,
      Pol2_0h_rep1 = 0,
      Pol2_0h_rep2 = 0,
      Pol2_1h_rep1 = 0,
      Pol2_1h_rep2 = 0,
      Pol2_6h_rep1 = 0,
      Pol2_6h_rep2 = 0
    ) 
  )

write_csv(all_DAS_events_005_bed_mark_final, "./analysis/all_DAS_events_005_bed_mark_final.csv", col_names = TRUE)




#generate the list of rep1 for linear regression (all events with FDR below 0.05)
all_DAS_events_005_bed_light <- all_DAS_events_005_bed_full %>%
  select(AS_type, Inclusion, comp, Incdf, event) %>%
  mutate(comp = str_remove_all(comp, "S")) %>%
  mutate(comp = str_replace(comp, ".vs.", "_"))

#generate the table with the intensity difference of marks of all peaks 
allAS_FDR005_markdf_IS_combined_r1 <- allAS_FDR005_mark_IS_combined_rep$rep1 %>%
  mutate(trt = str_c("trt_H", str_sub(trt, 1, 1))) %>%
  spread(key = trt, value = intensity) %>%
  replace_na(list(trt_H0=0, trt_H1=0, trt_H6=0)) %>%
  mutate(markdf_H10 = trt_H1 - trt_H0, markdf_H60 = trt_H6-trt_H0) %>%
  gather(key = comp, value = event, c("label_H1_H0", "label_H6_H0")) %>%
  mutate(comp = str_remove(comp, "label_")) %>%
  select(AS_type, Inclusion, comp, event, mark, trt_H0, trt_H1, trt_H6, markdf_H10, markdf_H60, rep)
  
#using the table of DAS (FDR < 0.05) with Incdf to connect the table of marks with mark peaks in the difference between treatments
#anchor: firstly use AS_type, Inclusion, comp, event to connect two tables
#then make sure Incdf and markdf come from the same comparison
df_mark_inc_lr_l <- all_DAS_events_005_bed_light %>%
  left_join(allAS_FDR005_markdf_IS_combined_r1) %>%
  drop_na() %>%
  gather(key = comp_mark, value = markdf, c("markdf_H10", "markdf_H60")) %>%
  mutate(comp_mark = str_remove(comp_mark, "markdf_")) %>%
  mutate(comp_mark = str_c(str_sub(comp_mark, 1, 2), "_H", str_sub(comp_mark, 3, 3))) %>%
  filter(comp == comp_mark) %>%
  group_by(comp, AS_type, Inclusion, mark) %>%
  nest()
  

mark_model <- function(df) {
  lm(Incdf ~ markdf, data = df)
}


df_mark_inc_lr_l <- df_mark_inc_lr_l %>%
  mutate(model = map(data, mark_model))

df_mark_inc_lr_l_result <- df_mark_inc_lr_l %>%
  mutate(glance = map(model, broom::glance)) %>%
  unnest(glance) 

df_mark_inc_lr_l_result_sig <- df_mark_inc_lr_l_result %>%
  arrange(p.value) %>%
  filter(p.value < 0.05) %>%
  arrange(p.value) 

df_mark_inc_lr_l_result_sig %>%
  select(AS_type, mark, comp, Inclusion, p.value) %>%
  filter(p.value == min(p.value))

df_mark_inc_lr_l$data 

summary(df_mark_inc_lr_l_result_sig$model[[17]])
plot(df_mark_inc_lr_l_result_sig$data[[1]]$Incdf, df_mark_inc_lr_l_result_sig$data[[1]]$markdf)

df_mark_inc_lr_l_result_sig$data[[1]] %>%
  ggplot() +
  geom_point(aes(Incdf, markdf))

summary(df_mark_inc_lr_l_result_sig$model[[1]])
