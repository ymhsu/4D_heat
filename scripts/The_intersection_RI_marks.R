library(pacman)
Packages <- c("tidyverse", "fc", "vroom", "gt", "plotly", "modelr")
p_load(Packages, character.only = TRUE)

#Since in plants, most alternative splicing (AS) is the intron retention, this analysis was only based this type of AS.
#another file called "RNAseq_fc_analysis was based on the workflow of this file, but extended for other AS types.
#The main idea of this analysis is to compare significant AS of RI with histone modification marks overlapping these RI events
#, and to check the trend of these histone modification marks.
#depending on the name of AS files to upload them together
#preliminary analysis on RI data which only focused on H6 vs H0 and H1 vs H0
#RI.MATS.JC.txt from the subfolders (H6.vs.H0/H1.vs.Ho) was renamed as RI_JC_(H6_H0/H1_H0).txt stored in the rMAT_out folder
AS_RI_files <- fs::dir_ls("./data/rMATS_out/", glob = "*H0.txt")

#based on this blog (https://heads0rtai1s.github.io/2022/02/24/r-python-read-separate/)
AS_RI_files_combined_raw <- vroom(AS_RI_files, col_types = cols(), id = "name", delim = "\t")

RI_IncLevel1_l <- str_split(AS_RI_files_combined_raw$IncLevel1, ",", simplify = TRUE)
RI_IncLevel2_l <- str_split(AS_RI_files_combined_raw$IncLevel2, ",", simplify = TRUE)

AS_RI_files_combined <- AS_RI_files_combined_raw %>%
  mutate(Inclvl1_r1 = as.double(RI_IncLevel1_l[,1]), Inclvl1_r2 = as.double(RI_IncLevel1_l[,2]),
         Inclvl2_r1 = as.double(RI_IncLevel2_l[,1]), Inclvl2_r2 = as.double(RI_IncLevel2_l[,2])) %>%
  replace_na(list(Inclvl1_r1 = 10^5, Inclvl1_r2 = 10^5, Inclvl2_r1 = 10^5, Inclvl2_r2 = 10^5)) %>%
  mutate(group = c(1:n())) %>%
  group_by(group) %>%
  mutate(Inclvl1_m = if_else(Inclvl1_r1 == 10^5 || Inclvl1_r2 == 10^5, Inclvl1_r1 + Inclvl1_r2 - 10^5, (Inclvl1_r1 + Inclvl1_r2)/2)) %>%
  mutate(Inclvl2_m = if_else(Inclvl2_r1 == 10^5 || Inclvl2_r2 == 10^5, Inclvl2_r1 + Inclvl2_r2 - 10^5, (Inclvl2_r1 + Inclvl2_r2)/2)) %>%
  filter(Inclvl1_m != Inclvl2_m) %>%
  mutate(Inc_max = if_else(Inclvl1_m > Inclvl2_m, Inclvl1_m, Inclvl2_m)) %>%
  ungroup() %>%
  filter(Inc_max >= 0.01)

#tidy the name of different files
AS_RI_files_l <- AS_RI_files_combined %>%
  mutate(name = str_replace(name, "[.]/data/rMATS_out/", "")) %>%
  mutate(name = str_replace(name, "[.]txt", "")) %>%
  split(.$name)

#produce the list of AS with significantly different (FDR < 0.05) AS events between treatments
AS_RI_files_sig_l <- AS_RI_files_l %>%
  map(. %>% filter(FDR < 0.05)) 


#produce bed files for IGV visualization 
RI_HS1_0_significant_bed <- tibble(
  chr = rep(AS_RI_files_sig_l[[1]]$chr, 2),
  str = c(AS_RI_files_sig_l[[1]]$upstreamES, AS_RI_files_sig_l[[1]]$downstreamES),
  end = c(AS_RI_files_sig_l[[1]]$upstreamEE, AS_RI_files_sig_l[[1]]$downstreamEE),
  FDR = rep(AS_RI_files_sig_l[[1]]$FDR, 2),
  ID = rep(AS_RI_files_sig_l[[1]]$ID...1, 2)
) %>%
  arrange(FDR, ID) %>%
  select(-ID)

write_delim(RI_HS1_0_significant_bed, "./analysis/RI_HS1_0_significant_bed", col_names = FALSE, delim = "\t")


#produce bed files for the intersection between AS and histone modification marks
#create the function for separating Inc (Inclusion) and Skip (Skipping) type of isoforms
produce_Inc_Sk_AS_l <- function(data){
  tibble(
    chr = data$chr,
    str = data$upstreamEE-1,
    end = data$downstreamES,
    FDR = data$FDR,
    Incdf = data$IncLevelDifference
  ) %>%
    filter(Incdf != 0) %>%
    mutate(Inclusion = if_else(Incdf > 0, "Inc", "Skip")) %>%
    arrange(FDR) %>%
    group_by(Inclusion) %>%
    mutate(event = c(1:n())) %>%
    arrange(chr, str) 
  
}  

#two lists for Inclusion or Skip type of H1/H0 and H6/H0
RI_JC_sig_bed_l <- vector("list", length = length(AS_RI_files_sig_l))

for (i in seq_along(AS_RI_files_sig_l)) {
  RI_JC_sig_bed_l[[i]] <- produce_Inc_Sk_AS_l(AS_RI_files_sig_l[[i]])
}


trt_name <-
  c("RI_JC_H1_H0_sig", "RI_JC_H6_H0_sig")

paths_RI_trt_table <-
  str_c(
    "data/AS_sig_FDR_bed/",
    trt_name,
    ".bed"
  )

pwalk(list(
  RI_JC_sig_bed_l,
  paths_RI_trt_table,
  delim = "\t",
  col_names = FALSE
),
write_delim)




#make violin plot
#depending on the name of AS files to upload them together
allAS_RI_marks_intersected_files <- fs::dir_ls("./data/Intersected_marks_AS/", regexp = "*rep[12]_all\\.bed")

#import Inc bed files together (IS means intersected) 
#produce the function to create files which can be used to plot

AS_RI_mark_IS_combined_rep_f <- function(data){
  a <- vroom(data, col_types = cols(), id = "name", delim = "\t", col_names = c("chr", "str", "end", "intensity", "Inclusion", "label_H1_H0", "label_H6_H0"))
  
  b <- str_split(a$name, "_", simplify = TRUE)
  
  a %>%
    mutate(trt = b[,6], mark = b[,7], rep = b[,8]) %>%
    select(-name) %>%
    split(.$rep)
}


#produce df reps of Inc/Skip 
allAS_RI_mark_IS_combined_rep <- AS_RI_mark_IS_combined_rep_f(allAS_RI_marks_intersected_files)


violin_plot_f <- function(data){
  data %>%
    ggplot(aes(trt, y=intensity)) +
    geom_violin(aes(fill = trt), trim = FALSE, draw_quantiles = TRUE) +
    scale_fill_manual(values=c("#7fc97f", "#beaed4", "#fdc086"))+
    geom_boxplot(width=0.1)+
    facet_wrap(~mark, scales = "free") + theme_grey()
}


#produce Inclusion/Skip isoform for rep1
RI_Inc_H610_r1 <- allAS_RI_mark_IS_combined_rep$rep1 %>%
  filter(Inclusion == "Inc") %>%
  violin_plot_f() 


RI_Skip_H610_r1 <- allAS_RI_mark_IS_combined_rep$rep1 %>%
  filter(Inclusion == "Skip") %>%
  violin_plot_f() 

ggsave("./analysis/RI_Inc_H610_r1.jpeg", RI_Inc_H610_r1, width = 300, height = 180, units = c("mm"), dpi = 600)
ggsave("./analysis/RI_Skip_H610_r1.jpeg", RI_Skip_H610_r1, width = 300, height = 180, units = c("mm"), dpi = 600)

#extract the difference of inclusion frequency at the top 10% of events
RI_H1_H0_q_90_l <- bind_rows(produce_Inc_Sk_AS_l(AS_RI_files_sig_l$RI_JC_H1_H0)) %>%
  mutate(name = names(AS_RI_files_sig_l)[[1]]) %>%
  group_by(Inclusion) %>%
  filter(Incdf >= quantile(Incdf, 0.9)) %>%
  ungroup() %>%
  select(label_H1_H0 = event, Inclusion) %>%
  mutate(q_90_H1_H0 = "yes")

RI_H6_H0_q_90_l <- bind_rows(produce_Inc_Sk_AS_l(AS_RI_files_sig_l$RI_JC_H6_H0)) %>%
  mutate(name = names(AS_RI_files_sig_l)[[2]]) %>%
  group_by(Inclusion) %>%
  filter(Incdf >= quantile(Incdf, 0.9)) %>%
  ungroup() %>%
  select(label_H6_H0 = event, Inclusion) %>%
  mutate(q_90_H6_H0 = "yes")


AS_RI_mark_plot_q_90_l <- bind_rows(allAS_RI_mark_IS_combined_rep) %>%
  left_join(RI_H1_H0_q_90_l) %>%
  left_join(RI_H6_H0_q_90_l) %>%
  replace_na(list(q_90_H1_H0 = "no", q_90_H6_H0 = "no")) %>%
  filter(q_90_H1_H0 != "no" | q_90_H6_H0 != "no") %>%
  split(.$Inclusion) %>%
  map(. %>% split(.$rep))

violin_plot_f(AS_mark_plot_q_90_l$Inc$rep1)