install.packages("ggpubr")
library(ggpubr)

#make histograms showing the distribution of PSI of all AS or only sig AS events
all_AS_hist_p <- all_AS_events_bed_TPM_q05 %>%
  mutate(event = if_else(FDR < 0.05, "sig", "non-sig")) %>%
  #filter(event == "sig") %>%
  mutate(Inc_former_trt = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_latter_trt = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  ggplot() +
  geom_histogram(aes(Inc_former_trt), bins = 50) +
  facet_grid(comp ~ AS_type)

sig_AS_hist_p <- all_AS_events_bed_TPM_q05 %>%
  mutate(event = if_else(FDR < 0.05, "sig", "non-sig")) %>%
  filter(event == "sig") %>%
  mutate(Inc_former_trt = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_latter_trt = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  ggplot() +
  geom_histogram(aes(Inc_former_trt), bins = 50) +
  facet_grid(comp ~ AS_type) 

plot_func_AS <- function(a){
  a +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))
  
}

all_AS_hist_p_f <- plot_func_AS(all_AS_hist_p)
sig_AS_hist_p_f <- plot_func_AS(sig_AS_hist_p)

ggsave(str_c("./analysis_output/all_AS_hist_p_f.jpeg"), all_AS_hist_p_f, width = 360, height = 270, units = c("mm"), dpi = 320)
ggsave(str_c("./analysis_output/sig_AS_hist_p_f.jpeg"), sig_AS_hist_p_f , width = 360, height = 270, units = c("mm"), dpi = 320)



#Use the data from "all_AS_events_TPM.R" to check whether PI_H and PI_L change in different ways in sig-AS events or non-sig AS events

all_AS_events_bed_TPM_q05_raw %>%
  filter(FDR < 0.05, "sig", "non-sig") %>%
  filter(comp == "HS1.vs.HS6") %>%
  filter(PI == "PI_H" | PI == "PI_L") %>%
  mutate(Inc_1H = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_6H = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  mutate(Inc_rise = if_else(Incdf > 0, "yes", if_else(Incdf < 0, "no", "equal"))) %>%
  group_by(AS_type, PI, Inc_rise) %>%
  summarise(count = n()) %>%
  View()

#making scatter plots for identifying the trend of changing PSI between former and latter trt

all_AS_events_bed_TPM_q05_v2 <- all_AS_events_bed_TPM_q05 %>%
  mutate(event = if_else(FDR < 0.05, "sig", "non-sig")) %>%
  mutate(Inc_former_trt = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_latter_trt = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  mutate(PSI_group = if_else(Inc_former_trt <= 0.5, "small", "large")) %>%
  group_by(comp, AS_type, PSI_group) %>%
  mutate(median_PSI_latter = median(Inc_latter_trt), mean_PSI_latter = mean(Inc_latter_trt)) %>%
  ungroup() %>%
  mutate(x_axis = if_else(PSI_group == "small", 0.25, 0.75))

temp_PSI_analysis <- all_AS_events_bed_TPM_q05_v2 %>%
  select(comp, AS_type, PSI_group, median_PSI_latter, mean_PSI_latter) %>%
  distinct()

temp_PSI_analysis_non <- all_AS_events_bed_TPM_q05_v2 %>%
  filter(event == "non-sig") %>%
  group_by(comp, AS_type, PSI_group) %>%
  mutate(median_PSI_latter_nons = median(Inc_latter_trt), mean_PSI_latter_nons = mean(Inc_latter_trt)) %>%
  ungroup() %>%
  select(comp, AS_type, PSI_group, median_PSI_latter_nons, mean_PSI_latter_nons) %>%
  distinct()

##----####
temp_PSI_analysis %>%
  left_join(temp_PSI_analysis_non) %>%
  mutate(df_median = median_PSI_latter - median_PSI_latter_nons, df_mean = mean_PSI_latter-mean_PSI_latter_nons) %>%
  View()
##---####

all_AS_events_bed_TPM_q05_v2_l <- all_AS_events_bed_TPM_q05_v2 %>%
  left_join(temp_PSI_analysis_non) %>%
  split(.$comp) 


for (i in seq_along(comp_title)) {
  a <- all_AS_events_bed_TPM_q05_v2_l[[i]] %>%
    ggplot() +
    geom_point(aes(x = Inc_former_trt, y = Inc_latter_trt, color = event), alpha = 0.2) +
    #geom_segment(aes(x = x_axis - 0.25, y = median_PSI_latter, xend = x_axis + 0.25, yend = median_PSI_latter), color = "yellow") +
    geom_segment(aes(x = x_axis - 0.25, y = mean_PSI_latter, xend = x_axis + 0.25, yend = mean_PSI_latter), color = "green") +
    #geom_segment(aes(x = x_axis - 0.25, y = median_PSI_latter_nons, xend = x_axis + 0.25, yend = median_PSI_latter_nons), color = "#CCCC00") +
    geom_segment(aes(x = x_axis - 0.25, y = mean_PSI_latter_nons, xend = x_axis + 0.25, yend = mean_PSI_latter_nons), color = "darkgreen") +
    facet_wrap( ~ AS_type, nrow = 2) +
    ggtitle(comp_title[[i]]) +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          plot.title = element_text(colour = "black", face = "bold", size = 18, hjust = 0.5), axis.text.y = element_text(size = 18, face = "bold"), 
          strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))
  ggsave(str_c("./analysis_output/all_AS_TPM_q05_sig_trend_", str_replace(comp_title [[i]], ".vs.", "_"), "_p.jpeg"), a, width = 360, height = 270, units = c("mm"), dpi = 320)
  
}





all_AS_events_bed_TPM_q05_v2 %>%
  filter(event == "sig" & comp == "HS1.vs.HS6" & AS_type == "RI") %>%
  ggplot() +
  geom_point(aes(x = Inc_former_trt, y = Inc_latter_trt, color = event), alpha = 0.2) +
  geom_segment(aes(x = x_axis - 0.25, y = median_PSI_latter, xend = x_axis + 0.25, yend = median_PSI_latter), color = "#CCCC00") +
  geom_segment(aes(x = x_axis - 0.25, y = mean_PSI_latter, xend = x_axis + 0.25, yend = mean_PSI_latter), color = "darkgreen") +
  facet_grid(comp ~ AS_type)



all_test_plot_AS_f <- plot_func_AS(all_test_plot_AS) 


sig_test_plot_AS_f <- plot_func_AS(sig_test_plot_AS)

ggsave("./analysis_output/sig_test_plot_AS_f.jpeg", sig_test_plot_AS_f, width = 400, height = 300, units = c("mm"), dpi = 320)
ggsave("./analysis_output/all_test_plot_AS_f.jpeg", all_test_plot_AS_f, width = 400, height = 300, units = c("mm"), dpi = 320)


all_AS_events_bed_TPM_q05 

test <- all_DAS_events_bed_TPM_q05_raw %>%
  filter(comp == "HS1.vs.HS6") %>%
  filter(PI == "PI_H" | PI == "PI_L") %>%
  mutate(Inc_1H = (Inclvl1_r1_n + Inclvl1_r2_n)/2, Inc_6H = (Inclvl2_r1_n + Inclvl2_r2_n)/2) %>%
  ggplot() +
  geom_point(aes(x = Inc_1H, y = Inc_6H)) +
  facet_wrap(~ AS_type)


all_DAS_events_bed_TPM_q05_raw %>%
  filter(comp == "HS1.vs.HS6") %>%
  filter(PI == "PI_H" | PI == "PI_L") %>%
  mutate(Inc_1H = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_6H = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  mutate(Inc_rise = if_else(Incdf > 0, "yes", if_else(Incdf < 0, "no", "equal"))) %>%
  group_by(comp, AS_type, PI, Inc_rise) %>%
  summarise(count = n()) %>%
  View()

all_DAS_events_bed_TPM_q05_raw %>%
  filter(comp == "HS6.vs.HS0") %>%
  filter(PI == "PI_H" | PI == "PI_L") %>%
  mutate(Inc_6H = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_0H = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  mutate(Inc_rise = if_else(Incdf < 0, "yes", if_else(Incdf > 0, "no", "equal"))) %>%
  group_by(comp, AS_type, PI, Inc_rise) %>%
  summarise(count = n()) %>%
  View()


all_AS_events_bed_TPM_raw %>%
  filter(IncLevel1 == "0.334,0.206")

all_DAS_events_bed_TPM_q05_raw %>%
  filter(comp == "HS1.vs.HS6") %>%
  View()

all_DAS_events_bed_TPM_q05_raw %>%
  mutate(Incdf_new = if_else(comp == "HS1.vs.HS6", -Incdf, Incdf)) %>%
  mutate(comp_new = if_else(comp == "HS1.vs.HS6", "HS6.vs.HS1", comp)) %>%
  filter(PI == "PI_H" | PI == "PI_L") %>%
  ggplot() +
  geom_point(aes(AS_type, Incdf_new, color = PI), alpha = 0.5) +
  geom_violin(aes(AS_type, Incdf_new, fill = PI)) +
  facet_wrap(~comp_new) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))


?geom_

all_DAS_events_bed_TPM_q05_raw %>%
  group_by(comp) %>%
  summarise()
