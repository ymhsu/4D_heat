all_AS_intersected_feature_raw %>% 
  left_join(M82_rMATs_anno_all_gene_light) %>%
  mutate(gene_size = gene_size %/% 1000) %>%
  group_by(gene_size, AS, event, feature, order) %>%
  summarise(count = n()) %>%
  View()

all_AS_intersected_feature_raw %>%
  group_by(AS, event, feature, order) %>%
  summarise(count = n()) %>%
  View()

all_AS_intersected_feature_raw %>%
  filter(anno == "Solyc_chr1_NPC00337.MC1")
gene_having_AS_info_l_f$exon %>%
  filter(anno == "Solyc_chr1_NPC00337.MC1")
M82_rMATs_anno_all_exon_order %>%
  group_by(anno) %>%
  mutate(pre_end = lag(end)) %>%
  drop_na() %>%
  mutate(dist_exon = str - pre_end) %>%
  ungroup() %>%
  summarise(q25 = quantile(dist_exon, 0.25), q50 = median(dist_exon), q75 = quantile(dist_exon, 0.75))


exon_intron_data_for_ML_list$exon$Solyc_chr2_NPC00216.MC2

make_AS_ML_data_f(exon_intron_data_for_ML_list$exon$Solyc_chr2_NPC00216.MC2, "HS6_HS0", "PI_H", "RI")

exon_intron_data_for_ML_list$exon$Solyc_chr2_NPC00216.MC2 %>%
  filter(event == "no") %>%
  filter(order != 4-1 & order != 4+1)

make_AS_ML_data_f <- function(data, a, b, c){ 
data_a <- data %>%
  ungroup()
AS_sig_m <- data_a %>%
  filter(comp == a[1] & PI == b[1] &  AS == c[1] & event == "sig") 

AS_non_sig_m <- data_a %>%
  filter(comp == a[1] & PI == b[1] &  AS == c[1] & event == "non-sig") 

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
length(exon_intron_data_for_ML_list$exon)
ML_data_test_head10 <- list(head(exon_intron_data_for_ML_list$exon, n = 10), head(exon_intron_data_for_ML_list$intron, n = 10))
system.time(ML_data_test_head10 %>%
  map(. %>% map(. %>% make_AS_ML_data_f(., comb_comp_PI_AS_list[[1]][[1]]))) %>%
  map(. %>% bind_rows()))

exon_intron_data_for_ML_comb_l[[1]][[1]][[4]][[4]] <- ML_data_test_head10[[1]] %>%
  map(. %>% make_AS_ML_data_f(., comp[[1]], PI_type[[4]], AS_type[[4]])) %>%
  bind_rows()

exon_intron_data_for_ML_list

registerDoParallel(cores = 10)
getDoParWorkers()


