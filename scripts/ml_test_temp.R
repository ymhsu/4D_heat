#machine learning by Boruta

install.packages("Boruta")
install.packages("mlbench")
library(tidyverse)
library(Boruta)
library(mlbench)
data(Sonar)
version

Sonar
as.tibble(Sonar) %>%
  select(Class) %>%
  arrange(Class) %>%
  group_by(Class) %>%
  summarise(count = n())
?Boruta

test_ML_df <- as.data.frame(test_ML)
test_ML_2_df <- as.data.frame(test_ML_2)
test_ML_3_df <- as.data.frame(test_ML_3)
test_ML_4_df <- as.data.frame(test_ML_4)

boruta.HS6_HS0_PI_H_exon <- Boruta(type~.,data=test_ML_df,doTrace = 2)
boruta.HS6_HS0_PI_L_exon <- Boruta(type~.,data=test_ML_2_df,doTrace = 2)
boruta.HS6_HS0_PI_H_L_ctrl_exon <- Boruta(type~.,data=test_ML_3_df,doTrace = 2)
boruta.HS6_HS0_PI_H_L_exon_four <- Boruta(type~.,data=test_ML_4_df,doTrace = 2)



Bor.son <- Boruta(Class~.,data=Sonar,doTrace=2)

plot(boruta.HS6_HS0_PI_H_exon)
plot(boruta.HS6_HS0_PI_L_exon)
plot(boruta.HS6_HS0_PI_H_L_ctrl_exon)
plot(boruta.HS6_HS0_PI_H_L_exon_four)
attStats(boruta.HS6_HS0_PI_H_exon)
attStats(boruta.HS6_HS0_PI_L_exon)
attStats(boruta.HS6_HS0_PI_H_L_ctrl_exon)
attStats(boruta.HS6_HS0_PI_H_L_exon_four)
attStats(Bor.son)

boruta.HS6_HS0_PI_H_L_ctrl_exon$finalDecision

tibble(Sonar)

#PCA
install.packages("ggfortify")
library(ggfortify)

ML_l_f_pca_df <- vector("list", length = length(ML_l_f))
ML_l_f_pca_r <- vector("list", length = length(ML_l_f))


for (i in seq_along(ML_l_f)) {
  ML_l_f_pca_df[[i]] <- ML_l_f[[i]] %>%
    select(-type)
  ML_l_f_pca_r[[i]] <- prcomp(ML_l_f_pca_df[[i]], scale. = TRUE)
  
}
autoplot(ML_l_f_pca_r[[20]], data = ML_l_f[[20]], colour = 'type')

install.packages("randomForest")

samp

#reference modify Boruta's result
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


feature = c("exon", "intron")
comp = c("HS6_HS0", "HS1_HS0", "HS1_HS6")
PI_type = c("PI_H", "PI_L")
title_name = list()

for (i in seq_along(feature)) {
  for (j in seq_along(comp)) {
    for (k in seq_along(PI_type)) {
      title_name <- append(title_name, list(str_c(feature[[i]], "_", comp[[j]], "_", PI_type[[k]])))
    }
  }
}

levels(ML_l_f[[3]]$type)

plot_list_test <- vector("list", length = length(title_name))

for (i in seq_along(plot_list_test)) {
  plot_list_test[[i]] <- process_the_Boruta_data(ML_l_result[[i]]) %>%
    pivot_longer(everything()) %>%
    mutate(mark = if_else(str_detect(name, "shadow")==TRUE, "shadowmarks", "marks")) %>%
    ggplot(aes(x = fct_reorder(name, value, median), y = value)) +
    geom_boxplot(aes(fill = mark)) +
    ggtitle(title_name[[i]]) +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 12, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
          axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, vjust = 0.5)) 
  
  
}
?ggarrange
library(ggpubr)
plot_list_test_v1 <- plot_list_test[1:6]
plot_list_test_v2 <- plot_list_test[7:12]
plot_list_test_v2[[1]]
ML_pairwise_exon_sig_ctrl <- ggarrange(plotlist=plot_list_test_v1, ncol = 2, nrow = 3)
ML_pairwise_intron_sig_ctrl <- ggarrange(plotlist=plot_list_test_v2, ncol = 2, nrow = 3)
ggsave("./analysis_output/ML_pairwise_exon_sig_ctrl.jpeg", ML_pairwise_exon_sig_ctrl, width = 400, height = 320, units = c("mm"), dpi = 320)
ggsave("./analysis_output/ML_pairwise_intron_sig_ctrl.jpeg", ML_pairwise_intron_sig_ctrl, width = 400, height = 320, units = c("mm"), dpi = 320)


plot_list_test[[4]]
process_the_Boruta_data(ML_l_f_combined_res[[6]]) %>%
  pivot_longer(everything()) %>%
  mutate(mark = if_else(str_detect(name, "shadow")==TRUE, "shadowmarks", "marks")) %>%
  ggplot(aes(x = fct_reorder(name, value, median), y = value)) +
  geom_boxplot(aes(fill = mark)) +
  #ggtitle(title_name[[i]]) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 12, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, vjust = 0.5))
ML_l_f_1_combined_res

ML_l_f_combined <- vector("list", length = length(c(1:6)))
ML_l_f_combined_res <- vector("list", length = length(c(1:6)))
for (i in seq_along(ML_l_f_combined)) {
  ML_l_f_combined[[i]] <- bind_rows(as_tibble(ML_l_f[[2*(i-1)+1]]), as_tibble(ML_l_f[[2*(i-1)+2]])) %>%
    mutate(type = str_replace(type, "PI_H", "major")) %>%
    mutate(type = str_replace(type, "PI_L", "major")) %>%
    mutate(type = as.factor(type)) %>%
    as.data.frame()
  ML_l_f_combined_res[[i]] <- Boruta(type~.,data=ML_l_f_combined[[i]], doTrace = 2)
}

ML_l_f_combined[[i]] <- bind_rows(as_tibble(ML_l_f[[2*(i-1)+1]]), as_tibble(ML_l_f[[2*(i-1)+2]])) %>%
  mutate(type = str_replace(type, "PI_H", "major")) %>%
  mutate(type = str_replace(type, "PI_L", "major")) %>%
  mutate(type = as.factor(type)) %>%
  as.data.frame()