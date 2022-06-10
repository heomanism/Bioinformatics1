########################################
#     Bioinformatics Free Project      #
#       Implemented by Min Heo         #
#     2022-05-23 ~  2022-06-10         # 
########################################

## Set Env.
library(data.table)
library(tidyverse)
library(Biostrings)
library(huxtable)
library(stringi)
library(caret)
library(doParallel)

WORK_DIRECTORY <- "/home2/mheo/Bioinformatics"
setwd(WORK_DIRECTORY)

## Save & Load R data 
#save.image(file="Free_Proj.Rdata")
#load(file="Free_Proj.Rdata")

## Data loading
# toy_data <- fread(file="toy_filtered.csv")
Data <- fread(file="Last_Final_strand_pileup35L33G_filtered.csv")

## Read filtering & CRES filtering
IDX_p_matches_filter <- which((nchar(Data$p_matches) > 50) & (Data$p_entropy > 0.8)) # 
IDX_m_matches_filter <- which((nchar(Data$m_matches) > 50) & (Data$m_entropy > 0.8)) # 

## unique(Data$chrom) # chr 1 ~ 19, X, Y, M     :: should remove chromosome M
Plus_Entropy_high_position <- Data[IDX_p_matches_filter,] %>% 
  dplyr::select(c("chrom","pos")) %>% 
  filter(chrom != "chrM") %>% 
  as.data.frame()

Minus_Entropy_high_position <- Data[IDX_m_matches_filter,] %>% 
  dplyr::select(c("chrom","pos")) %>% 
  filter(chrom != "chrM") %>% 
  as.data.frame()

############################################
## Get Motif each strand's -10 ~ +10      ##
############################################
reference_seq <-  readDNAStringSet("reference/mm39.fa")

Plus_motifs <- c()
for(i in 1:nrow(Plus_Entropy_high_position)){
  Plus_motifs[i] <- as.character(reference_seq[Plus_Entropy_high_position$chrom[i],][[1]][(Plus_Entropy_high_position$pos[i]-10):(Plus_Entropy_high_position$pos[i]+10)])
  print(i)
}

Minus_motifs <- c()
for(i in 1:nrow(Minus_Entropy_high_position)){
  Minus_motifs[i] <- as.character(reverseComplement(reference_seq[Minus_Entropy_high_position$chrom[i],][[1]][(Minus_Entropy_high_position$pos[i]-10):(Minus_Entropy_high_position$pos[i]+10)]))
  print(i)
}

## Convert thymine -> uracil
Plus_motifs <- str_replace_all(Plus_motifs,"T","U")
Minus_motifs <- str_replace_all(Minus_motifs,"T","U")

## Get -10 ~ 10
# Plus -10 ~ 10 motifs 
Plus_motifs_table <- table(Plus_motifs)
Plus_Positive_motifs_table <- sort(Plus_motifs_table,decreasing = T)

# Minus -10 ~ 10 motifs
Minus_motifs_table <- table(Minus_motifs)
Sorted_Minus_motifs_table <- sort(Minus_motifs_table,decreasing = T)

# Plus + Minus -10 ~ 10 motifs
Plus_Minus_motifs <- c(Plus_motifs , Minus_motifs)
Plus_Minus_motifs_table <- table(Plus_Minus_motifs)
Sorted_Plus_Minus_motifs_table <- sort(Plus_Minus_motifs_table,decreasing = T)

# save(Plus_Minus_motifs,file="Motifs_35L33G.Rdata")

## Get Hexamer
# Plus Strand
Plus_Hexamers <- substr(Plus_motifs,9,14)

Plus_Hexamers_table <- table(Plus_Hexamers)
Sorted_Plus_Hexamers_table <- sort(Plus_Hexamers_table,decreasing = T)

# Minus Strand
Minus_Hexamers <- substr(Minus_motifs,9,14)

Minus_Hexamers_table <- table(Minus_Hexamers)
Sorted_Minus_Hexamers_table <- sort(Minus_Hexamers_table,decreasing = T)

# Plus + Minus 
Hexamers <- c(Plus_Hexamers,Minus_Hexamers)

Hexamers_table <- table(Hexamers)
Sorted_Hexamers_table <- sort(Hexamers_table,decreasing = T)

# save(Hexamers,file="Hexamers_35L33G.Rdata")
Hexamers_table_frequency <- Sorted_Hexamers_table/sum(Hexamers_table)

Hexamers_table_frequency_df <- data.frame("Sequences"=names(Hexamers_table_frequency),
                                          "Frequency"=as.vector(round(Hexamers_table_frequency * 100,2)))

publication_table <- Hexamers_table_frequency_df[1:20,] %>% 
  as_huxtable() %>%
  theme_basic() %>% 
  set_top_border(1,everywhere) %>% 
  set_bold(1,everywhere,T) %>% 
  set_bold(everywhere,1,T) %>% 
  set_bottom_border(final(1),everywhere)
publication_table

table(substr(Hexamers,3,3))/sum(table(substr(Hexamers,3,3)))

## Make WebLogo's Input -> Process in Weblogo site 
# write.table(Plus_Minus_motifs,file="Motifs.txt",quote=F,sep='\n',row.names = F,col.names = F)

################################################
#              Motif Predictor                 #
################################################
## Step1. Data Filtering
IDX_p_matches_count_filter <- which(nchar(Data$p_matches) > 50) 
IDX_m_matches_count_filter <- which(nchar(Data$m_matches) > 50)

length(unique(IDX_p_matches_count_filter,IDX_m_matches_count_filter)) ## 928,288



Data_for_Motif_Prediction <- Data[unique(IDX_p_matches_count_filter,IDX_m_matches_count_filter),] %>% 
  as.data.frame() %>%
  filter(chrom != "chrM") %>%
  mutate("p_motif"=as.factor(ifelse(p_entropy > 0.8,"1","0")),
         "m_motif"=as.factor(ifelse(m_entropy > 0.8,"1","0")))

## Step2. Making Data for train and test  : Unbalanced Data 
# Non_motif Sampling : 모티프가 아닌 샘플이 너무 많다. 
set.seed(1234)
Idx_Non_motif <- sample(which(Data_for_Motif_Prediction$p_entropy < 0.8 | Data_for_Motif_Prediction$m_entropy < 0.8), 33000 , replace = FALSE, prob = NULL)
  
Sampled_Non_motif_Data <- Data_for_Motif_Prediction[Idx_Non_motif,] %>% 
  dplyr::select(c("chrom","pos")) %>% 
  mutate("motif"=0)

Entropy_high_position <- rbind(Plus_Entropy_high_position,Minus_Entropy_high_position)

Sampled_motif_Data <- Entropy_high_position %>% 
  mutate("motif"=1)

Final_Data_for_prediction <- rbind(Sampled_Non_motif_Data,Sampled_motif_Data)
Final_Data_for_prediction$motif <- as.factor(Final_Data_for_prediction$motif)

ExtractSequence <- function(chromosome_position,base_range){
  prediction_motifs_sequence <- c()
  for(i in 1:nrow(Final_Data_for_prediction)){
    prediction_motifs_sequence[i] <- as.character(reference_seq[Final_Data_for_prediction$chrom[i],][[1]][(Final_Data_for_prediction$pos[i]-base_range-1):(Final_Data_for_prediction$pos[i]+base_range+2)])
    print(i)
  }
  Final_Data_for_prediction$sequence <- prediction_motifs_sequence
  
  return(Final_Data_for_prediction)
}
Final_Data_for_prediction_1 <- ExtractSequence(chromosome_position,1)
#fwrite(Final_Data_for_prediction_1,file="Final_Data_for_prediction_1.csv",col.names = T)
# Final_Data_for_prediction_1 <- as.data.frame(fread(file="Final_Data_for_prediction_1.csv"))

# Final_Data_for_prediction_3 <- ExtractSequence(chromosome_position,3)
# fwrite(Final_Data_for_prediction_3,file="Final_Data_for_prediction_3.csv",col.names = T)
# Final_Data_for_prediction_5 <- ExtractSequence(chromosome_position,5)
# fwrite(Final_Data_for_prediction_5,file="Final_Data_for_prediction_5.csv",col.names = T)
# Final_Data_for_prediction_7 <- ExtractSequence(chromosome_position,7)
# fwrite(Final_Data_for_prediction_7,file="Final_Data_for_prediction_7.csv",col.names = T)

Sequence_to_data_frame <- function(Final_Data_for_prediction){
  sequence_to_matrix <- sapply(Final_Data_for_prediction$sequence,function(x){
    str_split(x,'')
  })
  sequence_to_df <- as.data.frame(do.call(rbind,sequence_to_matrix)) %>% 
    mutate("Motif"=Final_Data_for_prediction$motif)
  return(sequence_to_df)
}

## Step3. Making Model
Sequence_Predict_Data_1 <- Sequence_to_data_frame(Final_Data_for_prediction_1)
# Sequence_Predict_Data_3 <- Sequence_to_data_frame(Final_Data_for_prediction_3)
# Sequence_Predict_Data_5 <- Sequence_to_data_frame(Final_Data_for_prediction_5)
# Sequence_Predict_Data_7 <- Sequence_to_data_frame(Final_Data_for_prediction_7)

# fwrite(Seuquence_Predict_Data,file="Sequence_Prediction_Data.csv",quote = F,row.names = F,col.names = F)
# Seuquence_Predict_Data <- fread(file="Sequence_Prediction_Data.csv") %>% 
#   as.data.frame()
# Seuquence_Predict_Data$V7 <- as.factor(Seuquence_Predict_Data$V7)
# colnames(Seuquence_Predict_Data) <- c("V1","V2","V3","V4","V5","V6","Motif")

## Step4. 10-fold Cross Validation  
trControl <- trainControl(method  = "cv",
                          number  = 10,
                          savePredictions = T,
                          classProbs=T) #summaryFunction=twoClassSummary

Sequence_Predict_Data_1$Motif <- as.factor(ifelse(Sequence_Predict_Data_1$Motif == "1","Motif","Not"))

# KNN
cl <- makePSOCKcluster(32)
registerDoParallel(cl)

knn_fit <- train(Motif ~ .,
             method     = "knn",
             tuneGrid   = expand.grid(k = 4),
             trControl  = trControl,
             metric     = "Accuracy",
             data       = Sequence_Predict_Data_1)

stopCluster(cl)

# Random Forest
cl <- makePSOCKcluster(32)
registerDoParallel(cl)

mtry <- sqrt(ncol(Seuquence_Predict_Data))
tunegrid <- expand.grid(.mtry=mtry)
rf_default_fit <- train(Motif~., 
                    data=Sequence_Predict_Data_1, 
                    method='rf', 
                    metric='Accuracy', 
                    tuneGrid=tunegrid, 
                    trControl=trControl)

stopCluster(cl)

# SVM Linear
cl <- makePSOCKcluster(32)
registerDoParallel(cl)

svm_linear_fit <- train(Motif~.,
                        data=Sequence_Predict_Data_1, 
                        method = "svmLinear",
                        trControl = trControl)

stopCluster(cl)

# Logistic Regression
cl <- makePSOCKcluster(32)
registerDoParallel(cl)

logistic_fit <- train(Motif~.,
                      data = Sequence_Predict_Data_1,
                      trControl = trControl,
                      method = "glm",
                      family=binomial())

stopCluster(cl)

# SVM Radial
cl <- makePSOCKcluster(32)
registerDoParallel(cl)

svm_radial_fit <- train(Motif~.,
                        data=Sequence_Predict_Data_1,
                        method = "svmRadial",
                        trControl = trControl)

stopCluster(cl)
# save(svm_radial_fit,file="svm_radial_fit.Rdata")
# load(file="svm_radial_fit.Rdata")

# gbm : stochastic gradient boosting
cl <- makePSOCKcluster(32)
registerDoParallel(cl)

gbm_fit <- train(Motif~.,
                 data=Sequence_Predict_Data_1,
                 distribution="bernoulli",
                 method="gbm",
                 trControl=trControl)

stopCluster(cl)
# save(gbm_fit,file="gbm_fit.Rdata")
# load(file="gbm_fit.Rdata")

#################################################
#              Model Evaluation                 #
#################################################

## Accuracy 
All_Accuracy_Results <- data.frame("Accuracy"=c(knn_fit$results$Accuracy,
                        rf_default_fit$results$Accuracy,
                        svm_linear_fit$results$Accuracy,
                        logistic_fit$results$Accuracy,
                        svm_radial_fit$results$Accuracy[3],
                        gbm_fit$results$Accuracy[9]),"methods"=c("knn","Randomforest","SVM_Linear","Logistic","SVM_Radial","GBM"))


ggplot(All_Accuracy_Results,aes(x=reorder(methods,-Accuracy),y=Accuracy,col=methods,fill=methods))+
  geom_bar(stat='identity')+
  theme_classic()+
  geom_text(aes(label=round(Accuracy,3)),vjust=0)+
  labs(title= "Accuracy with multiple methods", x = "Methods", y = "Accuracy")

## ROC
Make_Sense_Specifi_DF <- function(model_df){
  tmp <- model_df
  
  tmp_sense <- sapply(seq(0, 1, by=0.05),function(x){
    predicted_values <-ifelse(tmp$Motif > x,1,0)
    actual_values<-ifelse(tmp$obs=="Motif",1,0)
    u <- union(predicted_values, actual_values)
    conf_matrix<-table(factor(predicted_values, u), factor(actual_values, u))
    sensitivity(conf_matrix)
  })
  
  tmp_specifi <- sapply(seq(0, 1, by=0.05),function(x){
    predicted_values <-ifelse(tmp$Motif > x,1,0)
    actual_values<-ifelse(tmp$obs=="Motif",1,0)
    u <- union(predicted_values, actual_values)
    conf_matrix<-table(factor(predicted_values, u), factor(actual_values, u))
    specificity(conf_matrix)
  })
  
  tmp_df <- data.frame("Sensitivity"=tmp_sense,"Specificity"=tmp_specifi)
  return(tmp_df)
}

# KNN
KNN_ROC <- Make_Sense_Specifi_DF(knn_fit$pred)
KNN_ROC$method <- "knn"

# Randomforest
RF_DF <- rf_default_fit$pred[rf_default_fit$pred$mtry==rf_default_fit$bestTune$mtry,]
RF_ROC <- Make_Sense_Specifi_DF(RF_DF)
RF_ROC$method <- "Randomforest"

# SVM Linear
SVM_L_DF <- svm_linear_fit$pred[svm_linear_fit$pred$C==1,] 
SVM_L_ROC <- Make_Sense_Specifi_DF(SVM_L_DF)
SVM_L_ROC$method <- "SVM_Linear"

# Logisitc Regression
Logistic_DF <- logistic_fit$pred
Logistic_ROC <- Make_Sense_Specifi_DF(Logistic_DF)
Logistic_ROC$method <- "Logistic"

# SVM Radial
SVM_R_DF <- svm_radial_fit$pred[svm_radial_fit$pred$C==1 & svm_radial_fit$pred$sigma==svm_radial_fit$bestTune$sigma,] 
SVM_R_ROC <- Make_Sense_Specifi_DF(SVM_R_DF)
SVM_R_ROC$method <- "SVM_Radial"

# gradient boosting model
GBM_DF <- gbm_fit$pred[gbm_fit$pred$n.trees==gbm_fit$bestTune$n.trees & gbm_fit$pred$interaction.depth==gbm_fit$bestTune$interaction.depth &
                         gbm_fit$pred$shrinkage==gbm_fit$bestTune$shrinkage & gbm_fit$pred$n.minobsinnode==gbm_fit$bestTune$n.minobsinnode,] 
GBM_ROC <- Make_Sense_Specifi_DF(GBM_DF)
GBM_ROC$method <- "GBM"

All_ROC_results <- rbind(KNN_ROC,RF_ROC,SVM_L_ROC,Logistic_ROC,SVM_R_ROC,GBM_ROC)

ggplot(All_ROC_results,aes(x=1-Specificity,y=Sensitivity,col=method,linetype=method))+
  geom_line(size = 1, alpha = 0.7)+
  theme_classic() +
  labs(title= "ROC curve with multiple methods", x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)")
  


##################################################################
#                       Trouble Shooting                         #
##################################################################

######################################
#   Bam file로 motif counting을 함   #
######################################
## Get Hexamer 
# library(Rsamtools)
# library(GenomicRanges)
# 
# .unlist <- function (x){
#   
#   ## do.call(c, ...) coerces factor to integer, which is undesired
#   x1 <- x[[1L]]
#   if (is.factor(x1)) {
#     structure(unlist(x), class = "factor", levels = levels(x1))
#   }else{
#     
#     do.call(c, x)
#   }
# }
# 
# GetHexamers <- function(bamfile,Entropy_high_position){
#   bampoint <- BamFile(bamfile)
#   Variation_names <- paste0(Entropy_high_position$chrom,"_",Entropy_high_position$pos)
#   print(paste0("# of Total Variation : ",nrow(Entropy_high_position)))
#   Hexamers_list <- list()
#   for(i in 1:nrow(Entropy_high_position)){
#     param <- ScanBamParam(which=GRanges(Entropy_high_position$chrom[i], 
#                                         IRanges(Entropy_high_position$pos[i]-2,Entropy_high_position$pos[i]+3,width=6)),
#                           what=c("strand","pos","qwidth","seq"))
#     bam <- scanBam(bampoint, param=param)
#     bam <- unname(bam) # names not useful in unlisted result
#     elts <- setNames(bamWhat(param), bamWhat(param))
#     lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))
#     bam_df <- as.data.frame(lst)
#     
#     Hexamers_per_variation <- sapply(1:nrow(bam_df),function(x){
#       substr(bam_df$seq[x],(Entropy_high_position$pos[i]-2-bam_df$pos[x]+1),(Entropy_high_position$pos[i]+3-bam_df$pos[x]+1))
#     })
#     Hexamers_per_variation <- gusb(' ','',Hexamers_per_variation)
#     Hexamers_per_variation <- Hexamers_per_variation[which(nchar(Hexamers_per_variation) == 6)]
#     Hexamers_list[[i]] <- Hexamers_per_variation
#     print(i)
#   }
#   names(Hexamers_list) <- Variation_names
#   return(Hexamers_list)
# }
# 
# Results <- GetHexamers(bamfile='CLIP-35L33G.bam',
#                        Entropy_high_position=Entropy_high_position)
# 
# #save(Results,file="Hexamers_list.Rdata")
# # load(file="Hexamers_list.Rdata")
# 
# All_hexamers <- unlist(Results)
# All_hexamers_toy <- gsub(' ','',All_hexamers)
# Real <- All_hexamers_toy[which(nchar(All_hexamers_toy) == 6)]
# tmp <- table(Real)
# sort(tmp,decreasing=T)

