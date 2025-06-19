hierarchical_predict <- function(model_list, data_ori, level_num, op_thr_list, rep_years, watch_Tstart, watch_Dt, N_ratio = 0.6, decay = 0.8){
  source("D:/Code/R/NPC_screen/data_trans.R")
  colnames <- c("id", "name", "sex", "age", "status", "year","VCA_IgA", "EBNA1_IgA", "event_year")
  data_anay <- data_trans(data_ori, rep_years, colnames)
  s <- sprintf("Finish data preprocessing, %i samples for evaluation", length(unique(data_anay$id)))
  print(s)
  
  library(JMbayes2)
  library(tidyverse)
  library(survival)
  source("D:/Code/R/NPC_screen/calculate_metrics.R")
  data_ori <- data_ori[which(data_ori$event_year > watch_Tstart), ]
  preds_score <- matrix(data = NA, nrow = nrow(data_ori), ncol = level_num)
  w_matrix <- matrix(data = NA, nrow = nrow(data_ori), ncol = level_num)
  mask_N <- matrix(data = NA, nrow = nrow(data_ori), ncol = level_num)
  for (i in seq_along(1:level_num)){
    tic(paste("Finished prediction in level", as.character(i)))
    jointFit <- model_list[i * 4]$jointFit
    tr_roc <- calculate_metrics(jointFit, newdata = data_anay, Tstart = watch_Tstart, Dt = watch_Dt)
    preds_score[, i] <- tr_roc$preds
    # generate weight matrix
    w_matrix[, i]  <- as.numeric(outer(1 - tr_roc$preds, op_thr_list[i], ">="))
    toc()
  }
  
  if (level_num > 1){
    for (j in seq_along(1:nrow(data_ori))){
      w <- w_matrix[j, ]
      ini_N <- which(w == 1)[1]
      if (!is.na(ini_N)){
        ini_N <- as.numeric(ini_N)
        N_ratio <- N_ratio * decay^(ini_N - 1) 
        w <- rep((1 - N_ratio) / (level_num - 1), level_num)
        w[ini_N] <- N_ratio
      }else{
        w <- rep(1 / level_num, level_num)
      }
      w_matrix[j, ] <- w
    }
    preds_score_w <- rowSums(preds_score * w_matrix)
  }else{
    preds_score_w <- preds_score
  }
  
  Thoriz <- watch_Tstart + watch_Dt + 1e-06
  ind <- data_ori$event_year < Thoriz & data_ori$status == 1
  thrs <- seq(0, 1, length = 101)
  Check_mean <- outer(1 - preds_score_w, thrs, "<")
  nTP <- colSums(Check_mean * c(ind))
  nFN <- sum(ind) - nTP
  nFP <- colSums(Check_mean * c(1 - ind))
  nTN <- sum(1 - ind) - nFP
  TP <- nTP / sum(ind)
  FP <- nFP / sum(1 - ind)
  youden <- TP - FP
  Youden <- median(thrs[youden == max(youden)])
  
  op_thr_v <- thrs[which.min(abs(thrs - Youden))]
  check <- outer(1 - preds_score_w, op_thr_v, "<")
  
  out <- list(preds_score = preds_score, w_matrix = w_matrix, 
              preds_score_w = preds_score_w, check = check, nTP = nTP, nFN = nFN, nFP = nFP,
              nTN = nTN, thrs = thrs, op_thr_v = op_thr_v, Youden = Youden)
}