hierarchical_train <- function(train_data_ori, exd_ratio, level, NM_r = 4, watch_Tstart = 1,  watch_Dt = 6){
  # 对训练集进行数据增强
  library(tidymodels)
  library(themis)
  xiaolan_1 <- train_data_ori %>%
    select(id, sex, age, status, VCA_IgA_1,VCA_IgA_2,VCA_IgA_3,VCA_IgA_4,VCA_IgA_5,VCA_IgA_6,EBNA1_IgA_1,EBNA1_IgA_2,EBNA1_IgA_3,EBNA1_IgA_4,EBNA1_IgA_5,EBNA1_IgA_6,event_year)
  xiaolan_1$sex <- as.numeric(xiaolan_1$sex)
  xiaolan_1$status <- as.factor(xiaolan_1$status)
  res <- nearmiss(xiaolan_1, var = "status", under_ratio = NM_r)
  res$status <- as.numeric(res$status) - 1
  genData = ANS(subset(res, select = -id), res$status, dupSize = 1)
  genData <- genData$syn_data
  genData <- genData[which(!is.na(genData$sex)), ]
  # 将性别、年龄和确诊年份转换成整数
  genData$sex <- round(genData$sex)
  genData$age <- round(genData$age)
  genData$event_year <- round(genData$event_year)
  # 去掉增强样本中多余的class类
  genData <- subset(genData, select = -class)
  # 过滤到可能不重要的增强样本
  genData <- neater(subset(res, select = -id), genData, k = 6, iterations = 100, smoothFactor = 1, classAttr = "status")
  # 参考原始数据格式，补齐增强样本中缺失的单元
  genData_2 <- data.frame(matrix(NA, nrow = nrow(genData), ncol = ncol(train_data_ori)))
  colnames(genData_2) <- colnames(train_data_ori)
  genData_2$id <- seq(max(train_data_ori$id) + 1, max(train_data_ori$id) + nrow(genData))
  genData_2$name <- "genID"
  genData_2[, 3:11] <- genData[, 1:9]
  genData_2[, 17:22] <- genData[, 10:15]
  genData_2[, 28] <- genData[, 16]
  # 拼接原始数据与增强样本，用于训练
  train_all <- rbind.data.frame(train_data_ori, genData_2)
  # 转换数据格式
  source("D:/Code/R/NPC_screen/data_trans.R")
  colnames <- c("id", "name", "sex", "age", "status", "year","VCA_IgA", "EBNA1_IgA", "event_year")
  train_all_anay <- data_trans(train_all, 11, colnames)
  
  id_train <- rbind(as.matrix(res$id), as.matrix(genData_2$id))
  train_data <- train_all_anay[which(train_all_anay$id %in% id_train), ]
  train_data_id <- train_data[which(train_data$year == 1),]
  watch_data <- train_all_anay[which(!train_all_anay$name == "genID"), ]
  s <- sprintf("Finish data preprocessing, %i samples for training, %i samples for watch", length(id_train), length(unique(watch_data$id)))
  print(s)
  
  library(JMbayes2)
  library(tidyverse)
  library(survival)
  
  fm1 <- lme(VCA_IgA ~ ns(year) * sex, data = train_data,
             random = ~ ns(year) | seq_id, control = lmeControl(opt = 'optim'))
  
  fm2 <- lme(EBNA1_IgA ~ ns(year) * sex, data = train_data, 
             random = ~ ns(year) | seq_id, control = lmeControl(opt = 'optim'))
  
  fForms <- list("VCA_IgA" = ~ value(VCA_IgA),
                 "EBNA1_IgA" = ~ value(EBNA1_IgA))
  
  train_data_id$status <- as.numeric(train_data_id$status)
  CoxFit <- coxph(Surv(event_year,status)~age+sex, data=train_data_id, method="breslow")
  
  # 拟合联合模型
  jointFit <- jm(CoxFit,list(fm1,fm2), time_var = 'year',  functional_forms = fForms)  # 设定时间起点
  
  # 将这一阶段的模型打包并保存
  model_list <- list(fm1 = fm1, fm2 = fm2, CoxFit = CoxFit, jointFit = jointFit)
  
  # 测试watch_data
  source("D:/Code/R/NPC_screen/calculate_metrics.R")
  tr_roc <- calculate_metrics(jointFit, newdata = watch_data, Tstart = watch_Tstart, Dt = watch_Dt)
  op_thr <- which.min(abs(tr_roc$thrs - tr_roc$Youden))
  ppv <- tr_roc$nTP[op_thr] / (tr_roc$nTP[op_thr] + tr_roc$nFP[op_thr])
  acc <- (tr_roc$nTP[op_thr]+tr_roc$nTN[op_thr])/ (tr_roc$nTP[op_thr] + tr_roc$nFP[op_thr] + tr_roc$nTN[op_thr] + tr_roc$nFN[op_thr])
  sen <- tr_roc$nTP[op_thr]/ (tr_roc$nTP[op_thr] + tr_roc$nFN[op_thr])
  spe <- tr_roc$nTN[op_thr]/ (tr_roc$nTN[op_thr] + tr_roc$nFP[op_thr])
  op_thr_v <- tr_roc$thrs[op_thr]
  s1 <- sprintf("Prediction performance on the watch data in level %i:", level)
  s2 <- sprintf("op_thr: %f; PPV: %f; acc: %f; sensitivity: %f; specificity: %f", op_thr_v, ppv, acc, sen, spe)
  s3 <- sprintf("%i NPC cases remained, nTP: %i; nFP: %i", tr_roc$nTP[op_thr] + tr_roc$nFN[op_thr], tr_roc$nTP[op_thr], tr_roc$nFP[op_thr])
  print(s1)
  print(s2)
  print(s3)
  
  # 1% TOP的TP和FP
  Thoriz <- watch_Tstart + watch_Dt + 1e-06
  train_data_ori <- train_data_ori[which(train_data_ori$event_year > watch_Tstart), ]
  score_all <- tr_roc$preds
  Check <- outer(1 - score_all, op_thr_v, "<")
  or <- order(score_all, decreasing = TRUE)
  ind_t1 <- or[1:round(length(or) * 0.01)]
  ind <- train_data_ori$event_year < Thoriz & train_data_ori$status == 1 & Check
  ind_N <- sum(ind[ind_t1])
  s3 <- sprintf("TOP 1/100 nTP: %i; TOP 1/100 nFP: %i", ind_N, round(length(or) * 0.01) - ind_N)
  print(s3)
  # 根据测试性能过滤当前层的负样本,强制保留标签为1的病人
  watch_data_id <- watch_data[which(watch_data$year == 1), ]
  watch_data_id <- watch_data_id[which(watch_data_id$event_year > watch_Tstart), ]
  status <-  watch_data_id$status == 1
  Check <- Check | status
  score_all_r <- 1- score_all
  score_N <- score_all_r[!Check]
  id_N <- watch_data_id$id[!Check]
  # 按照exd_ratio选出被过滤负样本
  or <- order(score_N, decreasing = TRUE)
  ind_rm <- or[1:round(length(or) * exd_ratio)]
  id_rm <- id_N[ind_rm]
  
  # 把留下样本的score传输回去， 用于下一轮训练
  m <- round(length(or) * exd_ratio) + 1
  ind_leave <- or[m:length(or)]
  score_leave <- rbind(as.matrix(score_all[Check]), as.matrix(score_N[ind_leave]))
  id_leave <- rbind(as.matrix(watch_data_id$id[Check]), as.matrix(id_N[ind_leave]))
  
  out <- list(object = model_list, id_rm = id_rm, score_leave = score_leave, id_leave = id_leave, op_thr = op_thr_v)
}