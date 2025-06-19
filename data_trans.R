data_trans <- function(data_ori, rep_years, colnames){
  rep_years <- as.numeric(rep_years)
  data_trans <- data.frame(matrix(NA, nrow = nrow(data_ori) * rep_years, ncol = length(colnames)))
  colnames(data_trans) <- colnames
  
  # 直接复制转换常规变量
  data_trans[, 1] <- rep(data_ori$id, each = rep_years)
  data_trans[, 2] <- rep(data_ori$name, each = rep_years)
 
  # 性别转换成因子变量后复制
  sex <- data_ori$sex
  sex[sex == 1] <- "male"
  sex[sex == 0] <- "female"
  data_ori$sex <- factor(sex, levels = c("male", "female"))
  data_trans[, 3] <- rep(data_ori$sex, each = rep_years)
  data_trans[, 4] <- rep(data_ori$age, each = rep_years)
  data_trans[, 5] <- rep(data_ori$status, each = rep_years)
  data_trans[, 9] <- rep(data_ori$event_year, each = rep_years)
  
  # 转换年份和两个血清指标
  h1 <- 6 + rep_years - 1
  vca <- data_ori[, 6:h1]
  vca <- vca %>% pivot_longer(cols = 1:rep_years , names_to = "year" , values_to = "VCA-IgA")%>%
    mutate(year = rep(1:rep_years, length.out = n()))
  data_trans[, 6:7] <- vca
  h1 <- h1 + 1
  h2 <- 6+2*rep_years-1
  ebn <- data_ori[, h1:h2]
  ebn <- ebn %>% pivot_longer(cols = 1:rep_years , names_to = "year" , values_to = "EBNA1-IgA")
  data_trans[,8] <- ebn$`EBNA1-IgA`
  
  # 最后一项数据预处理
  data_trans <- data_trans %>%
    mutate(seq_id = dense_rank(id))
  
  # 删除血清抗体为空的变量
  data_trans <- data_trans[!is.na(data_trans$VCA_IgA), ]
  library(dplyr)
  data_trans <- data_trans %>%
    mutate(status_label = ifelse(status == 1, "yes", "no"))
  
  data_trans$seq_id <- as.factor(data_trans$seq_id)
  data_trans$status_label <- as.factor(data_trans$status_label)
  data_trans$sex <- as.factor(data_trans$sex)
  data_trans$year<-as.numeric(data_trans$year)
  
  out <- data_trans
}