#!/usr/bin/env Rscript

library(vegan)
library(RColorBrewer)

options(bitmapType='cairo') #关闭服务器与界面的互动响应

get_color <- function(samples){
  
  #根据样本数目产生对应得颜色列表
  sam_num = length(samples)+1
  if (sam_num <= 12){
    colors <- brewer.pal(n=sam_num, name="Paired")
  }else{
    colors <- brewer.pal(n=12, name="Paired")
    pal <- colorRampPalette(colors)
    colors <- pal(sam_num)
  }
  
  return(colors)
}


plot_anosim<- function(Abundance,Group){  #anosim画图

  otu <- read.delim(Abundance, sep = '\t', header=TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  otu <- data.frame(t(otu), stringsAsFactors = FALSE) #将丰度表倒置
  otu <- otu[rowSums(otu[])>0,]
  group <- read.delim(Group, sep = '\t', stringsAsFactors = FALSE,header=TRUE)
  sample_id <- row.names(otu)
  part_group  <- data.frame(group,stringsAsFactors = FALSE)
  print(group)
  group <- part_group[which (part_group[,1] %in% sample_id),]
  colors<-get_color(unique(group$group))
  anosim_result_otu <- anosim(otu, group$group, permutations = 999,distance = 'bray') #跑anosim结果
  png(paste('anosim.all.png', sep = ''), width = 800, height = 400) #画所有分组图
  plot(anosim_result_otu, main="Anosim", xlab = "Sample grouping", ylab = "Bray-curtis distance", col = colors)
  dev.off()
  group_name <- unique(group$group)
  anosim_result_two <- NULL
  dir.create('anosim_group_by_group', recursive = TRUE)  #创建目录
  for (i in 1:(length(group_name) - 1)) {  #两两对比画图
    for (j in (i + 1):length(group_name)) {
      group_ij <- subset(group, group %in% c(group_name[i], group_name[j]))
      otu_ij <- otu[group_ij$sample, ]
      anosim_result_otu_ij <- anosim(otu_ij, group_ij$group, permutations = 999, distance = 'bray')	#随机置换检验 999 次
    
      if (anosim_result_otu_ij$signif <= 0.001) Sig <- '***'
      else if (anosim_result_otu_ij$signif <= 0.01) Sig <- '**'
      else if (anosim_result_otu_ij$signif <= 0.05) Sig <- '*'
      else Sig <- NA
      anosim_result_two <- rbind(anosim_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', anosim_result_otu_ij$statistic, anosim_result_otu_ij$signif, Sig))
    
    #每次循环输出图片
      png(paste('anosim_group_by_group/anosim.', group_name[i], '_', group_name[j], '.png', sep = ''), width = 600, height = 300)
      plot(anosim_result_otu_ij, main="Anosim", xlab = "Sample grouping", ylab = "Bray-curtis distance", col = c('gray', 'red', 'blue'))
      dev.off()
    }
    return(anosim_result_two)
  }
  
}

  anosim_result <- function(Abundance,Group){
    anosim_result_two <- plot_anosim(Abundance,Group)
    anosim_result_two <- data.frame(anosim_result_two, stringsAsFactors = FALSE)
    names(anosim_result_two) <- c('group', 'distance', 'R', 'P_value', 'Sig')
    write.table(anosim_result_two, 'anosim_group_by_group/ANOSIM.result_two.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
  }


  add_help_args <- function(args){
    
    if(length(args) != 2) {
      cat("Version: v1.0.0\n")
      cat("Author:Boya Xu, Xingguo Zhang\n")
      cat("Email:invicoun@foxmail.com\n")
      cat("Function:anosim analysis.\n")
      cat("Example:anosim_group.r abundance_species.tsv group.list\n")
      quit()
    }
  }
  
  
args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
anosim_result(args[1], args[2])
