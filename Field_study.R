setwd("/Users/leahsamuels/Box/EIS/Field-manuscript/Code")

################################################################################
###
### Leah Samuels 2/19/2025
### Script for Comparison of three monitoring methods in situ for threatened
### species detection
################################################################################

library('stringr') ##analysis

qPCR_all <- read.csv("qPCR_mastersheet_field.csv")

##### chose higher quants for repeat samples
list <- which(duplicated(qPCR_all$id_tag))
for (i in 1:length(list)){
  list <- which(duplicated(qPCR_all$id_tag))
  list2 <- which(qPCR_all$id_tag[] == qPCR_all$id_tag[list[1]])
  if (qPCR_all$mean_quant[list2[1]] <= qPCR_all$mean_quant[list2[2]]){
    qPCR_all <- qPCR_all[-list2[1],]
  }else{
    qPCR_all <- qPCR_all[-list2[2],]
  }
}

####make dataframe with qPCR data organized
qPCR <- data.frame(id = qPCR_all$id_tag, quant = qPCR_all$mean_quant)
qPCR[c('place','type','date')] <- str_split_fixed(qPCR$id, '_', 3)
qPCR['amp'] <- 0
qPCR$amp[(which(qPCR$quant > 0,))] <- 1
qPCR$quant <- qPCR$quant/2

###############################################################################
### Make a dataframe with all known-positive samples

positives <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(positives) <- colnames(qPCR)
for (i in 1:nrow(qPCR)){
  if ((grepl("EIS",qPCR$type[i])) | (grepl("IT",qPCR$type[i]))){
    positives[nrow(positives) + 1,] <- qPCR[i,]
  }
}
##add known positives from burrows and fences
short <- data.frame(matrix(ncol = 6, nrow = 0))
for (i in 1:nrow(qPCR)){
  if((grepl("d",qPCR$type[i])) & !(grepl("EIS",qPCR$type[i]))){
    short[nrow(short) + 1,] <- qPCR[i,]
  }
}

colnames(short)<- colnames(qPCR)
short[c('type2','time')] <- str_split_fixed(short$type, '\\+', 2)
short['id2'] <- paste(short$place,short$type2,short$date, sep = "_")

for (i in 1:nrow(short)){
  positives[nrow(positives)+1,]<- qPCR[which(qPCR$id[] == short$id[i]),]
  positives[nrow(positives) + 1,] <- qPCR[which(qPCR$id[] == short$id2[i]),]
  
}

###############################################################################
### Make dataframes for fences, burrows, and negatives

fence <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(fence) <- colnames(qPCR)
for (i in 1:nrow(qPCR)){
if (grepl("F",qPCR$type[i])){
  fence[nrow(fence) + 1,] <- qPCR[i,]
}
}


burrow <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(burrow) <- colnames(qPCR)
for (i in 1:nrow(qPCR)){
  if (grepl("B",qPCR$type[i])){
    burrow[nrow(burrow) + 1,] <- qPCR[i,]
  }
}

negatives <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(negatives) <- colnames(qPCR)
for (i in 1:nrow(qPCR)){
  if (grepl("NC",qPCR$type[i])){
    negatives[nrow(negatives) + 1,] <- qPCR[i,]
  }
}


sum(positives$amp) # positive 60/98
sum(burrow$amp) # positive 19/160
sum(fence$amp) # positive 3/59


