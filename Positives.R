setwd("/Users/leahsamuels/Box/EIS/Field-manuscript/Code")
source("Field_study.R")
source("../../Enclosure-manuscript/Code/Enclosure_study.R")

################################################################################
###
### Leah Samuels 2/19/2025
### Script for Comparison of three monitoring methods in situ for threatened
### species detection
###
### known-positives
###
################################################################################

library('lmerTest')

camera_positives <- read.csv("camera_positives.csv")

################################################################################
### add eDNA samples that were confirmed positive via camera

for (i in 1:nrow(camera_positives)){
  positives[nrow(positives) +1,] <- qPCR[which(qPCR$id[] == camera_positives$id[i]),]
}

for(i in 1:nrow(camera_positives)){
  positives$type[which(positives$id == 
                         camera_positives$id[i])] <- camera_positives$id_date[i]
}

##############################################################################
#### Adjusting format, adding location of eDNA, change date format, add hrs

pos_details <- read.csv("Positive Sample Details_TOS.csv")

###split by date
positives[c('pos','time')] <- str_split_fixed(positives$type, '\\+', 2)
positives <- positives[,-which(colnames(positives) == 'pos')]

for (i in 1:nrow(positives)){
  if (positives$time[i] == ""){
    positives$time[i] <- 0
  }
  else{
    positives$time[i] <- gsub("d", "", positives$time[i])
  }}

positives$pos[] <- gsub('([0-9])','',positives$pos[])
positives$time <- as.numeric(positives$time)

positives['int_sample'] <- positives$time * 24

positives['location'] <- NA
positives['location2'] <- NA 


##G - ground
## B - burrow
for (i in 1:nrow(positives)){
  if(grepl("ABP",positives$id[i])){
    positives$location2[i] <- "G"
  }
}

list <- c('ABP_IT1_1211','ABP_IT3_1211','ABP_EIS1+2d_1211','ABP_IT1_0104'
          ,'ABP_EIS1+1d_0104','ABP_EIS1_0116','ABP_EIS2_0221', 'ABP_EIS2+1d_0221',
          'ABP_IT2_1130','ABP_EIS1+1d_0103','ABP_EIS2+1d_1130',"ABP_B12_0224")

for(i in 1:length(list)){
  positives$location2[which(positives$id[] == list[i])] <- "B"
}


pos_details$Sample.ID[] <- gsub("\\-", "\\_", pos_details$Sample.ID[])
pos_details$Sample.ID[] <- gsub("\\/", "", pos_details$Sample.ID[])
pos_details$Repeat.Sample.ID[] <- gsub("\\-", "\\_", pos_details$Repeat.Sample.ID[])
pos_details$Repeat.Sample.ID[] <- gsub("\\/", "", pos_details$Repeat.Sample.ID[])


positives$date <- as.character(positives$date)
num <- strsplit(positives$date,"")
positives['one'] <- NA
positives['Yr'] <- NA
for (i in 1: nrow(positives)){
  positives$one[i] <- num[[i]][1]
  if(positives$one[i] == '0'){
    positives$Yr[i] <- '2024'
  }
  else if (positives$one[i] == '1'){
    positives$Yr[i] <- '2023'
  }
}

positives['date2'] <- paste(positives$date,positives$Yr, sep = "")
positives <- positives[,-which(colnames(positives) == 'one')]
positives <- positives[,-which(colnames(positives) == 'Yr')]

for (i in 1:nrow(pos_details)){
  positives$location[which(positives$id[] == pos_details$Sample.ID[i])
                     ] <- pos_details$Sample.Location[i]
}

for (i in 1:nrow(pos_details)){
  if(!(pos_details$Repeat.Sample.ID[i] == "")){
    positives$location[which(positives$id[] == pos_details$Repeat.Sample.ID[i])
                       ] <- pos_details$Repeat.Sample.Location[i]
  }
}

positives$location[which(positives$id[] == "CNF_IT1_1221")] <- "burrow opening"


for (i in 1:nrow(positives)){
  if(!(is.na(positives$location[i]))){
    if((grepl("Burrow",positives$location[i])) | (grepl("burrow", 
                                                        positives$location[i]))){
      positives$location2[i] <- "B"
    }
    else if (grepl("ground", positives$location[i])){
      positives$location2[i] <- "G"
      
    }
  }
}

##adding type (burrow or ground) to the four samples from camera data CNF
positives$location2[which(positives$id[] == 'CNF_CB2_1221')] <- 'B'
positives$location2[which(positives$id[] == 'CNF_CB5_0228')] <- 'B'
positives$location2[which(positives$id[] == 'CNF_RF1_0228')] <- 'G'
positives$location2[which(positives$id[] == 'CNF_RF1_1221')] <- 'G'

###adding type (burrow or ground) to the FS camera data samples
for (i in 1:nrow(positives)){
  if(is.na(positives$location2[i])){
    if(grepl('F',positives$type[i])){
      positives$location2[i] <- 'G'
    } else if (grepl('B', positives$type[i])){
      positives$location2[i] <- 'B'
    }
  }
}

###ADD a sample number so that repeat samples 
positives['typeX'] <- positives$type
positives$typeX <- gsub('\\+[0-9]d', '', positives$typeX)

positives['id2'] <- paste(positives$place, positives$typeX, positives$date, sep = "_")
positives <- positives[,-which(colnames(positives) == 'typeX')]

unique <- unique(positives$id2) 
j <- 1
positives['sample_id'] <- NA
for (i in 1:length(unique)){
  positives$sample_id[which(positives$id2[] == unique[i])] <- j
  j <- j +1
  
}

##add weather
weather <- read.csv("Field_weather.csv")
colnames(weather)[3] <- "Precip"

weather$Date <- as.character(weather$Date)
for(i in 1:nrow(weather)){
  if(nchar(weather$Date[i]) == 7){
    weather$Date[i] <- paste0("0", weather$Date[i])
  }
}

weather['mean_temp'] <- NA
for (i in 1:nrow(weather)){
  weather$mean_temp[i] <- mean(c(weather$Temp1[i],weather$Temp2[i]))
}

#Positives set new column 
positives['mean_temp'] <- NA
positives['precip'] <- NA
for(i in 1:nrow(positives)){
  positives$mean_temp[i] <- weather$mean_temp[which((grepl(positives$place[i], 
                                                           weather$Location)) &
  (positives$date2[i] == weather$Date))]
  positives$precip[i] <- weather$Precip[which((grepl(positives$place[i],
                    weather$Location)) & (positives$date2[i] == weather$Date))]
}

################################################################################
### Figure 2

ggplot(positives, aes(x=time, fill = as.factor(amp))
    ) + geom_histogram(position = 'dodge') + scale_fill_manual(values = c
    ( '#1a80bb','#f2c45f'), labels = c('No amplification', 'Amplification')) + 
    theme_classic() + labs(x = "Days since EIS exposure",y = "Number of samples"
    ) +  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth
    =0.5), axis.line = element_blank(),axis.text = element_text(colour = "black"
    ),  axis.title.y = element_text(size = 12,margin = margin(r=20)
    ), axis.title.x = element_text(size = 12, margin = margin(t=20)
    ),legend.position = c(0.86,0.88), plot.margin = unit(c(0.6,0.6,0.6,0.6),"cm"
    ), legend.key.size = unit(0.5, 'cm')) + guides(fill = guide_legend(
    title = "Samples with")) + scale_x_continuous(breaks = seq(0,11, by = 1)
      ) + scale_y_continuous(breaks = seq(0,40, by=5))


###############################################################################

positives2 <- positives[which(positives$quant[] > 0),]
dotchart(positives2$quant, main = "QUANT")
hist(positives2$quant)

###generalized linear modeling
positives$location2 <- as.factor(positives$location2)
positives$place <- as.factor(positives$place)
positives$sample_id <- as.factor(positives$sample_id)

#pairs(subset(positives[,c(8,15,14)]), lower.panel = panel.cor)

##Tests
###random effects structure ###
glmer1_RE <-  glmer(amp ~ int_sample + location2 + mean_temp + precip + (1|sample_id
      ), family = binomial(link = "logit"), data = positives)  
glmer2_RE <- glmer(amp ~ int_sample + location2 + mean_temp + precip + (1|place
      ), family = binomial(link = "logit"), data = positives) 
glm3_RE <- glm(amp ~ int_sample + location2 + mean_temp + precip,
      family = binomial(link = "logit"), data = positives) 
summary(glmer1_RE)
summary(glmer2_RE)
summary(glm3_RE)
AICcmodavg::AICc(glmer1_RE)
AICcmodavg::AICc(glmer2_RE)
AICcmodavg::AICc(glm3_RE)
anova(glmer2_RE, glm3_RE, test = "Chi")
anova(glmer1_RE, glm3_RE, test = "Chi")
###glmer1_RE best model for random effects = (1|sample_id)

###fixed effect structure
glmer1_FE <- glmer(amp ~ int_sample + (1|sample_id), family = binomial(
  link = "logit"), data = positives)
glmer2_FE <- glmer(amp ~ int_sample + location2 + (1|sample_id), family =
                     binomial(link = "logit"), data = positives) 
glmer3_FE <- glmer(amp ~ int_sample + mean_temp + (1|sample_id), family = 
                     binomial(link = "logit"), data = positives) 
glmer4_FE <- glmer(amp ~ int_sample + precip + (1|sample_id), family =
                     binomial(link = "logit"), data = positives) 
glmer5_FE <- glmer(amp ~ int_sample + location2 + precip + (1|sample_id), family =
                     binomial(link = "logit"), data = positives)  
glmer6_FE <- glmer(amp ~ int_sample + location2 + mean_temp + (1|sample_id
                  ), family = binomial(link = "logit"), data = positives)  
glmer7_FE <- glmer(amp ~ int_sample + precip + mean_temp + (1|sample_id), family =
                     binomial(link = "logit"), data = positives)
glmer8_FE <- glmer(amp ~ int_sample + location2 + mean_temp + precip + (1|sample_id
                    ), family = binomial(link = "logit"), data = positives)

AIC(glmer1_FE,glmer2_FE,glmer3_FE,glmer4_FE,glmer5_FE,glmer6_FE,glmer7_FE,glmer8_FE)
drop1(glmer8_FE, test = "Chi")
AICcmodavg::AICc(glmer2_FE)
AICcmodavg::AICc(glmer5_FE)
AICcmodavg::AICc(glmer6_FE)
#confint(glmer2_FE, level = 0.95, method = "boot", nsim = 500)


set.seed(500)
#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(positives), replace=TRUE, prob=c(0.7,0.3))
train <- positives[sample, ]
test <- positives[!sample, ]
test_glmer <- glmer(amp ~ int_sample + location2 + (1|sample_id), family = 
                      binomial(link = "logit"), data = train)
predicted <- predict(test_glmer, test, type="response", allow.new.levels = TRUE)

#calculate AUC
auc(test$amp, predicted)
#plot ROC
plot(roc(test$amp, predicted))

################################################################################
### Figure 3


glmer1 <- glmer(amp ~ int_sample + location2 + (1|sample_id), 
                family = binomial(link = "logit"), data = positives)
#re.G <- data.frame(location2 = rep('G', 400),
#                       int_sample = rep(seq(min(positives$int_sample), max(
#     positives$int_sample), length.out = 100), length(levels(positives$place))),
#                       place = factor(rep(levels(positives$place), each = 100)))
# #re.B <- data.frame(location2 = rep('B', 400),
#                        int_sample = rep(seq(min(positives$int_sample),
# max(positives$int_sample), length.out = 100), length(levels(positives$place))),
#                        place = factor(rep(levels(positives$place), each = 100)))
# re.G$pred <- predict(glmer1,re.G, type = "response", re.form = ~(1|place))
# re.B$pred <- predict(glmer1,re.B, type = "response", re.form = ~(1|place))

G <- data.frame(int_sample = c(seq(1,max(positives$int_sample),1)), 
                location2 = c(rep('G',max(positives$int_sample))))
B <- data.frame(int_sample = c(seq(1,max(positives$int_sample),1)),
                location2 = c(rep('B',max(positives$int_sample))))
G$pred <- predict(glmer1,G,type = "response", re.form = NA)
B$pred <- predict(glmer1,B,type = "response",  re.form = NA)


glmer1 <- glmer(amp ~ int_sample + (1|sample_id), family = binomial(link =
                                                  "logit"), data = positives)
sample <- data.frame(int_sample = c(seq(1,max(positives$int_sample),1)))
sample$pred <- predict(glmer1,sample,type = "response", re.form = NA)

data2$snake_id <- as.factor(data2$snake_id)
glmtest <- glmer(amp ~ int_sample + int_exposure + (1|snake_id), family = 
                   binomial(link = "logit"), data = data2)

hour2 <- data.frame(int_sample = c(seq(1,max(data2$int_sample),1)),
                    int_exposure = c(rep(120,max(data2$int_sample))))
hour1 <- data.frame(int_sample = c(seq(1,max(data2$int_sample),1)),
                    int_exposure = c(rep(60,max(data2$int_sample))))
min15 <- data.frame(int_sample = c(seq(1,max(data2$int_sample),1)),
                    int_exposure = c(rep(15,max(data2$int_sample))))
min2 <- data.frame(int_sample = c(seq(1,max(data2$int_sample),1)), 
                   int_exposure = c(rep(1.66,max(data2$int_sample))))
hour2$pred <- predict(glmtest,hour2,type = "response", re.form = NA)
hour1$pred <- predict(glmtest,hour1,type = "response",  re.form = NA)
min15$pred <- predict(glmtest,min15,type = "response",  re.form = NA)
min2$pred <- predict(glmtest,min2,type = "response",  re.form = NA)



plot(NULL, NULL, xlim = c(0, 270), ylim = c(0,1), las = 1, xlab =
       "Soil collection time (hrs)", ylab = "Pr (detection)")

xcoords <- min(hour2$int_sample):max(hour2$int_sample)
rangecolor <- rgb(180,180,180, alpha=80, maxColorValue=255)
polygon(x=c(xcoords,rev(xcoords)),y=c(hour2$pred, rev(min2$pred)),col=rangecolor,border=NA)

lines(hour2$int_sample,hour2$pred, col = "antiquewhite4", lty = 2)
lines(hour1$int_sample, hour1$pred, col = "antiquewhite4", lty = 2)
lines(min15$int_sample, min15$pred, col = "antiquewhite4", lty = 2)
lines(min2$int_sample, min2$pred, col = "antiquewhite4", lty = 2)

lines(G$int_sample,G$pred, col = "coral2", lwd = 2)
lines(B$int_sample, B$pred, col = "deepskyblue4", lwd = 2)

legend("topright", cex = 0.65, bty = "n",lwd = 2, lty = c(1,1,2), col = 
         c( "deepskyblue4", "coral2", "antiquewhite4"
),legend = c("Burrow", "Non-burrow", "Experimental"))

