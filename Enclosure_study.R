library('stringr') ##analysis
library('graphics') ## graphing
library('lattice') ##xyplot
library('lme4') ##glmer
library('pROC') ##ROC/AUC
library('car') ##vif
library('lubridate') ##data
library('dplyr') ##data
library('AICcmodavg') ##AICc
library('ggplot2')
library('lmerTest')
library('patchwork')


data <- read.csv("Enclosure_qPCR.csv")
schedule <- read.csv("Enclosure_schedule.csv")
weather <- read.csv("Enclosure_weather.csv")

###first and last day of the weather loggers doesn't have full weather and needs to be adjusted

### nothing on last day and can be taken out
weather <- weather[-nrow(weather),]
###first day corrections (effects 4 data points)
### changing lows and averages to the noaa station
weather$min_temp[which(weather$date == "08/18/2023")] <- 24.5
weather$average_temp[which(weather$data == "08/18/2023")] <- 28.88
### changing light to averages across whole experiment
weather$average_light[which(weather$date == "08/18/2023")] <- mean(weather$average_light[2:nrow(weather)])
weather$max_light[which(weather$date == "08/18/2023")] <- mean(weather$max_light[2:nrow(weather)])
###changing inches to mm
weather$mm_rain <- weather$inches_rain * 25.4

data2 <- data.frame( id = data$id_tag, quant = data$mean_quant)

### getting rid of "filter" in title
data2[c('snake_id', 'quadrat', 'sample_time', 'method')] <- str_split_fixed(data2$id, '_', 4)
data2 <- data2[,-which(colnames(data2) == 'method')]
data2 <- data2[,-which(colnames(data2) == 'id')]
data2['id'] <- paste(data2$snake_id,data2$quadrat,data2$sample_time, sep = "_")
data2 <- data2[,c(5,1,2,3,4)]

###get rid of duplicates in data set
list <- which(duplicated(data2$id))
for (i in 1:length(list)){
  list <- which(duplicated(data2$id))
  list2 <- which(data2$id[] == data2$id[list[1]])
  if (data2$quant[list2[1]] <= data2$quant[list2[2]]){
    data2 <- data2[-list2[1],]
  }else{
    data2 <- data2[-list2[2],]
  }
}


### organizing and adding exposure time to data frame
data2['exposure_time'] <- NA
data2['time_taken'] <- NA
colnames(schedule) <- c("quadrat","exposure_time","snake", "date", "time", "1 HR"
                        ,"4 HR", "24 HR", "48 HR", "72 HR", "6 DAY", "7 DAY", 
                        "9 DAY", "10 DAY", "NEG")

for (i in 1:length(data2$quadrat)){
  for (k in 1:length(schedule$quadrat)){
    if ((data2$quadrat[i] == schedule$quadrat[k])&(data2$sample_time[i] != "NEG")){
      data2$exposure_time[i] <- schedule$exposure_time[k]
    }
  }
}

###adding the date soil sample was taken to data frame and amp column (Y, N)

data2$sample_time[which(data2$sample_time == "10 DAYS")] <- "10 DAY"
data2$sample_time[which(data2$sample_time == "6 DAYS")] <- "6 DAY"
for (i in 1:length(data2$quadrat)){
  data2$time_taken[i] <- schedule[
    which(schedule$quadrat[] == data2$quadrat[i]), which(colnames(schedule) == data2$sample_time[i])]
}


data2['amp'] <- 0
data2$amp[(which(data2$quant > 0,))] <- 1

###adding integers for sample time and exposure time
data2['int_sample'] <- NA
data2['int_exposure'] <- NA
data2$int_exposure[which(data2$exposure_time[] == "2 hour(s)")] <- 120
data2$int_exposure[which(data2$exposure_time[] == "1 hour(s)")] <- 60
data2$int_exposure[which(data2$exposure_time[] == "15 minutes")] <- 15
data2$int_exposure[which(data2$exposure_time[] == "100 seconds")] <- 1.66

timing <- data.frame(int = c(0,1,4,24,48,72,6*24,7*24,9*24,10*24), exp = c(
  "NEG","1 HR", "4 HR","24 HR","48 HR","72 HR","6 DAY", "7 DAY", "9 DAY", "10 DAY"))
for (i in 1:length(data2$sample_time)){
  data2$int_sample[i] <- timing$int[which(timing$exp[] == data2$sample_time[i])]
}

##correcting for curve dilution issues with dsDNA and getting rid of negatives
data2 <- data2[-which(data2$sample_time == "NEG"),]

data2$quant <- data2$quant/2

##adding environmental covariates
data2['max_light'] <- NA
data2['max_temp'] <- NA
data2['average_light'] <- NA
data2['average_temp'] <- NA
data2['mm_rain'] <- NA
data2['average_rain'] <- NA


data2$time_taken <- mdy(data2$time_taken)
weather$date <- mdy(weather$date)

for (i in 1:length(data2$id)){
  enddate <- data2$time_taken[i] - days(as.integer(data2$int_sample[i]/24))
  list <- seq(enddate, data2$time_taken[i], by="days")
  total <- weather %>% filter(date %in% list)
  data2$max_temp[i] <- max(total$max_temp) ##mean or max
  data2$max_light[i] <- max(total$max_light) ##mean or max
  data2$average_light[i] <- mean(total$average_light)
  data2$average_temp[i] <- mean(total$average_temp)
  data2$mm_rain[i] <- sum(total$mm_rain)
  data2$average_rain[i] <- mean(total$mm_rain)
}

data2['day_rain'] <- NA
for (i in 1:length(data2$day_rain)){
  for (j in 1:length(weather$mm_rain)){
    if (data2$time_taken[i] == weather$date[j]){
      data2$day_rain[i] <- weather$mm_rain[j]
    }
  }
}

list <- c(0,0.29,3)
data2['rain_ranked'] <- NA
for (i in 1:length(data2$mm_rain)){
  for (j in 1:(length(list)-1)){
   if (between(data2$mm_rain[i], list[j], list[j+1])){
      data2$rain_ranked[i] <- j - 1
    }
  }
}


###amplification data - including information for amp table
amp <- data2[(which(data2$quant > 0,)),]
amp2 <- amp
amp2$sample_time[(which(amp2$sample_time[] == "6 DAY"))] <- "7 DAY"
my_table <- table(amp2$sample_time, amp2$exposure_time)
matrix_table <- as.matrix(my_table)
matrix_table <- matrix_table[c(1,4,3,5,7,6,2), c(2,3,1,4)]


###outliers
matrix_table
dotchart(data2$quant, main = "QUANT")
data3 <- data2[which(data2$quant < 20),]
dotchart(data3$quant,groups = data3$int_sample, main = "QUANT")


##z scores for outliers
z_scores <- (data2$quant-mean(data2$quant))/sd(data2$quant)
z_scores[which(z_scores[] > 3 | z_scores < -3)] ## z-scores of outliers
data2$quant[which(z_scores[] > 3 | z_scores < -3)] ##quant of outliers


##snakes comparison of quant
snake_df <- data.frame(table(amp$snake_id))
snake_df2 <- aggregate(cbind(quant) ~ snake_id, data= data2,FUN = median)
snake_df2 <- aggregate(cbind(quant) ~ snake_id, data= data3,FUN = mean)
snake_df ['average_quant'] <- snake_df2$quant

### Distribution - not normal
histogram( ~ quant | snake_id, type ="count", data = data3, layout = c(4,2))
hist(data3$quant)

kruskal.test(quant ~ snake_id, data = data2)

#visualize snake comparison
xyplot(quant ~ int_sample | factor(snake_id), type=c("p","r"), pch = 19,
       xlab = "Sample time taken (hr)", col = 1, ylab = "Quant", layout = c(4,2), data = data3)
xyplot(quant ~ int_exposure | factor(snake_id), type=c("p","r"), pch = 19,
       xlab = "Exposure time(min)", col = 1, ylab = "Quant", data = data3)



###quant distrubtion without outliers - FIGURE 3
plot1 <- ggplot(data3, aes(x=factor(exposure_time, level = c("100 seconds", "15 minutes",
         "1 hour(s)", "2 hour(s)")), y=quant)) + geom_violin(fill = "darkgray") + xlab("Time in enclosure")+ ylab("EIS mtDNA copies/rxn") +  
          theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_blank(),
          axis.text = element_text(colour = "black")) + geom_jitter(size = 0.9, position = position_jitter(seed = 1, width = 0.2))
                                                                                                                                                                                                                            
plot2 <- ggplot(data3, aes(x=factor(sample2, level = c("1 HR", "4 HR", "24 HR", 
          "48 HR", "72 HR", "7 DAY", "10 DAY")), y=quant)) + geom_violin(fill = "darkgray") + xlab("Soil collection"
          )+ theme_classic() + theme( axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
          axis.title.y=element_blank(), axis.text = element_text(colour = "black"), axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) + geom_jitter(size = 0.9, position = position_jitter(
          seed = 1, width = 0.2))
plot1 + plot2



###collinearity check
panel.cor <- function(x, y){
  usr <- par("usr")
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  text(0.5, 0.5, txt) 
}

pairs(subset(data2[,c(2,9,10)]), lower.panel = panel.cor)
pairs(subset(data2[,c(8,9,10,11,12,15)]),lower.panel = panel.cor)
pairs(subset(data2[,c(8,9,10,13,14,15,16,17,18)]), lower.panel = panel.cor)

cor(data2[,8:16], method = "pearson")


##correcting for scale. glmer_1 had issues. scaled data had same AIC and AUC as unscaled
scale_data <- data2
scale_data[,c(9:18)] <- scale(subset(scale_data[,c(9:18)]))


###modeling
###base model first
glmer_1 <- glmer(amp ~ int_exposure + int_sample + int_exposure:int_sample + (1 | snake_id), family = binomial(link = "logit"),
                 data = scale_data)
glmer_2 <- glmer(amp ~ int_exposure + int_sample + (1 | snake_id), family = binomial(link = "logit"), data = data2)
glmer_3 <- glmer(amp ~ int_exposure + (1| snake_id) , family = binomial(link = "logit"), data = scale_data)
glmer_4 <- glmer(amp ~  int_sample + (1|snake_id), family = binomial(link = "logit"), data = scale_data)
##confidence intervals
#(glmer_2, level = 0.95, method = "boot", nsim = 500)


## light and temp
glmer_1 <- glmer(amp ~ int_exposure + int_sample + int_exposure*int_sample + average_light + average_temp + (1 | snake_id),
                 family = binomial(link = "logit"), data = scale_data)
glmer_2 <- glmer(amp ~ int_exposure + int_sample + average_temp + average_light + (1|snake_id), family = binomial(link = "logit"), 
                 data = scale_data)
glmer_3 <- glmer(amp ~ int_exposure + int_sample + int_exposure*int_sample + average_temp + (1|snake_id), family = binomial(link = "logit"),
                 data = scale_data)
glmer_4 <- glmer(amp ~ int_exposure + int_sample + int_exposure*int_sample + average_light + (1|snake_id), family = binomial(link = "logit"),
                 data = scale_data)
##top model
glmer_5 <- glmer(amp ~ int_exposure + int_sample  + average_temp + average_light + (1|snake_id), family = binomial(link = "logit"),
                 data = scale_data)
glmer_6 <- glmer(amp ~ int_exposure + int_sample + int_exposure*int_sample + (1|snake_id) , family = binomial(link = "logit"),
                 data = scale_data)
glmer_7 <- glmer(amp ~ int_exposure + int_sample +  average_temp + (1|snake_id), family = binomial(link = "logit"), data = scale_data)
glmer_8 <- glmer(amp ~ int_exposure + int_sample + average_light + (1|snake_id), family = binomial(link = "logit"), data = scale_data)
glmer_9 <- glmer(amp ~ int_exposure + int_sample + (1|snake_id), family = binomial(link = "logit"), data = scale_data)

### rain and light
glmer_1 <- glmer(amp ~ int_exposure + int_sample + int_exposure*int_sample + average_light + average_rain + (1 | snake_id),
                 family = binomial(link = "logit"), data = scale_data)
glm_1 <- glm(amp ~ int_exposure + int_sample + int_exposure*int_sample + average_light + average_rain, family = binomial(link = "logit"), 
             data = scale_data)
glmer_2 <- glmer(amp ~ int_exposure + int_sample + average_light + average_rain + (1 | snake_id), family = binomial(link = "logit"), 
                 data = scale_data)
glmer_3 <- glmer(amp ~ int_exposure + int_sample + int_exposure*int_sample + average_rain + (1 | snake_id), family = binomial(link = "logit"), 
                 data = scale_data)
glmer_4 <- glmer(amp ~ int_exposure + int_sample + average_light + int_exposure*int_sample + (1 | snake_id), family = binomial(link = "logit"),
                 data = scale_data)
##top
glmer_5 <- glmer(amp ~ int_exposure + int_sample + average_rain + (1 | snake_id), family = binomial(link = "logit"), data = scale_data)
##top
glmer_6 <- glmer(amp ~ int_exposure + int_sample + average_light + (1 | snake_id), family = binomial(link = "logit"), data = scale_data)
glmer_7 <- glmer(amp ~ int_exposure + int_sample + int_exposure*int_sample + (1 | snake_id), family = binomial(link = "logit"),
                 data = scale_data)
glmer_8 <- glmer(amp ~ int_exposure + int_sample + (1 | snake_id), family = binomial(link = "logit"), data = scale_data)

##just covariates
glmer_3 <- glmer(amp ~ int_exposure + int_sample + average_rain + (1| snake_id), family = binomial(link = "logit"), data = data2)
glmer_1<- glmer(amp ~ int_exposure + int_sample + average_light + (1| snake_id), family = binomial(link = "logit"), data = data2)
glmer_2<- glmer(amp ~ int_exposure + int_sample + average_temp + (1| snake_id), family = binomial(link = "logit"), data = data2)
## confidence intervals
#confint(glmer_1, devtol = Inf)
# confint(glmer_2, level = 0.95, method = "boot", nsim = 500)
# confint(glmer_3, level = 0.95, method = "boot", nsim = 500)


AIC(glm_1,glmer_1)
anova(glm_1,glm_2,glm_3,glm_4, test  = "LRT")
drop1(glm_22, test = "Chi")
vif(glm_3)

avPlots(glm_1)

#Calculate AUC

set.seed(1)
#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(data2), replace=TRUE, prob=c(0.7,0.3))
train <- data2[sample, ]
test <- data2[!sample, ]
test_glmer <- glmer(amp ~ int_exposure + (1|snake_id) , family = binomial(link = "logit"), data = train)
predicted <- predict(test_glmer, test, type="response")

#calculate AUC
auc(test$amp, predicted)
#plot ROC
plot(roc(test$amp, test$prob))



##looking at homogeneity of variance
Fitted <- fitted(glmer_1)
Resid <- resid(glmer_1)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.2)
plot(x = Fitted, y = Resid,
     xlab = "Fitted values", 
     ylab = "Pearson Residuals")
abline(h = 0, lty = 2)

boxplot(Resid ~ data2$int_exposure)
boxplot(Resid ~ data2$int_sample)
plot(Resid ~ data2$average_temp)
plot(Resid ~ data2$average_rain)
plot(Resid ~ data2$average_light)

plot(data2$prob,Fitted)


##graph of detection probabilities - FIGURE 4 
data2$snake_id <- as.factor(data2$snake_id)
glmtest <- glmer(amp ~ int_sample + int_exposure + (1|snake_id), family = binomial(link = "logit"), data = data2)

re.hour2 <- data.frame(int_exposure = rep(120, 800),
                    int_sample = rep(seq(min(data2$int_sample), max(data2$int_sample), length.out = 100), length(levels(data2$snake_id))),
                    snake_id = factor(rep(levels(data2$snake_id), each = 100)))
re.hour1 <- data.frame(int_exposure = rep(60, 800),
                       int_sample = rep(seq(min(data2$int_sample), max(data2$int_sample), length.out = 100), length(levels(data2$snake_id))),
                       snake_id = factor(rep(levels(data2$snake_id), each = 100)))
re.min15 <- data.frame(int_exposure = rep(15, 800),
                       int_sample = rep(seq(min(data2$int_sample), max(data2$int_sample), length.out = 100), length(levels(data2$snake_id))),
                       snake_id = factor(rep(levels(data2$snake_id), each = 100)))
re.min2 <- data.frame(int_exposure = rep(1.66, 800),
                       int_sample = rep(seq(min(data2$int_sample), max(data2$int_sample), length.out = 100), length(levels(data2$snake_id))),
                       snake_id = factor(rep(levels(data2$snake_id), each = 100)))
re.hour2$pred <- predict(glmtest,re.hour2, type = "response", re.form = ~(1|snake_id))
re.hour1$pred <- predict(glmtest,re.hour1, type = "response", re.form = ~(1|snake_id))
re.min15$pred <- predict(glmtest,re.min15, type = "response", re.form = ~(1|snake_id))
re.min2$pred <- predict(glmtest,re.min2, type = "response", re.form = ~(1|snake_id))

hour2 <- data.frame(int_sample = c(seq(1,max(data2$int_sample),1)), int_exposure = c(rep(120,max(data2$int_sample))))
hour1 <- data.frame(int_sample = c(seq(1,max(data2$int_sample),1)), int_exposure = c(rep(60,max(data2$int_sample))))
min15 <- data.frame(int_sample = c(seq(1,max(data2$int_sample),1)), int_exposure = c(rep(15,max(data2$int_sample))))
min2 <- data.frame(int_sample = c(seq(1,max(data2$int_sample),1)), int_exposure = c(rep(1.66,max(data2$int_sample))))
hour2$pred <- predict(glmtest,hour2,type = "response", re.form = NA)
hour1$pred <- predict(glmtest,hour1,type = "response",  re.form = NA)
min15$pred <- predict(glmtest,min15,type = "response",  re.form = NA)
min2$pred <- predict(glmtest,min2,type = "response",  re.form = NA)

plot(NULL, NULL, xlim = c(0, 240), ylim = c(0,1), las = 1, xlab = "Soil collection time", ylab = "Pr (detection)")

list <- unique(data2$snake_id)
for (i in 1:length(list)){
  lines(pred ~ int_sample, data = subset(re.hour2, snake_id == list[i]), lwd = 1, col = rgb(1,0,0, alpha=0.15))
}
for (i in 1:length(list)){
  lines(pred ~ int_sample, data = subset(re.hour1, snake_id == list[i]), lwd = 1, col = rgb(0,1,1, alpha=0.15))
}
for (i in 1:length(list)){
  lines(pred ~ int_sample, data = subset(re.min15, snake_id == list[i]), lwd = 1, col = rgb(1,0.5,0, alpha=0.15))
}
for (i in 1:length(list)){
  lines(pred ~ int_sample, data = subset(re.min2, snake_id == list[i]), lwd = 1, col = rgb(0,0,1, alpha=0.15))
}

lines(hour2$int_sample,hour2$pred, col = "red", lwd = 2)
lines(hour1$int_sample, hour1$pred, col = "darkturquoise", lwd = 2)
lines(min15$int_sample, min15$pred, col = "orange", lwd = 2)
lines(min2$int_sample, min2$pred, col = "darkblue", lwd = 2)
legend("topright", cex = 0.65, bty = "n", lwd = 2, col = c("darkblue","orange", "darkturquoise", "red"
  ),legend = c("100 s", "15 min", "1 hr", "2 hr"), title = "Time in enclosure",)



###environmental data - FIGURE 2
par(mfrow = c(3,1))
par(mar=c(1,4.1,1,1))
par(oma = c(1,0,0,0))
xcoords <- min(weather$date):max(weather$date)
rangecolor <- rgb(30,144,255, alpha=80, maxColorValue=255)

# Plot for temperature
plot(y=c(min(0),max(weather$max_temp)),x=c(min(weather$date),max(weather$date)),
     type="n", ylab="Temp (˚C)", las = 2, xaxt = 'n')
polygon(x=c(xcoords,rev(xcoords)),y=c(weather$max_temp,weather$average_temp),col=rangecolor,border=NA)
polygon(x=c(xcoords,rev(xcoords)),y=c(weather$min_temp,weather$average_temp),col=rangecolor,border=NA)
lines(x=xcoords,y=weather$average_temp,col="black")
title (font.main = 1,"Temperature", line = -1.5, adj = 0.03, cex = 3)

# plot for light
plot(y=c(min(0),max(weather$max_light)),x=c(min(weather$date),max(weather$date)),
     type="n", ylab="Light (lux)", xaxt = 'n', yaxt = "n")
title (font.main = 1,"Light", line = -1.5, adj = 0.03)
axis(2, at = c(0,10000,20000,30000, 40000, 50000, 60000), labels = c("0", "10k", "20k", "30k","40k", "50k","60k"), las =2)
polygon(x=c(xcoords,rev(xcoords)),y=c(weather$max_light,weather$average_light),col=rangecolor,border=NA)
polygon(x=c(xcoords,rev(xcoords)),y=c(weather$min_light,weather$average_light),col=rangecolor,border=NA)
lines(x=xcoords,y=weather$average_light,col="black")

#par(mar=c(2,4.1,1,1))
# plot for rain
plot(y=c(min(0),max(weather$mm_rain)),x=c(min(weather$date),max(weather$date)),
     type="n", xaxt = "n", xlab="Date",ylab="Rainfall (mm)", las = 1)
title (font.main = 1,"Rain", line = -1.5, adj = 0.03)
axis(side = 1, at = seq(min(weather$date),max(weather$date), by =1), labels = NA)
axis(side = 1, at = seq(min(weather$date),max(weather$date), by =4), cex = 0.5,labels = c("Aug 18","Aug 22",
      "Aug 26", "Aug 30", "Sep 3", "Sep 7","Sep 11"))
lines(x=xcoords,y=weather$mm_rain,col="black")

