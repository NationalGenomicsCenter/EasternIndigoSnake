################################################################################
###
### Leah Samuels 2/19/2025
### Script for Comparison of three monitoring methods in situ for threatened
### species detection
###
### occupancy
###
################################################################################

#install.packages('RPresence', repo='https://www.mbr-pwrc.usgs.gov/mbrCRAN')
require("RPresence")
require('ggplot2')
require('dplyr')
require('emojifont') ##easiest way to include theta in exported pdf graph
require('patchwork')

data <- read.csv("presence_input.csv")

rownames(data) <- data$X
data <- data[,-which(colnames(data) == 'X')]

###create a covariate data frame
cov <- data.frame(row.names = rownames(data), site = rep(NA,nrow(data)), 
                  location = rep(NA,nrow(data)))

covariates <- function(df){
for (i in 1:nrow(df)){
  if(grepl('FSE',rownames(df)[i])){
    df$site[i] <- 'FSE'
  } else if (grepl('FSW',rownames(df)[i])){
    df$site[i] <- 'FSW'
  } else if (grepl('ABP',rownames(df)[i])){
    df$site[i] <- 'ABP'
  }else if (grepl('CNF',rownames(df)[i])){
    df$site[i] <- 'CNF'
  }
  }

for (i in 1:nrow(df)){
  if(grepl('B',rownames(df)[i])){
    df$location[i] <- 'B'
    
  } else if (grepl('F',rownames(df)[i])){
    df$location[i] <- 'G'
  }
}
  df
}
cov = covariates(cov)


### create occupancy Pao (data structure)
pao <- createPao(
  data = data, 
  frq = rep(1,nrow(data)), 
  unitcov = cov,
  unitnames = row.names(data),
  methods = 2
)

fixedvals = data.frame(param = c('psi'), value = c(1))

### occupancy models
mmmod_1 <- occMod(    
  data = pao,          
  type = 'so.mm',        #  so.mm : static occupancy, multi-method model
  model = list(       
    psi ~ 1,             #  occupancy constant for all sites
    theta ~ 1,    
    p ~ METHOD * location),
  fixed = fixedvals

)   

mmmod_2 <- occMod(    
  data = pao,          
  type = 'so.mm',        #  so.mm : static occupancy, multi-method model
  model = list(       
    psi ~ 1,             #  occupancy constant for all sites
    theta ~ location,      
    p ~ METHOD * location),
  fixed = fixedvals
  
)   

mmmod_3 <- occMod(    
  data = pao,          
  type = 'so.mm',        #  so.mm : static occupancy, multi-method model
  model = list(       
    psi ~ 1,             #  occupancy constant for all sites
    theta ~ location,   
    p ~ METHOD),
  fixed = fixedvals
)   

mmmod_4 <- occMod(    
  data = pao,          
  type = 'so.mm',        #  so.mm : static occupancy, multi-method model
  model = list(       
    psi ~ 1,             #  occupancy constant for all sites
    theta ~ 1,   
    p ~ METHOD),
  fixed = fixedvals
)   


### occupancy model analyses
summary(mmmod_1)
summary(mmmod_2)
summary(mmmod_3)
summary(mmmod_4)

mmmod_data <- function(df){
  
ggdata <- rbind(
  unique(df$real$psi),
  unique(df$real$p),
  unique(df$real$theta)
)
return(ggdata)
}

###make data frames for each of the models with parameter values
m1 <- mmmod_data(mmmod_1)
m2 <- mmmod_data(mmmod_2) 
m3 <- mmmod_data(mmmod_3) 
m4 <- mmmod_data(mmmod_4)

theta <- coef(mmmod_1, param = 'theta', prob = 0.95)
p <- coef(mmmod_1, param = 'p', prob = 0.95)
m1$parameter <- factor (x = c('psi','eDNA burrow','eDNA ground','Camera burrow',
      'Camera ground','theta'), levels = c('psi','eDNA burrow','eDNA ground',
      'Camera burrow','Camera ground','theta'))
m1$Beta <- c('NA',p$Beta_est[1],p$Beta_est[3],p$Beta_est[2],p$Beta_est[4],theta$Beta_est)

theta <- coef(mmmod_2, param = 'theta', prob = 0.95)
p <- coef(mmmod_2, param = 'p', prob = 0.95)
m2$parameter <- factor (x = c('psi', 'eDNA burrow_Camera burrow','eDNA ground',
      'Camera ground', 'Burrow','Ground'),levels = c('psi', 
      'eDNA burrow_Camera burrow','eDNA ground','Camera ground', 'Burrow','Ground'))
m2$Beta <- c('NA',paste(p$Beta_est[1],p$Beta_est[2], sep = '_'),
      p$Beta_est[3],p$Beta_est[4],theta$Beta_est[1],theta$Beta_est[2])
fitted(object = mmmod_2, param = 'theta')


theta <- coef(mmmod_3, param = 'theta', prob = 0.95)
p <- coef(mmmod_3, param = 'p', prob = 0.95)
m3$parameter <- factor (x = c('psi', 'eDNA','Camera','Burrow','Ground'),
                           levels = c('psi', 'eDNA','Camera','Burrow','Ground'))
m3$Beta <- c('NA',p$Beta_est,theta$Beta_est)

theta <- coef(mmmod_4, param = 'theta', prob = 0.95)
p <- coef(mmmod_4, param = 'p', prob = 0.95)
m4$parameter <- factor (x = c('psi', 'eDNA','Camera','Theta'),
                        levels = c('psi', 'eDNA','Camera','Theta'))
m4$Beta <- c('NA',p$Beta_est[1],p$Beta_est[2],theta$Beta_est)



###model comparison
mod.list <- list ( mmmod_1, mmmod_2, mmmod_3)
aictable <- createAicTable(mod.list = mod.list, use.aicc = TRUE)
aictable$table
mmmod_3



####goodness of fit###########################################################

### First, let's setup the expectation against which we'll compare our observed detection histories
### Model estimated parameters
thetaB <- m3$est[which(m3$parameter== 'Burrow')]
thetaG <- m3$est[which(m3$parameter== 'Ground')]
p_camera <- m3$est[which(m3$parameter== 'Camera')]
p_edna<- m3$est[which(m3$parameter== 'eDNA')]

### List of all 2^6 (64) possible detection histories
possible_detection_histories <- expand.grid(replicate(6, 0:1, simplify = FALSE)) # stolen from StackOverflow https://stackoverflow.com/questions/48978077/generate-all-possible-binary-vectors-of-length-n-in-r
head(possible_detection_histories)

g_histories <- possible_detection_histories   
collapse_history_g <- function(x){paste(c(x, "G"), collapse = "")}              # ground detection history names, function
rownames(g_histories) <- apply(possible_detection_histories, 
                               1, collapse_history_g)

b_histories <- possible_detection_histories 
collapse_history_b <- function(x){paste(c(x, "B"), collapse = "")}              # burrow detection history names, function
rownames(b_histories) <- apply(possible_detection_histories, 
                               1, collapse_history_b)

### List of detection probabilities
detection_probabilities <- c(rep(p_edna, 3), rep(p_camera, 3))

### Matrix of detection probabilities, across all possible detection histories, predicated on availability
matrix_detection_probabilities <- matrix(, nrow = 64, ncol = 6) # empty matrix to feed histories

for(i in 1:64){
  for(j in 1:6){
    matrix_detection_probabilities[i,j] <- ifelse(possible_detection_histories[i,j] == 1,
                                                  detection_probabilities[j],        # detection
                                                  1 - detection_probabilities[j])     # non-detection
  }
}

### Matrix of possibility that the snake wasn't available for detection anyways; GROUND
g_matrix_theta_histories <- matrix(, nrow = 64, ncol = 3) # empty matrix to feed unobserved theta histories

for(i in 1:64){
  ifelse(sum(possible_detection_histories[i,c(1,4)]) < 1, # if both eDNA and camera are negative on this visit, then it might just not be available
         g_matrix_theta_histories[i,1] <- 1 - thetaG,
         g_matrix_theta_histories[i,1] <- 0)
  
  ifelse(sum(possible_detection_histories[i,c(2,5)]) < 1, # if both eDNA and camera are negative on this visit, then it might just not be available
         g_matrix_theta_histories[i,2] <- 1 - thetaG,
         g_matrix_theta_histories[i,2] <- 0)
  
  ifelse(sum(possible_detection_histories[i,c(3,6)]) < 1, # if both eDNA and camera are negative on this visit, then it might just not be available
         g_matrix_theta_histories[i,3] <- 1 - thetaG,
         g_matrix_theta_histories[i,3] <- 0)
}

### Cumulate probability
g_cum_prob <- (matrix_detection_probabilities[,1] *     # probability of observed eDNA on first visit multiplied by...
                 matrix_detection_probabilities[,4] *   # probability of observed camera on first visit multiplied by...
                 thetaG +                               # probability that the snake was available plus
                 g_matrix_theta_histories[,1]) *          # probability that the snake wasn't available, all multiplied by...
  
  (matrix_detection_probabilities[,2] *                 # pr of observed eDNA on second visit multiplied by...
     matrix_detection_probabilities[,5] *               # pr of observed camera on second visit multiplied by...
     thetaG +                                           # pr that the snake was avail plus...
     g_matrix_theta_histories[,2]) *                      # pr that the snake wasn't avail, all multiplied by...
  
  (matrix_detection_probabilities[,3] *                 # third visit...
     matrix_detection_probabilities[,6] * 
     thetaG + 
     g_matrix_theta_histories[,3])

sum(g_cum_prob) # 1!!!

### Matrix of possibility that the snake wasn't available for detection anyways; BURROW
b_matrix_theta_histories <- matrix(, nrow = 64, ncol = 3) # empty matrix to feed unobserved theta histories

for(i in 1:64){
  ifelse(sum(possible_detection_histories[i,c(1,4)]) < 1, # if both eDNA and camera are negative on this visit, then it might just not be available
         b_matrix_theta_histories[i,1] <- 1 - thetaB,
         b_matrix_theta_histories[i,1] <- 0)
  
  ifelse(sum(possible_detection_histories[i,c(2,5)]) < 1, # if both eDNA and camera are negative on this visit, then it might just not be available
         b_matrix_theta_histories[i,2] <- 1 - thetaB,
         b_matrix_theta_histories[i,2] <- 0)
  
  ifelse(sum(possible_detection_histories[i,c(3,6)]) < 1, # if both eDNA and camera are negative on this visit, then it might just not be available
         b_matrix_theta_histories[i,3] <- 1 - thetaB,
         b_matrix_theta_histories[i,3] <- 0)
}

### Cumulate probability
b_cum_prob <- (matrix_detection_probabilities[,1] *     # probability of observed eDNA on first visit multiplied by...
                 matrix_detection_probabilities[,4] *   # probability of observed camera on first visit multiplied by...
                 thetaB +                               # probability that the snake was available plus
                 b_matrix_theta_histories[,1]) *          # probability that the snake wasn't available, all multiplied by...
  
  (matrix_detection_probabilities[,2] *                 # pr of observed eDNA on second visit multiplied by...
     matrix_detection_probabilities[,5] *               # pr of observed camera on second visit multiplied by...
     thetaB +                                           # pr that the snake was avail plus...
     b_matrix_theta_histories[,2]) *                      # pr that the snake wasn't avail, all multiplied by...
  
  (matrix_detection_probabilities[,3] *                 # third visit...
     matrix_detection_probabilities[,6] * 
     thetaB + 
     b_matrix_theta_histories[,3])

sum(b_cum_prob) # 1!!!

### Compare with observed data!################################################

gof <- data[c(1,3,5,2,4,6)]
gof$location <- cov$location
gof2 <-
  gof %>%
  group_by(across(everything())) %>%
  tally(name = "Frequency")
gof2 <- as.data.frame(gof2)
gof2$Hist <- apply(gof2[1:7], 1, function(row) paste(row, collapse = ""))

g_histories['Prob'] <- g_cum_prob
b_histories['Prob'] <- b_cum_prob
histories <- rbind(g_histories, b_histories)
histories$Prob <- format(histories$Prob, scientific = FALSE)

gof2$expected_prob <- NA
for (i in 1:nrow(gof2)){
  gof2$expected_prob[i] <- 
    histories$Prob[which(rownames(histories) == gof2$Hist[i])]
}


g_other_sum <- 0
b_other_sum <- 0
for (i in 1:nrow(g_histories)){
  if (!(rownames(g_histories)[i] %in% gof2$Hist)){
    g_other_sum <- g_other_sum + g_histories$Prob[i]
  }}
for (i in 1:nrow(b_histories)){
  if (!(rownames(b_histories)[i] %in% gof2$Hist)){
    b_other_sum <- b_other_sum + b_histories$Prob[i]
  }}

gof2[nrow(gof2) + 1,] <- c(NA,NA,NA,NA,NA,NA,'B',0,'Everything else B',b_other_sum)
gof2[nrow(gof2) + 1,] <- c(NA,NA,NA,NA,NA,NA,'G',0,'Everything else G',g_other_sum)


gof2$Frequency <- as.numeric(gof2$Frequency)
b_sum <- sum(gof2$Frequency[gof2$location == 'B'])
g_sum <- sum(gof2$Frequency[gof2$location == 'G'])

###check
gof2$expected_prob <- as.numeric(gof2$expected_prob)
sum(gof2$expected_prob[gof2$location == 'B'])
sum(gof2$expected_prob[gof2$location == 'G'])

gof2$Expected <- NA
for(i in 1:nrow(gof2)){
  if(gof2$location[i] == 'B'){
    gof2$Expected[i] <- as.numeric(gof2$expected_prob)[i]*b_sum
  } else if (gof2$location[i] == "G"){
    gof2$Expected[i] <- as.numeric(gof2$expected_prob)[i]*g_sum
  }
}

gof <- data.frame(gof2$Frequency,gof2$Expected)
fisher.test(gof, simulate.p.value = TRUE)

###############################################################################
### further analysis of best fitting model

print(mmmod_3$real) ##estimates of real 

coef(mmmod_3, param = 'psi', prob = 0.95)
coef(mmmod_3, param = 'theta', prob = 0.95)
coef(mmmod_3, param = 'p', prob = 0.95)
fitted(object = mmmod_3, param = 'theta')


results1 <- mmmod_3$real$p %>% 
  distinct(.keep_all = TRUE) %>%
  mutate(method = as.factor(c("eDNA","Camera")))
results1

################################################################################
### Figure 4

plot1 <- ggplot(results1, aes(x = method, y = est)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(
    aes(ymin = ciLow, ymax = ciHi), 
    width = 0.2,
    position = "dodge") +  theme_classic() + theme(
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
    axis.line = element_blank(),axis.text = element_text(colour = "black", size = 12),
    axis.title.x = element_text(margin = margin(t=20), size = 14),
    axis.title.y = element_text(margin = margin(r = 20), size =14),
    ) + labs(x = "Method",y = "Predicted probability of detection") +
  geom_text(aes(label = "*", y = ciHi), position = position_dodge(width = .9), 
    vjust = -.1, size = 20 / .pt) + lims(y = c(0, 1))

results2 <- mmmod_3$real$theta %>% 
  distinct(.keep_all = TRUE) %>%
  mutate(region = as.factor(c("Burrow","Drift fence")))
results2

plot2 <- ggplot(results2, aes(x = region, y = est)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(
    aes(ymin = ciLow, ymax = ciHi), 
    width = 0.2,
    position = "dodge") + theme_classic(
    ) + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
      axis.line = element_blank(),axis.text = element_text(colour = "black", size = 12),
      axis.title.x = element_text(margin = margin(t=20), size = 14),
      axis.title.y = element_text(margin = margin(r = 20), size =14),
    ) + labs(x = "Sampling location", y = paste0("Predicted availability (","\U03B8",")")) +
  geom_text(aes(label = "*", y=ciHi), position = position_dodge(
  width = .9), vjust = -.1, size = 20 / .pt) +
  lims(y = c(0, 1))

(plot_spacer() + plot1 + plot_spacer() + plot2 + plot_spacer()) + 
  plot_layout(widths = c(0.05, 1, 0.05, 1, 0.05)) 

################################################################################
################################################################################
### cost benefit analysis


################################################################################
###visual search

mean <- mean(126,178,147,178,129,216) ## mean of sample studies
mean1 <- 90 # average across
mean2 <- 40 ## lower end but common
mean3 <- 150 ## higher end

x <- c((2/126),(1/178),(1/147),(3/178),(3/129),(3/216))
p_search <- mean(x)
# search <- function(n){
#   return(1-(1-snakes)^(n))}
##per site kind of, with thetaB, needs help probably
search_h <- function(n){
  return(1-(1-p_search)^(n*mean3))}
search_l <- function(n){
  return(1-(1-p_search)^(n*mean2))}
#95 = (1-(1-0.01583333))^x

cost_search_low <- function(n){
  return(n*180.50)
}
cost_search_high <- function(n){
  return(n*447)
}

search_high <- as.data.frame(matrix(nrow = 10,ncol = 2))
colnames(search_high) <- c('cost','det_prob')
for (i in 1:10){
  search_high$cost[i] <- cost_search_high(i)
  search_high$det_prob[i] <- search_l(i)
}
search_low <- as.data.frame(matrix(nrow = 10,ncol = 2))
colnames(search_low) <- c('cost','det_prob')
for (i in 1:10){
  search_low$cost[i] <- cost_search_low(i)
  search_low$det_prob[i] <- search_h(i)
}

################################################################################
### eDNA

#burrow eDNA
burrow_edna <- function(n) {
  return(1-(1-thetaB*p_edna)^n)}
#ground eDNA
ground_edna <- function(n) {
  return(1-(1-thetaG*p_edna)^n)}

cost_edna_low <- function(n){
  return(60+33.50+(n * 95))}
cost_edna_high <- function(n){
  return(90+134+(n * 95))}

##using averages
cost_edna <- function(n){
  return(158.75+(n * 95))}

###eDNA dataframes

burrow_edna_df <- as.data.frame(matrix(nrow = 40,ncol = 2))
colnames(burrow_edna_df) <- c('cost','det_prob')
for (i in 1:40){
  burrow_edna_df$cost[i] <- cost_edna(i)
  burrow_edna_df$det_prob[i] <- burrow_edna(i)
}
ground_edna_df <- as.data.frame(matrix(nrow = 40,ncol = 2))
colnames(ground_edna_df) <- c('cost','det_prob')
for (i in 1:40){
  ground_edna_df$cost[i] <- cost_edna(i)
  ground_edna_df$det_prob[i] <- ground_edna(i)
}

################################################################################
### camera

#1-(1-p_day)^14 = p_camera(2 weeks)
#pday = 1.941278 or 0.05872178
p_camera_day = 0.05872178
#approximatedly 120 days in the survey season
days = 120
#ground camera
ground_camera <- function(n){
  return( 1-(1-thetaG*p_camera_day)^(n*days))}
#burrow camera
burrow_camera <- function(n){
  return(1-(1-thetaB*p_camera_day)^(n*days))}

bcost1 <- 1995.15
fcost1 <- 958.03
cost_bcamera_avg <- function(n){
  return(n * bcost1)}

cost_fcamera_avg <- function(n){
  return(n * fcost1)}


burrow_camera_df <- as.data.frame(matrix(nrow = 4,ncol = 2))
colnames(burrow_camera_df) <- c('cost','det_prob')
for (i in 1:4){
  burrow_camera_df$cost[i] <- cost_bcamera_avg(i)
  burrow_camera_df$det_prob[i] <- burrow_camera(i)
}

ground_camera_df <- as.data.frame(matrix(nrow = 4,ncol = 2))
colnames(ground_camera_df) <- c('cost','det_prob')
for (i in 1:4){
  ground_camera_df$cost[i] <- cost_fcamera_avg(i)
  ground_camera_df$det_prob[i] <- ground_camera(i)
}

###############################################################################
### combo eDNA and VES

p_edna1 <- (1-(1-thetaB*p_edna)^20)#  probability of an eDNA "effort"
p_search1 <- 1-(1-p_search)^(90)# 90 = average burrows at a site, probability of an active search "effort
p_multi <- p_search1 + (1 - p_search1) * p_edna1 # combined; 0.9

edna_cost <- 60 + (95 * 10)# cost for eDNA sampling effort
search_cost <- (240 + 120)/2 + (27+73)/2 #  cost for a search
travel_cost <- (33.50 + 134)/2 # cost for travel

number_visits <- 1:10 # number of visits to consider

eDNA_dollars <- (edna_cost + travel_cost) * number_visits # eDNA total cost
eDNA_pr <- 1 - (1 - p_edna1)^number_visits # eDNA total detection pr

search_dollars <- (search_cost + travel_cost) * number_visits # eDNA total cost
search_pr <- 1 - (1 - p_search1)^(number_visits)# search total detection pr

multi_dollars_low <- search_dollars
multi_dollars_hi <- (search_cost + edna_cost + travel_cost) * number_visits
multi_dollars_mean <- (search_cost + edna_cost * (1 - p_search1) + travel_cost) * number_visits
multi_pr <- 1 - (1 - p_multi)^number_visits

################################################################################
### Figure 7

cost_barplot <- data.frame (matrix(,nrow = 6,ncol = 4))
colnames(cost_barplot) <- c('type', 'cost_95', 'location','amount')

#search
#0.95 =  (1-(1-p_search)^(n*90))
n = (log(0.05))/(90*log(1-p_search))
n = ceiling(n)
cost_barplot$type[1] <- 'VES'
cost_barplot$cost_95[1] <- n *((180.50+447)/2)
cost_barplot$amount[1] <- paste(as.integer(n),'surveys')
cost_barplot$location[1] <- 'Full site'

##burrow edna
#0.95 = (1-(1-0.1134109)^n)
n = (log(0.05))/(log(1-(thetaB*p_edna)))
n = ceiling(n)
cost_barplot$type[2] <- 'Burrow eDNA'
cost_barplot$cost_95[2] <- cost_edna(n)
cost_barplot$amount[2] <- paste(as.integer(n),'samples')
cost_barplot$location[2] <- 'Burrow'

#ground eDNA
#0.95 = (1-(1-thetaG*p_edna)^n)
n = (log(0.05))/(log(1-(thetaG*p_edna)))
n = ceiling(n)
cost_barplot$type[3] <- 'Fence eDNA'
cost_barplot$cost_95[3] <- cost_edna(n)
cost_barplot$amount[3] <- paste(as.integer(n),'samples')
cost_barplot$location[3] <- 'Fence'

#burrow camera
#0.95 = (1-(1-0.thetaB*p_camera_day)^(n*120))
n = (log(0.05))/(120*log(1-(thetaB*p_camera_day)))
n = ceiling(n)
cost_barplot$type[4] <- 'Burrow camera'
cost_barplot$cost_95[4] <- cost_bcamera_avg(n)
cost_barplot$amount[4] <- paste(as.integer(n),'cameras')
cost_barplot$location[4] <- 'Burrow'

#ground camera
#0.95 = 1-(1-thetaG*p_camera_day)^(n*120)
n = (log(0.05))/(120*log(1-(thetaG*p_camera_day)))
n = ceiling(n)
cost_barplot$type[5] <- 'Drift fence'
cost_barplot$cost_95[5] <- cost_fcamera_avg(n)
cost_barplot$amount[5] <- paste(as.integer(n),'cameras')
cost_barplot$location[5] <- 'Fence'


#multi
#0.95 <- 1 - (1 - p_multi)^n
n = (log(0.05))/(log(1-(p_multi)))
n = ceiling(n)
cost_barplot$type[6] <- 'VES and eDNA'
cost_barplot$cost_95[6] <-(search_cost + edna_cost * (1 - p_search1) + travel_cost) * n
cost_barplot$amount[6] <- paste(as.integer(n),'surveys')
cost_barplot$location[6] <- 'Full site'

cost_barplot$type <- factor(cost_barplot$type, levels = unique(cost_barplot$type), 
                            ordered = TRUE)

ggplot(data = cost_barplot, aes(x = type, y = cost_95, fill = location)) + ylim(0,4300) +
  geom_bar(stat = 'identity') + labs(x = "",y = "Cost($)", fill = '') +
  theme_classic() + theme(panel.border = element_rect(
    colour = "black", fill=NA, linewidth=0.5), axis.line = element_blank(
    ),axis.text = element_text(colour = "black"), plot.title = element_text(
      hjust = 0.05, vjust = -9, face = "bold", size = 15), axis.title.y = 
      element_text(margin = margin(r=20)),legend.text = element_text(size=1)
  ) + geom_text(aes(label = amount), vjust = 1.5, size = 4, color = 'black') +
  theme(axis.text.x = element_text (size = 14),
        axis.text.y = element_text(size=14), axis.title = element_text (size = 14, 
  face = 'bold'))+ scale_fill_manual(values = c( '#1a80bb','#ea801c','#b8b8b8')) +
  ggtitle('Cost for 95% detection rate of eastern indigo snake') 

################################################################################       
### Figure 5

###eDNA BURROW
burrow_edna_thetaf <- function(n) {
  return(1-(1-n*p_edna)^25)}

burrow_edna_theta <- as.data.frame(matrix(nrow = 26, ncol = 2))
colnames(burrow_edna_theta) <- c('theta','Pr')
B <- 0.27
for(i in 1:nrow(burrow_edna_theta)){
  burrow_edna_theta$theta[i] <- B
  burrow_edna_theta$Pr[i] <- burrow_edna_thetaf(B)
  B <- B - 0.01
}
####eDNA Ground
ground_edna_thetaf <- function(n) {
  return(1-(1-n*p_edna)^30)}

ground_edna_theta <- as.data.frame(matrix(nrow = 22, ncol = 2))
colnames(ground_edna_theta) <- c('theta','Pr')
G <- thetaG
for(i in 1:nrow(ground_edna_theta)){
  ground_edna_theta$theta[i] <- G
  ground_edna_theta$Pr[i] <- ground_edna_thetaf(G)
  G <- G - 0.01
}
####CAMERA burrow
burrow_camera_thetaf <- function(n){
  return(1-(1-n*p_camera_day)^(2*120))
}
burrow_camera_theta <- as.data.frame(matrix(nrow = 26, ncol = 2))
colnames(burrow_camera_theta) <- c('theta','Pr')
B <- 0.27
for(i in 1:nrow(burrow_camera_theta)){
  burrow_camera_theta$theta[i] <- B
  burrow_camera_theta$Pr[i] <- burrow_camera_thetaf(B)
  B <- B - 0.01
}
####CAMERA ground
ground_camera_thetaf <- function(n){
  return(1-(1-n*p_camera_day)^(3 *120))
}
ground_camera_theta <- as.data.frame(matrix(nrow = 22, ncol = 2))
colnames(ground_camera_theta) <- c('theta','Pr')
G <- thetaG
for(i in 1:nrow(ground_camera_theta)){
  ground_camera_theta$theta[i] <- G
  ground_camera_theta$Pr[i] <- ground_camera_thetaf(G)
  G <- G - 0.01
}


plot(NULL, NULL, xlim = c(0,0.27), ylim = c(0,1), las = 1, xlab = 
       paste("\U03B8", "(burrow)"), ylab = "Pr (detection)")
lines(burrow_edna_theta$theta,burrow_edna_theta$Pr, col = "chocolate", lwd = 2, lty = 6)
lines(burrow_camera_theta$theta, burrow_camera_theta$Pr, col = "cadetblue", lwd = 2, lty = 6)
lines(c(0.26,0.26), c(1.0,0.9), col = "black", lwd = 3)                                                                                                          


legend("topleft", cex = 0.8, bty = "n", lwd = c(2,2,3), lty = c(6,6,1),
       col = c("chocolate","cadetblue","black"),
       legend = c("Burrow eDNA", "Burrow camera", "Our theta value"))

################################################################################
### Figure 6

detection_rate <- read.csv("rate_of_detection.csv")

rownames(detection_rate) <- detection_rate$X
detection_rate <- detection_rate[,-which(colnames(detection_rate) == 'X')]
rownames(detection_rate)[4:5] <- c('p_eDNA * θ_eDNA','p_camera * θ_camera')

rod_df <- as.data.frame(matrix(nrow=25,ncol=3))
colnames(rod_df) <- c('Site','Method','Probability')

x <-1
for(i in 1:nrow(detection_rate)){
  for(j in 1:ncol(detection_rate)){
    rod_df$Method[x] <- rownames(detection_rate)[i]
    rod_df$Site[x] <- colnames(detection_rate)[j]
    rod_df$Probability[x] <-detection_rate[i,j]
    x=x+1
  }
}
rod_df<-na.omit(rod_df)


rod_df2 <- rod_df
rod_df2['Prob2'] <- NA
for(i in 1:nrow(rod_df2)){
  if(grepl('eDNA',rod_df2$Method[i])){
    rod_df2$Prob2[i] <- 1-(1- rod_df2$Probability[i])^(20*3)
    
  } else if (grepl('amera',rod_df2$Method[i])){
    rod_df2$Prob2[i] <- 1-(1- rod_df2$Probability[i])^(120)
    
  } else if (grepl('VES',rod_df2$Method[i])){
    rod_df2$Prob2[i] <- 1-(1- rod_df2$Probability[i])^((3*90))
}}

rod_df2$Method = factor(rod_df2$Method, levels = c("Camera","eDNA","VES",
"p_camera * θ_camera","p_eDNA * θ_eDNA"))

ggplot(data = rod_df2, aes(x = Site, y = Prob2, fill = Method)
    ) + geom_bar(position = 'dodge', stat = 'identity',width = 0.7) + theme_classic() +
   theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
  axis.title.y = element_text(size = 14, margin = margin(r = 20)), 
  axis.title.x = element_text(size = 14, margin = margin(t = 20)
  ), axis.text.x = element_text (size = 13),
  axis.text.y = element_text(size = 13), legend.text = element_text(size =11),
  legend.title = element_text(size = 12), plot.margin = margin(1,1.5,1,1, "cm")
  ) + scale_fill_manual(breaks = c("Camera","eDNA","VES","p_camera * θ_camera",
  "p_eDNA * θ_eDNA"),values = c("#648FFF", "#FFB000","ivory3","slateblue4","chocolate"
                                    )) + labs(y = "Rate of detection")
