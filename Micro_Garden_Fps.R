rm(list = ls())

##################################

# Load librarys
library(reshape2)
library(matrixTests)
library(ggpubr)
library(tidyverse)
library(vegan)

##################################

# Load o2 respiration data
prelim <- read.csv('/path/to/file/O2_FP.csv', header=T)

##################################

# Split data by experiment (# 4) and run linear regression between O2 and time


exp1 <- prelim[prelim$Exp==1,]
exp2 <- prelim[prelim$Exp==2,]
exp3 <- prelim[prelim$Exp==3,]
exp4 <- prelim[prelim$Exp==4,]

all <- list(exp1, exp2, exp3, exp4)


coef <- c()
type <- c()
pv_i <- c()
pv_c <- c()
exp.all <- c()
r2 <- c()
for (j in 1:length(all)){
  exp.no <- rep(j,6)
  exp.all <- c(exp.all, exp.no)
  for (i in 3:8){
    exp <- as.data.frame(all[j])
    x <- exp$Time
    y <- exp[,i]
    r <- lm(y~x)
    print(names(exp[i]))
    type <- c(type, names(exp[i]))
    #print(summary(r))
    coef <- c(coef, coef(r)[2])
    #f <- summary(r)$fstatistic 
    #print(as.numeric(summary(r)$coefficients[,4][1]))
    #print(as.numeric(summary(r)$coefficients[,4][2]))
    pv_i <- as.numeric(c(pv_i, summary(r)$coefficients[,4][1]))
    pv_c <- as.numeric(c(pv_c, summary(r)$coefficients[,4][2]))
    r2 <- c(r2, summary(r)$r.squared)
  }
}
length <- rep(c('control', 'control', 'long', 'long', 'short', 'short'), 4)
col <- rep(c('control', 'control', 'brown', 'green', 'brown', 'green'), 4)

exp.coef <- data.frame(exp.all, type, coef, length, col, pv_i, pv_c, r2)

# make O2 uptake umol h-1 positive
exp.coef$coef <- exp.coef$coef * -1

##########################

# Remove control O2 uptake slopes from experimental vials slopes
exp.coef$coef.blnkd <- NA

for (i in 3:6){
  exp.coef$coef.blnkd[i] <- exp.coef$coef[i] - mean(exp.coef$coef[1:2])
}

for (i in 9:12){
  exp.coef$coef.blnkd[i] <- exp.coef$coef[i] - mean(exp.coef$coef[7:8])
}

for (i in 15:18){
  exp.coef$coef.blnkd[i] <- exp.coef$coef[i] - mean(exp.coef$coef[13:14])
}

for (i in 21:24){
  exp.coef$coef.blnkd[i] <- exp.coef$coef[i] - mean(exp.coef$coef[19:20])
}

# remove control data dfrom data set
exp.coef <- exp.coef[exp.coef$length!='control',]

# set length to be a factor
exp.coef$length <- factor(exp.coef$length)
exp.coef$type <- droplevels(exp.coef$type)

# Correct labelling of long green pellets
library(plyr)
exp.coef$type = revalue(exp.coef$type, c("Gl"="GL"))


################################################

# Add POC data to calculate microbial turnover of POC

poc <- read.csv('/path/to/file/POC_FP.csv', header=T)

# merge POC to O2 uptake results
exp.coef$poc <- poc$Cumol_L

# compute uptake (/h)
exp.coef$k_h <- exp.coef$coef.blnkd/exp.coef$poc

# compute uptake (/d)
exp.coef$k_d <- exp.coef$k_h * 24


################################################

# Analyse FP size data

# Find all .csv in folder with files BL1.csv.....GS4.csv
files <- dir("/path/to/files/folder-with-perimeter-files", recursive=TRUE, full.names=TRUE, pattern="\\.csv$")

# Add columns of pellet type and experiment number
poop = c(rep(c('BL', 'BS', 'GL', 'GS'), 4))
col = rep(c('brown', 'brown', 'green', 'green'),4)
exp = c(rep(1:4, each=4))

fp_length = data.frame()
for (i in 1:length(files)){
  df <- read.csv(files[i])
  df$exp = rep(exp[i], nrow(df))
  df$poop = rep(poop[i], nrow(df))
  df$col = rep(col[i], nrow(df))
  fp_length = rbind(fp_length, df)
}

# remove very low perimeter values
fp_length = fp_length[fp_length$Perim. > 0.05,]

# Calculate mean and sd per experiment/poo type
detach(package:plyr)
Perim_data = fp_length %>%
  group_by( poop) %>%
  summarise(mean=mean(Perim.), sd=sd(Perim.))


# Match with o2/poc data frame
exp.coef$perim<-Perim_data$mean

#######################################################

# Krill fragmentation experiment (for figure 3)

# read in file
frag = read.csv("/path/to/file/In_Situ_fragmentation.csv", header=TRUE)

# Find mean FP perimeter length at each sampling point
frag_mean = frag %>%
  group_by(Label) %>% # Label = sampling event
  summarise(mean=mean(Perim.), se=sd(Perim.)/sqrt(n()), n = n())

# Add in time (as decimals)
frag_mean$time = c(10.25, 11.25, 12.50,  14.50,16.75) 


#######################################################

# Figures for paper

# 1 images of 4 FP types 

# 2

# 2a Perim (size) 
a = ggplot(fp_length, aes(x=poop, y=Perim., fill=col)) +
  geom_boxplot(outlier.shape = NA)  +
  xlab('') + ylab(expression(paste('Faecal pellet perimeter (mm)')))+
  scale_fill_manual(values = c('brown', 'chartreuse4')) +
  theme_classic() + theme(legend.position="none", axis.title.y = element_text(size = 16),axis.text=element_text(size=10)) +
  scale_x_discrete(labels = c("Brown Long", "Brown Short", "Green Long", "Green Short"))+ 
  ggtitle("a) Faecal pellet perimeter length") +
  stat_summary(fun.y=mean, aes(group=as.factor(poop)), geom="point", shape=20, size=5, color="black", fill="black")+
  coord_cartesian(ylim=c(0, 21))

# 2b O2 uptake
b = ggplot(exp.coef, aes(x=type, y=coef.blnkd, fill=col)) +
  geom_boxplot(outlier.shape = NA)  +
  xlab('') + ylab(expression(paste(O[2], ' uptake (',mu,'mol L'^'-1'*' h'^'-1'*')')))+
  scale_fill_manual(values = c('brown', 'chartreuse4')) +
  theme_classic() + theme(legend.position="none", axis.title.y = element_text(size = 16),axis.text=element_text(size=10)) +
  scale_x_discrete(labels = c("Brown Long", "Brown Short", "Green Long", "Green Short"))+ 
  ggtitle("b) Microbial oxygen uptake") +
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="black")


# 2c POC conc
c = ggplot(exp.coef, aes(x=type, y=poc/1000, fill=col)) +
  geom_boxplot()  +
  xlab('') +   ylab(expression(paste('POC concentration (mmol L'^'-1'*')')))+
  scale_fill_manual(values = c('brown', 'chartreuse4')) +
  theme_classic() + theme(legend.position="none", axis.title.y = element_text(size = 16),axis.text=element_text(size=10)) +
  scale_x_discrete(labels = c("Brown Long", "Brown Short", "Green Long", "Green Short"))+ 
  ggtitle("c) Faecal pellet POC content") +
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="black")

# 2d POC turnover
d = ggplot(exp.coef, aes(x=type, y=k_d, fill=col)) +
  geom_boxplot()  +
  xlab('') + ylab(expression(paste('Microbial POC turnover (d'^'-1'*')')))+
  scale_fill_manual(values = c('brown', 'chartreuse4')) +
  theme_classic() + theme(legend.position="none", axis.title.y = element_text(size = 16),axis.text=element_text(size=10)) +
  scale_x_discrete(labels = c("Brown Long", "Brown Short", "Green Long", "Green Short"))+ 
  ggtitle("d) Microbial POC turnover") +
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="black")




quartz(width=9,height=9) # code for mac
par(mar=c(5,5,4,1), mfrow=c(2,2), bty='o', cex=1)

ggarrange(a, b, c, d, ncol=2, nrow=2)


# Fig 3 - schematic with fragmentation data using base plot

quartz(width=6,height=4) # code for mac
par(mar=c(4,4,1,5), mfrow=c(1,1), bty='o', cex=1)

plot(frag_mean$time, frag_mean$mean, xlim=c(9,17), ylim=c(0,6), pch=16, cex=frag_mean$n/max(frag_mean$n)*3,
     ylab='Mean faecal pellet perimeter (mm)', xlab='Time of day',xpd=NA,
     col=c('chartreuse2', 'chartreuse3', 'chartreuse4', 'chocolate3' , 'chocolate4'), type='n')

xy.error.bars<- function (x,y,ybar){
  plot(frag_mean$time, frag_mean$mean, xlim=c(9,17), ylim=c(0,7), pch=16, cex=frag_mean$n/max(frag_mean$n)*3,
       ylab='Mean faecal pellet perimeter (mm)', xlab='Time of day',xpd=NA,
       col=c('chartreuse2', 'chartreuse3', 'chartreuse4', 'chocolate3' , 'chocolate4'), type='n')
  arrows(x, y-ybar, x, y+ybar, code=3, angle=90, length=0.05, lwd=1.8, xpd=FALSE,
         col=c('chartreuse2', 'chartreuse3', 'chartreuse4', 'chocolate3' , 'chocolate4')) 
}
x<-frag_mean$time
y<-frag_mean$mean
yb<-frag_mean$se
xy.error.bars(x,y,yb)

time = c(9,12.3,14)
food = c(6, 4, 0)
lines(time, food, col='chartreuse2', lwd=2)
time = c(9, 17)
poop = c(0, 6)
lines(time, poop, col='chocolate4', lwd=2)
points(frag_mean$time, frag_mean$mean,  pch=16, cex=frag_mean$n/max(frag_mean$n)*3,
       col=c('chartreuse2', 'chartreuse3', 'chartreuse4', 'chocolate3' , 'chocolate4'))
mtext(c('Food', 'Faecal pellets'), side=4, line=1.5, at=c(1,4) , col=c("chartreuse3", 'chocolate4'))
arrows(x0=12.3, y0=7, x1=12.3, y1=5.5, col='grey', lwd=2)
arrows(x0=17.8, y0=7, x1=17.8, y1=0, col='chartreuse3', lwd=2, xpd=NA)
arrows(x0=18.6, y0=0, x1=18.6, y1=7, col='chocolate4', lwd=2, xpd=NA)


#############################################



# Supplementary material


######################################################

# Statistical tests between k_d of different faecal fractions
# t.test across different groups

poop_tmp = c('BL', 'BS', 'GL', 'GS')
comb = as.matrix(combn(poop_tmp,2))

matrix_colnames       <- c("Group.1", "Group.2", "POC_turn.1","POC_turn.2","p.POC_turn")
pair           <- matrix(ncol=length(matrix_colnames),nrow=6) # need to specify nrow. e.g. length(Dates)
colnames(pair) <- matrix_colnames
counter       <- 1 


for (i in 1:6){
  t = (col_t_welch(exp.coef[exp.coef$type==comb[1,i],12], exp.coef[exp.coef$type==comb[2,i],12]))

  pair[counter,1] <- comb[1,i]
  pair[counter,2] <- comb[2,i]
  pair[counter,3] <- t[1,4]            
  pair[counter,4] <- t[1,5]   
  pair[counter,5] <- t[1,12]   
  #pair[counter,6] <- perm$tab[1,4]   # F-value
  counter <- counter+1
}

pair = as.data.frame(pair)


# Is non-significant t-test between BL and both GL and GS because of difference in variance?
p = exp.coef %>% group_by(type) %>% summarise(var=var(k_d))


######################################################

