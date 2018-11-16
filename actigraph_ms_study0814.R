## MS, Actigraph study 
setwd('C:/Users/fanti/Documents/Dolores Document/University/Biostatistics Project/Summer Research/Data cleaned/MS_CrossSec_study')

########################### 
#       Libararies
###########################
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(reshape2)
library(plyr)
library(zoo)
library(tibble)

# For plotting
library(RColorBrewer)
library(gridExtra)
library(grid)
library(ggpubr)
library(lemon)
library(ggfortify)

# For FOSR
library(mgcv)
library(refund)

############################
#      Functions
############################
# 1.
### Function to create a dataframe contain 1400 log(activity counts) and date info. for each day
# 1) Starts from 00:00:00 on each day 
# 2) If start time is different from 00:00am, assume 0 counts before the time started.
# 3) The returned dataframe has first 1400 columns of log-activity counts, and the 1441th column is date
data_completedays <- function(dataframe){
  
  dataframe$timedate <- as.POSIXct(dataframe$timedate)
  unique.date <- unique(date(dataframe$timedate))
  
  ac.matrix <- data.frame(matrix(0, nrow=length(unique.date), ncol = 1441)) # Initiate an 0 dataframe
  
  ac.matrix[,1441] <- format(date(unique.date), '%Y-%m-%d') # Last column is date
  ac.matrix[,1441] <- as.POSIXct(ac.matrix[,1441])
  names(ac.matrix)[1441] <- 'date'
  
  for(i in 1:dim(ac.matrix)[1]){
    date.ind <- date(dataframe$timedate)==ac.matrix$date[i] # Date indicator
    data.len <- length(dataframe[date.ind,]$timedate) # Data length on that day
    if(data.len != 1440){ 
      sub.start <- 1440-data.len+1 # If not from 00:00am, starts on the exact start time
      ac.matrix[i,sub.start:1440] <- log(dataframe[date.ind,]$vectormagnitude+1)
    }
    else{
      ac.matrix[i,1:1440] <- log(dataframe[date.ind,]$vectormagnitude+1)
    }
  }
  return(ac.matrix)
}

# 2. 
### Funtion to generate valid days
# Non-valid days are defined as:  
# 1) Days with more than a consecutive of 2.5 hours non-wearing time 
# 2) Days with less than 480 active minutes (8 hrs) counts

valid_full_days <- function(data, missing_hrs=2.5, active_threshold=8){
  # input data has last column of dates
  # work on first 1:1440 columns of log-activity-counts only
  
  ## Define thresholds for 1) missing hours and 2) non-wear hours
  consec_missing_minutes = missing_hrs*60
  active_minute_threshold = active_threshold*60
  
  # Remove days that have active minutes below threshold
  below_active <- (which(rowSums(data[,1:1440] > 0) < active_minute_threshold))
  data <- data %>% filter(!(row_number() %in% below_active))
  
  # Remove days that have consecutive non-wearing time
  data <- data[!apply(data[,1:1440], 1, FUN = function(x) any(with(rle(x ==0), lengths[values]) > consec_missing_minutes)),]
  return(data)
}

# 3.
### Functions to generate 10 summaries (Vadim's code)
stat.subject.generate <- function(data.subject){
  
  sub.date <- data.subject$date
  ###log-transform data ################################################
  ## data already log-transformed
  data.subject = t(data.subject[,1:1440])
  
  ###calculation of summary statistics#################################
  mean.stat <- apply(data.subject,2,mean)
  sd.stat <- apply(data.subject,2,sd)
  coefvar.stat <- apply(data.subject,2,function(a){sd(a)/mean(a)})
  num.actmins <- apply(data.subject,2,function(a){sum(a>0)})
  act.mean <- 1440*mean.stat/num.actmins
  
  max6min <- function(day.counts){
    max6min.window <- rbind(day.counts[1:1435],day.counts[2:1436],day.counts[3:1437],
                            day.counts[4:1438],day.counts[5:1439],day.counts[6:1440])
    return(max(apply(max6min.window,2,sum)))
  }
  max6.sum <- apply(data.subject,2,max6min)
  total.variablity <- function(day.counts){
    head.count <- head(day.counts,-1)
    tail.count <- tail(day.counts,-1)
    zero.switch <- sum(head.count==0&tail.count==0)
    return(sum(abs(diff(day.counts)))/(1440-zero.switch))
    #return(sum(sign(diff(day.counts))==1)/length(sign(diff(day.counts))))
  }
  total.var <- apply(data.subject,2,total.variablity)
  fragmentation <- function(day.counts){
    indicator.vec <- ifelse(day.counts>0,1,0)
    return(sum(abs(diff(indicator.vec))))
  }
  frag.stat <- apply(data.subject,2,fragmentation)
  ###Relative Amplitude###
  roll <- function(day.counts,k){
    kvec <- rollapplyr(day.counts, k, function(x) mean(x), fill = NA)
    kvec <- kvec[!is.na(kvec)]
  }
  M10.comp <- function(day.counts){
    return(max(roll(day.counts,600)))
  }
  L5.comp <- function(day.counts){
    return(min(roll(day.counts,300)))
  }
  M10 <- apply(data.subject, 2, M10.comp)
  L5 <- apply(data.subject, 2, L5.comp)
  RA <- (M10 - L5)/(M10 + L5)
  ###Intradaily Variability###
  intradaily.variability <- function(day.counts){
    mean.counts <- mean(day.counts)
    numerator <- sum(diff(day.counts)^2)
    denominator <- sum((day.counts-mean.counts)^2)
    return(numerator/denominator)
  }
  IV <- apply(data.subject, 2, intradaily.variability)
  subj.characteristic.data <- data.frame(days=sub.date, mean.stat=mean.stat,sd.stat=sd.stat,
                                         coefvar.stat=coefvar.stat,num.actmins=num.actmins,act.mean=act.mean,max6.sum=max6.sum,total.var=total.var,frag.stat=frag.stat, RA=RA, IV=IV)
  return(subj.characteristic.data = subj.characteristic.data)
}


### Other Functions used in plotting #####

# 1. For plotting activity profile based on day of the week
### Function to return a list of lists that splits weekdays to different weeks
# Input data fomat is the dataframe of size 1441xn (1440 log-ac + 1 date )
# n=number of days
# return is a list of lists, such that each sublist is ordered by week1, week2 etc.
# each sublist contains dataframe of log-ac + date 
Sub.weeks <- function(subject.df){
  
  n <- dim(subject.df)[1]
  day.of.week <- weekdays(subject.df$date)
  Sun.sub <- grep("Sunday",day.of.week)
  if(Sun.sub[1]==1){
    if(length(Sun.sub)==1){
      Week1 <- subject.df[Sun.sub[1]:n,]
      Weeks <- Week1
    }else{
      if(length(Sun.sub)==2){
        Week1 <- subject.df[Sun.sub[1]:Sun.sub[2],]
        Week2 <- subject.df[(Sun.sub[2]+1):n,]
        Weeks <- list(Week1, Week2)
      }else{
        if(length(Sun.sub)==3){
          Week1 <- subject.df[Sun.sub[1]:Sun.sub[2],]
          Week2 <- subject.df[(Sun.sub[2]+1):Sun.sub[3],]
          Week3 <- subject.df[(Sun.sub[3]+1):n,]
          Weeks <- list(Week1, Week2, Week3)
        }}}
  }else{
    if(length(Sun.sub)==1){
      Week1 <- subject.df[1:Sun.sub[1],]
      Week2 <- subject.df[(Sun.sub[1]+1):n,]
      Weeks <- list(Week1, Week2)
    }else{
      if(length(Sun.sub)==2){
        Week1 <- subject.df[1:Sun.sub[1],]
        Week2 <- subject.df[(Sun.sub[1]+1):Sun.sub[2],]
        Week3 <- subject.df[(Sun.sub[2]+1):n,]
        Weeks <- list(Week1, Week2, Week3)
      }else{
        if(length(Sun.sub==3)){
          Week1 <- subject.df[1:Sun.sub[1],]
          Week2 <- subject.df[(Sun.sub[1]+1):Sun.sub[2],]
          Week3 <- subject.df[(Sun.sub[2]+1):Sun.sub[3],]
          Week4 <- subject.df[(Sun.sub[3]+1):n,]
          Weeks <- list(Week1, Week2, Week3, Week4)
        }}}
  }
  return(Weeks)
}

# 2. For plotting activity profile based on day of the week
### Function to generate layout matrix for plotting
# The input varible plot.list is the list of lists based on weeks: week1, week2, etc.
plt.layout <- function(plot.list){
  n <- length(plot.list)
  lay.mat <- matrix(NA,nrow=7,ncol=n)
  
  plt.week1 <- plot.list[[1]]
  lay.mat[(7-length(plt.week1)+1):7,1] <- rep(1:length(plt.week1))
  
  plt.week2 <- plot.list[[2]]
  if(n>2){
    plt.week3 <- plot.list[[3]]
    lay.mat[,2] <- rep((length(plt.week1)+1):(length(plt.week1)+7))
    if(n==3){
      lay.mat[1:length(plt.week3),3] <- rep((length(plt.week1)+8):length(unlist(plt.weeks,recursive=FALSE)))
    }else{
      if(n==4){
      plt.week4 <- plot.list[[4]]
      lay.mat[,3] <- rep((length(plt.week1)+8):(length(plt.week1)+14))
      lay.mat[1:length(plt.week4),4] <- rep((length(plt.week1)+15):length(unlist(plt.weeks,recursive=FALSE)))
      }
    }
  }else{
    lay.mat[1:length(plt.week2),2] <- rep((length(plt.week1)+1):length(unlist(plt.weeks,recursive=FALSE)))
  }
  return(lay.mat)
}

# 3. For plotting ten summaries
### Function to average ten-summaries statistics by weekdays
# Input is the ten-summaries for all days
ave.by.day <- function(sub.ten.sum){
  df.week.ave <- aggregate(. ~ weekdays(sub.ten.sum$days), sub.ten.sum, mean)
  df.week.ave <- subset(df.week.ave, select = -c(days))
  names(df.week.ave)[1] <- 'day'
  df.week.ave$day <- ordered(df.week.ave$day, levels=c("Monday", "Tuesday", "Wednesday", "Thursday", 
                                                       "Friday", "Saturday", "Sunday"))
  df.week.ave <- df.week.ave[order(df.week.ave$day),]
  return(df.week.ave)
}

# 4. For plotting linear model
### Function to plot Regression Models
# Input is regression result
ggplotRegression <- function (fit) {
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

###############################
#       Read Data Files
###############################


###### Read Kathy study dataset #######
K.datafolder <- 'kathleen_data/actigraph_data' #Kathy data folder path
K.subinfo <- 'kathleen_data/K_MS Accelerometer Data.csv' #Kathy clinical file name

K.files <- list.files(path=K.datafolder) #read file names from folder
K.names <- sapply(K.files,word) #extract datafile name without datetime.csv extension => this gives subject ID

K.sub.data <- read.csv(K.subinfo,header = TRUE, stringsAsFactors = FALSE) #read clinical data

# Read actigraph data from Kathy's study
for(i in 1:length(K.names)){
  assign(K.names[[i]],read.csv(file.path(K.datafolder,K.files[i]),header = TRUE, stringsAsFactors=FALSE))
}


###### Read Mowry Baseline study dataset  ######
M.datafolder <- 'mowry_data/actigraph_data baseline' #Mowry data folder path
M.subinfo <- 'mowry_data/M_MS Accelerometer Data Baseline.csv' #Mowry clinical file name

M.files <- list.files(path=M.datafolder) #read file names from folder
M.names <- sapply(M.files,word) #extract datafile name without extension

M.sub.data <- read.csv(M.subinfo,header = TRUE, stringsAsFactors = FALSE) #read clinical data

# Read actigraph data from Mowry's study
for(i in 1:length(M.names)){
  assign(M.names[[i]],read.csv(file.path(M.datafolder,M.files[i]),header = TRUE, stringsAsFactors=FALSE))
}

# Concatenate raw dataframes into a list
datalist <- mget(c(ls(pattern="^ACC"),ls(pattern="^TDPM"))) # list of raw data

######################################
#        Data Processing (Stage 1)
######################################
# Stage 1:
# We will process raw data for each subject by days start from 00:00am to 23:59pm
# We will not consider invalid days and outliers at this stage
# We will also process clinical info at this stage, to include concerned predictor variables
# Full data will be used for plotting activity profiles for each subject
# It will also be used to plot TEN-SUMMARIES for each subject

#### 1. Process Actigraph Data
# Convert raw data to log-activity counts for each day
# Re-assign dataframes
for(i in 1:length(datalist)){
  assign(names(datalist)[i],data_completedays(datalist[[i]]))
}

# Concatenate new dataframes into another list
data.full<- mget(c(ls(pattern="^ACC"),ls(pattern="^TDPM"))) ## total 64 subjects

# data.full is a list of dataframes named by Subject ID, and contains log-ac on all dates recorded


#### 2. Process Subject clinical information data
# Create a dataframe to include subject name, EDSS score, T25fw score, and TUG(t2) score
clinical_info <- data.frame(matrix(nrow=length(datalist),ncol=7)) # Initiate dataframe
names(clinical_info) <- c('Subject.ID','EDSS','T25fw','TUG','SEX','AGE','BMI') 

two_names <- unlist(c(K.names,M.names))
clinical_info$Subject.ID <- two_names

# EDSS score for each subject
for(i in 1:length(datalist)){
  for(j in 1:dim(K.sub.data)[1]){
    if(K.sub.data$Subject.ID[j] == clinical_info$Subject.ID[i]){
      clinical_info$EDSS[i] = K.sub.data$edss_score[j]
      clinical_info$T25fw[i] = K.sub.data$t25.Average[j]
      clinical_info$TUG[i] = K.sub.data$TUG.t2[j]
      clinical_info$SEX[i] = K.sub.data$SEX[j]
      clinical_info$AGE[i] = K.sub.data$AGE[j]
      clinical_info$BMI[i] = K.sub.data$BMI[j]
      
    }
  }
  for(j in 1:dim(M.sub.data)[1]){
    if(M.sub.data$Subject.ID[j] == clinical_info$Subject.ID[i]){
      clinical_info$EDSS[i] = M.sub.data$edss_score[j]
      clinical_info$T25fw[i] = M.sub.data$t25.Average[j]
      clinical_info$TUG[i] = M.sub.data$TUG.t2[j]
      clinical_info$SEX[i] = M.sub.data$SEX[j]
      clinical_info$AGE[i] = M.sub.data$AGE[j]
      clinical_info$BMI[i] = M.sub.data$BMI[j]
    }
  }
}


# Classify EDSS, T25fw, TUG groups
# EDSS intervals: [1,2];[3,4];[6,7]
# T25fw intervals: (0,6.35);(6.35,9.53),(9.53,19.05),(>19.05)
# TUG intervals identified by T25fw intervals

edss.int1 <- (clinical_info$EDSS >= 1 ) & (clinical_info$EDSS <= 2)
edss.int2 <- (clinical_info$EDSS >= 2.5 ) & (clinical_info$EDSS <= 3.5)
edss.int3 <- (clinical_info$EDSS >= 4 ) & (clinical_info$EDSS <= 7)

clinical_info$EDSS.level <- rep(NA,dim(clinical_info)[1])
clinical_info$EDSS.level[edss.int1] <- "EDSS-(1, 2)"
clinical_info$EDSS.level[edss.int2] <- "EDSS-(2.5, 3.5)"
clinical_info$EDSS.level[edss.int3] <- "EDSS-(4, 7)"
clinical_info$EDSS.level <- factor(clinical_info$EDSS.level, 
                               levels = c("EDSS-(1, 2)", "EDSS-(2.5, 3.5)","EDSS-(4, 7)"))


t25.int1 <- (clinical_info$T25fw > 0 ) & (clinical_info$T25fw < 6.35)
t25.int2 <- (clinical_info$T25fw >= 6.35 ) & (clinical_info$T25fw < 9.53)
t25.int3 <- (clinical_info$T25fw >= 9.53 ) & (clinical_info$T25fw < 19.05)
t25.int4 <- clinical_info$T25fw >= 19.05 

clinical_info$T25.level <- rep(NA,dim(clinical_info)[1])
clinical_info$T25.level[t25.int1] <- "T25fw-( < 6.35)"
clinical_info$T25.level[t25.int2] <- "T25fw-(6.35, 9.53)"
clinical_info$T25.level[t25.int3] <- "T25fw-(9.53, 19.05)"
clinical_info$T25.level[t25.int4] <- "T25fw-( >19.05)"
clinical_info$T25.level <- factor(clinical_info$T25.level,
                              levels =c("T25fw-( < 6.35)","T25fw-(6.35, 9.53)",
                                        "T25fw-(9.53, 19.05)","T25fw-( >19.05)"))

######### End of Stage 1 Data Process ############

##########################################
##             Plots (Stage 1)          ##
##########################################

######### Figure 0: Histograms
# of EDSS, T25fw and TUG scores
# Plot histogram of EDSS, T25fw, TUG based on EDSS intervals
co.EDSS = c('royalblue1','paleturquoise3','red')
p1 <- ggplot(clinical_info,aes(EDSS,fill=EDSS.level))+geom_bar(alpha=0.5)+
  ggtitle('Histogram of EDSS')+
  scale_y_discrete("Number of subjects", limits=seq(0,12,1))+
  scale_x_discrete("EDSS Score", limits = seq(0,7.5,0.5)) +
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.title.x = element_text(colour="grey20",size=14),
        axis.title.y = element_text(colour="grey20",size=14))+
  scale_fill_manual(name="EDSS intervals",values=co.EDSS)

p2<-ggplot(clinical_info, aes(T25fw, fill=EDSS.level))+
  geom_histogram(breaks=seq(0, 80, by = 2),col='gray', alpha=0.5)+
  ggtitle('Histogram of T25fw')+
  scale_y_discrete("Number of subjects", limits=seq(0,22,2))+
  scale_x_discrete("T25fw Seconds", limits = seq(0,80, 5))+
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.title.x = element_text(colour="grey20",size=14),
        axis.title.y = element_text(colour="grey20",size=14))+
  scale_fill_manual(name="EDSS intervals",values=co.EDSS)

p3<-ggplot(clinical_info, aes(TUG, fill=EDSS.level))+
  geom_histogram(breaks=seq(0, 80, by = 2),col='gray', alpha=0.5)+
  ggtitle('Histogram of TUG')+
  scale_y_discrete("Number of subjects", limits=seq(0,20,2))+
  scale_x_discrete("TUG Seconds", limits = seq(0,80, 5))+
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.title.x = element_text(colour="grey20",size=14),
        axis.title.y = element_text(colour="grey20",size=14))+
  scale_fill_manual(name="EDSS intervals",values=co.EDSS)

png('plot/hist_1.png',  units='px', width = 1300, height =300)
ggarrange(p1,p2,p3,nrow=1,ncol=3,common.legend = TRUE, legend="right")
dev.off()

# three people in High EDSS group with T25fw below 4 seconds
clinical_info[clinical_info$EDSS >=4 & clinical_info$T25fw < 4,1]
#"TDPMSA029" "TDPMSA046" "TDPMSA048"

# Plot histogram of EDSS, T25fw, TUG based on T25fw intervals
co.t25 = c('steelblue','lightseagreen','plum','orangered')
p1<-ggplot(clinical_info, aes(EDSS, fill=T25.level))+geom_bar(alpha=0.5)+
  ggtitle('Histogram of EDSS')+
  scale_y_discrete("Number of subjects", limits=seq(0,12,1))+
  scale_x_discrete("EDSS Score", limits = seq(0,7.5,0.5)) +
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.title.x = element_text(colour="grey20",size=14),
        axis.title.y = element_text(colour="grey20",size=14))+
  scale_fill_manual(name="T25fw intervals",values=co.t25)

p2<-ggplot(clinical_info, aes(T25fw, fill=T25.level))+ 
  geom_histogram(breaks=seq(0, 80, by = 2),col='gray', alpha=0.6)+
  ggtitle('Histogram of T25fw')+
  scale_y_discrete("Number of subjects", limits=seq(0,22,2))+
  scale_x_discrete("T25fw Seconds", limits = seq(0,80, 5))+
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.title.x = element_text(colour="grey20",size=14),
        axis.title.y = element_text(colour="grey20",size=14))+
  scale_fill_manual(name="T25fw intervals",values=co.t25)


p3<-ggplot(clinical_info, aes(TUG, fill=T25.level))+ 
  geom_histogram(breaks=seq(0, 80, by = 2),col='gray', alpha=0.6)+
  ggtitle('Histogram of TUG')+
  scale_y_discrete("Number of subjects", limits=seq(0,20,2))+
  scale_x_discrete("TUG Seconds", limits = seq(0,80, 5))+
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.title.x = element_text(colour="grey20",size=14),
        axis.title.y = element_text(colour="grey20",size=14))+
  scale_fill_manual(name="T25fw intervals",values=co.t25)

png('plot/hist_2.png',  units='px', width = 1300, height =300)
ggarrange(p1,p2,p3,nrow=1,ncol=3,common.legend = TRUE, legend="right")
dev.off()


######### Figure 1 - Activity Profile plot: plot by weeks
# Plot on minute-level subject log-ac profile
# Colored based on EDSS intervals: Low/Medium/Hight <-> DeepBlue/LightBlue/Red
# Plot by week1, week2 etc.

co.EDSS = c('royalblue1','paleturquoise3','red')
for(i in 1:length(data.full)){
  df.name <- names(data.full)[i] # Subject ID
  
  name.ind <- clinical_info$Subject.ID == df.name
  
  # Color by EDSS score
  sub.edss <- clinical_info[name.ind,]$EDSS 
  if((sub.edss >= 1) && (sub.edss <=2))
    co.i <- co.EDSS[1]
  else{
    if((sub.edss >= 2.5) && (sub.edss <=3.5)){
      co.i <- co.EDSS[2]
    }else{co.i <- co.EDSS[3]}
  }
  
  # Tile: ID, EDSS score, T25fw score, TUG score
  maintitle <- paste("Activity Profile. ID: ", df.name,
                     ", Age: ", clinical_info$AGE[name.ind],
                     ", Gender: ", clinical_info$SEX[name.ind],
                     ", EDSS=", clinical_info$EDSS[name.ind], 
                     ", T25fw=", clinical_info$T25fw[name.ind], 
                     ", TUG=", clinical_info$TUG[name.ind],
                     sep='')
  
  # Multiple plots for n days for Subject.ID
  sub.by.week <- Sub.weeks(data.full[[i]])
  plt.weeks <- c()
  for(j in 1:length(sub.by.week)){
    plt.week.j <- c()
    df <- sub.by.week[[j]]
    
    for(k in 1:dim(df)[1]){
      t.range <- seq.POSIXt(df$date[k], by="min", length.out = 1440)
      d <- data.frame(x=t.range,y=t(df[k,1:1440]))
      names(d) <- c('time','lac')
      no.counts <- sum(d$lac>0)
      plt.week.j[[k]] <- ggplot(d, aes(x=time, ymin=0, ymax=lac)) + geom_linerange(col=co.i) + 
        ylab(paste(weekdays(df$date[k]),"LAC")) + theme(axis.title.x=element_blank())+
        geom_text(aes(time[100], 10), label=paste('number of counts:', no.counts))
    }
    plt.weeks[[j]] <- plt.week.j
  }
  
  n <- length(plt.weeks)
  lay.matrix <- plt.layout(plt.weeks)
  # Save as png file
  file.name.png=paste("plot/profile/Figure_1_", df.name, ".png", sep = "")
  png(file.name.png, width = n*700, height = 700)
  grid.arrange(grobs=unlist(plt.weeks,recursive=FALSE),layout_matrix=lay.matrix,top=textGrob(maintitle, gp=gpar(fontsize=24)))
  
  dev.off()
}


######## Figure 2 - ten summaries 

### 1. grouping individual 10 summaries by EDSS group
# Each subject has 10-summaries averaged by weekdays
# Individual 10-summaries are grouped to low/medium/high EDSS groups
# Note: this is before removing invalid days and outlier subjects
co.EDSS = c('royalblue1','paleturquoise3','red')
y.label = c("Mean of LAC",
            "Standard deviation of LAC",
            "Coefficient of variation",
            "Number of active minutes",
            "Mean LAC during Active Time",
            "Max6",
            "Total variability",
            "Fragmentation",
            "Relative amplitude",
            "Intra-daily variability")

low_edss <- c()
med_edss <- c()
high_edss <- c()
for(i in 1:length(data.full)){
  df <- stat.subject.generate(data.full[[i]])
  ten.stat <- ave.by.day(df)
  
  df.name <- names(data.full)[i]# Subject ID
  name.ind <- clinical_info$Subject.ID == df.name
  
  # Group ten.sum dataframes by EDSS score
  sub.edss <- clinical_info[name.ind,]$EDSS 
  if((sub.edss >= 1) && (sub.edss <=2)){
    low_edss[[df.name]] <- ten.stat
  }else{
    if((sub.edss >= 2.5) && (sub.edss <=3.5)){
      med_edss[[df.name]] <- ten.stat
    }else{
      high_edss[[df.name]] <- ten.stat}
  }
}

### 2. Average 10-sum in each edss group
ten.sum.low <- rbind.fill(low_edss)
ten.sum.low_mean <- aggregate(ten.sum.low[,2:11], by=list(ten.sum.low$day),mean)
ten.sum.low_mean$e.group <- "Low EDSS"
names(ten.sum.low_mean)[1] <- "days"

ten.sum.med <- rbind.fill(med_edss)
ten.sum.med_mean <- aggregate(ten.sum.med[,2:11], by=list(ten.sum.med$day),mean)
ten.sum.med_mean$e.group <- "Medium EDSS"
names(ten.sum.med_mean)[1] <- "days"

ten.sum.high <- rbind.fill(high_edss)
ten.sum.high_mean <- aggregate(ten.sum.high[,2:11], by=list(ten.sum.high$day),mean)
ten.sum.high_mean$e.group <- "High EDSS"
names(ten.sum.high_mean)[1] <- "days"

ten.sum.ave <- rbind.fill(ten.sum.low_mean,ten.sum.med_mean,ten.sum.high_mean)
ten.sum.ave$e.group <- factor(ten.sum.ave$e.group, levels = c('Low EDSS', 'Medium EDSS', 'High EDSS'))

# Plot low/medium/high Ten-summaries on the same plot
ten.plt <- c()
for(k in 1:10){
  d <- data.frame('days'=ten.sum.ave$days,'stat'=ten.sum.ave[,k+1], 'e.group'=ten.sum.ave$e.group )
  ten.plt[[k]] <- ggplot(d, aes(x=days, y=stat, group=e.group, colour=e.group))+geom_line(size=1)+
    scale_color_manual(name = "", values =co.EDSS)+
    ggtitle(y.label[k])+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'hidden',legend.text=element_text(size=20))
}
tensum.file.name <-'plot/Figure_2_full_10sum_ave.png'
png(tensum.file.name,width = 1600, height = 1200)
legend <- g_legend(ten.plt[[1]]+ theme(legend.position='bottom'))
ten.sum.title <- textGrob("Ten summaries by average EDSS group (before remove invalid days and outlier subject)", gp=gpar(fontsize=24))
grid.arrange(grobs=ten.plt, nrow=5, ncol=2, top=ten.sum.title, bottom=legend$grobs[[1]])
dev.off()

### 3. 10-summaries plot for each individual 
co.EDSS = c('royalblue1','paleturquoise3','red')
y.label = c("Mean of LAC",
            "Standard deviation of LAC",
            "Coefficient of variation",
            "Number of active minutes",
            "Mean LAC during Active Time",
            "Max6",
            "Total variability",
            "Fragmentation",
            "Relative amplitude",
            "Intra-daily variability")

ten.sum.list <- c(low_edss,med_edss,high_edss)
for(i in 1:length(ten.sum.list)){
  df <- ten.sum.list[[i]]
  df.name <- names(ten.sum.list)[i]# Subject ID
  
  name.ind <- clinical_info$Subject.ID == df.name
  
  # Color by EDSS score
  sub.edss <- clinical_info[name.ind,]$EDSS 
  if((sub.edss >= 1) && (sub.edss <=2))
    co.i <- co.EDSS[1]
  else{
    if((sub.edss >= 2.5) && (sub.edss <=3.5)){
      co.i <- co.EDSS[2]
    }else{co.i <- co.EDSS[3]}
  }
  
  
  maintitle <- paste("10 Summaries. ID: ", df.name,
                     ", Age: ", clinical_info$AGE[name.ind],
                     ", Gender: ", clinical_info$SEX[name.ind],
                     ", EDSS=", clinical_info$EDSS[name.ind],
                     ", T25fw=", clinical_info$T25fw[name.ind], 
                     ", TUG=", clinical_info$TUG[name.ind],
                     sep='')
  
  n <- 10
  ten.plt <- c()
  for(k in 1:n){
    d <- data.frame('days'=df$day, 'stat'=df[,k+1])
    ten.plt[[k]] <- ggplot(d, aes(x=days,y=stat, group = 1))+geom_line(size=1, col=co.i)+
      ggtitle(y.label[k])+
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text=element_text(size=14),
            axis.title=element_text(size=16,face="bold"),
            plot.title = element_text(hjust = 0.5))
  }
  tensum.file.name <- paste('plot/ten_sum/Figure_2_10sum_',df.name, '.png', sep='')
  png(tensum.file.name,width = 1600, height = 1200)
  ten.sum.title <- textGrob(maintitle, gp=gpar(fontsize=24))
  grid.arrange(grobs=ten.plt, nrow=5, ncol=2, top=ten.sum.title )
  dev.off()
}

######################################
#        Data Processing (Stage 2)
######################################
# Stage 2:
# Since we don't have non-valid subject (less than 3 days of data)
# At stage 2, we will remove invalid days for each subject 
# Non-valid days are defined as:  
# 1) Days with more than a consecutive of 2.5 hours non-wearing time 
# 2) Days with less than 480 active minutes (8 hrs) counts


# Remove invalid days, and re-assign valid data to each subject
for(i in 1:length(data.full)){
  df <- data.full[[i]]
  id <- names(data.full)[i]
  valid.df <- valid_full_days(df)
  assign(id, valid.df)
}

## Combine new valid dataframes into a list
data.valid <- mget(c(ls(pattern="^ACC"),ls(pattern="^TDPM"))) 
# We now have 64 subjects with data of valid days

######### End of Stage 2 Data Process ############

##########################################
##             Plots (Stage 2)          ##
##########################################
######### Figure 0 
## Distribution of valid/invalid days
full.counts <- c()
full.ID <- c()
valid_list <- c()
all_dates <- c()
for(i in 1:length(data.full)){
  df <- data.full[[i]]
  full.ID[[i]] <- rep(names(data.full)[i],dim(df)[1])
  sub.counts <- c()
  for(k in 1:dim(df)[1]){
    sub.counts[k] <- sum(df[k,1:1440]>0)
  }
  full.counts[[i]] <- sub.counts
  valid_list[[i]] <- df$date %in% data.valid[[i]]$date
  all_dates <- append(all_dates, df$date)
}

days.count <- data.frame(matrix(nrow=length(unlist(full.ID)),ncol = 3))
names(days.count) <- c("subject.ID","counts","status")
days.count$subject.ID <- unlist(full.ID)
days.count$counts <- unlist(full.counts)
days.count$dates <- all_dates

valid_days_list <- unlist(valid_list)
length(valid_days_list)
days.count$status[valid_days_list] <- "valid"
days.count$status[!valid_days_list] <- "invalid"

days.count$status <- factor(days.count$status, levels = c("valid","invalid"))

# Plot histogram on valid and invalid days minute counts for all subjects
col <- c("royalblue1","orangered")
png("plot/hist_counts.png",width = 1000, height = 600)
ggplot(data = days.count,aes(x=counts, fill=status))+geom_histogram(bins=30,alpha=0.5,colour="snow")+
  scale_x_continuous("Daily active-minute counts for all subjects", breaks = seq(0,1300,100))+
  scale_y_discrete("", limits=seq(0,160,20))+
  ggtitle("Histogram of Daily Active-Minute Counts")+
  theme(text = element_text(size=14),
        legend.text=element_text(size=14),
        plot.title = element_text(size = 22, hjust = 0.5))+   
  scale_fill_manual(name="",values=c("royalblue1","orangered"))
dev.off()

### Write invalid days into a text file
invalid_dates <- cbind(days.count[days.count$status == 'invalid',1],
            as.character(days.count[days.count$status == 'invalid',4]))
write.csv(invalid_dates , "invalid dates record.csv")

### Record number of valid days for each subject#
no_valid <- matrix(nrow=length(data.valid),ncol=2)
names(no_valid) <- c('subject.ID','no_days')
for(i in 1:length(data.valid)){
  no_valid[i,1] <- names(data.valid)[i]
  no_valid[i,2] <- dim(data.valid[[i]])[1]
}
write.csv(no_valid , "number valid.csv")

######## Figure 2. Replot 10-sum for valid data 

### 10 summaries plot for each individual of valid data
low_edss.valid <- c()
med_edss.valid <- c()
high_edss.valid <- c()
for(i in 1:length(data.valid)){
    df <- stat.subject.generate(data.valid[[i]])
    ten.stat <- ave.by.day(df)
    
    df.name <- names(data.valid)[i]# Subject ID
    name.ind <- clinical_info$Subject.ID == df.name
    
    # Group ten.sum dataframes by EDSS score
    sub.edss <- clinical_info[name.ind,]$EDSS 
    if((sub.edss >= 1) && (sub.edss <=2)){
      low_edss.valid[[df.name]] <- ten.stat
    }else{
      if((sub.edss >= 2.5) && (sub.edss <=3.5)){
        med_edss.valid[[df.name]] <- ten.stat
      }else{
        high_edss.valid[[df.name]] <- ten.stat}
    }
  }

# Plot 10-summaries for valid days data for each subject
co.EDSS = c('royalblue1','paleturquoise3','red')
y.label = c("Mean of LAC",
            "Standard deviation of LAC",
            "Coefficient of variation",
            "Number of active minutes",
            "Mean LAC during Active Time",
            "Max6",
            "Total variability",
            "Fragmentation",
            "Relative amplitude",
            "Intra-daily variability")

valid_ten.sum <- c(low_edss.valid,med_edss.valid,high_edss.valid)
for(i in 1:length(valid_ten.sum)){
  df <- valid_ten.sum[[i]]
  df.name <- names(valid_ten.sum)[i]# Subject ID
  
  name.ind <- clinical_info$Subject.ID == df.name
  
  # Color by EDSS score
  sub.edss <- clinical_info[name.ind,]$EDSS 
  if((sub.edss >= 1) && (sub.edss <=2))
    co.i <- co.EDSS[1]
  else{
    if((sub.edss >= 2.5) && (sub.edss <=3.5)){
      co.i <- co.EDSS[2]
    }else{co.i <- co.EDSS[3]}
  }
  
  
  maintitle <- paste("10 Summaries. ID: ", df.name,
                     ", Age: ", clinical_info$AGE[name.ind],
                     ", Gender: ", clinical_info$SEX[name.ind],
                     ", EDSS=", clinical_info$EDSS[name.ind],
                     ", T25fw=", clinical_info$T25fw[name.ind], 
                     ", TUG=", clinical_info$TUG[name.ind],
                     sep='')
  
  n <- 10
  ten.plt <- c()
  for(k in 1:n){
    d <- data.frame('days'=df$day, 'stat'=df[,k+1])
    ten.plt[[k]] <- ggplot(d, aes(x=days,y=stat, group = 1))+geom_line(size=1, col=co.i)+
      ggtitle(y.label[k])+
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text=element_text(size=14),
            axis.title=element_text(size=16,face="bold"),
            plot.title = element_text(hjust = 0.5))
  }
  tensum.file.name <- paste('plot/ten_sum_valid/Figure_2_10sum_valid_',df.name, '.png', sep='')
  png(tensum.file.name,width = 1600, height = 1200)
  ten.sum.title <- textGrob(maintitle, gp=gpar(fontsize=24))
  grid.arrange(grobs=ten.plt, nrow=5, ncol=2, top=ten.sum.title )
  dev.off()
}


### 2. Average ten summaries 
# Average of ten.sum by EDSS group: ONLY VALID SUBJECTS with valid dates

# Ten.sum by average of edss group
ten.sum.low <- rbind.fill(low_edss.valid)
ten.sum.low_mean <- aggregate(ten.sum.low[,2:11], by=list(ten.sum.low$day),mean)
ten.sum.low_mean$e.group <- "Low EDSS"
names(ten.sum.low_mean)[1] <- "days"

ten.sum.med <- rbind.fill(med_edss.valid)
ten.sum.med_mean <- aggregate(ten.sum.med[,2:11], by=list(ten.sum.med$day),mean)
ten.sum.med_mean$e.group <- "Medium EDSS"
names(ten.sum.med_mean)[1] <- "days"

ten.sum.high <- rbind.fill(high_edss.valid)
ten.sum.high_mean <- aggregate(ten.sum.high[,2:11], by=list(ten.sum.high$day),mean)
ten.sum.high_mean$e.group <- "High EDSS"
names(ten.sum.high_mean)[1] <- "days"

ten.sum.ave <- rbind.fill(ten.sum.low_mean,ten.sum.med_mean,ten.sum.high_mean)
ten.sum.ave$e.group <- factor(ten.sum.ave$e.group, levels = c('Low EDSS', 'Medium EDSS', 'High EDSS'))

# Plot Ten.sum by average of edss group
co.EDSS = c('royalblue1','paleturquoise3','red')
y.label = c("Mean of LAC",
            "Standard deviation of LAC",
            "Coefficient of variation",
            "Number of active minutes",
            "Mean LAC during Active Time",
            "Max6",
            "Total variability",
            "Fragmentation",
            "Relative amplitude",
            "Intra-daily variability")
ten.plt.valid <- c()
for(k in 1:10){
  d <- data.frame('days'=ten.sum.ave$days,'stat'=ten.sum.ave[,k+1], 'e.group'=ten.sum.ave$e.group )
  ten.plt.valid[[k]] <- ggplot(d, aes(x=days, y=stat, group=e.group, colour=e.group))+geom_line(size=1)+
    scale_color_manual(name = "", values =co.EDSS)+
    ggtitle(y.label[k])+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'hidden',legend.text=element_text(size=20))
}
tensum.file.name <-'plot/Figure_2_valid_10sum_ave.png'
png(tensum.file.name,width = 1600, height = 1200)
legend <- g_legend(ten.plt.valid[[1]]+ theme(legend.position='bottom'))
ten.sum.title <- textGrob("Ten summaries by average EDSS group (valid days, outlier removed)", gp=gpar(fontsize=24))
grid.arrange(grobs=ten.plt.valid, nrow=5, ncol=2, top=ten.sum.title, bottom=legend$grobs[[1]])
dev.off()

### 3. Plot only low & high EDSS group 
ten.sum.ave <- rbind.fill(ten.sum.low_mean, ten.sum.high_mean)
ten.sum.ave$e.group <- factor(ten.sum.ave$e.group, levels = c('Low EDSS', 'High EDSS'))

co.EDSS = c('royalblue1','red')
ten.plt.valid <- c()
for(k in 1:10){
  d <- data.frame('days'=ten.sum.ave$days,'stat'=ten.sum.ave[,k+1], 'e.group'=ten.sum.ave$e.group )
  ten.plt.valid[[k]] <- ggplot(d, aes(x=days, y=stat, group=e.group, colour=e.group))+geom_line(size=1)+
    scale_color_manual(name = "", values =co.EDSS)+
    ggtitle(y.label[k])+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'hidden',legend.text=element_text(size=20))
}
tensum.file.name <-'plot/Figure_2_valid_10sum_ave_binary.png'
png(tensum.file.name,width = 1600, height = 1200)
legend <- g_legend(ten.plt.valid[[1]]+ theme(legend.position='bottom'))
ten.sum.title <- textGrob("Ten summaries by average EDSS group (valid days, outlier removed)", gp=gpar(fontsize=24))
grid.arrange(grobs=ten.plt.valid, nrow=5, ncol=2, top=ten.sum.title, bottom=legend$grobs[[1]])
dev.off()

######################################
#        Data Processing (Stage 3)
######################################
# Stage 3:
# At stage 3, we will remove outlier subjects
# One subject (ACC003) with T25fw>=75 seconds, and detected by Cook's distance
# One subject (TDPMSA042) detected by Cook's distance

# Drop outlier subject from list: ACC003
data.valid$ACC003 <- NULL

# Drop additional outlier subject from list: TDPMSA042 
data.valid$TDPMSA042 <- NULL

length(data.valid) # We now have 62 subjects with data of valid days


#########  Linear regression models ##########

#  Process activity data as TLAC10 (total log-ac at 10min interval)
df.mlac <- data.frame(matrix(nrow=length(data.valid), ncol=145)) 
for(i in 1:dim(df.mlac)[1]){
  if(dim(data.valid[[i]])[1] > 3){  # at least three valid days
    df.mlac[i,1] <- names(data.valid)[i]
    mean.over.days <- apply(data.valid[[i]][,1:1440],2,mean,na.rm=TRUE)
    df.mlac[i,2:145] <- rollapply(mean.over.days,width=10,by=10,FUN=sum) # sum over 10 minutes intervals
  }
}
names(df.mlac)[1] <- "Subject.ID"
df.mlac <- na.omit(df.mlac)
dim(df.mlac) 
# We have 64 subjects before remove outlier; Remove ACC003 and TDPMSA042, 62 subjects

# Work with a clean dataframe of linear-model results and predictor variables for each subject
ac <- data.frame(matrix(nrow = dim(df.mlac)[1],ncol=14))
colnames(ac) <- c("Subject.ID","Mean.lac_24hrs","Mean.lac_interval1","Mean.lac_interval2",
                  "Mean.lac_interval3","Mean.lac_interval4","Mean.lac_interval5",
                  "Mean.lac_interval6","Age","Sex","BMI","EDSS","T25fw","TUG.t2")

ac$Subject.ID <- df.mlac$Subject.ID
ac$Mean.lac_24hrs <- as.numeric(rowMeans(df.mlac[,2:145], na.rm = FALSE))

for(i in 1:dim(df.mlac)[1]){
  ac[i,3:8] <- rollapply(as.numeric(df.mlac[i,2:145]),width=24,by=24,FUN=mean)
}

for(i in 1:dim(clinical_info)[1]){
  for(j in 1:dim(ac)[1]){
    if(ac$Subject.ID[j] == clinical_info$Subject.ID[i]){
      ac$Age[j] <- clinical_info$AGE[i]
      ac$Sex[j] <- clinical_info$SEX[i]
      ac$BMI[j] <- clinical_info$BMI[i]
      ac$EDSS[j] <- clinical_info$EDSS[i]
      ac$T25fw[j] <- clinical_info$T25fw[i]
      ac$TUG.t2[j] <- clinical_info$TUG[i]
    }
  }
}
# Normalize age and bmi
ave.age <- mean(ac$Age)
ave.bmi <- mean(ac$BMI)

ac$Age <- ac$Age-ave.age
ac$BMI <- ac$BMI - ave.bmi

# Convert SEX to factors
ac$Sex <- factor(ac$Sex, levels = c('Female','Male'))

######### End of Stage 3 Data Process ############

##########################################
##             Plots (Stage 3)          ##
##########################################

#### Figure 3. Linear model plot on 24-hr intervals
# Linear models
fit.EDSS <- lm(Mean.lac_24hrs ~ EDSS + Age+Sex + BMI, data=ac)
fit.T25 <- lm(Mean.lac_24hrs ~ T25fw + Age+Sex + BMI, data=ac)
fit.TUG <- lm(Mean.lac_24hrs ~ TUG.t2 + Age+Sex + BMI, data=ac)

# 1. Plot linear model diagnostics
png(file = "plot/Figure_24hrEDSS_diag.png", width = 850, height = 700)
autoplot(fit.EDSS, which = 1:6, ncol = 2, nrow=3, label.size = 3)
dev.off()

png(file = "plot/Figure_24hrT25_diag.png", width = 850, height = 700)
autoplot(fit.T25, which = 1:6, ncol = 2, nrow=3, label.size = 3)
dev.off()

png(file = "plot/Figure_24hrTUG_diag.png", width = 850, height = 700)
autoplot(fit.TUG, which = 1:6, ncol = 2, nrow=3, label.size = 3)
dev.off()

# 2. Summary linear model results
sink("lm_24hr_2OUTLIERS.txt")
print(summary(fit.EDSS))
print(summary(fit.T25))
print(summary(fit.TUG))
sink()

# 3. Plot 24hr model by time of the day
# Plot raw log-ac data background 
# Plot linear model estimations for EDSS/T25fw/TUG 

# Subjects' profile data
df.mlac$Subject.ID <- as.factor(df.mlac$Subject.ID)
df.sub <- as.data.frame(t(df.mlac[,2:145]))
df.sub$time <- seq.POSIXt(as.POSIXct("2010-01-01 00:00:00"),by="10 min", length.out = 144)
names(df.sub)[1:62] <- df.mlac$Subject.ID

# Melt TLAC10 data for plotting
meltdf <- melt(df.sub,id="time")
colnames(meltdf) <- c("time","Subject.ID","Mean.log.ac")

# Colors
co <- colorRampPalette(brewer.pal(9,"Blues"))(9)[2:9]
co[1] <- "gray0"
col.sub <- grDevices::rainbow(length(data.valid))
df.color <- c()
for(i in 1:length(data.valid)){
  df.color[[i]] <- rep(col.sub[i],144)
}
meltdf$color <- unlist(df.color)

## Initiate dataframe to record estimated values
y.est <- data.frame(matrix(nrow=8,ncol=2))
colnames(y.est) <- c("EDSS","Y.EDSS")

# 1) Plot EDSS on 24hr interval
y.est$EDSS <- seq(0,7,1)
y.est$Y.EDSS <- coef(fit.EDSS)['(Intercept)']+coef(fit.EDSS)['EDSS']*y.est$EDSS
a <- signif(coef(fit.EDSS)['(Intercept)'], digits = 2)
b <- signif(coef(fit.EDSS)['EDSS'], digits = 2)
R2 <- signif(summary(fit.EDSS)$adj.r.squared,3)
lm_pvalue <- signif(summary(fit.EDSS)$coefficients[,4]['EDSS'],3)
textlab <- paste("Tlac10 = ",a,"+(", b,")*EDSS", ", Adj.R^2=",R2, ", p.value=",lm_pvalue, sep="")

colour <- c("Baseline")
for(i in 1:7){
  colour[i+1] <- paste("EDSS=",i,sep="")
}

plot.df.edss <- data.frame(
  x = rep(df.sub$time[1],8),
  xend = rep(df.sub$time[144],8),
  y = y.est$Y.EDSS,
  yend=y.est$Y.EDSS,
  colour=colour
)


plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #colour=meltdf$color
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("TLAC10 Linear Model Over 24hr - EDSS")+
  theme(axis.text.x= element_text(angle=90),
        axis.text=element_text(size=12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

plt2 <- plt+geom_segment(data = plot.df.edss, aes(x=x, xend=xend, 
                              y=y,yend=yend, colour=colour), size=1.2,inherit.aes =FALSE)

plt3 <- plt2+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)+
   annotate("text",x=df.sub$time[20], y=75, label = textlab, color="black", size=4, parse=FALSE)

png(file = "plot/Figure_3_24hrEDSS_2OUTLIER.png", width = 1000, height = 450)
plt3
dev.off()
 
# 2) Plot T25fw on 24hr interval
coef(fit.T25)
y.est$T25fw <- NA
y.est$T25fw[1:6] <- seq(0,25,by=5)
y.est$Y.T25fw<- coef(fit.T25)['(Intercept)']+coef(fit.T25)['T25fw']*y.est$T25fw
a <- signif(coef(fit.T25)['(Intercept)'], digits = 2)
b <- signif(coef(fit.T25)['T25fw'], digits = 2)
R2 <- signif(summary(fit.T25)$adj.r.squared,3)
lm_pvalue <- signif(summary(fit.T25)$coefficients[,4]['T25fw'],3)
textlab <- paste("Tlac10 = ",a,"+(", b,")*T25fw", ", Adj.R^2=",R2, ", p.value=",lm_pvalue, sep="")

colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("T25fw=",y.est$T25fw[i+1],sep="")
}

plot.df.t25 <- data.frame(
  x = rep(df.sub$time[1],6),
  xend = rep(df.sub$time[144],6),
  y = y.est$Y.T25fw[complete.cases(y.est$Y.T25fw)],
  yend=y.est$Y.T25fw[complete.cases(y.est$Y.T25fw)],
  colour=colour
)
plot.df.t25$colour <- factor(plot.df.t25$colour ,levels = colour)

plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #colour=meltdf$color
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("TLAC10 Linear Model Over 24hr - T25fw")+
  theme(axis.text.x= element_text(angle=90),
        axis.text=element_text(size=12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

plt2 <- plt+geom_segment(data = plot.df.t25, aes(x=x, xend=xend, 
                                                 y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt3 <- plt2+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)+
  annotate("text",x=df.sub$time[20], y=75, label = textlab, color="black", size=4, parse=FALSE)

png(file = "plot/Figure_3_24hrT25fw_2OUTLIER.png", width = 1000, height = 450)
plt3
dev.off()


# 3) Plot TUG on 24hr interval
coef(fit.TUG)
y.est$TUG <- NA
y.est$TUG[1:6]<- seq(0,25,5)
y.est$Y.TUG<- coef(fit.TUG)['(Intercept)']+coef(fit.TUG)['TUG.t2']*y.est$TUG
a <- signif(coef(fit.TUG)['(Intercept)'], digits = 2)
b <- signif(coef(fit.TUG)['TUG.t2'], digits = 2)
R2 <- signif(summary(fit.TUG)$adj.r.squared,3)
lm_pvalue <- signif(summary(fit.TUG)$coefficients[,4]['TUG.t2'],3)
textlab <- paste("Tlac10 = ",a, "+(", b,")*TUG", ", Adj.R^2=",R2, ", p.value=",lm_pvalue, sep="")

colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("TUG=",y.est$TUG[i+1],sep="")
}

plot.df.tug <- data.frame(
  x = rep(df.sub$time[1],6),
  xend = rep(df.sub$time[144],6),
  y = y.est$Y.TUG[complete.cases(y.est$Y.TUG)],
  yend=y.est$Y.TUG[complete.cases(y.est$Y.TUG)],
  colour=colour
)
plot.df.tug$colour <- factor(plot.df.tug$colour ,levels = colour)

plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts") +  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("TLAC10 Linear Model Over 24hr - TUG")+
  theme(axis.text.x = element_text(angle=90),
        axis.text= element_text(size=12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

plt2 <- plt+geom_segment(data = plot.df.tug, aes(x=x, xend=xend, 
                                                y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt3 <- plt2+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)+
  annotate("text",x=df.sub$time[20], y=75, label = textlab, color="black", size=4, parse=FALSE)

png(file = "plot/Figure_3_24hrTug_2OUTLIER.png", width = 1000, height = 450)
plt3
dev.off()

####### Figure 4
# Plot Linear Regression Model 
p1 <- ggplotRegression(fit.EDSS)+xlab('EDSS score')+ylab('MLAC over 24 hrs')
p2 <- ggplotRegression(fit.T25)+xlab('T25fw seconds')+ylab('MLAC over 24 hrs')
p3 <- ggplotRegression(fit.TUG)+xlab('TUG seconds')+ylab('MLAC over 24 hrs')
png('plot/Figure4_LM_24hr_2OUTLIER.png', width = 800, height = 800)
grid.arrange(p1,p2,p3,nrow=3,top=textGrob('Linear Model of TLAC10 over 24 hrs', gp=gpar(fontsize=16)))
dev.off()

####### Figure 5. Linear model plot on 4-hr intervals

######### EDSS ############
# Linear models for EDSS
for(i in 1:6){
  assign(paste("edss4hr.fit", i, sep = ""),lm(ac[,i+2]~EDSS + Age + Sex + BMI, data=ac))
}

edss4hr.fit.models <- c(list(edss4hr.fit1),list(edss4hr.fit2),list(edss4hr.fit3),
                        list(edss4hr.fit4),list(edss4hr.fit5),list(edss4hr.fit6))

# Summary linear models
sink("edss_lm_4hr.txt")
for(i in 1:6){
  print(summary(edss4hr.fit.models[[i]]))
}
sink()

# 1) plot EDSS over evenly spaced 4hr interval
### Estimated EDSS values over 4hr intervals
edss.est_4hr <- data.frame(matrix(nrow=8,ncol=7))
names(edss.est_4hr) <- c("EDSS", "lm1","lm2","lm3","lm4","lm5","lm6")
edss.est_4hr$EDSS <- seq(0,7)
for(i in 1:6){
  edss.int <- coef(edss4hr.fit.models[[i]])['(Intercept)']
  edss.beta <- coef(edss4hr.fit.models[[i]])['EDSS']
  edss.est_4hr[,i+1] <- edss.int+edss.beta*edss.est_4hr$EDSS
}

colour <- c("Baseline")
for(i in 1:7){
  colour[i+1] <- paste("EDSS=",i,sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,6)

x0 <- c(rep(df.sub$time[1],8),rep(df.sub$time[4*6+1],8), rep(df.sub$time[8*6+1],8),
        rep(df.sub$time[12*6+1],8),rep(df.sub$time[16*6+1],8), rep(df.sub$time[20*6+1],8))
x1 <- c(rep(df.sub$time[4*6+1],8),rep(df.sub$time[8*6+1],8), rep(df.sub$time[12*6+1],8),
        rep(df.sub$time[16*6+1],8),rep(df.sub$time[20*6+1],8), rep(df.sub$time[24*6],8))

y0 <- c()
for(i in 1:6){
  y0[[i]] <- edss.est_4hr[,i+1]
}
y0 <- unlist(y0)

plot.edss.4hr <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y0,
  colour=colour
)
co <- colorRampPalette(brewer.pal(9,"Blues"))(9)[2:9]
co[1] <- "gray0"
plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ # colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("TLAC10 Linear Model Over 4hr Interval - EDSS")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.edss.4hr, aes(x=x, xend=xend, 
                                                 y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)

plt3 <- plt2+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[4*6+1],df.sub$time[8*6+1],df.sub$time[12*6+1],
                                       df.sub$time[16*6+1],df.sub$time[20*6+1],df.sub$time[24*6]), alpha=0.5)

png(file = "plot/Figure_5_4hrEDSS.png", width = 1000, height = 450)
plt3
dev.off()

# 2) plot EDSS over connected 4hr interval
colour <- c("Baseline")
for(i in 1:7){
  colour[i+1] <- paste("EDSS=",i,sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,7)

x0 <- c(rep(df.sub$time[1],8),rep(df.sub$time[2*6+1],8), rep(df.sub$time[6*6+1],8),
        rep(df.sub$time[10*6+1],8),rep(df.sub$time[14*6+1],8), rep(df.sub$time[18*6+1],8),
        rep(df.sub$time[22*6+1],8))
x1 <- c(rep(df.sub$time[2*6+1],8),rep(df.sub$time[6*6+1],8), rep(df.sub$time[10*6+1],8),
        rep(df.sub$time[14*6+1],8),rep(df.sub$time[18*6+1],8), rep(df.sub$time[22*6+1],8),
        rep(df.sub$time[24*6],8))

y0 <- c(edss.est_4hr$lm1,edss.est_4hr$lm1,edss.est_4hr$lm2,edss.est_4hr$lm3,
        edss.est_4hr$lm4,edss.est_4hr$lm5,edss.est_4hr$lm6)
y1 <- c(edss.est_4hr$lm1,edss.est_4hr$lm2,edss.est_4hr$lm3,
        edss.est_4hr$lm4,edss.est_4hr$lm5,edss.est_4hr$lm6,edss.est_4hr$lm6)


plot.edss.4hr_con <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y1,
  colour=colour
)

plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ # colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("Connected Plot: TLAC10 Linear Model Over 4hr Interval - EDSS")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.edss.4hr_con, aes(x=x, xend=xend, 
                                                  y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)+
  geom_vline(xintercept = c(df.sub$time[1],df.sub$time[4*6+1],
                            df.sub$time[8*6+1],df.sub$time[12*6+1],
                            df.sub$time[16*6+1],df.sub$time[20*6+1],df.sub$time[24*6]), alpha=0.5)


EDSS.pvalue <- data.frame(matrix(nrow=length(edss4hr.fit.models),ncol=2))
names(EDSS.pvalue) <- c('Interval','p.value')
for(i in 1:length(edss4hr.fit.models)){
  EDSS.pvalue[i,1] <- paste('lm',i,sep="")
  EDSS.pvalue[i,2] <- summary(edss4hr.fit.models[[i]])$coef[2,4]
}
EDSS.pvalue[EDSS.pvalue$p.value <0.05,1]

int3 <- seq.POSIXt(df.sub$time[12*6+1], by='10 min', length.out = 24)
int4 <- seq.POSIXt(df.sub$time[16*6+1], by='10 min', length.out = 24)
significant.range <- c(int3,int4)

plt3 <- plt2+geom_vline(xintercept=significant.range, colour='slategray2',size=4, alpha=0.25)

png(file = "plot/Figure_6_4hrEDSS_con.png", width = 1000, height = 450)
plt3
dev.off()


############### T25fw ###################
# Linear models
for(i in 1:6){
  assign(paste("t25fw4hr.fit", i, sep = ""),lm(ac[,i+2]~T25fw + Age + Sex + BMI, data=ac))
}

t25fw4hr.fit.models <- mget(ls(pattern="t25fw4hr.fit"))

# Summary linear regression
sink("t25fw_lm_4hr.txt")
for(i in 1:6){
  print(summary(t25fw4hr.fit.models[[i]]))
}
sink()

# 1) plot T25fw over evenly spaced 4hr interval
### Estimated T25fw values over 4hr intervals
t25fw.est_4hr <- data.frame(matrix(nrow=6,ncol=7))
names(t25fw.est_4hr) <- c("T25fw", "lm1","lm2","lm3","lm4","lm5","lm6")
t25fw.est_4hr$T25fw <- seq(0,25,5)
for(i in 1:6){
  t25fw.int <- coef(t25fw4hr.fit.models[[i]])['(Intercept)']
  t25fw.beta <- coef(t25fw4hr.fit.models[[i]])['T25fw']
  t25fw.est_4hr[,i+1] <- t25fw.int+t25fw.beta*t25fw.est_4hr$T25fw
}

colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("T25fw=",t25fw.est_4hr$T25fw[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,6)

x0 <- c(rep(df.sub$time[1],6),rep(df.sub$time[4*6+1],6), rep(df.sub$time[8*6+1],6),
        rep(df.sub$time[12*6+1],6),rep(df.sub$time[16*6+1],6), rep(df.sub$time[20*6+1],6))
x1 <- c(rep(df.sub$time[4*6+1],6),rep(df.sub$time[8*6+1],6), rep(df.sub$time[12*6+1],6),
        rep(df.sub$time[16*6+1],6),rep(df.sub$time[20*6+1],6), rep(df.sub$time[24*6],6))

y0 <- c()
for(i in 1:6){
  y0[[i]] <- t25fw.est_4hr[,i+1]
}
y0 <- unlist(y0)

plot.t25fw.4hr <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y0,
  colour=colour
)
co =c("gray0", "#9ECAE1" , "#6BAED6", "#4292C6", "#2171B5", "#08519C")
plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #, colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M") +ggtitle ("TLAC10 Linear Model Over 4hr Interval - T25fw")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.t25fw.4hr, aes(x=x, xend=xend, 
                                                  y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)

plt3 <- plt2+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[4*6+1],df.sub$time[8*6+1],df.sub$time[12*6+1],
                                       df.sub$time[16*6+1],df.sub$time[20*6+1],df.sub$time[24*6]), alpha=0.5)
png(file = "plot/Figure_5_4hrT25fw.png", width = 1000, height = 450)
plt3
dev.off()

# 2) plot T25fw over connected 4hr interval
T25.pvalue <- data.frame(matrix(nrow=length(t25fw4hr.fit.models),ncol=2))
names(T25.pvalue) <- c('Interval','p.value')
for(i in 1:length(t25fw4hr.fit.models)){
  T25.pvalue[i,1] <- paste('lm',i,sep="")
  T25.pvalue[i,2] <- summary(t25fw4hr.fit.models[[i]])$coef[2,4]
}
T25.pvalue[T25.pvalue$p.value <0.05,1]

int1 <- seq.POSIXt(df.sub$time[1], by='10 min', length.out = 24)
int3 <- seq.POSIXt(df.sub$time[8*6+1], by='10 min', length.out = 24)
significant.range <- c(int1,int3)

colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("T25fw=",t25fw.est_4hr$T25fw[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,7)

x0 <- c(rep(df.sub$time[1],6),rep(df.sub$time[2*6+1],6), rep(df.sub$time[6*6+1],6),
        rep(df.sub$time[10*6+1],6),rep(df.sub$time[14*6+1],6), rep(df.sub$time[18*6+1],6),
        rep(df.sub$time[22*6+1],6))
x1 <- c(rep(df.sub$time[2*6+1],6),rep(df.sub$time[6*6+1],6), rep(df.sub$time[10*6+1],6),
        rep(df.sub$time[14*6+1],6),rep(df.sub$time[18*6+1],6), rep(df.sub$time[22*6+1],6),
        rep(df.sub$time[24*6],6))

y0 <- c(t25fw.est_4hr$lm1,t25fw.est_4hr$lm1,t25fw.est_4hr$lm2,t25fw.est_4hr$lm3,
        t25fw.est_4hr$lm4,t25fw.est_4hr$lm5,t25fw.est_4hr$lm6)
y1 <- c(t25fw.est_4hr$lm1,t25fw.est_4hr$lm2,t25fw.est_4hr$lm3,
        t25fw.est_4hr$lm4,t25fw.est_4hr$lm5,t25fw.est_4hr$lm6,t25fw.est_4hr$lm6)


plot.t25fw.4hr_con <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y1,
  colour=colour
)
co =c("gray0", "#9ECAE1" , "#6BAED6", "#4292C6", "#2171B5", "#08519C")
plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #, colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M") +
  ggtitle ("Connected Plot: TLAC10 Linear Model Over 4hr Interval - T25fw")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.t25fw.4hr_con, aes(x=x, xend=xend, 
                                                      y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)+
  geom_vline(xintercept = c(df.sub$time[1],df.sub$time[4*6+1],
                            df.sub$time[8*6+1],df.sub$time[12*6+1],
                            df.sub$time[16*6+1],df.sub$time[20*6+1],df.sub$time[24*6]), alpha=0.5)

plt3 <- plt2+geom_vline(xintercept=significant.range, colour='slategray2',size=4, alpha=0.25)
png(file = "plot/Figure_6_4hrT25fw_con.png", width = 1000, height = 450)
plt3
dev.off()

############### TUG ###################
# Linear models
for(i in 1:6){
  assign(paste("tug4hr.fit", i, sep = ""),lm(ac[,i+2]~TUG.t2 + Age + Sex + BMI, data=ac))
}

tug4hr.fit.models <- mget(ls(pattern="tug4hr.fit"))

# Summary linear regression
sink("tug_lm_4hr.txt")
for(i in 1:6){
  print(summary(tug4hr.fit.models[[i]]))
}
sink()


# 1) plot TUG over evenly spaced 4hr interval
### Estimated TUG values over 4hr intervals
tug.est_4hr <- data.frame(matrix(nrow=6,ncol=7))
names(tug.est_4hr) <- c("tug", "lm1","lm2","lm3","lm4","lm5","lm6")
tug.est_4hr$tug <- seq(0,25,5)
for(i in 1:6){
  tug.int <- coef(tug4hr.fit.models[[i]])['(Intercept)']
  tug.beta <- coef(tug4hr.fit.models[[i]])['TUG.t2']
  tug.est_4hr[,i+1] <- tug.int+tug.beta*tug.est_4hr$tug
}

colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("TUG=",tug.est_4hr$tug[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,6)

x0 <- c(rep(df.sub$time[1],6),rep(df.sub$time[4*6+1],6), rep(df.sub$time[8*6+1],6),
        rep(df.sub$time[12*6+1],6),rep(df.sub$time[16*6+1],6), rep(df.sub$time[20*6+1],6))
x1 <- c(rep(df.sub$time[4*6+1],6),rep(df.sub$time[8*6+1],6), rep(df.sub$time[12*6+1],6),
        rep(df.sub$time[16*6+1],6),rep(df.sub$time[20*6+1],6), rep(df.sub$time[24*6],6))

y0 <- c()
for(i in 1:6){
  y0[[i]] <- tug.est_4hr[,i+1]
}
y0 <- unlist(y0)

plot.tug.4hr <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y0,
  colour=colour
)
co =c("gray0", "#9ECAE1" , "#6BAED6", "#4292C6", "#2171B5", "#08519C")
plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #, colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("TLAC10 Linear Model Over 4hr Interval - TUG")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.tug.4hr, aes(x=x, xend=xend, 
                                                   y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)

plt3 <- plt2+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[4*6+1],df.sub$time[8*6+1],df.sub$time[12*6+1],
                                       df.sub$time[16*6+1],df.sub$time[20*6+1],df.sub$time[24*6]),alpha=0.5)
png(file = "plot/Figure_5_4hrTUG.png", width = 1000, height = 450)
plt3
dev.off()

# 2) plot TUG over connected 4hr interval
TUG.pvalue <- data.frame(matrix(nrow=length(tug4hr.fit.models),ncol=2))
names(TUG.pvalue) <- c('Interval','p.value')
for(i in 1:length(tug4hr.fit.models)){
  TUG.pvalue[i,1] <- paste('lm',i,sep="")
  TUG.pvalue[i,2] <- summary(tug4hr.fit.models[[i]])$coef[2,4]
}
TUG.pvalue[TUG.pvalue$p.value <0.05,1]

int1 <- seq.POSIXt(df.sub$time[1], by='10 min', length.out = 24)
int3 <- seq.POSIXt(df.sub$time[8*6+1], by='10 min', length.out = 24)
significant.range <- c(int1,int3)

colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("TUG=",tug.est_4hr$tug[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,7)

x0 <- c(rep(df.sub$time[1],6),rep(df.sub$time[2*6+1],6), rep(df.sub$time[6*6+1],6),
        rep(df.sub$time[10*6+1],6),rep(df.sub$time[14*6+1],6), rep(df.sub$time[18*6+1],6),
        rep(df.sub$time[22*6+1],6))
x1 <- c(rep(df.sub$time[2*6+1],6),rep(df.sub$time[6*6+1],6), rep(df.sub$time[10*6+1],6),
        rep(df.sub$time[14*6+1],6),rep(df.sub$time[18*6+1],6), rep(df.sub$time[22*6+1],6),
        rep(df.sub$time[24*6],6))

y0 <- c(tug.est_4hr$lm1,tug.est_4hr$lm1,tug.est_4hr$lm2,tug.est_4hr$lm3,
        tug.est_4hr$lm4,tug.est_4hr$lm5,tug.est_4hr$lm6)
y1 <- c(tug.est_4hr$lm1,tug.est_4hr$lm2,tug.est_4hr$lm3,
        tug.est_4hr$lm4,tug.est_4hr$lm5,tug.est_4hr$lm6,tug.est_4hr$lm6)


plot.t25fw.4hr_con <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y1,
  colour=colour
)
co =c("gray0", "#9ECAE1" , "#6BAED6", "#4292C6", "#2171B5", "#08519C")
plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #, colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("Connected Plot: TLAC10 Linear Model Over 4hr Interval - TUG")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.t25fw.4hr_con, aes(x=x, xend=xend, 
                                                       y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)+
  geom_vline(xintercept = c(df.sub$time[1],df.sub$time[4*6+1],
                            df.sub$time[8*6+1],df.sub$time[12*6+1],
                            df.sub$time[16*6+1],df.sub$time[20*6+1],df.sub$time[24*6]), alpha=0.5)

plt3 <- plt2+geom_vline(xintercept=significant.range, colour='slategray2',size=4, alpha=0.25)

png(file = "plot/Figure_6_4hrTUG_con.png", width = 1000, height = 450)
plt3
dev.off()

####### Figure 7&8. Linear model plot on 2-hr intervals

# Work on a clean dataframe of 2-hr estimates and predictor varibales for each subject
ac.2hr<- data.frame(matrix(nrow = dim(df.mlac)[1],ncol=19))

MLAC.2hr <- c()
for(i in 1:12){
  MLAC.2hr[i] <- paste("MLAC_2hr_int",i,sep="")
}

colnames(ac.2hr) <- c("Subject.ID",MLAC.2hr,"Age","Sex","BMI","EDSS","T25fw","TUG.t2")
ac.2hr$Subject.ID <- df.mlac$Subject.ID


for(i in 1:dim(df.mlac)[1]){
  ac.2hr[i,2:13] <- rollapply(as.numeric(df.mlac[i,2:145]),width=12,by=12,FUN=mean)
}

for(i in 1:dim(clinical_info)[1]){
  for(j in 1:dim(ac.2hr)[1]){
    if(ac.2hr$Subject.ID[j] == clinical_info$Subject.ID[i]){
      ac.2hr$Age[j] <- clinical_info$AGE[i]
      ac.2hr$Sex[j] <- clinical_info$SEX[i]
      ac.2hr$BMI[j] <- clinical_info$BMI[i]
      ac.2hr$EDSS[j] <- clinical_info$EDSS[i]
      ac.2hr$T25fw[j] <- clinical_info$T25fw[i]
      ac.2hr$TUG.t2[j] <- clinical_info$TUG[i]
    }
  }
}


# Normalize age and bmi
ave.age <- mean(ac.2hr$Age)
ave.bmi <- mean(ac.2hr$BMI)

ac.2hr$Age <- ac.2hr$Age-ave.age
ac.2hr$BMI <- ac.2hr$BMI - ave.bmi

# Convert SEX to factors
ac.2hr$Sex <- as.factor(ac.2hr$Sex)

# Plotting
############### EDSS ###################
# Linear models
for(i in 1:12){
  assign(paste("edss2hr.fit", i, sep = ""),lm(ac.2hr[,i+1]~EDSS + Age + Sex + BMI, data=ac.2hr))
}

edss2hr.fit.models <- c(list(edss2hr.fit1),list(edss2hr.fit2),list(edss2hr.fit3),list(edss2hr.fit4),list(edss2hr.fit5),list(edss2hr.fit6),
                        list(edss2hr.fit7),list(edss2hr.fit8),list(edss2hr.fit9),list(edss2hr.fit10),list(edss2hr.fit11),list(edss2hr.fit12))

# Summary linear regression
sink("edss_lm_2hr.txt")
for(i in 1:12){
  print(summary(edss2hr.fit.models[[i]]))
}
sink()

# 1) plot EDSS over evenly spaced 2hr interval
### Estimated EDSS values over 2hr intervals
edss.est_2hr <- data.frame(matrix(nrow=8,ncol=13))
names(edss.est_2hr) <- c("EDSS", "lm1","lm2","lm3","lm4","lm5","lm6",
                         "lm7","lm8","lm9","lm10","lm11","lm12")
edss.est_2hr$EDSS <- seq(0,7)
for(i in 1:12){
  edss.int <- coef(edss2hr.fit.models[[i]])['(Intercept)']
  edss.beta <- coef(edss2hr.fit.models[[i]])['EDSS']
  edss.est_2hr[,i+1] <- edss.int+edss.beta*edss.est_2hr$EDSS
}

colour <- c("Baseline")
for(i in 1:7){
  colour[i+1] <- paste("EDSS=",i,sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,12)

x0 <- c(rep(df.sub$time[1],8))
a <- seq(2,22,2)
for(i in 1:11){
  x0 <- append(x0, rep(df.sub$time[a[i]*6+1],8))
}

x1 <- c()
for(i in 1:11){
  x1 <- append(x1, rep(df.sub$time[a[i]*6+1],8))
}
x1 <- append(x1, rep(df.sub$time[24*6],8))

y0 <- c()
for(i in 1:12){
  y0[[i]] <- edss.est_2hr[,i+1]
}
y0 <- unlist(y0)

plot.edss.2hr <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y0,
  colour=colour
)

co <- colorRampPalette(brewer.pal(9,"Blues"))(9)[2:9]
co[1] <- "gray0"
plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #, colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("TLAC10 Linear Model Over 2hr Interval - EDSS")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.edss.2hr, aes(x=x, xend=xend, 
                                                  y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)

plt3 <- plt2+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[2*6+1],df.sub$time[4*6+1],df.sub$time[6*6+1],
                                       df.sub$time[8*6+1],df.sub$time[10*6+1],df.sub$time[12*6+1],
                                       df.sub$time[14*6+1],df.sub$time[16*6+1],df.sub$time[18*6+1],
                                       df.sub$time[20*6+1],df.sub$time[22*6+1],df.sub$time[24*6]),alpha=0.5)
png(file = "plot/Figure_7_2hrEDSS.png", width = 1000, height = 450)
plt3
dev.off()

# 2) plot EDSS over connected 2hr interval
EDSS2hr.pvalue <- data.frame(matrix(nrow=length(edss2hr.fit.models),ncol=2))
names(EDSS2hr.pvalue) <- c('Interval','p.value')
for(i in 1:length(edss2hr.fit.models)){
  EDSS2hr.pvalue[i,1] <- paste('lm',i,sep="")
  EDSS2hr.pvalue[i,2] <- summary(edss2hr.fit.models[[i]])$coef[2,4]
}
EDSS2hr.pvalue[EDSS2hr.pvalue$p.value <0.05,1]

int7 <- seq.POSIXt(df.sub$time[12*6+1], by='10 min', length.out = 12)
int8 <- seq.POSIXt(df.sub$time[14*6+1], by='10 min', length.out = 12)
int9 <- seq.POSIXt(df.sub$time[16*6+1], by='10 min', length.out = 12)
int10 <- seq.POSIXt(df.sub$time[18*6+1], by='10 min', length.out = 12)
int11 <- seq.POSIXt(df.sub$time[20*6+1], by='10 min', length.out = 12)

significant.range <- c(int7,int8,int9,int10,int11)


colour <- c("Baseline")
for(i in 1:7){
  colour[i+1] <- paste("EDSS=",i,sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,13)


x0 <- c(rep(df.sub$time[1],8))
b <- seq(1,23,2)
for(i in 1:12){
  x0 <- append(x0, rep(df.sub$time[b[i]*6+1],8))
}

x1 <- c()
for(i in 1:12){
  x1 <- append(x1, rep(df.sub$time[b[i]*6+1],8))
}
x1 <- append(x1, rep(df.sub$time[24*6],8))

y0<-c(list(edss.est_2hr$lm1))
for(i in 2:13){
  y0[[i]] <- edss.est_2hr[,i]
}
y0 <- unlist(y0)

y1<-c()
for(i in 1:12){
  y1[[i]] <- edss.est_2hr[,i+1]
}
y1 <- unlist(y1)
y1 <- append(y1,edss.est_2hr$lm12)

plot.edss.2hr_con <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y1,
  colour=colour
)

co <- colorRampPalette(brewer.pal(9,"Blues"))(9)[2:9]
co[1] <- "gray0"
plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1) + #, colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts") +  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("Connected Plot: TLAC10 Linear Model Over 2hr Interval - EDSS")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.edss.2hr_con, aes(x=x, xend=xend, 
                                                      y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)

plt3 <- plt2+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[2*6+1],df.sub$time[4*6+1],df.sub$time[6*6+1],
                                       df.sub$time[8*6+1],df.sub$time[10*6+1],df.sub$time[12*6+1],
                                       df.sub$time[14*6+1],df.sub$time[16*6+1],df.sub$time[18*6+1],
                                       df.sub$time[20*6+1],df.sub$time[22*6+1],df.sub$time[24*6]), alpha=0.5)
plt4 <- plt3+geom_vline(xintercept=significant.range, colour='slategray2',size=4, alpha=0.25)


png(file = "plot/Figure_8_2hrEDSS_con.png", width = 1000, height = 450)
plt4
dev.off()

############### T25fw ###################
# Linear models
for(i in 1:12){
  assign(paste("t25fw2hr.fit", i, sep = ""),lm(ac.2hr[,i+1]~T25fw + Age + Sex + BMI, data=ac.2hr))
}

t25fw2hr.fit.models <- c(list(t25fw2hr.fit1),list(t25fw2hr.fit2),list(t25fw2hr.fit3),
                         list(t25fw2hr.fit4),list(t25fw2hr.fit5),list(t25fw2hr.fit6),
                         list(t25fw2hr.fit7),list(t25fw2hr.fit8),list(t25fw2hr.fit9),
                         list(t25fw2hr.fit10),list(t25fw2hr.fit11),list(t25fw2hr.fit12))


# Summary linear regression
sink("t25fw_lm_2hr.txt")
for(i in 1:12){
  print(summary(t25fw2hr.fit.models[[i]]))
}
sink()



# 1) plot T25 over evenly spaced 2hr interval
# Estimated T25fw values over 2hr intervals
t25.est_2hr <- data.frame(matrix(nrow=6,ncol=13))
names(t25.est_2hr) <- c("T25", "lm1","lm2","lm3","lm4","lm5","lm6",
                         "lm7","lm8","lm9","lm10","lm11","lm12")
t25.est_2hr$T25 <- seq(0,25,5)
for(i in 1:12){
  int <- coef(t25fw2hr.fit.models[[i]])['(Intercept)']
  beta <- coef(t25fw2hr.fit.models[[i]])['T25fw']
  t25.est_2hr[,i+1] <- int+beta*t25.est_2hr$T25
}

colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("T25fw=",t25.est_2hr$T25[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,12)

x0 <- c(rep(df.sub$time[1],6))
a <- seq(2,22,2)
for(i in 1:11){
  x0 <- append(x0, rep(df.sub$time[a[i]*6+1],6))
}

x1 <- c()
for(i in 1:11){
  x1 <- append(x1, rep(df.sub$time[a[i]*6+1],6))
}
x1 <- append(x1, rep(df.sub$time[24*6],6))

y0 <- c()
for(i in 1:12){
  y0[[i]] <- t25.est_2hr[,i+1]
}
y0 <- unlist(y0)

plot.t25.2hr <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y0,
  colour=colour
)
co =c("gray0", "#9ECAE1" , "#6BAED6", "#4292C6", "#2171B5", "#08519C")
plt <-plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #, colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts") +  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("TLAC10 Linear Model Over 2hr Interval - T25fw")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data =plot.t25.2hr, aes(x=x, xend=xend, 
                                                  y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)

plt3 <- plt2+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[2*6+1],df.sub$time[4*6+1],df.sub$time[6*6+1],
                                       df.sub$time[8*6+1],df.sub$time[10*6+1],df.sub$time[12*6+1],
                                       df.sub$time[14*6+1],df.sub$time[16*6+1],df.sub$time[18*6+1],
                                       df.sub$time[20*6+1],df.sub$time[22*6+1],df.sub$time[24*6]),alpha=0.5)
png(file = "plot/Figure_7_2hrT25.png", width = 1000, height = 450)
plt3
dev.off()

# 2) plot T25fw over connected 2hr interval
T252hr.pvalue <- data.frame(matrix(nrow=length(t25fw2hr.fit.models),ncol=2))
names(T252hr.pvalue) <- c('Interval','p.value')
for(i in 1:length(t25fw2hr.fit.models)){
  T252hr.pvalue[i,1] <- paste('lm',i,sep="")
  T252hr.pvalue[i,2] <- summary(t25fw2hr.fit.models[[i]])$coef[2,4]
}
T252hr.pvalue[T252hr.pvalue$p.value <0.05,1]

int1 <- seq.POSIXt(df.sub$time[1], by='10 min', length.out = 12)
int2 <- seq.POSIXt(df.sub$time[2*6+1], by='10 min', length.out = 12)
int4 <- seq.POSIXt(df.sub$time[6*6+1], by='10 min', length.out = 12)
int5 <- seq.POSIXt(df.sub$time[8*6+1], by='10 min', length.out = 12)
int12 <- seq.POSIXt(df.sub$time[22*6+1], by='10 min', length.out = 12)

significant.range <- c(int1,int2,int4,int5,int12)

colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("T25fw=",t25.est_2hr$T25[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,13)

x0 <- c(rep(df.sub$time[1],6))
b <- seq(1,23,2)
for(i in 1:12){
  x0 <- append(x0, rep(df.sub$time[b[i]*6+1],6))
}

x1 <- c()
for(i in 1:12){
  x1 <- append(x1, rep(df.sub$time[b[i]*6+1],6))
}
x1 <- append(x1, rep(df.sub$time[24*6],6))

y0<-c(list(t25.est_2hr$lm1))
for(i in 2:13){
  y0[[i]] <- t25.est_2hr[,i]
}
y0 <- unlist(y0)

y1<-c()
for(i in 1:12){
  y1[[i]] <- t25.est_2hr[,i+1]
}
y1 <- unlist(y1)
y1 <- append(y1,t25.est_2hr$lm12)

plot.t25.2hr_con <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y1,
  colour=colour
)

co =c("gray0", "#9ECAE1" , "#6BAED6", "#4292C6", "#2171B5", "#08519C")
plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ # , colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts")+  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("Connected Plot: TLAC10 Linear Model Over 2hr Interval - T25fw")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.t25.2hr_con , aes(x=x, xend=xend, 
                                                      y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)

plt3 <- plt2+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[2*6+1],df.sub$time[4*6+1],df.sub$time[6*6+1],
                                       df.sub$time[8*6+1],df.sub$time[10*6+1],df.sub$time[12*6+1],
                                       df.sub$time[14*6+1],df.sub$time[16*6+1],df.sub$time[18*6+1],
                                       df.sub$time[20*6+1],df.sub$time[22*6+1],df.sub$time[24*6]),alpha=0.5)
plt4 <- plt3+geom_vline(xintercept=significant.range, colour='slategray2',size=4, alpha=0.25)

png(file = "plot/Figure_8_2hrT25_con.png", width = 1000, height = 450)
plt4
dev.off()

############### TUG ###################
# Linear models
for(i in 1:12){
  assign(paste("tug2hr.fit", i, sep = ""),lm(ac.2hr[,i+1]~TUG.t2 + Age + Sex + BMI, data=ac.2hr))
}

tug2hr.fit.models <- c(list(tug2hr.fit1),list(tug2hr.fit2),list(tug2hr.fit3),
                       list(tug2hr.fit4),list(tug2hr.fit5),list(tug2hr.fit6),
                       list(tug2hr.fit7),list(tug2hr.fit8),list(tug2hr.fit9),
                       list(tug2hr.fit10),list(tug2hr.fit11),list(tug2hr.fit12))

# Summary linear regression
sink("tug_lm_2hr.txt")
for(i in 1:12){
  print(summary(tug2hr.fit.models[[i]]))
}
sink()

# 1) plot T25 over evenly spaced 2hr interval
# Estimated T25fw values over 2hr intervals
tug.est_2hr <- data.frame(matrix(nrow=6,ncol=13))
names(tug.est_2hr) <- c("TUG", "lm1","lm2","lm3","lm4","lm5","lm6",
                        "lm7","lm8","lm9","lm10","lm11","lm12")
tug.est_2hr$TUG <- seq(0,25,5)
for(i in 1:12){
  int <- coef(tug2hr.fit.models[[i]])['(Intercept)']
  beta <- coef(tug2hr.fit.models[[i]])['TUG.t2']
  tug.est_2hr[,i+1] <- int+beta*tug.est_2hr$TUG
}

colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("TUG=",tug.est_2hr$TUG[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,12)

x0 <- c(rep(df.sub$time[1],6))
a <- seq(2,22,2)
for(i in 1:11){
  x0 <- append(x0, rep(df.sub$time[a[i]*6+1],6))
}

x1 <- c()
for(i in 1:11){
  x1 <- append(x1, rep(df.sub$time[a[i]*6+1],6))
}
x1 <- append(x1, rep(df.sub$time[24*6],6))

y0 <- c()
for(i in 1:12){
  y0[[i]] <- tug.est_2hr[,i+1]
}
y0 <- unlist(y0)

plot.tug.2hr <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y0,
  colour=colour
)
co =c("gray0", "#9ECAE1" , "#6BAED6", "#4292C6", "#2171B5", "#08519C")
plt <- plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #, colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts") +  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("TLAC10 Linear Model Over 2hr Interval - TUG")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data =plot.tug.2hr, aes(x=x, xend=xend, 
                                                y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)

plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)

plt3 <- plt2+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[2*6+1],df.sub$time[4*6+1],df.sub$time[6*6+1],
                                       df.sub$time[8*6+1],df.sub$time[10*6+1],df.sub$time[12*6+1],
                                       df.sub$time[14*6+1],df.sub$time[16*6+1],df.sub$time[18*6+1],
                                       df.sub$time[20*6+1],df.sub$time[22*6+1],df.sub$time[24*6]),alpha=0.5)
png(file = "plot/Figure_7_2hrTUG.png", width = 1000, height = 450)
plt3
dev.off()

# 2) plot TUG over connected 2hr interval
TUG2hr.pvalue <- data.frame(matrix(nrow=length(tug2hr.fit.models),ncol=2))
names(TUG2hr.pvalue) <- c('Interval','p.value')
for(i in 1:length(tug2hr.fit.models)){
  TUG2hr.pvalue[i,1] <- paste('lm',i,sep="")
  TUG2hr.pvalue[i,2] <- summary(tug2hr.fit.models[[i]])$coef[2,4]
}
TUG2hr.pvalue[TUG2hr.pvalue$p.value <0.05,1]

int1 <- seq.POSIXt(df.sub$time[1], by='10 min', length.out = 12)
int2 <- seq.POSIXt(df.sub$time[2*6+1], by='10 min', length.out = 12)
int4 <- seq.POSIXt(df.sub$time[6*6+1], by='10 min', length.out = 12)
int5 <- seq.POSIXt(df.sub$time[8*6+1], by='10 min', length.out = 12)
int12 <- seq.POSIXt(df.sub$time[22*6+1], by='10 min', length.out = 12)

significant.range <- c(int1,int2,int4,int5,int12)
colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("TUG=",tug.est_2hr$TUG[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,13)

x0 <- c(rep(df.sub$time[1],6))
b <- seq(1,23,2)
for(i in 1:12){
  x0 <- append(x0, rep(df.sub$time[b[i]*6+1],6))
}

x1 <- c()
for(i in 1:12){
  x1 <- append(x1, rep(df.sub$time[b[i]*6+1],6))
}
x1 <- append(x1, rep(df.sub$time[24*6],6))

y0<-c(list(tug.est_2hr$lm1))
for(i in 2:13){
  y0[[i]] <- tug.est_2hr[,i]
}
y0 <- unlist(y0)

y1<-c()
for(i in 1:12){
  y1[[i]] <- tug.est_2hr[,i+1]
}
y1 <- unlist(y1)
y1 <- append(y1,tug.est_2hr$lm12)

plot.tug.2hr_con <- data.frame(
  x = x0,
  xend = x1,
  y = y0,
  yend=y1,
  colour=colour
)

plt <-plt <- ggplot(meltdf,aes(x=time,y=Mean.log.ac,group=Subject.ID))+
  geom_point(alpha = 0.1)+ #, colour=meltdf$color)+
  xlab ("Time of the day") + ylab ("Total 10-min Log Activity Counts") +  
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
  ggtitle ("Connected Plot: TLAC10 Linear Model Over 2hr Interval - TUG")+
  theme (axis.text.x= element_text (angle=90),axis.text=element_text(size=12), plot.title = element_text(hjust = 0.5))

plt <- plt+geom_segment(data = plot.tug.2hr_con , aes(x=x, xend=xend, 
                                                      y=y,yend=yend,colour=colour),size=1.2,inherit.aes =FALSE)
co =c("gray0", "#9ECAE1" , "#6BAED6", "#4292C6", "#2171B5", "#08519C")
plt2 <- plt+theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)

plt3 <- plt2+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[2*6+1],df.sub$time[4*6+1],df.sub$time[6*6+1],
                                       df.sub$time[8*6+1],df.sub$time[10*6+1],df.sub$time[12*6+1],
                                       df.sub$time[14*6+1],df.sub$time[16*6+1],df.sub$time[18*6+1],
                                       df.sub$time[20*6+1],df.sub$time[22*6+1],df.sub$time[24*6]),alpha=0.5)
plt4 <- plt3+geom_vline(xintercept=significant.range, colour='slategray2',size=4, alpha=0.25)

png(file = "plot/Figure_8_2hrTUG_con.png", width = 1000, height = 450)
plt4
dev.off()

######## Figure 9. FOSR MODEL AND PLOTS ########
Y.model = as.matrix(t(df.sub[,1:62]))

####### FOSR: EDSS ###########
# create matrix for the covariates
X.model.edss = model.matrix ( ~ 1 + ac$EDSS +  ac$Age + ac$Sex + ac$BMI)
dim(X.model.edss)
tval = (1 : 144) / 144
fosr.edss = fosr (Y = Y.model, X = X.model.edss, argvals = tval, method = "OLS",
                   gam.method = "REML", nbasis = 12)
# estimated beta coefficients
SBeta = fosr.edss$est.func
VSBeta = fosr.edss$se.func

co <- colorRampPalette(brewer.pal(9,"Blues"))(9)[2:9]
co[1] <- "gray0"
hour<-seq.POSIXt(as.POSIXct("2010-01-01 00:00:00"),by="10 min", length.out = 144)
fosr.plt <- ggplot(as.data.frame(SBeta),aes(x=hour))+
  xlab("Time of the day") + ylab ("Total 10-min Log Activity Counts") +
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M") +ggtitle ("FOSR: EDSS")+
  theme(axis.text.x= element_text (angle=90),axis.text=element_text(size=12),plot.title = element_text(hjust = 0.5))
fosr.plt2 <- fosr.plt + geom_line(aes(x=hour,y=SBeta[,1], colour="Baseline"),size=1.2)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2], colour="EDSS=1"),size=1.2)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*2, colour="EDSS=2"),size=1.2)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*3, colour="EDSS=3"),size=1.2)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*4, colour="EDSS=4"),size=1.2)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*5, colour="EDSS=5"),size=1.2)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*6, colour="EDSS=6"),size=1.2)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*7, colour="EDSS=7"),size=1.2)
fosr.plt3 <- fosr.plt2+ 
  theme(legend.title=element_blank())+scale_colour_manual(name="", values=co)
#+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[4*6+1],df.sub$time[8*6+1],df.sub$time[12*6+1],
#                            df.sub$time[16*6+1],df.sub$time[20*6+1],df.sub$time[24*6]))
fosr.plt4 <- fosr.plt3+geom_point(data=meltdf,aes(x=time,y=Mean.log.ac, group=Subject.ID),
                                  alpha = 0.1,inherit.aes =FALSE)
significant.range = hour[which(((SBeta[, 2] - 2 * VSBeta[, 2]) * (SBeta[,2] + 2 * VSBeta[, 2]) > 0))]
fosr.plt5 <- fosr.plt4+geom_vline(xintercept=significant.range, colour='slategray2',size=4, alpha=0.25)

# color used to indicate significant regions

png(file = "plot/Figure_9_fosr_edss.png", width = 1000, height = 450)
fosr.plt5
dev.off()

####### FOSR: T25fw ###########
X.model.t25 = model.matrix ( ~ 1 + ac$T25fw +  ac$Age + ac$Sex + ac$BMI)
dim(X.model.t25)
tval = (1 : 144) / 144
fosr.t25 = fosr (Y = Y.model, X = X.model.t25, argvals = tval, method = "OLS",
                  gam.method = "REML", nbasis = 12)
# estimated beta coefficients
SBeta = fosr.t25$est.func
VSBeta = fosr.t25$se.func

hour<-seq.POSIXt(as.POSIXct("2010-01-01 00:00:00"),by="10 min", length.out = 144)
fosr.plt <- ggplot(as.data.frame(SBeta),aes(x=hour))+
  xlab("Time of the day") + ylab ("Total 10-min Log Activity Counts") +
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M") +ggtitle ("FOSR: T25fw")+
  theme(axis.text.x= element_text (angle=90),axis.text=element_text(size=12),plot.title = element_text(hjust = 0.5))

co =c("gray0" , "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#9ECAE1")
colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("T25fw=",y.est$T25fw[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,144)

fosr.plt2 <- fosr.plt + geom_line(aes(x=hour,y=SBeta[,1], colour="Baseline"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*5, colour="T25fw=5"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*10, colour="T25fw=10"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*15, colour="T25fw=15"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*20, colour="T25fw=20"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*25, colour="T25fw=25"),size=1.2,inherit.aes =FALSE)

fosr.plt3 <- fosr.plt2+ 
  theme(legend.title=element_blank())+scale_colour_manual(name="", breaks=colour, values=co)
#+geom_vline(xintercept = c(df.sub$time[1],df.sub$time[4*6+1],df.sub$time[8*6+1],df.sub$time[12*6+1],
#                            df.sub$time[16*6+1],df.sub$time[20*6+1],df.sub$time[24*6]))
fosr.plt4 <- fosr.plt3+geom_point(data=meltdf,aes(x=time,y=Mean.log.ac, group=Subject.ID),
                                 alpha = 0.1,inherit.aes =FALSE)
significant.range = hour[which(((SBeta[, 2] - 2 * VSBeta[, 2]) * (SBeta[,2] + 2 * VSBeta[, 2]) > 0))]
fosr.plt5 <- fosr.plt4+geom_vline(xintercept=significant.range, colour='slategray2',size=4, alpha=0.25)

png(file = "plot/Figure_9_fosr_t25.png", width = 1000, height = 450)
fosr.plt5 
dev.off()

####### FOSR: TUG ###########
X.model.tug = model.matrix ( ~ 1 + ac$TUG.t2+  ac$Age + ac$Sex + ac$BMI)
tval = (1 : 144) / 144
fosr.tug = fosr (Y = Y.model, X = X.model.tug, argvals = tval, method = "OLS",
                 gam.method = "REML", nbasis = 12)
# estimated beta coefficients
SBeta = fosr.tug$est.func
VSBeta = fosr.tug$se.func

hour<-seq.POSIXt(as.POSIXct("2010-01-01 00:00:00"),by="10 min", length.out = 144)
fosr.plt <- ggplot(as.data.frame(SBeta),aes(x=hour))+
  xlab("Time of the day") + ylab ("Total 10-min Log Activity Counts") +
  scale_x_datetime(date_breaks = "2 hour", date_labels ="%H:%M") +ggtitle ("FOSR: TUG")+
  theme(axis.text.x= element_text (angle=90),axis.text=element_text(size=12),plot.title = element_text(hjust = 0.5))

co =c("gray0" , "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#9ECAE1")
colour <- c("Baseline")
for(i in 1:5){
  colour[i+1] <- paste("TUG=",y.est$TUG[i+1],sep="")
}
colour <- factor(colour ,levels = colour)
colour <- rep(colour,144)

fosr.plt2 <- fosr.plt + geom_line(aes(x=hour,y=SBeta[,1], colour="Baseline"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*5, colour="TUG=5"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*10, colour="TUG=10"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*15, colour="TUG=15"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*20, colour="TUG=20"),size=1.2,inherit.aes =FALSE)+
  geom_line(aes(x=hour,y=SBeta[,1]+SBeta[,2]*25, colour="TUG=25"),size=1.2,inherit.aes =FALSE)

fosr.plt3 <- fosr.plt2+ 
  theme(legend.title=element_blank())+scale_colour_manual(name="", breaks=colour, values=co)
#+  geom_vline(xintercept = c(df.sub$time[1],df.sub$time[4*6+1],df.sub$time[8*6+1],df.sub$time[12*6+1],
#                            df.sub$time[16*6+1],df.sub$time[20*6+1],df.sub$time[24*6]))
fosr.plt4 <- fosr.plt3+geom_point(data=meltdf,aes(x=time,y=Mean.log.ac, group=Subject.ID), #colour=meltdf$color, 
                                 alpha = 0.1, inherit.aes =FALSE)
significant.range = hour[which(((SBeta[, 2] - 2 * VSBeta[, 2]) * (SBeta[,2] + 2 * VSBeta[, 2]) > 0))]
fosr.plt5 <- fosr.plt4+geom_vline(xintercept=significant.range, colour='slategray2',size=4, alpha=0.25)

png(file = "plot/Figure_9_fosr_tug.png", width = 1000, height = 450)
fosr.plt5
dev.off()