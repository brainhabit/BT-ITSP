#==========
# Script for data preparation for the study: Individual variability in sensory processing: etiological structure and developmental links to autism.
# 
# G. Bussu, October 2022
#==========

rm(list=ls())

### set a working directory (otherwise create a new R project):

setwd('C:/Users/your_username/Documents/BT-ITSP')
list.files() # show files in the working directory

### import and check data:

data_itsp <- read.csv('C:/Users/your_username/Documents/BT-ITSP/data/ITSP_items_clean.csv')
names(data_itsp);dim(data_itsp)

data_demo <- readxl::read_xlsx('C:/Users/your_username/Documents/BT-MSEL/data_batss_msel_vabs/Demographics_5m.xlsx')
names(data_demo);dim(data_demo)

data_exclusion<-readxl::read_xlsx('C:/Users/your_username/Documents/BT-MSEL/data_batss_msel_vabs/BabyTwins background summary and exclusion 20201102.xlsx')
names(data_exclusion); dim(data_exclusion)

data_background<- readxl::read_xlsx('C:/Users/your_username/Documents/BT-MSEL/data_batss_msel_vabs/Background_5m.xlsx')

data_pregnancy<- readxl::read_xlsx('C:/Users/your_username/Documents/BT-MSEL/data_batss_msel_vabs/pregnancy_time_data.xlsx')

### matching participants across files
indx<-match(data_itsp$id,data_demo$`Participant EASE Code`)
data_demo<-data_demo[indx,]

indx_ex<-match(data_itsp$id,data_exclusion$Code)
data_exclusion<-data_exclusion[indx_ex,]

length(which(data_itsp$id!=data_demo$`Participant EASE Code`| data_itsp$id!=data_exclusion$Code | data_demo$`Participant EASE Code`!=data_exclusion$Code))

# add dataset background info

backmatch<-match(data_itsp$id,data_background$Kod)
data_background<-data_background[backmatch,]

length(which(data_itsp$id!=data_background$Kod))

# add data on pregnancy term
indx_preg<-match(data_itsp$id,data_pregnancy$ID)

data_pregnancy<-data_pregnancy[indx_preg,]

### exclusion based on general criteria

excl_indx<-which(data_exclusion$`Exclusion version A`==1)

# clean data files from excluded participants
data_demo_clean<-data_demo[-excl_indx,]
data_itsp_clean<-data_itsp[-excl_indx,]
data_background_clean<-data_background[-excl_indx,]
data_pregnancy_clean<-data_pregnancy[-excl_indx,]

# check pregnancy
preg_day<-data_pregnancy_clean$Gdays
preg_day[which(is.na(preg_day))]<-0

pregnancy_term<-data_pregnancy_clean$Gweeks*7+preg_day

### dataset included in our analyses

data<-data.frame(cbind(data_itsp_clean[-c(38)],data_demo_clean$Gender,data_demo_clean$`Age at Date of Assessment (days)`,pregnancy_term,data_demo_clean$`Bio Mum Age`,data_demo_clean$`Bio Dad Age`,data_demo_clean$`A. Highest level of education`,data_demo_clean$`B. Highest level of education`,data_background_clean$`F16. Ungefär hur hög är familjens* gemensamma månadsinkomst innan skatt (lön och andra ersättningar/bidrag)?`,data_demo_clean$TWAB,data_demo_clean$`Twin pair no.`))
names(data)<-c(names(data_itsp_clean[-38]),'sex','age','term_age','Mum_age','Dad_age','A_edu','B_edu','family_income','TWAB','Twinpair')

# parental education level to MAX within family
data$A_edu_num<-as.numeric(as.factor(data$A_edu))
data$A_edu_num[which(data$A_edu_num==4)]<-6
data$A_edu_num[which(data$A_edu_num==5)]<-4
data$A_edu_num[which(data$A_edu_num==6)]<-5

data$B_edu_num<-as.numeric(as.factor(data$B_edu))
data$B_edu_num[which(data$B_edu_num==4)]<-6
data$B_edu_num[which(data$B_edu_num==5)]<-4
data$B_edu_num[which(data$B_edu_num==6)]<-5

#data$edu_max<-pmax(data$A_edu_num,data$B_edu_num)
data$edu_mean<-rowMeans(cbind(data$A_edu_num,data$B_edu_num),na.rm = T)


# parental age to mean across mum and dad (NA.RM=TRUE for 3 missing Dad age)
data$Mum_age<-as.numeric(data$Mum_age)
data$Dad_age<-as.numeric(data$Dad_age)
data$parental_age<-rowMeans(cbind(data$Mum_age,data$Dad_age),na.rm = T)

## add zygosity
data_exclusion_clean<-data_exclusion[-excl_indx,]
data$zygosity<-data_exclusion_clean$Zygosity

names(data);dim(data)

##### basic data tidying #####

data<-data[,-c(2:37)]
data<-data[,-c(17:20,24,25)]

table(data$zygosity)

## table demo and descriptives ##

library(tidyverse)
library(finalfit)

data$sex_factor<-as.factor(data$sex)
data$family_income[which(data$family_income=='vet ej')]<-'unknown'
data$income_factor<-as.factor(data$family_income)
data$zyg_factor<-as.factor(data$zygosity)
explanatory <- c('age','sex_factor','term_age','edu_factor','income_factor','parental_age')
data$edu_factor<-as.factor(round(data$edu_mean))
data <- data %>%
mutate(
income_factor = ff_label(income_factor, "Family income")
)
data <- data %>%
mutate(
edu_factor = ff_label(edu_factor, "Mean parental education")
)
data <- data %>%
mutate(
term_age = ff_label(term_age, "Gestation age (days)")
)
data <- data %>%
mutate(
age = ff_label(age, "Age (days)")
)
data <- data %>%
mutate(
sex_factor = ff_label(sex_factor, "Sex")
)
table1 <- data %>%
summary_factorlist('zyg_factor', explanatory,
p=TRUE, na_include=TRUE,
add_dependent_label=TRUE,
dependent_label_prefix = "Zygosity: "
)
write.table(table1, file = "table1.txt", sep = ",", quote = FALSE, row.names = F)
explanatory <- c('Auditory','Visual','Tactile','Vestibular','General','LowRegistration','SensationSeeking','SensationAvoiding','SensorySensitivity')
table2 <- data %>%
summary_factorlist('zyg_factor', explanatory,
p=TRUE, na_include=TRUE,
add_dependent_label=TRUE,
dependent_label_prefix = "Zygosity: "
)
write.table(table2, file = "table2.txt", sep = ",", quote = FALSE, row.names = F)
##########################

## load QCHAT
data_quest <- read.csv('imputed_data_BT.csv')
data$qchat<-data_quest$qchat_total_score_36m[match(data$id,data_quest$id)]

########################

# check how many pairs are incomplete by making a variable that counts 
#   the frequency of their pair ID:

data <- merge(data,data.frame(table(Twinpair=data$Twinpair)),by='Twinpair')
table(data$Freq) # check how many pair IDs appear only once 

# carry out the exclusion:

data$Twinpair[which(data$Freq==1)]
data_test<-rbind(data[1:89,],data[89:94,],data[94:179,],data[179:220,],data[220:281,],data[281:448,],data[448:453,],data[453:506,],data[506:519,],data[519:534,],data[534:length(data$id),])
data_test <- merge(data_test,data.frame(table(id=data_test$id)),by='id')
naindx<-matrix(which(data_test$Freq.y==2),nrow=2)[1,]

data_test$LowRegistration[naindx]<-NA
data_test$SensationSeeking[naindx]<-NA
data_test$SensorySensitivity[naindx]<-NA
data_test$SensationAvoiding[naindx]<-NA
data_test$Hyposensitivity[naindx]<-NA
data_test$Hypersensitivity[naindx]<-NA
data_test$SumScore[naindx]<-NA
data_test$General[naindx]<-NA
data_test$Auditory[naindx]<-NA
data_test$Visual[naindx]<-NA
data_test$Tactile[naindx]<-NA
data_test$Vestibular[naindx]<-NA
data_test$f1[naindx]<-NA
data_test$f2[naindx]<-NA
data_test$f3[naindx]<-NA
data_test$f4[naindx]<-NA
#data_test$f5[naindx]<-NA
data_test$TWAB[naindx[c(1,2,4,7,9,10)]]<-2
data_test$TWAB[naindx[c(3,5,6,8)]]<-1

data<-data_test

dim(data) 

# remove the frequency variable:

data <- data[-c(24,25)]

##### transform data for twin analysis #####

# binary sex: 0=Females; 1=Males.
data$sex[which(data$sex=='Male')]<-1
data$sex[which(data$sex=='Female')]<-0
data$sex<-as.numeric(data$sex)

# ordinal discrete income
data$income<-as.numeric(as.factor(data$family_income))
data$income[which(data$income==2)]<-12
data$income[which(data$income==1)]<-2
data$income[which(data$income==11)]<-1
data$income[which(data$income==12)]<-11

data<-data[,-c(18)]

# numeric variables
vars <- colnames(data)[c(2:20,22:27)] # before was 22, 23
data[vars]<-lapply(data[vars],as.numeric)

### check the distribution of the ITSP scales:

vars <- colnames(data)[c(22:26)]
lapply(data[vars],psych::describe)

# transform scales
data$SensationSeeking_t<-I(data$SensationSeeking)^3
data$SensationAvoiding_t<-log(data$SensationAvoiding)
data$SumScore_t<-log(data$SumScore)
data$Hypersensitivity_t<-log(data$Hypersensitivity)
data$SensorySensitivity_t<-log(data$SensorySensitivity)
data$Hyposensitivity_t<-log(data$Hyposensitivity)
data$LowRegistration_t<-sqrt(data$LowRegistration)
data$Tactile_t<-sqrt(data$Tactile)
data$General_t<-sqrt(data$General)

## 4 factors solution (March 2023)
data$f1_t<-log(data$f1)
data$f4_t<-log(data$f4)

# scale variables
vars <- colnames(data)[c(16,17,19,20,24,27:36)]
data[vars]<-lapply(data[vars],scale)

### control for the effects of covariates: sex, age, parental age, income, education. Factor example, run for quadrant and section scores too.

library(drgee)

summary(gee( formula = f3~income,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f2~income,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f4_t~income,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f1_t~income,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f3~edu_mean,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f2~edu_mean,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f4_t~edu_mean,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f1_t~edu_mean,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f3~sex,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f2~sex,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f4_t~sex,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f1_t~sex,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f3~parental_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f2~parental_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f4_t~parental_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f1_t~parental_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f3~term_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f2~term_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f4_t~term_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f1_t~term_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f3~age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f2~age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f4_t~age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = f1_t~age,data = data, clusterid = "Twinpair", cond = F))

# fit regression models that adjust for sex and then save the residuals:

data$resSensSeek <- resid(lm(SensationSeeking_t~parental_age+edu_mean,data=data,na.action=na.exclude))
data$resSensAvoid <- resid(lm(SensationAvoiding_t~parental_age+edu_mean,data=data,na.action=na.exclude))
data$resLowReg <- resid(lm(LowRegistration_t~parental_age+edu_mean,data=data,na.action=na.exclude))
data$resSensSens <- resid(lm(SensorySensitivity_t~parental_age+edu_mean,data=data,na.action=na.exclude))
data$resSum <- resid(lm(SumScore_t~parental_age+edu_mean,data=data,na.action=na.exclude))
data$resQCHAT <- resid(lm(qchat~parental_age+edu_mean,data=data,na.action=na.exclude))


data$resHypo <- resid(lm(Hyposensitivity_t~parental_age+edu_mean,data=data,na.action=na.exclude))
data$resHyper <- resid(lm(Hypersensitivity_t~parental_age+edu_mean,data=data,na.action=na.exclude))

# add income to domain score residuals since income linked to general domain, parental age and education to general, visual and vestibular, while tactile and auditory showed no associations to covariates
data$resGeneral <- resid(lm(General_t~parental_age+edu_mean+income,data=data,na.action=na.exclude))
data$resVisual <- resid(lm(Visual~parental_age+edu_mean,data=data,na.action=na.exclude))
data$resAuditory <- resid(lm(Auditory~parental_age+edu_mean,data=data,na.action=na.exclude))
data$resTactile <- resid(lm(Tactile_t~parental_age+edu_mean,data=data,na.action=na.exclude))
data$resVestibular <- resid(lm(Vestibular~parental_age+edu_mean,data=data,na.action=na.exclude))

# edu mean F2 est=-.12, err=.06, p=.03035 not boferroni corrected, therefore not corrected
#data$resF1 <-resid(lm(f1_t~edu_mean,data=data,na.action=na.exclude))
#data$resF2 <-resid(lm(f2_t~edu_mean,data=data,na.action=na.exclude))
#data$resF3 <-resid(lm(f3_t~edu_mean,data=data,na.action=na.exclude))
#data$resF4 <-resid(lm(f4_t~edu_mean,data=data,na.action=na.exclude))
#data$resF5 <-resid(lm(f5_t~sex+parental_age+edu_mean,data=data,na.action=na.exclude))

# check that the standardization worked:

stand <- colnames(data)[c(22:27,28:29)]
lapply(data[stand],psych::describe)

data[stand]<-lapply(data[stand],scale)

########################### Option 1: split twin1 data and twin2 data based on TWAB ####################################
#twin1 <- subset(data,TWAB==1)
#twin2 <- subset(data,TWAB==2)
#
#twin1labs <- paste(colnames(twin1),1,sep='')
#twin2labs <- paste(colnames(twin2),2,sep='')
#
#names(twin1) <- twin1labs
#names(twin2) <- twin2labs
#
############################ Option 2: try another randomization than TWAB ################################################

### give each person a random number (only use two possible numbers, and never 0 and 1):

set.seed(2022)
npairs <- nrow(data)/2
rand <- c(rep(4,npairs),rep(5,npairs))
data$rand <- sample(rand) # randomly choose 1 number for each person

### divide the twins into two subsets: twin 1 and twin 2

twin1 <- subset(data,TWAB==1)
twin2 <- subset(data,TWAB==2)

### now we need to flip twin 2's random number so that it is the opposite of twin 1:

# basically we are going to add twin 1's random number to the twin 2 dataframe and 
#  then recode it so it's the opposite of what it currently is. We have to start by
#   making a dataframe containing just the pair ID number and twin 1's random number:

randVar <- data.frame(cbind(twin1$Twinpair,twin1$rand))
colnames(randVar) <- c('Twinpair','rand')

# now take away twin 2's random number and replace it with the twin 1 random number:

twin2 <- twin2[-c(30)]
twin2 <- merge(twin2,randVar,by='Twinpair')

# now we check that the frequency of each number is the same in both files (i.e.
#   twin 1 and twin 2 should both have the same random number):

table(twin1$rand);table(twin2$rand)

# for twin 1, convert the random numbers to 0s and 1s:

twin1$rand[twin1$rand==4] <- 0
twin1$rand[twin1$rand==5] <- 1

# and then recode twin 2, but in the opposite direction to twin 1:

twin2$rand[twin2$rand==4] <- 1
twin2$rand[twin2$rand==5] <- 0

# check that it worked:

table(twin1$rand);table(twin2$rand) # twin 1 should have as many 0s as twin 2 has 1s

# now join twin 1 and twin 2 back together:

data <- data.frame(rbind(twin1,twin2))
names(data);dim(data)

# sort the data by pairnr and check the first few rows:

data <- data[order(data[,2]),]
head(data)

# split the data by random number:

twin1 <- subset(data,rand==1)
twin2 <- subset(data,rand==0)

# check that the dimensions are the same:

dim(twin1);dim(twin2)

# add '1' to the end of the variables for twin 1 and '2' for twin 2:

twin1labs <- paste(colnames(twin1),1,sep='')
twin2labs <- paste(colnames(twin2),2,sep='')

names(twin1) <- twin1labs
names(twin2) <- twin2labs

##########################################################################################################################

# combine data so that 1 pair each line
dataD <- data.frame(cbind(twin1,twin2))

# remove unused variables:

dataD <- dataD[-c(25:26,50:52,65:77,100)]

# relabel a few variables:

names(dataD)[names(dataD)=='Twinpair1'] <- 'Twinpair'
names(dataD)[names(dataD)=='id1'] <- 'id'
names(dataD)[names(dataD)=='zygosity1'] <- 'zygosity'
names(dataD)[names(dataD)=='sex1'] <- 'sex'
names(dataD)[names(dataD)=='age1'] <- 'age'
names(dataD)[names(dataD)=='term_age1'] <- 'term_age'
names(dataD)[names(dataD)=='parental_age1'] <- 'parental_age'
names(dataD)[names(dataD)=='income1'] <- 'income'
names(dataD)[names(dataD)=='TWAB1'] <- 'TWAB'
names(dataD)[names(dataD)=='edu_mean1'] <- 'edu_mean'

dataD <- dataD[-c(3:14,24:32,37,39:67,69:78,85:92)]
#dataD <- dataD[-c(14:18,23)]

#### done. Save two versions of the data: 

write.csv(data,'C:/Users/your_username/Documents/BT-ITSP/data/phenotypic_data.csv',row.names = F)
write.csv(dataD,'C:/Users/your_username/Documents/BT-ITSP/data/twin_data.csv',row.names = F)
