#=========================================================================================================================================
#
# Main for factor analysis based on the ITSP items at 5 months in BATSS
#
# run for the study: Individual variability in sensory processing: etiological structure and developmental links to autism.
# 
# 
# G. Bussu, January 2023
#
#=========================================================================================================================================

rm(list=ls())

library(RcmdrMisc)
library(psych)

#################################################################################################################################################
## prep data
#################################################################################################################################################


setwd('C:/Users/your_username/Documents/BT-ITSP')
list.files() # show files in the working directory

### import and check data:

data_itsp <- read.csv('C:/Users/your_username/Documents/BT-ITSP/data/ITSP_items_clean.csv')
data_itsp<-data_itsp[-1]
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

# select one twin in a pair for EFA and set aside other twin for later CFA
data_CFA <- data[-which(data$TWAB==2),]
data<-data[-which(data$TWAB==1),]

#################################################################################################################################################
## check data
#################################################################################################################################################

# Kaiser-Meyer-Olkin factor adequacy (>0.49 unacceptable, >0.6 mediocre, >0.8 meritorious, >0.9 marvelous)
mydata<-data[,2:37]

kmo_itsp<-KMO(mydata) ## NB items 1, 17, 18, 21, 23, 31 below 0.6 mediocre
excl_item<-which(kmo_itsp$MSAi<0.6)

round( kmo_itsp$MSA, 2 )

mydata<-mydata[-excl_item]
round(KMO(mydata)$MSA,2)

# Barlett's test for sphericity (should be able to reject H0=no correlation among variables)
cortest.bartlett(mydata)

#################################################################################################################################################
## determine number of factors
#################################################################################################################################################

eigenvalues<-eigen(cor(mydata))
eigenvalues$values

scree(mydata,pc=F) # change in the curve direction after f=4
fa.parallel(mydata,fa='fa') # suggests f=10 but definitely too many for meaningful, interpretable factors. Trust scree more.

#################################################################################################################################################
## extract and rotate factors
#################################################################################################################################################

nfacs<-4

library(GPArotation)

fit <- factanal(mydata, nfacs, rotation="oblimin")
print(fit, digits=2, cutoff=0.3, sort=TRUE)

## NB factor 4 has 2 items loading up on that only. A minimum of 3 items is necessary for identifiability through CFA.

fa.none <- fa(r=mydata,nfactors = 4,fm='pa', max.iter=100, rotate='oblimin')
print(fa.none)

fa.diagram(fa.none)


#################################################################################################################################################
## export factor loadings & visualize results
#################################################################################################################################################

loads <- fit$loadings
fa.diagram(loads)

load <- fit$loadings[,1:2]
plot(load,type="n") # set up plot
text(load,labels=names(mydata),cex=.7)


FactorLoadings <- round(fit$loadings[1:28, ], 3)

write.csv(FactorLoadings, file="FacLoads_28items_oblimin.csv")

#################################################################################################################################################
## factor evaluation (consistency)
#################################################################################################################################################

f1<-mydata[,c('itsp7','itsp9','itsp12','itsp15','itsp20','itsp22')]
f2<-mydata[,c('itsp24','itsp25','itsp26','itsp27','itsp28','itsp35')]
f3<-mydata[,c('itsp5','itsp6','itsp8','itsp13','itsp14')]
f4<-mydata[,c('itsp32','itsp34')]

alpha(f1,check.keys=TRUE)
alpha(f2,check.keys=TRUE)
alpha(f3,check.keys=TRUE)
alpha(f4,check.keys=TRUE)

data$f1<-rowSums(f1)
data$f2<-rowSums(f2)
data$f3<-rowSums(f3)
data$f4<-rowSums(f4)

## check correlation with quadrant scores
corr.test(data[,c(38:41,65:68)])

write.csv(data,'C:/Users/giobu365/Documents/BT-ITSP/data/data_ITSP_factor_long.csv',row.names = F)

################## CFA
library(lavaan)

data<-data_CFA

mydata<-data[,2:37]

f1<-mydata[,c('itsp7','itsp9','itsp12','itsp15','itsp20','itsp22')]
f2<-mydata[,c('itsp24','itsp25','itsp26','itsp27','itsp28','itsp35')]
f3<-mydata[,c('itsp5','itsp6','itsp8','itsp13','itsp14')]
f4<-mydata[,c('itsp32','itsp34')]

data$f1<-rowSums(f1)
data$f2<-rowSums(f2)
data$f3<-rowSums(f3)
data$f4<-rowSums(f4)

#correlated two factor solution, marker method
mfac <- 'f1 =~ itsp3 + itsp7 + itsp8 + itsp9 + itsp10 + itsp12 + itsp15 + itsp20 + itsp29
        f2 =~ itsp4 + itsp5 + itsp6 + itsp13 + itsp14 + itsp19 + itsp25 + itsp27 + itsp33 + itsp34 + itsp36
        f3 =~ itsp11 + itsp16 + itsp32
        f4 =~ itsp22 + itsp24 + itsp26 + itsp28 + itsp35
' 
mcfa <- cfa(mfac, data=mydata[which(data$TWAB==1),],std.lv=TRUE) 
summary(mcfa,fit.measures=TRUE,standardized=TRUE)
