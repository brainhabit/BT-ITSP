#==========
# Script for bivariate twin modelling for the study: Individual variability in sensory processing: etiological structure and developmental links to autism.
# 
# G. Bussu, October 2022
#==========

rm(list=ls())


setwd('C:/Users/your_username/Documents/BT-ITSP/data')
list.files()

require(OpenMx)
source('C:/Users/your_username/Documents/twin_modelling/tutorial/miFunctions.R')

##### prepare data

### import data

data <- read.csv(file='twin_data_ITSP_all_factor_wQCHAT.csv',header=T,sep=',')
names(data);dim(data)

# # eventually change
# data$qchat1<-scale(data$qchat1)
# data$qchat2<-scale(data$qchat2)

### select varianbles for analysis:

Vars <- c('f3','qchat')
nv <- 2 # number of phenotypes
ntv <- nv*2 # number of measures/variables
(selVars <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep=''))

### select subsets of data for analysis:

mz <- subset(data,zygosity=='MZ',c(selVars))
dz <- subset(data,zygosity=='DZ',c(selVars))

##### saturated model #####

### specify and fit a fully saturated model

# starting values for the covariances:

svCov <- c(1,rep(.5,3),1,rep(.5,2),1,.5,1)

# means:

expMeanMZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=.001,
                      labels=labFull('m_mz',1,ntv),name='ExpMeanMZ')
expMeanDZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=.001,
                      labels=labFull('m_dz',1,ntv),name='ExpMeanDZ')

# variances and covariances:

expCovMZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=svCov,
                     labels=labSymm('mz_cov',ntv),name='ExpCovMZ')
expCovDZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=svCov,
                     labels=labSymm('dz_cov',ntv),name='ExpCovDZ')

# convert the covariances to correlations:

matI <- mxMatrix(type='Iden',nrow=ntv,ncol=ntv,name='I')

expCorMZ <- mxAlgebra(solve(sqrt(I*ExpCovMZ))%&%ExpCovMZ,name='ExpCorMZ')  
expCorDZ <- mxAlgebra(solve(sqrt(I*ExpCovDZ))%&%ExpCovDZ,name='ExpCorDZ')

#expCorMZ <- mxAlgebra(solve(sqrt(I*ExpCovMZ))%*%ExpCovMZ%*%solve(sqrt(I*ExpCovMZ)),name='ExpCorMZ')  
#expCorDZ <- mxAlgebra(solve(sqrt(I*ExpCovDZ))%*%ExpCovDZ%*%solve(sqrt(I*ExpCovDZ)),name='ExpCorDZ')

  
# observed data:

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')  

# objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMeanMZ',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMeanDZ',dimnames=selVars)

# estimation method:

funcML <- mxFitFunctionML()
  
# specify the submodels for each zygosity group:

modelMZ <- mxModel('MZ',expMeanMZ,expCovMZ,matI,expCorMZ,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',expMeanDZ,expCovDZ,matI,expCorDZ,dataDZ,objDZ,funcML)

# combine submodels into one object to make sure they are evaluated together:

multi <- mxFitFunctionMultigroup(c('MZ','DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.ExpCorMZ','DZ.ExpCorDZ'))
  
# combine all model objects:

SatModel <- mxModel('Sat',modelMZ,modelDZ,multi,ci)

# fit the model:

SatFit <- mxTryHard(SatModel,intervals=T)
summary(SatFit)

### add constraints to estimate the twin correlations:

CorModel <- mxModel(SatFit,name='Cor')

# equal means across twins and zygosity:

CorModel <- omxSetParameters(CorModel,labels=c(labFull('m_mz',1,ntv),
                                               labFull('m_dz',1,ntv)),
                             free=T,values=.001,newlabels=labFull('m',1,nv))

# equal variances:

CorModel <- omxSetParameters(CorModel,labels=c(labDiag('mz_cov',ntv),
                                               labDiag('dz_cov',ntv)),
                             free=T,values=1,newlabels=labDiag('v',nv))

# equal phenotypic correlation in twin 1 and twin 2:

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_2_1','mz_cov_4_3',
                                              'dz_cov_2_1','dz_cov_4_3'),
                            free=T,values=.5,newlabel='rph')

# fix the CTCT to be the same between twins:

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_3_2','mz_cov_4_1'),
                            free=T,values=.5,newlabel='ctct_mz')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_3_2','dz_cov_4_1'),
                             free=T,values=.5,newlabel='ctct_dz')

# fit the model:

CorFit <- mxTryHard(CorModel,intervals=T)
summary(CorFit)

##### Cholesky decomposition #####

### full model:

# path coefficients:

pathA <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=c(1,.5,1),
                  labels=labLower('a',nv),name='a')
pathC <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=c(1,.5,1),
                  labels=labLower('c',nv),name='c')  
pathE <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=c(1,.5,1),
                  labels=labLower('e',nv),name='e')  

# calculate the variance components:

varA <- mxAlgebra(a%*%t(a),name='A')
varC <- mxAlgebra(c%*%t(c),name='C')
varE <- mxAlgebra(e%*%t(e),name='E')

# calculate the total variance and the standard deviation:

varP <- mxAlgebra(A+C+E,name='V')

matI <- mxMatrix(type='Iden',nrow=nv,ncol=nv,name='I')
isd <- mxAlgebra(solve(sqrt(I*V)),name='iSD')

# standardize the variance components and then square them:

cholA <- mxAlgebra((iSD%*%a)^2,name='CholA')
cholC <- mxAlgebra((iSD%*%c)^2,name='CholC')
cholE <- mxAlgebra((iSD%*%e)^2,name='CholE')

# calculate the phenotypic correlation:

corP <- mxAlgebra(solve(sqrt(I*V))%&%V,name='rPH')

# means:

expMean <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=.001,
                    labels=labFull('m',1,nv),name='ExpMean')

# variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

# observed data:

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',dimnames=selVars)

# method of estimation:

funcML <- mxFitFunctionML()

# specify submodels for each zygosity group:

pars <- list(pathA,pathC,pathE,varA,varC,varE,varP,matI,isd,cholA,cholC,cholE,corP,expMean)

modelMZ <- mxModel('MZ',pars,expCovMZ,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',pars,expCovDZ,dataDZ,objDZ,funcML)

# specify that we need to evaluate all submodels together:

multi <- mxFitFunctionMultigroup(c('MZ','DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.CholA','MZ.CholC','MZ.CholE'))

# combine all objects:

Cholesky <- mxModel('Chol',modelMZ,modelDZ,multi,ci)

# fit the model:

CholFit <- mxTryHard(Cholesky,intervals=T)
summary(CholFit)

# compare fit to saturated model:

mxCompare(SatFit,CholFit)

### nested models

# AE model

CholAE <- mxModel(CholFit,name='AE')
CholAE <- omxSetParameters(CholAE,labels=labLower('c',nv),free=F,values=0)
CholAEFit <- mxRun(CholAE,intervals=T)
summary(CholAEFit)

# CE model

CholCE <- mxModel(CholFit,name='CE')
CholCE <- omxSetParameters(CholCE,labels=labLower('a',nv),free=F,values=0)
CholCEFit <- mxRun(CholCE,intervals=T)
summary(CholCEFit)

# E model

CholE <- mxModel(CholFit,name='E')
CholE <- omxSetParameters(CholE,labels=c(labLower('a',nv),
                                         labLower('c',nv)),free=F,values=0)
CholEFit <- mxRun(CholE,intervals=T)
summary(CholEFit)

# compare models

mxCompare(CholFit,c(CholAEFit,CholCEFit,CholEFit))

save.image("~/BT-ITSP/cholesky_F3_QCHAT.RData")













































































