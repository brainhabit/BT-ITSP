#==========
# Script for multivariate twin modelling for the study: Individual variability in sensory processing: etiological structure and developmental links to autism.
# 
# of note: here it runs on ITSP sensory quadrants, but separate runs can be used to model ITSP section and factor scores.
#
# G. Bussu, October 2022
#==========

rm(list=ls())

setwd('C:/Users/your_username/Documents/BT-ITSP/data')
list.files()

### load OpenMx and functions:
require(OpenMx)
source('C:/Users/your_username/Documents/twin_modelling/tutorial/miFunctions.R')

mxOption(key="Number of Threads",
         value=parallel::detectCores())

### prepare some data:
data <- read.csv(file='twin_data_ITSP.csv',header=T,sep=',')
names(data);dim(data)

Vars <- c('resSensSens','resSensSeek','resSensAvoid','resLowReg')
nv <- 4 # number of phenotypes
ntv <- nv*2 # number of measures/variables
(selVars <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep=''))


mz <- subset(data,zygosity=='MZ',c(selVars))
dz <- subset(data,zygosity=='DZ',c(selVars))

############## Fully Saturated model ##################################################################################

svCov <- c(1,rep(.5,7),1,rep(.5,6),1,rep(.5,5),1,rep(.5,4),1,rep(.5,3),1,rep(.5,2),1,.5,1)


# Create Matrices for Covariates and linear Regression Coefficients
expMeanMZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.01,
                      labels=labFull('m_mz',1,ntv),name='ExpMeanMZ')
expMeanDZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.01,
                      labels=labFull('m_dz',1,ntv),name='ExpMeanDZ')


expCovMZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=svCov,
labels=labSymm('mz_cov',ntv),name='ExpCovMZ')

expCovDZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=svCov,
labels=labSymm('dz_cov',ntv),name='ExpCovDZ')

matI <- mxMatrix(type='Iden',nrow=ntv,ncol=ntv,name='I')

expCorMZ <- mxAlgebra(solve(sqrt(I*ExpCovMZ))%&%ExpCovMZ,name='ExpCorMZ')
expCorDZ <- mxAlgebra(solve(sqrt(I*ExpCovDZ))%&%ExpCovDZ,name='ExpCorDZ')

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:
objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMeanMZ',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMeanDZ',dimnames=selVars)

funcML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups

# specify the submodels for each zygosity group:
modelMZ <- mxModel('MZ',expMeanMZ,expCovMZ,matI,expCorMZ,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',expMeanDZ,expCovDZ,matI,expCorDZ,dataDZ,objDZ,funcML)

# combine submodels into one object to make sure they are evaluated together:
multi <- mxFitFunctionMultigroup(c('MZ','DZ'))
ci <- mxCI(c('MZ.ExpCorMZ','DZ.ExpCorDZ'))

# combine all model objects:
SatModel <- mxModel('Sat',modelMZ,modelDZ,multi,ci)

mxOption(model= SatModel, key="Number of Threads", value= (omxDetectCores() - 1))

SatFit <- mxTryHard(SatModel,intervals=T)
sumSatFit<-summary(SatFit)

####################################### Constrained Model (estimate twin correlation) ############################################################

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
CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_2_1','mz_cov_6_5',
'dz_cov_2_1','dz_cov_6_5'),
free=T,values=.5,newlabel='rph_SensSeek')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_3_1','mz_cov_7_5',
'dz_cov_3_1','dz_cov_7_5'),
free=T,values=.5,newlabel='rph_SensSA')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_4_1','mz_cov_8_5',
'dz_cov_4_1','dz_cov_8_5'),
free=T,values=.5,newlabel='rph_SensLR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_3_2','mz_cov_7_6',
                                               'dz_cov_3_2','dz_cov_7_6'),
                             free=T,values=.5,newlabel='rph_SeekSA')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_4_2','mz_cov_8_6',
                                               'dz_cov_4_2','dz_cov_8_6'),
                             free=T,values=.5,newlabel='rph_SeekLR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_4_3','mz_cov_8_7',
                                               'dz_cov_4_3','dz_cov_8_7'),
                             free=T,values=.5,newlabel='rph_SALR')


# fix the CTCT to be the same between twins:
CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_6_1','mz_cov_5_2'),
free=T,values=.5,newlabel='ctct_mz_SensSeek')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_6_1','dz_cov_5_2'),
free=T,values=.5,newlabel='ctct_dz_SensSeek')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_5_3','mz_cov_7_1'),
free=T,values=.5,newlabel='ctct_mz_SensSA')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_5_3','dz_cov_7_1'),
free=T,values=.5,newlabel='ctct_dz_SensSA')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_5_4','mz_cov_8_1'),
                             free=T,values=.5,newlabel='ctct_mz_SensLR')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_5_4','dz_cov_8_1'),
                             free=T,values=.5,newlabel='ctct_dz_SensLR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_6_3','mz_cov_7_2'),
                             free=T,values=.5,newlabel='ctct_mz_SeekSA')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_6_3','dz_cov_7_2'),
                             free=T,values=.5,newlabel='ctct_dz_SeekSA')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_6_4','mz_cov_8_2'),
free=T,values=.5,newlabel='ctct_mz_SeekLR')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_6_4','dz_cov_8_2'),
free=T,values=.5,newlabel='ctct_dz_SeekLR')

CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_7_4','mz_cov_8_3'),
                             free=T,values=.5,newlabel='ctct_mz_SALR')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_7_4','dz_cov_8_3'),
                             free=T,values=.5,newlabel='ctct_dz_SALR')


# fit the model
CorFit <- mxTryHard(CorModel,intervals=T)
sumCorFit<-summary(CorFit)

# compare with fully saturated to check assumptions
mxCompare(SatFit,CorFit)

########################################## correlated factors solution ####################################################

### full model:


# path coefficients:
coefpath<-c(1,rep(.5,3),1,rep(.5,2),1,.5,1)

pathA <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=coefpath,
                  labels=labLower('a',nv),name='a')
pathC <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=coefpath,
                  labels=labLower('c',nv),name='c')  
pathE <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=coefpath,
                  labels=labLower('e',nv),name='e')  

# calculate the variance components:

varA <- mxAlgebra(a%*%t(a),name='A')
varC <- mxAlgebra(c%*%t(c),name='C')
varE <- mxAlgebra(e%*%t(e),name='E')

# calculate the total variance and the standard deviation:

varP <- mxAlgebra(A+C+E,name='V')

matI <- mxMatrix(type='Iden',nrow=nv,ncol=nv,name='I')
isd <- mxAlgebra(solve(sqrt(I*V)),name='iSD')

# calculate phenotypic and etiological correlations:

corP <- mxAlgebra(solve(sqrt(I*V))%&%V,name='rPH')
corA <- mxAlgebra(solve(sqrt(I*A))%&%A,name='rA')
corC <- mxAlgebra(solve(sqrt(I*C))%&%C,name='rC')
corE <- mxAlgebra(solve(sqrt(I*E))%&%E,name='rE')

# means:

expMean <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.001,
                  labels=labFull('m',1,nv),name='ExpMean')


# variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

# convert the variance components to proportions:

estVC <- mxAlgebra(cbind(A/V,C/V,E/V),name='EstVC')

# observed data:

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',dimnames=selVars)

# method of estimation:

funcML <- mxFitFunctionML()

# specify submodels for each zygosity group:

pars <- list(pathA,pathC,pathE,varA,varC,varE,varP,matI,isd,corP,corA,corC,corE,expMean)

modelMZ <- mxModel('MZ',pars,estVC,expCovMZ,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',pars,estVC,expCovDZ,dataDZ,objDZ,funcML)

# specify that we need to evaluate all submodels together:

multi <- mxFitFunctionMultigroup(c('MZ','DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.EstVC','MZ.rPH','MZ.rA','MZ.rC','MZ.rE'))

# combine all objects:

CorFac<- mxModel('Cor',modelMZ,modelDZ,multi,ci)

# fit the model:

FacFit <- mxTryHard(CorFac,intervals=T)
sumCorFac<-summary(FacFit)

# compare fit to saturated model:

mxCompare(SatFit,FacFit)


######################################## Independent pathways #######################################################################################################################################33

nf        <- 1       # number of factors

# Create Matrices for Covariates and linear Regression Coefficients
expMean <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.001,
                  labels=labFull('m',1,nv),name='ExpMean')


# Matrices ac, cc, and ec to store a, d, and e path coefficients for common factors
pathAc    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE, values=.5, labels=labFull("ac",nv,nf), name="ac" )
pathCc    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE, values=.5, labels=labFull("cc",nv,nf), name="cc" )
pathEc    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE, values=.5, labels=labFull("ec",nv,nf), name="ec" )

# Matrices as, ds, and es to store a, d, and e path coefficients for specific factors
pathAs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("as",nv), name="as" )
pathCs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("cs",nv), name="cs" )
pathEs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("es",nv), name="es" )

# Matrices A, D, and E compute variance components
varA      <- mxAlgebra( expression=ac %*% t(ac) + as %*% t(as), name="A" )
varC      <- mxAlgebra( expression=cc %*% t(cc) + cs %*% t(cs), name="C" )
varE      <- mxAlgebra( expression=ec %*% t(ec) + es %*% t(es), name="E" )

varP      <- mxAlgebra( expression= A+C+E, name="V" )


### expected variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

# raw data
dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:
objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',dimnames=selVars)

funcML <- mxFitFunctionML()

# calculate proportions of variance explained by common and specific 
# variance components:

estVC <- mxAlgebra(rbind(cbind(A/V,C/V,E/V),
                         cbind((ac%*%t(ac))/V,(cc%*%t(cc))/V,(ec%*%t(ec))/V),
                         cbind((as%*%t(as))/V,(cs%*%t(cs))/V,(es%*%t(es))/V)),
                   name='EstVC')

# Create Model Objects for Multiple Groups
pars      <- list(pathAc, pathCc, pathEc, pathAs, pathCs, pathEs, varA, varC, varE, varP, expMean)


modelMZ   <- mxModel( name="MZ", pars, expCovMZ, dataMZ, objMZ, funcML, estVC )
modelDZ   <- mxModel( name="DZ", pars, expCovDZ, dataDZ, objDZ, funcML, estVC )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

### specify the confidence intervals:

ci <- mxCI(c('MZ.EstVC'))


# Build & Run Model 
modelIP   <- mxModel( "mulIPc", modelMZ, modelDZ, multi, ci )
fitIP     <- mxTryHard( modelIP, intervals=T)
sumIP     <- summary( fitIP )

mxCompare( FacFit, fitIP )



########################################### Common Pathway ACE Model ######################################################################################################################
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

nl        <- 1

# Create Matrices for Covariates and linear Regression Coefficients
expMean <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.001,
                  labels=labFull('m',1,nv),name='ExpMean')

# Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
pathAl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.5, labels=labLower("al",nl), name="al" )
pathCl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.5, labels=labLower("cl",nl), name="cl" )
pathEl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.5, labels=labLower("el",nl), name="el" )


# Matrix and Algebra for constraint on variance of latent phenotype
varLP   <- mxAlgebra( expression=al %*% t(al) + cl %*% t(cl) + el %*% t(el), name="VarLP" )

unit      <- mxMatrix( type="Unit", nrow=nl, ncol=1, name="Unit")
varLP1    <- mxConstraint(diag2vec(VarLP)==Unit, name="varLP1")

# Matrix f for factor loadings on latent phenotype
pathFl    <- mxMatrix( type="Full", nrow=nv, ncol=nl, free=TRUE, values=1, labels=labFull("fl",nv,nl), name="fl" )

# Matrices as, cs, and es to store a, d, and e path coefficients for specific factors
pathAs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("as",nv), name="as" )
pathCs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("cs",nv), name="cs" )
pathEs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, labels=labDiag("es",nv), name="es" )


# Matrices A, C, and E compute variance components
covA      <- mxAlgebra( expression=fl %&% (al %*% t(al)) + as %*% t(as), name="A" )
covC      <- mxAlgebra( expression=fl %&% (cl %*% t(cl)) + cs %*% t(cs), name="C" )
covE      <- mxAlgebra( expression=fl %&% (el %*% t(el)) + es %*% t(es), name="E" )

covP      <- mxAlgebra( expression= A+C+E, name="V" )

### expected variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

# # Create Algebra for Standardization
# matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
# invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")
# 
# # Calculate genetic and environmental correlations
# corA      <- mxAlgebra( expression=solve(sqrt(I*A))%&%A, name ="rA" ) #cov2cor()
# corC      <- mxAlgebra( expression=solve(sqrt(I*C))%&%C, name ="rC" )
# corE      <- mxAlgebra( expression=solve(sqrt(I*E))%&%E, name ="rE" )

# raw data
dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:
objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',dimnames=selVars)

funcML <- mxFitFunctionML()

# calculate proportions:
estVC <- mxAlgebra(cbind(A/V,C/V,E/V),name='EstVC')

estVCl <- mxAlgebra(cbind((al%*%t(al))/VarLP,(cl%*%t(cl))/VarLP,(el%*%t(el))/VarLP),
                    name='EstVCl')
estVCc <- mxAlgebra(cbind((fl%&%(al%*%t(al)))/V,(fl%&%(cl%*%t(cl)))/V,
                          (fl%&%(el%*%t(el)))/V),name='EstVCc')
estVCs <- mxAlgebra(cbind((as%*%t(as))/V,(cs%*%t(cs))/V,(es%*%t(es))/V),name='EstVCs')


# Create Model Objects for Multiple Groups
pars      <- list(expMean, pathAl, pathCl, pathEl, varLP, unit, pathFl, pathAs, pathCs, pathEs, covA, covC, covE, covP)

modelMZ   <- mxModel( name="MZ", pars, expCovMZ, dataMZ, objMZ, funcML,estVC,estVCl,estVCc,estVCs )
modelDZ   <- mxModel( name="DZ", pars, expCovDZ, dataDZ, objDZ, funcML,estVC,estVCl,estVCc,estVCs )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

### specify the confidence intervals:

ci <- mxCI(c('MZ.EstVC','MZ.EstVCl','MZ.EstVCc','MZ.EstVCs'))


# Build & Run Model 
modelCP   <- mxModel( "mulCPc", pars, varLP1, modelMZ, modelDZ, multi, ci )
#modelCP   <- mxModel( "mulCPc", modelMZ, modelDZ, multi, ci )
fitCP    <- mxTryHard(modelCP, intervals=T)

sumCP     <- summary( fitCP )
mxCompare( FacFit,c(fitIP, fitCP ))

mxCompare( SatFit,c(FacFit, fitIP, fitCP ))

save.image("~/multivariate_quadrants.RData")

########################################### Plots ######################################################################################################################
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Plot variance
IPvariance<-matrix(round(sumIP$CI[,2],digits=2),nrow=12)

Aunique<-IPvariance[rbind(c(9,1),c(10,2),c(11,3),c(12,4))]
Cunique<-IPvariance[rbind(c(9,5),c(10,6),c(11,7),c(12,8))]
Eunique<-IPvariance[rbind(c(9,9),c(10,10),c(11,11),c(12,12))]

Ashared<-IPvariance[rbind(c(5,1),c(6,2),c(7,3),c(8,4))]
Cshared<-IPvariance[rbind(c(5,5),c(6,6),c(7,7),c(8,8))]
Eshared<-IPvariance[rbind(c(5,9),c(6,10),c(7,11),c(8,12))]

Atot<-IPvariance[rbind(c(1,1),c(2,2),c(3,3),c(4,4))]
Ctot<-IPvariance[rbind(c(1,5),c(2,6),c(3,7),c(4,8))]
Etot<-IPvariance[rbind(c(1,9),c(2,10),c(3,11),c(4,12))]

IPvariance_low<-matrix(round(sumIP$CI[,1],digits=2),nrow=12)

Aunique_low<-IPvariance_low[rbind(c(9,1),c(10,2),c(11,3),c(12,4))]
Cunique_low<-IPvariance_low[rbind(c(9,5),c(10,6),c(11,7),c(12,8))]
Eunique_low<-IPvariance_low[rbind(c(9,9),c(10,10),c(11,11),c(12,12))]

Ashared_low<-IPvariance_low[rbind(c(5,1),c(6,2),c(7,3),c(8,4))]
Cshared_low<-IPvariance_low[rbind(c(5,5),c(6,6),c(7,7),c(8,8))]
Eshared_low<-IPvariance_low[rbind(c(5,9),c(6,10),c(7,11),c(8,12))]

Atot_low<-IPvariance_low[rbind(c(1,1),c(2,2),c(3,3),c(4,4))]
Ctot_low<-IPvariance_low[rbind(c(1,5),c(2,6),c(3,7),c(4,8))]
Etot_low<-IPvariance_low[rbind(c(1,9),c(2,10),c(3,11),c(4,12))]

IPvariance_high<-matrix(round(sumIP$CI[,3],digits=2),nrow=12)

Aunique_high<-IPvariance_high[rbind(c(9,1),c(10,2),c(11,3),c(12,4))]
Cunique_high<-IPvariance_high[rbind(c(9,5),c(10,6),c(11,7),c(12,8))]
Eunique_high<-IPvariance_high[rbind(c(9,9),c(10,10),c(11,11),c(12,12))]

Ashared_high<-IPvariance_high[rbind(c(5,1),c(6,2),c(7,3),c(8,4))]
Cshared_high<-IPvariance_high[rbind(c(5,5),c(6,6),c(7,7),c(8,8))]
Eshared_high<-IPvariance_high[rbind(c(5,9),c(6,10),c(7,11),c(8,12))]

Atot_high<-IPvariance_high[rbind(c(1,1),c(2,2),c(3,3),c(4,4))]
Ctot_high<-IPvariance_high[rbind(c(1,5),c(2,6),c(3,7),c(4,8))]
Etot_high<-IPvariance_high[rbind(c(1,9),c(2,10),c(3,11),c(4,12))]

vartype<-c(rep('Au',4),rep('Cu',4),rep('Eu',4),rep('Ac',4),rep('Cc',4),rep('Ec',4))
pheno<-rep(c('Sensory Sensitivity','Sensory Seeking','Sensory Avoidance','Low Registration'),6)
vartype_general<-c(rep('Unique variance',12),rep('Common variance',12))
plotdata<-data.frame(cbind(c(Aunique_low,Cunique_low,Eunique_low,Ashared_low,Cshared_low,Eshared_low),c(Aunique,Cunique,Eunique,Ashared,Cshared,Eshared),c(Aunique_high,Cunique_high,Eunique_high,Ashared_high,Cshared_high,Eshared_high),pheno,vartype,vartype_general))
plotdata$pheno<-ordered(pheno,level=c('Sensory Sensitivity','Sensory Seeking','Sensory Avoidance','Low Registration'))
plotdata$vartype<-as.factor(vartype)
plotdata$vartype_general<-as.factor(vartype_general)
names(plotdata)<-c('CI_low','Variance','CI_high','Phenotype','Variance_Components','Variance_Type')
plotdata$CI_low<-as.numeric(plotdata$CI_low)
plotdata$CI_high<-as.numeric(plotdata$CI_high)
plotdata$Variance<-as.numeric(plotdata$Variance)

library(ggplot2)
library(viridis)

png('variance.png',width=2000,height=2000,res=300)
ggplot(plotdata, aes(fill=Type, y=varlist, x=pheno)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T, option = "E") +
  facet_wrap(~vartype_general) +
  theme_bw() +
  theme(legend.position="right") +
  xlab("MSEL scale score") + ylab('Variance (%)')
dev.off()

# correlation barplot
rph<-sumCorFit$CI[c(2,3,4,11,12,20),2]
rph_low<-sumCorFit$CI[c(2,3,4,11,12,20),1]
rph_high<-sumCorFit$CI[c(2,3,4,11,12,20),3]
ctct_mz<-sumCorFit$CI[c(6,21,29,22,30,31),2]
ctct_mz_low<-sumCorFit$CI[c(6,21,29,22,30,31),1]
ctct_mz_high<-sumCorFit$CI[c(6,21,29,22,30,31),3]
ctct_dz<-sumCorFit$CI[c(70,85,93,86,94,95),2]
ctct_dz_low<-sumCorFit$CI[c(70,85,93,86,94,95),1]
ctct_dz_high<-sumCorFit$CI[c(70,85,93,86,94,95),3]
pheno<-rep(c('SensSens-Seek','SensSens-SA','SensSens-LR','Seek-SA','Seek-LR','SA-LR'),3)
 
corrtype<-c(rep('Rph',6),rep('CTCT-MZ',6),rep('CTCT-DZ',6))
corrplot<-data.frame(cbind(c(rph,ctct_mz,ctct_dz),c(rph_low,ctct_mz_low,ctct_dz_low),c(rph_high,ctct_mz_high,ctct_dz_high),pheno,corrtype))
names(corrplot)<-c('Correlation','Lower bound','Upper bound','Phenotype','Type')
corrplot$Phenotype<-ordered(corrplot$Phenotype,levels=c('SensSens-Seek','SensSens-SA','SensSens-LR','Seek-SA','Seek-LR','SA-LR'))
corrplot$Type<-ordered(corrplot$Type,levels=c('Rph','CTCT-MZ','CTCT-DZ'))
corrplot$Correlation<-as.numeric(corrplot$Correlation)
corrplot$`Lower bound`<-as.numeric(corrplot$`Lower bound`)
corrplot$`Upper bound`<-as.numeric(corrplot$`Upper bound`)

png('CTCT_quadrants.png',width=6000,height=6000,res=300)
ggplot(corrplot, aes(fill=Type, y=Correlation, x=Phenotype)) +
  geom_bar(position="dodge", stat="identity") + geom_errorbar(position=position_dodge(width=0.8),aes(ymin = `Lower bound`, ymax = `Upper bound`), width =.2, col = "black")+
  scale_fill_brewer(palette = "Pastel2", name = "Type",labels=c('Phenotypic','CTCT_MZ','CTCT-DZ')) +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(text=element_text(size=26),legend.position="right")+
  xlab("ITSP quadrant") + ylab('Correlation')
dev.off()

rDZ<-sumCorFit$CI[c(69,78,87,96),2]
rMZ<-sumCorFit$CI[c(5,14,23,32),2]
rDZlow<-sumCorFit$CI[c(69,78,87,96),1]
rDZhigh<-sumCorFit$CI[c(69,78,87,96),3]
rMZlow<-sumCorFit$CI[c(5,14,23,32),1]
rMZhigh<-sumCorFit$CI[c(5,14,23,32),3]
type<-c(rep('MZ',4),rep('DZ',4))
type<-ordered(type,levels=c('MZ','DZ'))
pheno<-rep(c('Sensory Sensitivity','Sensory Seeking','Sensory Avoidance','Low Registration'),2)
pheno<-ordered(pheno,levels=c('Sensory Sensitivity','Sensory Seeking','Sensory Avoidance','Low Registration'))
twinplot<-data.frame(cbind(c(rMZ,rDZ),c(rMZlow,rDZlow),c(rMZhigh,rDZhigh),pheno,type))
names(twinplot)<-c('Correlation','Low','High','Phenotype','Type')
twinplot$Correlation<-as.numeric(twinplot$Correlation)
twinplot$Low<-as.numeric(twinplot$Low)
twinplot$High<-as.numeric(twinplot$High)
twinplot$Type<-as.factor(twinplot$Type)
twinplot$Phenotype<-as.factor(twinplot$Phenotype)
levels(twinplot$Phenotype)<-c('Low Registration','Sensation Seeking','Sensory Sensitivity','Sensation Avoiding')
levels(twinplot$Type)<-c('MZ','DZ')


png('twincorr_quadrants.png',width=5000,height=5000,res=300)
ggplot(twinplot, aes(fill=Type, y=Correlation, x=Phenotype)) +
geom_bar(position="dodge", stat="identity") + geom_errorbar(position=position_dodge(width=0.8),aes(ymin = Low, ymax = High), width =.2, col = "black")+
  scale_fill_brewer(palette = "Set3", name = "Zygosity") +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(text=element_text(size=26),legend.position="right")+
xlab("ITSP quadrant score") + ylab('Correlation')
dev.off()

comptype<-c(rep('A',4),rep('C',4),rep('E',4))
pheno<-rep(c('Low Registration','Sensation Seeking','Sensory Sensitivity','Sensation Avoiding'),3)
corrplot<-data.frame(cbind(c(.29,.30,.44,.38,.62,.65,.52,.56,.09,.06,.04,.06),pheno,comptype))
names(corrplot)<-c('Value','Phenotype','Etiology')
corrplot$Phenotype<-ordered(corrplot$Phenotype,levels=c('Low Registration','Sensation Seeking','Sensory Sensitivity','Sensation Avoiding'))
corrplot$Etiology<-ordered(corrplot$Etiology,levels=c('A','C','E'))
corrplot$Value<-as.numeric(corrplot$Value)

png('unitot_variance_quadrant.png',width=5000,height=5000,res=300)
ggplot(corrplot, aes(fill=Etiology, y=Value, x=Phenotype)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Pastel2", name = "Etiology",labels=c('A','C','E')) +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(text=element_text(size=26),legend.position="right")+
  xlab("ITSP quadrant score") + ylab('Variance')
dev.off()
