###
###
### This code needs to be revised
###
###

library(survey)
library(Hmisc)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(stringr)
library(car)
library(data.table)

sampname = 'SampleData'
S=1000 #Number of samples drawn
B=500 #Number of bootstrap runs within each sample
idx=1:1000
seed = 20190414+idx[1]+3
cat('The seed is: ',seed)
set.seed(seed)

files = paste0('./RData/',sampname,"_",str_pad(1:S,nchar(S),pad=0),".RData")
m1.list.raking = list()
m2.list.raking = list()

calib = function(w1,w2,z){
    # This function implements equation 1 from https://www.stata.com/merror/call.pdf (page 83) or https://www.stata.com/merror/rcal.pdf
    # Best linear approximation = exact conditional expectation under joint normality
    N = nrow(z)
    
    wbar = rowMeans(cbind(w1,w2))
    muz = matrix(colMeans(z),nrow = N,ncol = ncol(z),byrow = T)
    nu = N*2 - (N*4)/(N*2)
    
    suu = sum((w1-wbar)^2+(w2-wbar)^2)/(N*(2-1))
    sxx = (sum(2*((wbar-mean(wbar))^2))-(N-1)*suu)/nu
    sxz = colSums(2*(wbar-mean(wbar))*(z-muz))/nu
    szz = cov(z)
    
    return(as.numeric(c(sxx,sxz)%*%solve(rbind(c(sxx+suu/2,sxz),cbind(sxz,szz)))%*%rbind(wbar-mean(wbar),t((z-muz)))))
}

for(i in idx){
    load(file=files[i])
    # Creating sampling probabilities
    samp <- as.data.table(samp) 
    samp[,c('P.Bg','P.hh','P.sub') := .(1/W.bg,1/W.hh,1/W.sub)]
    samp.solnas <- samp[(solnas==T),]

    # Step 1: fitting calibration model
    # The code below fits the calibration equation using one of the two available biomarker replicates
    lm.lsodi = glm(c_ln_na_bio1 ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,data=samp.solnas)
    samp[,c_ln_na_calib := predict(lm.lsodi,newdata=samp,'response')]
    
    # Step 2: Fit DDM using calib, get influence functions, and add to data.
    # In Lumley's book, the following models are fitted without sampling weights and design as a simple GLM
    ddm1 = glm(hypertension ~ c_age + c_bmi + c_ln_na_calib + high_chol + usborn + female + bkg_pr + bkg_o,data=samp,family=quasibinomial())
    ddm2 = glm(sbp ~ c_age + c_bmi + c_ln_na_calib + high_chol + usborn + female + bkg_pr + bkg_o,data=samp,family=gaussian())
    
    inf.ddm1 = dfbeta(ddm1);colnames(inf.ddm1)[1] = 'Int';colnames(inf.ddm1) = paste0('Inf.ddm1.',colnames(inf.ddm1));samp = cbind(samp,inf.ddm1)
    inf.ddm2 = dfbeta(ddm2);colnames(inf.ddm2)[1] = 'Int';colnames(inf.ddm2) = paste0('Inf.ddm2.',colnames(inf.ddm2));samp = cbind(samp,inf.ddm2)
    
    # Step 3: Create twophase design and calibrate with influence functions
    twophase.design = twophase(id=list(~BGid+hhid+subid,~1),subset=~solnas,strata=list(~strat,NULL),probs=list(~P.Bg+P.hh+P.sub,NULL),data=samp,method='full')
    twophase.design.ddm1 = calibrate(twophase.design,phase=2,calfun='raking',as.formula(paste0('~',paste(colnames(inf.ddm1),collapse='+'))))
    twophase.design.ddm2 = calibrate(twophase.design,phase=2,calfun='raking',as.formula(paste0('~',paste(colnames(inf.ddm2),collapse='+'))))
    
    # Step 4: Fit DDM models
    m1.raking = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_bio1 + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm1,family=quasibinomial())
    m2.raking = svyglm(sbp ~ c_age + c_bmi + c_ln_na_bio1 + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm2,family=gaussian())
    
    samp.design = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=samp)
    
    m1.true = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_true + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
    m2.true = svyglm(sbp ~ c_age + c_bmi + c_ln_na_true + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
    
    m1.avg = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
    m2.avg = svyglm(sbp ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
    
    # Saving the results
    m1.list.raking[[i]] = data.frame(Sim=i,Coeff=names(m1.true$coefficients),
                                     True.Est=as.numeric(m1.true$coefficients),True.SE = as.numeric(SE(m1.true)),
                                     Avg.Est=as.numeric(m1.avg$coefficients),Avg.SE = as.numeric(SE(m1.avg)),
                                     Raking.Est=as.numeric(m1.raking$coefficients),Raking.SE = as.numeric(SE(m1.raking)))
    
    m2.list.raking[[i]] = data.frame(Sim=i,Coeff=names(m2.true$coefficients),
                                     True.Est=as.numeric(m2.true$coefficients),True.SE = as.numeric(SE(m2.true)),
                                     Avg.Est=as.numeric(m2.avg$coefficients),Avg.SE = as.numeric(SE(m2.avg)),
                                     Raking.Est=as.numeric(m2.raking$coefficients),Raking.SE = as.numeric(SE(m2.raking)))
}
cat('Done!',"\n")

if(!dir.exists('./Output')){system('mkdir Output')}
save(m1.list.raking,m2.list.raking,
     file=paste0('./Output/raking_BetaVersion.RData'))