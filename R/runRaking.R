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
  
  m1.list.raking[[i]] = list()
  m2.list.raking[[i]] = list()
  
  # Step 0: Adjusting biomarker nutrient using the attenuation factor (R.J. Carroll & D. Ruppert, 2002) 
  # See https://www.stata.com/merror/call.pdf (beginning on page 82) and https://www.stata.com/merror/rcal.pdf
  # Because Biomarker = True + error, we can regress Biomarker ~ Naive + Adjustments to estimate the parameters in the regression of True ~ Naive + Adjustments
  
  samp.solnas <- samp[(solnas==T),]
  samp.solnas[,c_ln_na_bio_moments_adj := calib(w1 = samp.solnas$c_ln_na_bio1,w2 = samp.solnas$c_ln_na_bio2,z = as.matrix(samp.solnas[,c('c_age','c_bmi','high_chol','usborn','female','bkg_pr','bkg_o')]))]
  samp.solnas[,c_ln_na_bio_glm_adj := predict(glm(c_ln_na_bio2 ~ c_age + c_bmi + c_ln_na_bio1 + high_chol + usborn + female + bkg_pr + bkg_o,data = samp.solnas),newdata=samp.solnas,'response')]
  samp <- merge(samp,samp.solnas[,c('subid','c_ln_na_bio_moments_adj','c_ln_na_bio_glm_adj')],by = 'subid',all.x = T)
  
  # Step 1: fitting calibration model
  # The code below fits the calibration equation using one of the two available biomarker replicates
  # I have tried fitting this model using both c_ln_na_bio_moments_adj  and c_ln_na_bio_glm_adj, but the results did not look good
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
  
  # Step 4: Fit DDM without and with calibration
  m1.raking.moments = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_bio_moments_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm1,family=quasibinomial())
  m2.raking.moments = svyglm(sbp ~ c_age + c_bmi + c_ln_na_bio_moments_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm2,family=gaussian())
  
  m1.raking.glm = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_bio_glm_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm1,family=quasibinomial())
  m2.raking.glm = svyglm(sbp ~ c_age + c_bmi + c_ln_na_bio_glm_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm2,family=gaussian())
  
  # Extra: fitting 'true' (unobserved), 'naive' (2-day avg.) model, and 'biomarker' models for comparison purposes.
  # For the biomarker models, we fit the model with the two-phase design without raking.
  # At some point, we wanted to see if the estimates are at least unbiased in the subset
  
  samp.design = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=samp)
  
  m1.true = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_true + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
  m2.true = svyglm(sbp ~ c_age + c_bmi + c_ln_na_true + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
  
  m1.avg = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
  m2.avg = svyglm(sbp ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
  
  m1.bio.moments = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_bio_moments_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design,family=quasibinomial())
  m2.bio.moments = svyglm(sbp ~ c_age + c_bmi + c_ln_na_bio_moments_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design,family=gaussian())
  
  m1.bio.glm = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_bio_glm_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design,family=quasibinomial())
  m2.bio.glm = svyglm(sbp ~ c_age + c_bmi + c_ln_na_bio_glm_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design,family=gaussian())
  
  # Bootstrap Raking starts now
  # Not sure if vairance estimates using the raking approach needs to be corrected by bootstrapping, in any case I am doing it below
  samp[,paste0(colnames(inf.ddm1))] = NULL
  samp[,paste0(colnames(inf.ddm2))] = NULL
  samp$c_ln_na_bio_moments_adj = NULL
  samp$c_ln_na_bio_glm_adj = NULL
  samp$c_ln_na_calib = NULL
  
  idx.solnas = which(samp$solnas==T)
  idx.notsolnas = which(samp$solnas==F)
  
  for(j in 1:B){cat('Sample: ',i,". Boot: ",j,".\n")
    samp.boot <- samp[c(idx.notsolnas,sample(x = idx.solnas,size = length(idx.solnas),replace = T)),]
    samp.boot[,subid := str_replace(make.names(samp.boot$subid,unique=T),"X","")] # Renaming the duplicated IDs
    
    # Adjusting the biomarker
    samp.boot.solnas <- samp.boot[(solnas==T),]
    samp.boot.solnas[,c_ln_na_bio_moments_adj := calib(w1 = samp.boot.solnas$c_ln_na_bio1,w2 = samp.boot.solnas$c_ln_na_bio2,z = as.matrix(samp.boot.solnas[,c('c_age','c_bmi','high_chol','usborn','female','bkg_pr','bkg_o')]))]
    samp.boot.solnas[,c_ln_na_bio_glm_adj := predict(glm(c_ln_na_bio2 ~ c_age + c_bmi + c_ln_na_bio1 + high_chol + usborn + female + bkg_pr + bkg_o,data = samp.boot.solnas),newdata=samp.boot.solnas,'response')]
    samp.boot <- merge(samp.boot,samp.boot.solnas[,c('subid','c_ln_na_bio_moments_adj','c_ln_na_bio_glm_adj')],by = 'subid',all.x = T)
    
    # Now, fitting calibration model
    lm.lsodi = glm(c_ln_na_bio1 ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,data=samp.boot.solnas)
    samp.boot[,c_ln_na_calib := predict(lm.lsodi,newdata=samp.boot,'response')]
    
    # Step 2: Fit DDM using calib, get influence functions, and add to data.
    ddm1 = glm(hypertension ~ c_age + c_bmi + c_ln_na_calib + high_chol + usborn + female + bkg_pr + bkg_o,data = samp.boot,family=quasibinomial())
    ddm2 = glm(sbp ~ c_age + c_bmi + c_ln_na_calib + high_chol + usborn + female + bkg_pr + bkg_o,data = samp.boot,family=gaussian())
    
    inf.ddm1 = dfbeta(ddm1);colnames(inf.ddm1)[1] = 'Int';colnames(inf.ddm1) = paste0('Inf.ddm1.',colnames(inf.ddm1));samp.boot = cbind(samp.boot,inf.ddm1)
    inf.ddm2 = dfbeta(ddm2);colnames(inf.ddm2)[1] = 'Int';colnames(inf.ddm2) = paste0('Inf.ddm2.',colnames(inf.ddm2));samp.boot = cbind(samp.boot,inf.ddm2)
    
    # Step 3: Create twophase design, alibrate with influence functions
    twophase.design = twophase(id=list(~BGid+hhid+subid,~1),subset=~solnas,strata=list(~strat,NULL),probs=list(~P.Bg+P.hh+P.sub,NULL),data=samp.boot,method='full')
    twophase.design.ddm1 = calibrate(twophase.design,phase=2,calfun='raking',as.formula(paste0('~',paste(colnames(inf.ddm1),collapse='+'))))
    twophase.design.ddm2 = calibrate(twophase.design,phase=2,calfun='raking',as.formula(paste0('~',paste(colnames(inf.ddm2),collapse='+'))))
    
    # Step 4: Fit DDM without and with calibration
    m1.raking.moments.boot = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_bio_moments_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm1,family=quasibinomial())
    m2.raking.moments.boot = svyglm(sbp ~ c_age + c_bmi + c_ln_na_bio_moments_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm2,family=gaussian())
    
    m1.raking.glm.boot = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_bio_glm_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm1,family=quasibinomial())
    m2.raking.glm.boot = svyglm(sbp ~ c_age + c_bmi + c_ln_na_bio_glm_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design.ddm2,family=gaussian())
    
    # Extra: we wanted to see how the variance estimates from the usual model fitted with the biomarker in the subset improve with bootstrap
    m1.bio.moments.boot = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_bio_moments_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design,family=quasibinomial())
    m2.bio.moments.boot = svyglm(sbp ~ c_age + c_bmi + c_ln_na_bio_moments_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design,family=gaussian())
    
    m1.bio.glm.boot = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_bio_glm_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design,family=quasibinomial())
    m2.bio.glm.boot = svyglm(sbp ~ c_age + c_bmi + c_ln_na_bio_glm_adj + high_chol + usborn + female + bkg_pr + bkg_o,design=twophase.design,family=gaussian())
    
    m1.list.raking[[i]][[j]] = data.frame(Sim=i,Boot=j,Coeff=names(m1.true$coefficients),
                                          Bio.Moments.Est.Boot=as.numeric(m1.bio.moments.boot$coefficients), Bio.Moments.SE.Boot=as.numeric(SE(m1.bio.moments.boot)),
                                          Raking.Moments.Est.Boot=as.numeric(m1.raking.moments.boot$coefficients), Raking.Moments.SE.Boot=as.numeric(SE(m1.raking.moments.boot)),
                                          Bio.GLM.Est.Boot=as.numeric(m1.bio.glm.boot$coefficients), Bio.GLM.SE.Boot=as.numeric(SE(m1.bio.glm.boot)),
                                          Raking.GLM.Est.Boot=as.numeric(m1.raking.glm.boot$coefficients), Raking.GLM.SE.Boot=as.numeric(SE(m1.raking.glm.boot)),
                                          True.Est=as.numeric(m1.true$coefficients),True.SE = as.numeric(SE(m1.true)),
                                          Avg.Est=as.numeric(m1.avg$coefficients),Avg.SE = as.numeric(SE(m1.avg)),
                                          Bio.Moments.Est=as.numeric(m1.bio.moments$coefficients),Bio.Moments.SE = as.numeric(SE(m1.bio.moments)),
                                          Raking.Moments.Est=as.numeric(m1.raking.moments$coefficients),Raking.Moments.SE = as.numeric(SE(m1.raking.moments)),
                                          Bio.GLM.Est=as.numeric(m1.bio.glm$coefficients),Bio.GLM.SE = as.numeric(SE(m1.bio.glm)),
                                          Raking.GLM.Est=as.numeric(m1.raking.glm$coefficients),Raking.GLM.SE = as.numeric(SE(m1.raking.glm)))
    
    m2.list.raking[[i]][[j]] = data.frame(Sim=i,Boot=j,Coeff=names(m2.true$coefficients),
                                          Bio.Moments.Est.Boot=as.numeric(m2.bio.moments.boot$coefficients), Bio.Moments.SE.Boot=as.numeric(SE(m2.bio.moments.boot)),
                                          Raking.Moments.Est.Boot=as.numeric(m2.raking.moments.boot$coefficients), Raking.Moments.SE.Boot=as.numeric(SE(m2.raking.moments.boot)),
                                          Bio.GLM.Est.Boot=as.numeric(m2.bio.glm.boot$coefficients), Bio.GLM.SE.Boot=as.numeric(SE(m2.bio.glm.boot)),
                                          Raking.GLM.Est.Boot=as.numeric(m2.raking.glm.boot$coefficients), Raking.GLM.SE.Boot=as.numeric(SE(m2.raking.glm.boot)),
                                          True.Est=as.numeric(m2.true$coefficients),True.SE = as.numeric(SE(m2.true)),
                                          Avg.Est=as.numeric(m2.avg$coefficients),Avg.SE = as.numeric(SE(m2.avg)),
                                          Bio.Moments.Est=as.numeric(m2.bio.moments$coefficients),Bio.Moments.SE = as.numeric(SE(m2.bio.moments)),
                                          Raking.Moments.Est=as.numeric(m2.raking.moments$coefficients),Raking.Moments.SE = as.numeric(SE(m2.raking.moments)),
                                          Bio.GLM.Est=as.numeric(m2.bio.glm$coefficients),Bio.GLM.SE = as.numeric(SE(m2.bio.glm)),
                                          Raking.GLM.Est=as.numeric(m2.raking.glm$coefficients),Raking.GLM.SE = as.numeric(SE(m2.raking.glm))) 
  }
  m1.list.raking[[i]] = do.call(rbind,m1.list.raking[[i]])
  m2.list.raking[[i]] = do.call(rbind,m2.list.raking[[i]])
}
cat('Done!',"\n")

if(!dir.exists('./Output')){system('mkdir Output')}
save(m1.list.raking,m2.list.raking,
     file=paste0('./Output/raking.RData'))