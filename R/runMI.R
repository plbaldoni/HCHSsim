library(survey)
library(Hmisc)
library(dplyr)
library(MASS)
library(magrittr)
library(ggplot2)
library(stringr)
library(car)
library(data.table)

sampname = 'SampleData'
S=1000 #Number of samples drawn
B=500 #Number of bootstrap runs within each sample
idx=1:1000
seed = 20190414+idx[1]+2
cat('The seed is: ',seed)
set.seed(seed)

files = paste0('./RData/',sampname,"_",str_pad(1:S,nchar(S),pad=0),".RData")
m1.list.mi = list()
m2.list.mi = list()

m1.list.cov.mi = list()
m2.list.cov.mi = list()

for(i in idx){
  load(file=files[i])
  samp <- as.data.table(samp) 
  m1.list.mi[[i]] = list()
  m2.list.mi[[i]] = list()
  
  m1.list.cov.mi[[i]] = list()
  m2.list.cov.mi[[i]] = list()
  
  # Fitting calibration model
  samp.solnas <- samp[(solnas==T),]
  lm.lsodi = glm(c_ln_na_bio1 ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,data=samp.solnas)
  x.lsodi = model.matrix(lm.lsodi) #X
  xtx.lsodi = t(x.lsodi)%*%x.lsodi #X`X`
  ixtx.lsodi = solve(xtx.lsodi) #(X`X)^-1
  samp[,c_ln_na_calib := predict(lm.lsodi,newdata=samp,'response')]
  
  # Fitting 'true' (unobserved), 'naive' (2-day avg.) model, 'calibrated', biomarker SOLNAS, and biomarker entire cohort (unobserved)
  samp.design = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=samp)
  
  m1.true = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_true + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
  m2.true = svyglm(sbp ~ c_age + c_bmi + c_ln_na_true + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
  
  m1.avg = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
  m2.avg = svyglm(sbp ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
  
  m1.calib = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_calib + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
  m2.calib = svyglm(sbp ~ c_age + c_bmi + c_ln_na_calib + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
  
  # Multiple Imputation starts now
  for(j in 1:B){cat('Sample: ',i,". MI: ",j,".\n")
    
    # Sampling new variance-covariance matrix from inverse chi-squared distribution
    samp.solnas.boot <- samp.solnas[sample(1:nrow(samp.solnas),nrow(samp.solnas),replace = T),]
    vcov.star = (((sigma(lm.lsodi)^2)*(lm.lsodi$df.residual))/rchisq(n=1,df=lm.lsodi$df.residual))*ixtx.lsodi
    samp.design <- update(samp.design,c_ln_na_calib.mi = with(samp,cbind(1,c_age,c_bmi,c_ln_na_avg,high_chol,usborn,female,bkg_pr,bkg_o))%*%mvrnorm(n=1,mu=as.numeric(lm.lsodi$coefficients),Sigma=vcov.star))
    
    # Fitting diet-disease model
    m1.calib.mi = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_calib.mi + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
    m2.calib.mi = svyglm(sbp ~ c_age + c_bmi + c_ln_na_calib.mi + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
    
    # Saving parameters
    m1.list.mi[[i]][[j]] = data.frame(Sim=i,MI=j,Coeff=names(m1.calib.mi$coefficients),
                                      Est=as.numeric(m1.calib.mi$coefficients),SE=as.numeric(SE(m1.calib.mi)),
                                      True.Est=as.numeric(m1.true$coefficients),True.SE = as.numeric(SE(m1.true)),
                                      Avg.Est=as.numeric(m1.avg$coefficients),Avg.SE = as.numeric(SE(m1.avg)),
                                      Calib.Est=as.numeric(m1.calib$coefficients),Calib.SE = as.numeric(SE(m1.calib)))
    
    m2.list.mi[[i]][[j]] = data.frame(Sim=i,MI=j,Coeff=names(m2.calib.mi$coefficients),
                                      Est=as.numeric(m2.calib.mi$coefficients),SE=as.numeric(SE(m2.calib.mi)),
                                      True.Est=as.numeric(m2.true$coefficients),True.SE = as.numeric(SE(m2.true)),
                                      Avg.Est=as.numeric(m2.avg$coefficients),Avg.SE = as.numeric(SE(m2.avg)),
                                      Calib.Est=as.numeric(m2.calib$coefficients),Calib.SE = as.numeric(SE(m2.calib)))
    
    m1.list.cov.mi[[i]][[j]] = m1.calib.mi$cov.unscaled
    m2.list.cov.mi[[i]][[j]] = m2.calib.mi$cov.unscaled
    
  }
  m1.list.mi[[i]] = do.call(rbind,m1.list.mi[[i]])
  m2.list.mi[[i]] = do.call(rbind,m2.list.mi[[i]])
}
cat('Done!',"\n")

if(!dir.exists('./Output')){system('mkdir Output')}
save(m1.list.mi,m2.list.mi,
     m1.list.cov.mi,m2.list.cov.mi,
     file=paste0('./Output/MI.RData'))