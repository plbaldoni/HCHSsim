library(survey)
library(Hmisc)
library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)
library(car)
library(data.table)
library(microbenchmark)

sampname = 'SampleData'
S=1000 #Number of samples drawn
B=500 #Number of bootstrap runs within each sample
idx=1:1000
seed = 20190414+idx[1]+1
cat('The seed is: ',seed)
set.seed(seed)

files = paste0('./RData/',sampname,"_",str_pad(1:S,nchar(S),pad=0),".RData")
m1.list.boot = list()
m2.list.boot = list()

m1.list.cov.boot = list()
m2.list.cov.boot = list()

m1.list.mitime = list()
m2.list.mitime = list()

for(i in idx){
  load(file=files[i])
  samp <- as.data.table(samp) 
  m1.list.boot[[i]] = list()
  m2.list.boot[[i]] = list()
  
  m1.list.cov.boot[[i]] = list()
  m2.list.cov.boot[[i]] = list()
  
  # Fitting calibration model
  time.calibration <- microbenchmark({
    samp.solnas <- samp[(solnas==T),]
    lm.lsodi <- glm(c_ln_na_bio1 ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,data=samp.solnas)
    samp[,c_ln_na_calib := predict(lm.lsodi,newdata=samp,'response')]
  },times = 1)
  
  # Fitting 'true' (unobserved), 'naive' (2-day avg.) model, 'calibrated', biomarker SOLNAS, and biomarker entire cohort (unobserved)
  samp.design = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=samp)
  
  m1.true = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_true + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
  m2.true = svyglm(sbp ~ c_age + c_bmi + c_ln_na_true + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
  
  m1.avg = svyglm(hypertension ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
  m2.avg = svyglm(sbp ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
  
  time.model.m1 <- microbenchmark({
    m1.calib <- svyglm(hypertension ~ c_age + c_bmi + c_ln_na_calib + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())
  },times = 1)
  time.model.m2 <- microbenchmark({
    m2.calib <- svyglm(sbp ~ c_age + c_bmi + c_ln_na_calib + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
  },times = 1)
  
  # Bootstrap starts now
  time.bootstrap.m1 <- microbenchmark({
    for(j in 1:B){cat('Sample: ',i,". Boot: ",j,".\n")
      
      # Fit calibration equation
      samp.solnas.boot <- samp.solnas[sample(1:nrow(samp.solnas),nrow(samp.solnas),replace = T),]
      lm.lsodi.boot <- glm(c_ln_na_bio1 ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,data=samp.solnas.boot)
      samp.design <- update(samp.design,c_ln_na_calib.boot = predict(lm.lsodi.boot,newdata=samp,'response'))
      
      # Fitting diet-disease model
      m1.calib.boot <- svyglm(hypertension ~ c_age + c_bmi + c_ln_na_calib.boot + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=quasibinomial())

      # Saving parameters
      m1.list.boot[[i]][[j]] <- data.frame(Sim=i,Boot=j,Coeff=names(m1.calib.boot$coefficients),
                                          Est=as.numeric(m1.calib.boot$coefficients),SE=as.numeric(SE(m1.calib.boot)),
                                          True.Est=as.numeric(m1.true$coefficients),True.SE = as.numeric(SE(m1.true)),
                                          Avg.Est=as.numeric(m1.avg$coefficients),Avg.SE = as.numeric(SE(m1.avg)),
                                          Calib.Est=as.numeric(m1.calib$coefficients),Calib.SE = as.numeric(SE(m1.calib)))
    
      m1.list.cov.boot[[i]][[j]] <- m1.calib.boot$cov.unscaled
    }
  },times = 1)
  
  time.bootstrap.m2 <- microbenchmark({
    for(j in 1:B){cat('Sample: ',i,". Boot: ",j,".\n")
      
      # Fit calibration equation
      samp.solnas.boot <- samp.solnas[sample(1:nrow(samp.solnas),nrow(samp.solnas),replace = T),]
      lm.lsodi.boot <- glm(c_ln_na_bio1 ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,data=samp.solnas.boot)
      samp.design <- update(samp.design,c_ln_na_calib.boot = predict(lm.lsodi.boot,newdata=samp,'response'))
      
      # Fitting diet-disease model
      m2.calib.boot <- svyglm(sbp ~ c_age + c_bmi + c_ln_na_calib.boot + high_chol + usborn + female + bkg_pr + bkg_o,design=samp.design,family=gaussian())
      
      # Saving parameters
      m2.list.boot[[i]][[j]] <- data.frame(Sim=i,Boot=j,Coeff=names(m2.calib.boot$coefficients),
                                          Est=as.numeric(m2.calib.boot$coefficients),SE=as.numeric(SE(m2.calib.boot)),
                                          True.Est=as.numeric(m2.true$coefficients),True.SE = as.numeric(SE(m2.true)),
                                          Avg.Est=as.numeric(m2.avg$coefficients),Avg.SE = as.numeric(SE(m2.avg)),
                                          Calib.Est=as.numeric(m2.calib$coefficients),Calib.SE = as.numeric(SE(m2.calib)))
      
      m2.list.cov.boot[[i]][[j]] <- m2.calib.boot$cov.unscaled
    }
  },times = 1)
  
  m1.list.boot[[i]] <- do.call(rbind,m1.list.boot[[i]])
  m2.list.boot[[i]] <- do.call(rbind,m2.list.boot[[i]])
  
  m1.list.mitime[[i]] <- data.table(Sim = i, Time.Calibration = summary(time.calibration,'s')$mean, Time.Model = summary(time.model.m1,'s')$mean, Time.Resampling = summary(time.bootstrap.m1,'s')$mean)
  m2.list.mitime[[i]] <- data.table(Sim = i, Time.Calibration = summary(time.calibration,'s')$mean, Time.Model = summary(time.model.m2,'s')$mean, Time.Resampling = summary(time.bootstrap.m2,'s')$mean)
}
cat('Done!',"\n")

if(!dir.exists('./Output')){system('mkdir Output')}
save(m1.list.boot,m2.list.boot,
     m1.list.cov.boot,m2.list.cov.boot,
     m1.list.mitime,m2.list.mitime,
     file=paste0('./Output/bootstrap.RData'))