### Generate data for TargetPopulation.RData ###

### Loading libraries
lapply(c('data.table','mvtnorm','magrittr','survey','dplyr','tidyr','ggplot2','ggpubr','tidyverse','kableExtra','ggjoy'),library,character.only = T)

### Expit function
expit = function(xb){return(exp(xb)/(1+exp(xb)))}

### Setting seed
set.seed(04142019)

### Loading data for Multivariate Gaussian ###
load('./data/multiNorm.RData')

### Loading parameter estimates to simulate binary variables
load('./data/logiReg.RData')

### Loading Target Population Data ###
load(paste0('./RData/TargetPopulation.RData'))
pop = pop[v.num==1,] 

### sim_sexbkg: simulate multinomial distribution of sex*background ###
sim_sexbkg = function(dat,p,sex=c('F','M'),bkg=c('D','PR','O')){
  tb =  expand.grid(sex,bkg);colnames(tb) = c('sex','bkg')
  tb$idx = 1:nrow(tb)
  tb$prob = p
  
  dat$idx = apply(rmultinom(nrow(dat),size=1,prob=tb$prob),2,which.max)
  dat = merge(dat,subset(tb,select=-prob),by='idx',all.x=T)
  dat = dat[order(dat$strat,dat$BGid,dat$hhid,dat$subid),]
  return(dat)
}

### sim_covar: simulate covariates ###
sim_covar = function(dat,ranges,means,covs,dat.par,count.seed=1,
                     covars=c('AGE','BMI','LN_NA_AVG','LN_K_AVG','LN_KCAL_AVG','LN_PROTEIN_AVG')){
  dat <- as.data.frame(dat)
  ncovars = length(covars)
  df = data.frame()
  groups = names(table(dat$idx))
  namevar = c('age','bmi','ln_na_avg','ln_k_avg','ln_kcal_avg','ln_protein_avg')
  
  m = list()
  v = list()
  
  # F=Female, M=Male
  # D= Dominican, PR=Puerto Rican, O=other
  m[[1]] = as.numeric(means[1,]) #F,D
  m[[2]] = as.numeric(means[2,]) #M,D
  m[[3]] = as.numeric(means[3,]) #F,PR
  m[[4]] = as.numeric(means[4,]) #M,PR
  m[[5]] = as.numeric(means[5,]) #F,O
  m[[6]] = as.numeric(means[6,]) #M,O
  
  v[[1]] = matrix(as.numeric(covs[1,]),ncovars,ncovars,byrow = T) #F,D
  v[[2]] = matrix(as.numeric(covs[2,]),ncovars,ncovars,byrow = T) #M,D
  v[[3]] = matrix(as.numeric(covs[3,]),ncovars,ncovars,byrow = T) #F,PR
  v[[4]] = matrix(as.numeric(covs[4,]),ncovars,ncovars,byrow = T) #M,PR
  v[[5]] = matrix(as.numeric(covs[5,]),ncovars,ncovars,byrow = T) #F,O
  v[[6]] = matrix(as.numeric(covs[6,]),ncovars,ncovars,byrow = T) #M,O
  
  for(i in groups){cat('Group: ',i,".\n")
    subdat = subset(dat,idx==i)
    
    m.topass = m[[as.numeric(i)]]
    v.topass = v[[as.numeric(i)]]
    
    ### Simulating Age, BMI, and log-transformed nutrients
    subdat = cbind(subdat,rmvnorm(nrow(subdat),mean=m.topass,sigma=v.topass))
    colnames(subdat) = c(colnames(dat),namevar)
    
    idx.age = with(subdat,which(!between(age,ranges[1,1],ranges[1,2]) |
                                  !between(bmi,ranges[2,1],ranges[2,2]) |
                                  !between(ln_na_avg,ranges[3,1],ranges[3,2]) |
                                  !between(ln_k_avg,ranges[4,1],ranges[4,2]) |
                                  !between(ln_kcal_avg,ranges[5,1],ranges[5,2]) |
                                  !between(ln_protein_avg,ranges[6,1],ranges[6,2])))
    ### Checking for values out of the range in the population data
    for(j in idx.age){cat(count.seed,"\r")
      flag = T
      while(flag){
        subdat[j,namevar] <- rmvnorm(1,mean=m.topass,sigma=v.topass)
        count.seed = count.seed + 1
        flag = with(subdat,!between(age[j],ranges[1,1],ranges[1,2]) |
                      !between(bmi[j],ranges[2,1],ranges[2,2]) |
                      !between(ln_na_avg[j],ranges[3,1],ranges[3,2]) |
                      !between(ln_k_avg[j],ranges[4,1],ranges[4,2]) |
                      !between(ln_kcal_avg[j],ranges[5,1],ranges[5,2]) |
                      !between(ln_protein_avg[j],ranges[6,1],ranges[6,2]))
      }
    }
    ### Now, simulating binary variables ###
    subdat.par = subset(dat.par,idx==i & MODEL =='USBORN',select=c(idx,Variable,Estimate))
    subdat.par = subdat.par[match(c("Intercept","AGE","BMI","LN_NA_AVG","LN_K_AVG","LN_KCAL_AVG","LN_PROTEIN_AVG"),subdat.par$Variable),]
    p.par = expit(as.matrix(cbind('Intercept'=1,subdat[,namevar]))%*%subdat.par$Estimate)
    subdat$usborn = as.numeric(1*(runif(nrow(p.par))<=p.par))
    
    subdat.par = subset(dat.par,idx==i & MODEL =='CHOL',select=c(idx,Variable,Estimate))
    subdat.par = subdat.par[match(c("Intercept","AGE","BMI","LN_NA_AVG","LN_K_AVG","LN_KCAL_AVG","LN_PROTEIN_AVG","US_BORN"),subdat.par$Variable),]
    p.par = expit(as.matrix(cbind('Intercept'=1,subdat[,c(namevar,'usborn')]))%*%subdat.par$Estimate)
    subdat$high_chol = as.numeric(1*(runif(nrow(p.par))<=p.par))
    
    ### Saving data
    df = rbind(df,subdat)
  }
  df1 = df[order(df$strat,df$BGid,df$hhid,df$subid),];rm(df)
  df1$age.strat = (df1$age>=45)
  df1$v.num=1
  
  #df2 = df1
  #df2$v.num = 2
  
  #df = rbind(df1,df2)
  df = df1
  df = df[order(df$strat,df$BGid,df$hhid,df$subid,df$v.num),]
  
  cat('Done!',"\n")
  df = as.data.table(df)
  
  return(list('Data'=df,'HCHS.Mean'=m,'HCHS.Cov'=v))
}

ls.pop = sim_sexbkg(dat=pop,p=c(0.1953,0.1300,0.2065,0.1984,0.1421,0.1277))
ls.pop = sim_covar(dat=ls.pop,ranges = ranges,means = means,covs = covs,
                   dat.par = dat.par)

### Checking the distribution of the simulated data

ls.pop$Data[,.(usborn = mean(usborn),high_chol = mean(high_chol)),by=c('sex','bkg')]
ls.pop$Data[,.(age = mean(age),bmi = mean(bmi),
               ln_na_avg = mean(ln_na_avg),ln_k_avg = mean(ln_k_avg),
               ln_kcal_avg = mean(ln_kcal_avg),ln_protein_avg = mean(ln_protein_avg)),by=c('sex','bkg')]
do.call(rbind,ls.pop$HCHS.Mean)

### Centering continous naive nutrients
ls.pop$Data$c_ln_na_avg = as.numeric(scale(ls.pop$Data$ln_na_avg,scale = F))
ls.pop$Data$c_ln_k_avg = as.numeric(scale(ls.pop$Data$ln_k_avg,scale = F))
ls.pop$Data$c_ln_kcal_avg = as.numeric(scale(ls.pop$Data$ln_kcal_avg,scale = F))
ls.pop$Data$c_ln_protein_avg = as.numeric(scale(ls.pop$Data$ln_protein_avg,scale = F))

ls.pop$Data$c_age = as.numeric(scale(ls.pop$Data$age,scale = F))
ls.pop$Data$c_bmi = as.numeric(scale(ls.pop$Data$bmi,scale = F))
ls.pop$Data$female = 1*(ls.pop$Data$sex=='F')
ls.pop$Data$bkg_pr = 1*(ls.pop$Data$bkg=='PR')
ls.pop$Data$bkg_o = 1*(ls.pop$Data$bkg=='O')

### Saving the data
dt.pop <- ls.pop$Data

### Now loading model coefficients and variance parameter for simulation of true and biomarker intake ###
load('./data/calibrCoeff.RData')

dt.pop$ln_na_true = as.numeric(as.matrix(cbind('Int'=1,subset(dt.pop,select=c('c_ln_na_avg','c_age','c_bmi','female','usborn','high_chol','bkg_pr','bkg_o'))))%*%sodicoeff$Estimate+rnorm(nrow(dt.pop),0,sqrt(var.df$Var.True[1])))
dt.pop$ln_k_true = as.numeric(as.matrix(cbind('Int'=1,subset(dt.pop,select=c('c_ln_k_avg','c_age','c_bmi','female','usborn','high_chol','bkg_pr','bkg_o'))))%*%potacoeff$Estimate+rnorm(nrow(dt.pop),0,sqrt(var.df$Var.True[2])))
dt.pop$ln_kcal_true = as.numeric(as.matrix(cbind('Int'=1,subset(dt.pop,select=c('c_ln_kcal_avg','c_age','c_bmi','female','usborn','high_chol','bkg_pr','bkg_o'))))%*%kcalcoeff$Estimate+rnorm(nrow(dt.pop),0,sqrt(var.df$Var.True[2])))
dt.pop$ln_protein_true = as.numeric(as.matrix(cbind('Int'=1,subset(dt.pop,select=c('c_ln_protein_avg','c_age','c_bmi','female','usborn','high_chol','bkg_pr','bkg_o'))))%*%protcoeff$Estimate+rnorm(nrow(dt.pop),0,sqrt(var.df$Var.True[2])))

# Centering the true nutrients

dt.pop$c_ln_na_true = as.numeric(scale(dt.pop$ln_na_true,scale=F))
dt.pop$c_ln_k_true = as.numeric(scale(dt.pop$ln_k_true,scale=F))
dt.pop$c_ln_kcal_true = as.numeric(scale(dt.pop$ln_kcal_true,scale=F))
dt.pop$c_ln_protein_true = as.numeric(scale(dt.pop$ln_protein_true,scale=F))

#Checking correlations:
with(dt.pop,cor(c_ln_na_true,c_ln_na_avg))
with(dt.pop,cor(c_ln_k_true,c_ln_k_avg))
with(dt.pop,cor(c_ln_kcal_true,c_ln_kcal_avg))
with(dt.pop,cor(c_ln_protein_true,c_ln_protein_avg))

### Now, simulating biomarker
dt.pop$ln_na_bio1 = dt.pop$ln_na_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[1]))
dt.pop$ln_na_bio2 = dt.pop$ln_na_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[1]))

dt.pop$ln_k_bio1 = dt.pop$ln_k_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[2]))
dt.pop$ln_k_bio2 = dt.pop$ln_k_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[2]))

dt.pop$ln_kcal_bio1 = dt.pop$ln_kcal_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[3]))
dt.pop$ln_kcal_bio2 = dt.pop$ln_kcal_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[3]))

dt.pop$ln_protein_bio1 = dt.pop$ln_protein_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[4]))
dt.pop$ln_protein_bio2 = dt.pop$ln_protein_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[4]))

### Centering the biomarkers
dt.pop$c_ln_na_bio1 = as.numeric(scale(dt.pop$ln_na_bio1,scale = F))
dt.pop$c_ln_na_bio2 = as.numeric(scale(dt.pop$ln_na_bio2,scale = F))

dt.pop$c_ln_k_bio1 =  as.numeric(scale(dt.pop$ln_k_bio1,scale = F))
dt.pop$c_ln_k_bio2 = as.numeric(scale(dt.pop$ln_k_bio2,scale = F))

dt.pop$c_ln_kcal_bio1 = as.numeric(scale(dt.pop$ln_kcal_bio1,scale = F))
dt.pop$c_ln_kcal_bio2 = as.numeric(scale(dt.pop$ln_kcal_bio2,scale = F))

dt.pop$c_ln_protein_bio1 = as.numeric(scale(dt.pop$ln_protein_bio1,scale = F))
dt.pop$c_ln_protein_bio2 = as.numeric(scale(dt.pop$ln_protein_bio2,scale = F))

#Checking correlations:
with(dt.pop,cor(c_ln_na_bio1,c_ln_na_true)) #0.81
with(dt.pop,cor(c_ln_k_bio1,c_ln_k_true)) #0.80
with(dt.pop,cor(c_ln_kcal_bio1,c_ln_kcal_true)) #0.98
with(dt.pop,cor(c_ln_protein_bio1,c_ln_protein_true)) #0.85

with(dt.pop,cor(c_ln_na_bio1,c_ln_na_avg)) #0.39
with(dt.pop,cor(c_ln_k_bio1,c_ln_k_avg)) # 0.43
with(dt.pop,cor(c_ln_kcal_bio1,c_ln_kcal_avg)) #0.26
with(dt.pop,cor(c_ln_protein_bio1,c_ln_protein_avg)) #0.31

### Now, I will simulate two outcomes (continuous sbp and hypertension2_indicator)
load('./data/outPar.RData')

### Now, simulating hypertension outcome ###
suboutcome.par = data.frame(Variable = c('Intercept',all.vars(m1.formula[-2])), Estimate = m1.coefficients)
# Changing the main effect size: 1.3 OR for a 20% increase (1.2) in Sodium. From Prentice (2017)
suboutcome.par$Estimate[suboutcome.par$Variable=='C_LN_NA_CALIBR'] = log(1.3)/log(1.2)
# Adjusting the intercept for an overall prevalence of 0.08 (from real data)
suboutcome.par$Estimate[suboutcome.par$Variable=='Intercept'] = log(0.08/(1-0.08))-sum(suboutcome.par$Estimate[-1]*apply(subset(dt.pop,select=c('c_age','c_bmi','c_ln_na_true','high_chol','usborn','female','bkg_pr','bkg_o')),2,mean))

suboutcome.par = suboutcome.par[match(c("Intercept","C_AGE","C_BMI","C_LN_NA_CALIBR",'HIGH_TOTAL_CHOL','US_BORN','FEMALE','BKGRD_PR','BKGRD_OTHER'),suboutcome.par$Variable),]
p.par = expit(as.matrix(cbind('Intercept'=1,dt.pop[,c('c_age','c_bmi','c_ln_na_true','high_chol','usborn','female','bkg_pr','bkg_o')]))%*%suboutcome.par$Estimate)
dt.pop$hypertension = as.numeric(1*(runif(nrow(p.par))<=p.par))

htn_parameters <- suboutcome.par
if(!dir.exists('./Output')){system('mkdir Output')}
save(htn_parameters,file = './Output/htn_parameters.RData')
rm(htn_parameters)

### Now, simulating SBP outcome ###
# Based on real data, 1 unit increase in log_NA (or exp(1)~2.7 times increase in sodium)
# leads to an increment of 1.4058215 units in SBP. Very small effect.
# I will make the effect of sodium larger.
# I will do: 
# leads to an increment of 5 units in SBP. Large effect

suboutcome.par = data.frame(Variable = c('Intercept',all.vars(m2.formula[-2])), Estimate = m2.coefficients)
# Increasing the main effect size: 0.182 unit increase in log_NA (or exp(0.182)~1.2 times increase in sodium)
suboutcome.par$Estimate[suboutcome.par$Variable=='C_LN_NA_CALIBR'] = 5/0.182
# Adjusting the intercept for an overall average of 120 SBP
suboutcome.par$Estimate[suboutcome.par$Variable=='Intercept'] = 120-sum(suboutcome.par$Estimate[-1]*apply(subset(dt.pop,select=c('c_age','c_bmi','c_ln_na_true','high_chol','usborn','female','bkg_pr','bkg_o')),2,mean))

suboutcome.par = suboutcome.par[match(c("Intercept","C_AGE","C_BMI","C_LN_NA_CALIBR",'HIGH_TOTAL_CHOL','US_BORN','FEMALE','BKGRD_PR','BKGRD_OTHER'),suboutcome.par$Variable),]
xb.par = as.matrix(cbind('Intercept'=1,dt.pop[,c('c_age','c_bmi','c_ln_na_true','high_chol','usborn','female','bkg_pr','bkg_o')]))%*%suboutcome.par$Estimate
dt.pop$sbp = as.numeric((xb.par+rnorm(nrow(xb.par),0,m2.sigma)))

sbp_parameters <- list()
sbp_parameters[['coefficients']] <- suboutcome.par
sbp_parameters[['variance']] <- m2.sigma^2
if(!dir.exists('./Output')){system('mkdir Output')}
save(sbp_parameters,file = './Output/sbp_parameters.RData')
rm(sbp_parameters)

### Organizing, Centering and Saving data
dt.pop$bkg %<>% as.character()

pop1 = dt.pop
pop1$v.num=1

pop2 = dt.pop
pop2$v.num=2

pop = rbind(pop1,pop2)

pop = pop[order(pop$strat,pop$BGid,pop$hhid,pop$subid,pop$v.num),]
save(pop,file=paste0('./RData/TargetPopulationData.RData'),compress = 'xz')

