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

### Organizing, Centering and Saving data
dt.pop$bkg %<>% as.character()

pop1 = dt.pop
pop1$v.num=1

pop2 = dt.pop
pop2$v.num=2

pop = rbind(pop1,pop2)

pop = pop[order(pop$strat,pop$BGid,pop$hhid,pop$subid,pop$v.num),]
save(pop,file=paste0('./RData/TargetPopulationData.RData'),compress = 'xz')

###########################################
############# Creating Plots ##############
###########################################

pt.title = 10
pt.text = 10

sdnutr = round(apply(pop[pop$v.num==1,c('ln_na_bio1','ln_na_avg','ln_na_true')],2,sd),2)
pop.na <- as_tibble(pop[pop$v.num==1,c('ln_na_true','ln_na_bio1','ln_na_avg')]) %>%
  gather(Type,Sodium,ln_na_true:ln_na_avg) %>%
  mutate(Type = recode(Type,ln_na_true = 'True', ln_na_bio1 = 'Biomarker', ln_na_avg='Self-reported'))

naplot = ggplot(pop.na,aes(Sodium,color=Type,fill=Type))+
  geom_density(alpha=0.15)+
  theme_bw()+ylab('Density')+xlab('Log-Sodium')+
  scale_fill_manual(labels=as.character(apply(cbind(c('Biomarker','Self-reported','True'),sdnutr),1,FUN = function(x){paste0(x[1],' (SD=',x[2],')')})),values=c('darkgreen','darkred','darkblue'))+
  scale_color_manual(labels=as.character(apply(cbind(c('Biomarker','Self-reported','True'),sdnutr),1,FUN = function(x){paste0(x[1],' (SD=',x[2],')')})),values=c('darkgreen','darkred','darkblue'))+
  theme(legend.title=element_blank(),legend.position = c(0.275,0.85),legend.text = element_text(size=pt.title),legend.background =  element_rect(fill = "transparent", colour = "transparent"))+
  guides(color = guide_legend(override.aes = list(size = rel(0.3))))+ylim(0,1.35)

sdnutr = round(apply(pop[pop$v.num==1,c('ln_k_bio1','ln_k_avg','ln_k_true')],2,sd),2)
pop.k <- as_tibble(pop[pop$v.num==1,c('ln_k_true','ln_k_bio1','ln_k_avg')]) %>%
  gather(Type,Potassium,ln_k_true:ln_k_avg) %>%
  mutate(Type = recode(Type,ln_k_true = 'True', ln_k_bio1 = 'Biomarker', ln_k_avg='Self-reported'))

kplot = ggplot(pop.k,aes(Potassium,color=Type,fill=Type))+
  geom_density(alpha=0.15)+
  theme_bw()+ylab('Density')+xlab('Log-Potassium')+
  scale_fill_manual(labels=as.character(apply(cbind(c('Biomarker','Self-reported','True'),sdnutr),1,FUN = function(x){paste0(x[1],' (SD=',x[2],')')})),values=c('darkgreen','darkred','darkblue'))+
  scale_color_manual(labels=as.character(apply(cbind(c('Biomarker','Self-reported','True'),sdnutr),1,FUN = function(x){paste0(x[1],' (SD=',x[2],')')})),values=c('darkgreen','darkred','darkblue'))+
  theme(legend.title=element_blank(),legend.position = c(0.275,0.85),legend.text = element_text(size=pt.title),legend.background =  element_rect(fill = "transparent", colour = "transparent"))+
  guides(color = guide_legend(override.aes = list(size = rel(0.3))))+ylim(0,1.35)


sdnutr = round(apply(pop[pop$v.num==1,c('ln_kcal_bio1','ln_kcal_avg','ln_kcal_true')],2,sd),2)
pop.kcal <- as_tibble(pop[pop$v.num==1,c('ln_kcal_true','ln_kcal_bio1','ln_kcal_avg')]) %>%
  gather(Type,Kcal,ln_kcal_true:ln_kcal_avg) %>%
  mutate(Type = recode(Type,ln_kcal_true = 'True', ln_kcal_bio1 = 'Biomarker', ln_k_avg='Self-reported'))

kcalplot = ggplot(pop.kcal,aes(Kcal,color=Type,fill=Type))+
  geom_density(alpha=0.15)+
  theme_bw()+ylab('Density')+xlab('Log-Energy')+
  scale_fill_manual(labels=as.character(apply(cbind(c('Biomarker','Self-reported','True'),sdnutr),1,FUN = function(x){paste0(x[1],' (SD=',x[2],')')})),values=c('darkgreen','darkred','darkblue'))+
  scale_color_manual(labels=as.character(apply(cbind(c('Biomarker','Self-reported','True'),sdnutr),1,FUN = function(x){paste0(x[1],' (SD=',x[2],')')})),values=c('darkgreen','darkred','darkblue'))+
  theme(legend.title=element_blank(),legend.position = c(0.275,0.85),legend.text = element_text(size=pt.title),legend.background =  element_rect(fill = "transparent", colour = "transparent"))+
  guides(color = guide_legend(override.aes = list(size = rel(0.3))))+ylim(0,1.35)

sdnutr = round(apply(pop[pop$v.num==1,c('ln_protein_bio1','ln_protein_avg','ln_protein_true')],2,sd),2)
pop.protein <- as_tibble(pop[pop$v.num==1,c('ln_protein_true','ln_protein_bio1','ln_protein_avg')]) %>%
  gather(Type,Protein,ln_protein_true:ln_protein_avg) %>%
  mutate(Type = recode(Type,ln_protein_true = 'True', ln_protein_bio1 = 'Biomarker', ln_k_avg='Self-reported'))

proteinplot = ggplot(pop.protein,aes(Protein,color=Type,fill=Type))+
  geom_density(alpha=0.15)+
  theme_bw()+ylab('Density')+xlab('Log-Protein')+
  scale_fill_manual(labels=as.character(apply(cbind(c('Biomarker','Self-reported','True'),sdnutr),1,FUN = function(x){paste0(x[1],' (SD=',x[2],')')})),values=c('darkgreen','darkred','darkblue'))+
  scale_color_manual(labels=as.character(apply(cbind(c('Biomarker','Self-reported','True'),sdnutr),1,FUN = function(x){paste0(x[1],' (SD=',x[2],')')})),values=c('darkgreen','darkred','darkblue'))+
  theme(legend.title=element_blank(),legend.position = c(0.275,0.85),legend.text = element_text(size=pt.title),legend.background =  element_rect(fill = "transparent", colour = "transparent"))+
  guides(color = guide_legend(override.aes = list(size = rel(0.3))))+ylim(0,1.35)

figure1 <- ggarrange(naplot+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                     kplot+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                     kcalplot+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                     proteinplot+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                     nrow = 2,ncol = 2,labels = list('A','B','C','D'))

# Organizing plot
if(!dir.exists('./Output')){system('mkdir Output')}
ggsave(figure1,filename = './Output/Figure1.eps',
       dpi = 'retina',width = 8.5,height = 7,device = cairo_ps,fallback_resolution = 300)

## Sodium-only figure

sdnutr = round(apply(pop[pop$v.num==1,c('ln_na_true','ln_na_bio1','ln_na_avg')],2,sd),2)
pop.na <- as_tibble(pop[pop$v.num==1,c('ln_na_true','ln_na_bio1','ln_na_avg')]) %>%
  gather(Type,Sodium,ln_na_true:ln_na_avg) %>%
  mutate(Type = recode(Type,ln_na_true = 'True', ln_na_bio1 = 'Biomarker', ln_na_avg='Self-reported'))
pop.na$Type %<>% factor(levels = c('True','Biomarker','Self-reported'))

figure1.sodiumonly.naplot = ggplot(pop.na,aes(x=Sodium,linetype = Type))+
  stat_density(geom="line", position="identity")+
  theme_bw()+ylab('Density')+xlab('Log-Sodium')+
  theme_bw() + 
  scale_linetype_manual(values = c('solid','longdash','dotted'),
                        labels=as.character(apply(cbind(c('True','Biomarker','Self-reported'),sdnutr),1,FUN = function(x){paste0(x[1],' (SD=',x[2],')')})))+
  guides(linetype = guide_legend(override.aes = list(size = rel(0.3))))+ylim(0,1)+
  theme(legend.title=element_blank(),legend.position = 'top',legend.direction = 'vertical',legend.text = element_text(size=pt.title),legend.background =  element_rect(fill = "transparent", colour = "transparent"))

dt1 <- melt(pop[pop$v.num==1,c('subid','sex','bkg','age','bmi','ln_na_true','ln_na_bio1','ln_na_avg')],
            id.vars = c('subid','sex','bkg'),measure.vars = c('ln_na_true','ln_na_bio1','ln_na_avg'),variable.name = 'Nutrient',value.name = 'Value')
dt1$sex %<>% plyr::mapvalues(from = c('M','F'),to = c('Male','Female')) %<>% factor(levels = c('Female','Male'))
dt1$bkg %<>% plyr::mapvalues(from = c('D','PR','O'),to = c('Dominican','Puerto Rican','Other')) %<>% factor(levels = c('Dominican','Puerto Rican','Other'))
dt1$Nutrient %<>% plyr::mapvalues(from = c('ln_na_true','ln_na_bio1','ln_na_avg'),to = c('True','Biomarker','Self-\nreported')) %<>% factor(levels = c('True','Biomarker','Self-\nreported'))

figure1.sodiumonly.boxplot <-  ggplot(dt1,aes(y = Value,x = Nutrient))+
  facet_grid(cols = vars(bkg),rows = vars(sex))+
  geom_boxplot(outlier.alpha = 0.25,outlier.size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_blank(),legend.title = element_blank(),
        legend.position = 'bottom',legend.direction = 'horizontal',
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.text = element_text(size=pt.title),axis.text.x = element_text(angle = 30, hjust = 1))+
  ylab('Log-Sodium') + xlab('Nutrient')

# figure1.sodiumonly.reg <- ggplot(merge(dt1,pop[pop$v.num==1,c('subid','sbp')],by = 'subid',all.x = T))+
#   facet_grid(cols = vars(bkg),rows = vars(sex))+
#   geom_smooth(method = lm,aes(x = Value,y = sbp,linetype = Nutrient),color = 'black',se = T,size = 0.5) +
#   scale_linetype_manual(values = c('solid','longdash','dotted')) +
#   theme_bw()+
#   guides(linetype = guide_legend(override.aes = list(size = rel(0.3))))+
#   xlab('Log-Sodium')+ylab('Systolic Blood Pressure')+
#   theme(legend.title=element_blank(),legend.position = 'none')

# figure1.sodiumonly <- ggarrange(figure1.sodiumonly.naplot+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
#                                 figure1.sodiumonly.bottom,
#                                 ncol = 1,nrow = 2,labels = list('A'),heights = c(0.35,0.65))

figure1.sodiumonly <- ggarrange(figure1.sodiumonly.naplot+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                                       figure1.sodiumonly.boxplot+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                                       ncol = 2,nrow = 1,labels = list('A','B'))

# Organizing plot
if(!dir.exists('./Output')){system('mkdir Output')}
ggsave(figure1.sodiumonly,filename = './Output/Figure1_SodiumOnly.eps',
       dpi = 'retina',width = 8,height = 4,device = cairo_ps,fallback_resolution = 300)

###########################################
############# Creating Table ##############
###########################################

pop.v1 <- pop[(v.num==1),]

dt1 <- pop.v1[,.(N = format(.N,big.mark=',',trim=T),Age = round(100*mean(age.strat),2)),by=c('strat','hisp.strat','bkg','sex')]

dt2 <- pop.v1[,.(NBG = length(unique(BGid)), NHH = format(length(unique(hhid)),big.mark=',',trim=T)),by = c('strat')]

dt3 <- pop.v1[,.(N = length(unique(hhid))),by=c('strat','hisp.strat')]
dt3[,HispStrat := paste0(ifelse(hisp.strat,'Hisp.','Non-hisp.'),' (',format(N,big.mark = ',',trim = T),')')][,N := NULL]

dt <- merge(dt2,dt3,by='strat')
dt <- merge(dt,dt1,by = c('strat','hisp.strat'))

dt[,hisp.strat := NULL]
dt$bkg %<>% plyr::mapvalues(from = c('PR','D','O'),to = c('Puerto Rican','Dominican','Other'))
dt$sex %<>% plyr::mapvalues(from = c('M','F'),to = c('Male','Female'))
setorder(dt,strat,HispStrat,bkg,sex)

sink('./Output/Table1.tex')
kable(dt,booktab = T,caption = '\\label{svysamp}Characteristics of the simulated target population',
      col.names=c('','(N)','(N)','(N)','','','(N)','(\\%)'),escape=F,format = "latex",
      align = c('c','c','c','c','l','l','r','r')) %>%
  add_header_above(c('Stratum'=1,'Block Groups'=1,'Households'=1,'Household Type'=1,'Background'=1,'Sex'=1,'Individuals'=1,'45+ years old'=1), line = F) %>%
  collapse_rows(columns = 1:5,latex_hline = 'major',valign = 'top') %>%
  row_spec(0,align = 'c') %>%
  kable_styling(latex_options = "scale_down")
sink()

