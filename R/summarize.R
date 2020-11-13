library(plyr)
library(magrittr)
library(ggplot2)
library(haven)
library(data.table)
library(readr)
library(ggpubr)
library(kableExtra)
library(tidyr)
library(dplyr)

# Defining plot theme

convertPt <- function(x){
    return(x/(ggplot2::.pt*72.27/96))
}

theme_my <- function(base_size = 10, base_family = "sans") {
    
    base_size <- base_size/.pt # See https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#font-size why
    
    txt <- element_text(size = base_size, colour = "black", face = "plain",family = base_family)
    
    theme_bw(base_size = base_size, base_family = base_family) +
        theme(strip.text = txt,
              axis.ticks = element_line(colour = 'black', size = convertPt(1)),
              panel.border = element_rect(linetype = "solid", colour = "black", size=convertPt(1)),
              
              panel.grid = element_blank(),
              
              legend.key = element_blank(), 
              strip.background = element_blank(), 
              
              text = txt, 
              plot.title = txt, 
              
              axis.title = txt, 
              axis.text = txt, 
              
              legend.title = txt, 
              legend.text = txt) 
}

equal_breaks <- function(n = 3, s = 0.05, ...){
    function(x){
        # rescaling
        d <- s * diff(range(x)) / (1+2*s)
        seq(min(x)+d, max(x)-d, length=n)
    }
}

# Loading output from simulation study

load('./Output/bootstrap.RData')

df.m1.boot <- m1.list.boot %<>% rbindlist()
df.m2.boot <- m2.list.boot %<>% rbindlist()
df.m1.boot.time <- m1.list.mitime %<>% rbindlist()
df.m2.boot.time <- m2.list.mitime %<>% rbindlist()

load('./Output/MI.RData')

df.m1.mi <- m1.list.mi %<>% rbindlist()
df.m2.mi <- m2.list.mi %<>% rbindlist()
df.m1.mi.time <- m1.list.mitime %<>% rbindlist()
df.m2.mi.time <- m2.list.mitime %<>% rbindlist()

# Removing unecessary data
rm(m1.list.boot,m2.list.boot,m1.list.mitime,m2.list.mitime,m1.list.mi,m2.list.mi)

# Color of plots
Coeff = c('Intercept','Age','BMI','Log-Sodium','High Cholesterol','US Born','Sex: Female','Background: Puerto Rico','Background: Other')
coeffcol = c('#000000','#000000','#000000',"#000000",'#000000','#000000','#000000','#000000','#000000')
names(coeffcol) = Coeff

### Getting Population Coefficients
load('./RData/TargetPopulationData.RData')
pop <- pop[(v.num==1),-1]

pop.m1 = glm(hypertension~c_age+c_bmi+c_ln_na_true+high_chol+usborn+female+bkg_pr+bkg_o,data=pop,family=binomial())
pop.m2 = glm(sbp~c_age+c_bmi+c_ln_na_true+high_chol+usborn+female+bkg_pr+bkg_o,data=pop,family=gaussian())

pop.coeff = rbind(data.frame(Coeff = names(pop.m1$coefficients), Estimate = pop.m1$coefficients, Model = 'Hypertension'),
                  data.frame(Coeff = names(pop.m2$coefficients), Estimate = pop.m2$coefficients, Model = 'SBP'))
pop.coeff$Coeff %<>% mapvalues(from=unique(pop.coeff$Coeff),to=Coeff) %<>% as.character()
pop.coeff %<>% as.data.table()

rm(pop.m1,pop.m2,pop)

### Summary function
summ = function(df,method,model,raking,nboot=500,
                func.avg = 'mean',func.sd='sd',
                func.var.avg = 'mean',func.var.sd = 'sd',
                ci=95,digits=3,plots=F){
    df$Coeff %<>% mapvalues(from=unique(df$Coeff),to=Coeff) %<>% as.character()
    
    funclist = list(func.avg,func.sd)
    f.avg = match.fun(funclist[[1]]);f.sd = match.fun(funclist[[2]])
    
    funcvarlist = list(func.var.avg,func.var.sd)
    f.var.avg = match.fun(funcvarlist[[1]]);f.var.sd = match.fun(funcvarlist[[2]])
    
    # True Intake
    df.true.est = df[,.(True.Est = unique(True.Est)),by=c('Coeff','Sim')]
    df.true.se = df[,.(True.SE = unique(True.SE)),by=c('Coeff','Sim')]
    df.true.est.avg = df.true.est[,.(True.Est = f.avg(True.Est)),by=c('Coeff')]
    df.true.est.ese = df.true.est[,.(True.Est = f.sd(True.Est)),by=c('Coeff')]
    df.true.se.avg = df.true.se[,.(True.SE = f.avg(True.SE)),by=c('Coeff')]
    df.true.cov = merge(df.true.est,df.true.se,by=c('Coeff','Sim'),all.x=T)
    df.true.cov = merge(df.true.cov,subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),by='Coeff',all.x=T)
    df.true.cov[,True.Cov := 1*(Estimate>=(True.Est-qnorm(1-(1-ci/100)/2)*True.SE) & Estimate<=(True.Est+qnorm(1-(1-ci/100)/2)*True.SE))]
    df.true.cov = df.true.cov[,.(True.Cov = mean(True.Cov)),by=c('Coeff')]
    df.true.out = merge(df.true.est.avg,df.true.est.ese,by='Coeff',all.x=T)
    df.true.out = merge(df.true.out,df.true.se.avg,by='Coeff',all.x=T)
    df.true.out = merge(df.true.out,df.true.cov,by='Coeff',all.x=T);colnames(df.true.out) = c('Coeff','True.Est','True.ESE','True.SE','True.Cov')
    rm(df.true.est,df.true.se,df.true.est.avg,df.true.est.ese,df.true.se.avg,df.true.cov)
    
    # Naive Intake
    df.avg.est = df[,.(Avg.Est = unique(Avg.Est)),by=c('Coeff','Sim')]
    df.avg.se = df[,.(Avg.SE = unique(Avg.SE)),by=c('Coeff','Sim')]
    df.avg.est.avg = df.avg.est[,.(Avg.Est = f.avg(Avg.Est)),by=c('Coeff')]
    df.avg.est.ese = df.avg.est[,.(Avg.Est = f.sd(Avg.Est)),by=c('Coeff')]
    df.avg.se.avg = df.avg.se[,.(Avg.SE = f.avg(Avg.SE)),by=c('Coeff')]
    df.avg.cov = merge(df.avg.est,df.avg.se,by=c('Coeff','Sim'),all.x=T)
    df.avg.cov = merge(df.avg.cov,subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),by='Coeff',all.x=T)
    df.avg.cov[,Avg.Cov := 1*(Estimate>=(Avg.Est-qnorm(1-(1-ci/100)/2)*Avg.SE) & Estimate<=(Avg.Est+qnorm(1-(1-ci/100)/2)*Avg.SE))]
    df.avg.cov = df.avg.cov[,.(Avg.Cov = mean(Avg.Cov)),by=c('Coeff')]
    df.avg.out = merge(df.avg.est.avg,df.avg.est.ese,by='Coeff',all.x=T)
    df.avg.out = merge(df.avg.out,df.avg.se.avg,by='Coeff',all.x=T)
    df.avg.out = merge(df.avg.out,df.avg.cov,by='Coeff',all.x=T);colnames(df.avg.out) = c('Coeff','Naive.Est','Naive.ESE','Naive.SE','Naive.Cov')
    rm(df.avg.est,df.avg.se,df.avg.est.avg,df.avg.est.ese,df.avg.se.avg,df.avg.cov)
    
    if(raking==F){
        # Calibrated intake
        df.calib.est = df[,.(Calib.Est = unique(Calib.Est)),by=c('Coeff','Sim')]
        df.calib.se = df[,.(Calib.SE = unique(Calib.SE)),by=c('Coeff','Sim')]
        df.calib.est.avg = df.calib.est[,.(Calib.Est = f.avg(Calib.Est)),by=c('Coeff')]
        df.calib.est.ese = df.calib.est[,.(Calib.Est = f.sd(Calib.Est)),by=c('Coeff')]
        df.calib.se.avg = df.calib.se[!is.na(Calib.SE),.(Calib.SE = f.avg(Calib.SE)),by=c('Coeff')]
        df.calib.cov = merge(df.calib.est,df.calib.se,by=c('Coeff','Sim'),all.x=T)
        df.calib.cov = merge(df.calib.cov,subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),by='Coeff',all.x=T)
        df.calib.cov[,Calib.Cov := 1*(Estimate>=(Calib.Est-qnorm(1-(1-ci/100)/2)*Calib.SE) & Estimate<=(Calib.Est+qnorm(1-(1-ci/100)/2)*Calib.SE))]
        df.calib.cov = df.calib.cov[!is.na(Calib.SE),.(Calib.Cov = mean(Calib.Cov)),by=c('Coeff')]
        df.calib.out = merge(df.calib.est.avg,df.calib.est.ese,by='Coeff',all.x=T)
        df.calib.out = merge(df.calib.out,df.calib.se.avg,by='Coeff',all.x=T)
        df.calib.out = merge(df.calib.out,df.calib.cov,by='Coeff',all.x=T);colnames(df.calib.out) = c('Coeff','Calib.Est','Calib.ESE','Calib.SE','Calib.Cov')
        rm(df.calib.se,df.calib.est.avg,df.calib.est.ese,df.calib.se.avg,df.calib.cov)
        
        # Bootstrap/MI correction for SE
        df.var.within = df[,.(withinVar = f.var.avg(SE^2)),by=c('Coeff','Sim')]
        df.var.between = df[,.(betweenSE = f.var.sd(Est)),by=c('Coeff','Sim')]
        df.var = merge(x=df.var.within,y=df.var.between,by=c('Coeff','Sim'))
        
        # Applying Rubin's formula (equation 12.1 from Handbook of Missing Data)
        df.var[,Corrected.SE := sqrt(withinVar + betweenSE^2)]
        
        df.var.mean = df.var[,.(Corrected.SE = f.avg(Corrected.SE)),by=c('Coeff')]
        df.var.cov = merge(df,df.var,by = c('Sim','Coeff'),all.x = T)
        df.var.cov = df.var.cov[,.(t1 = quantile((Est-Calib.Est)/(SE+betweenSE),probs=(1-ci/100)/2),
                                   t2 = quantile((Est-Calib.Est)/(SE+betweenSE),probs=1-(1-ci/100)/2),
                                   q1 = quantile(Est,probs=(1-ci/100)/2),
                                   q2 = quantile(Est,probs=1-(1-ci/100)/2)),by = c('Sim','Coeff')]
        df.var.cov = merge(df.var.cov,subset(df.var,select=c(Coeff,Sim,Corrected.SE)),by=c('Coeff','Sim'),all.x=T)
        df.var.cov = merge(df.var.cov,df.calib.est,by=c('Coeff','Sim'),all.x=T)
        df.var.cov = merge(df.var.cov,subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),by='Coeff',all.x=T)
        df.var.cov[,c('Corrected.Cov.Normal','Corrected.Cov.Effron','Corrected.Cov.Pivot','Corrected.Cov.Quantile') :=
                       .(1*(Estimate>=(Calib.Est-qnorm(1-(1-ci/100)/2)*Corrected.SE) & Estimate<=(Calib.Est+qnorm(1-(1-ci/100)/2)*Corrected.SE)),
                         1*(Estimate>=(Calib.Est-t2*Corrected.SE) & Estimate<=(Calib.Est-t1*Corrected.SE)),
                         1*(Estimate>=(2*Calib.Est-q2) & Estimate<=(2*Calib.Est-q1)),
                         1*(Estimate>=q1 & Estimate<=q2))]
        df.var.cov = df.var.cov[,.(Corrected.Cov.Normal = mean(Corrected.Cov.Normal),
                                   Corrected.Cov.Effron = mean(Corrected.Cov.Effron),
                                   Corrected.Cov.Pivot = mean(Corrected.Cov.Pivot),
                                   Corrected.Cov.Quantile = mean(Corrected.Cov.Quantile)),by='Coeff']
        
        #Output table
        df.out = merge(df.true.out,df.calib.out,by='Coeff',all.x=T)
        df.out = merge(df.out,df.avg.out,by='Coeff',all.x=T)
        df.out = merge(df.out,df.var.mean,by='Coeff',all.x=T)
        df.out = merge(df.out,df.var.cov,by='Coeff',all.x=T)
        
        df.out = df.out[match(df[1:nrow(df.out),]$Coeff,df.out$Coeff),
                        c('Coeff','True.Est','True.ESE','True.SE','True.Cov', 
                          'Naive.Est','Naive.ESE','Naive.SE','Naive.Cov',
                          'Calib.Est','Calib.ESE','Calib.SE','Calib.Cov',
                          'Corrected.SE','Corrected.Cov.Normal','Corrected.Cov.Effron','Corrected.Cov.Pivot','Corrected.Cov.Quantile')]
        df.out[,2:ncol(df.out)] %<>% round(digits)
        return(df.out)
    }
    if(raking==T){
        # Raking models
        df.raking.est = df[,.(Raking.Est = unique(Raking.Est)),by=c('Coeff','Sim')]
        df.raking.se = df[,.(Raking.SE = unique(Raking.SE)),by=c('Coeff','Sim')]
        df.raking.est.avg = df.raking.est[,.(Raking.Est = f.avg(Raking.Est)),by=c('Coeff')]
        df.raking.est.ese = df.raking.est[,.(Raking.Est = f.sd(Raking.Est)),by=c('Coeff')]
        df.raking.se.avg = df.raking.se[,.(Raking.SE = f.avg(na.omit(Raking.SE))),by=c('Coeff')]
        df.raking.cov = merge(df.raking.est,df.raking.se,by=c('Coeff','Sim'),all.x=T)
        df.raking.cov = merge(df.raking.cov,subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),by='Coeff',all.x=T)
        df.raking.cov[,Raking.Cov := 1*(Estimate>=(Raking.Est-qnorm(1-(1-ci/100)/2)*Raking.SE) & Estimate<=(Raking.Est+qnorm(1-(1-ci/100)/2)*Raking.SE))]
        df.raking.cov = df.raking.cov[,.(Raking.Cov = mean(Raking.Cov,na.rm = T)),by=c('Coeff')] #Removing NaN from a single case 
        df.raking.out = merge(df.raking.est.avg,df.raking.est.ese,by='Coeff',all.x=T)
        df.raking.out = merge(df.raking.out,df.raking.se.avg,by='Coeff',all.x=T)
        df.raking.out = merge(df.raking.out,df.raking.cov,by='Coeff',all.x=T);colnames(df.raking.out) = c('Coeff','Raking.Est','Raking.ESE','Raking.SE','Raking.Cov')
        rm(df.raking.est,df.raking.se,df.raking.est.avg,df.raking.est.ese,df.raking.se.avg,df.raking.cov)
        
        # This needs to be revised if we use raking estimators
        # # Biomarker Intake
        # df.bio.est = df[,.(Bio.Moments.Est = unique(Bio.Moments.Est)),by=c('Coeff','Sim')]
        # df.bio.se = df[,.(Bio.Moments.SE = unique(Bio.Moments.SE)),by=c('Coeff','Sim')]
        # df.bio.est.avg = df.bio.est[,.(Bio.Moments.Est = f.avg(Bio.Moments.Est)),by=c('Coeff')]
        # df.bio.est.ese = df.bio.est[,.(Bio.Moments.Est = f.sd(Bio.Moments.Est)),by=c('Coeff')]
        # df.bio.se.avg = df.bio.se[,.(Bio.Moments.SE = f.avg(Bio.Moments.SE)),by=c('Coeff')]
        # df.bio.cov = merge(df.bio.est,df.bio.se,by=c('Coeff','Sim'),all.x=T)
        # df.bio.cov = merge(df.bio.cov,subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),by='Coeff',all.x=T)
        # df.bio.cov[,Bio.Moments.Cov := 1*(Estimate>=(Bio.Moments.Est-qnorm(1-(1-ci/100)/2)*Bio.Moments.SE) & Estimate<=(Bio.Moments.Est+qnorm(1-(1-ci/100)/2)*Bio.Moments.SE))]
        # df.bio.cov = df.bio.cov[,.(Bio.Moments.Cov = mean(Bio.Moments.Cov)),by=c('Coeff')]
        # df.bio.out = merge(df.bio.est.avg,df.bio.est.ese,by='Coeff',all.x=T)
        # df.bio.out = merge(df.bio.out,df.bio.se.avg,by='Coeff',all.x=T)
        # df.bio.out = merge(df.bio.out,df.bio.cov,by='Coeff',all.x=T);colnames(df.bio.out) = c('Coeff','Bio.Moments.Est','Bio.Moments.ESE','Bio.Moments.SE','Bio.Moments.Cov')
        # rm(df.bio.est,df.bio.se,df.bio.est.avg,df.bio.est.ese,df.bio.se.avg,df.bio.cov)
        # 
        # # Bootstrap correction for SE
        # df.var.within = df[,.(withinVar = f.var.avg(Bio.Moments.SE.Boot^2)),by=c('Coeff','Sim')]
        # df.var.between = df[,.(betweenSE = f.var.sd(Bio.Moments.Est.Boot)),by=c('Coeff','Sim')]
        # df.var = merge(x=df.var.within,y=df.var.between,by=c('Coeff','Sim'))
        # df.var[,Bio.Moments.Corrected.SE := betweenSE]
        # 
        # df.var.mean = df.var[,.(Bio.Moments.Corrected.SE = f.avg(Bio.Moments.Corrected.SE)),by=c('Coeff')]
        # df.var.cov = merge(df,df.var,by = c('Sim','Coeff'),all.x = T)
        # df.var.cov = merge(df.var.cov,subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),by='Coeff',all.x=T)
        # df.var.cov[,Bio.Moments.Corrected.Cov:= 1*(Estimate>=(Bio.Moments.Est-qnorm(1-(1-ci/100)/2)*Bio.Moments.Corrected.SE) & Estimate<=(Bio.Moments.Est+qnorm(1-(1-ci/100)/2)*Bio.Moments.Corrected.SE))]
        # df.var.cov = df.var.cov[,.(Bio.Moments.Corrected.Cov = mean(Bio.Moments.Corrected.Cov)),by='Coeff']
        # 
        # # Raking Intake
        # df.raking.est = df[,.(Raking.Moments.Est = unique(Raking.Moments.Est)),by=c('Coeff','Sim')]
        # df.raking.se = df[,.(Raking.Moments.SE = unique(Raking.Moments.SE)),by=c('Coeff','Sim')] 
        # df.raking.est.avg = df.raking.est[,.(Raking.Moments.Est = f.avg(Raking.Moments.Est)),by=c('Coeff')]
        # df.raking.est.ese = df.raking.est[,.(Raking.Moments.Est = f.sd(Raking.Moments.Est)),by=c('Coeff')]
        # df.raking.se.avg = df.raking.se[,.(Raking.Moments.SE = mean(Raking.Moments.SE,na.rm=T)),by=c('Coeff')] #Fix is to f.avg because there was an NA in the data
        # df.raking.cov = merge(df.raking.est,df.raking.se,by=c('Coeff','Sim'),all.x=T)
        # df.raking.cov = merge(df.raking.cov,subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),by='Coeff',all.x=T)
        # df.raking.cov[,Raking.Moments.Cov := 1*(Estimate>=(Raking.Moments.Est-qnorm(1-(1-ci/100)/2)*Raking.Moments.SE) & Estimate<=(Raking.Moments.Est+qnorm(1-(1-ci/100)/2)*Raking.Moments.SE))]
        # df.raking.cov = df.raking.cov[,.(Raking.Moments.Cov = mean(Raking.Moments.Cov,na.rm=T)),by=c('Coeff')] #Fix is to f.avg because there was an NA in the data
        # df.raking.out = merge(df.raking.est.avg,df.raking.est.ese,by='Coeff',all.x=T)
        # df.raking.out = merge(df.raking.out,df.raking.se.avg,by='Coeff',all.x=T)
        # df.raking.out = merge(df.raking.out,df.raking.cov,by='Coeff',all.x=T);colnames(df.raking.out) = c('Coeff','Raking.Moments.Est','Raking.Moments.ESE','Raking.Moments.SE','Raking.Moments.Cov')
        # rm(df.raking.est,df.raking.se,df.raking.est.avg,df.raking.est.ese,df.raking.se.avg,df.raking.cov)
        
        #Output table
        df.out = merge(df.true.out,df.avg.out,by='Coeff',all.x=T)
        df.out = merge(df.out,df.raking.out,by='Coeff',all.x=T)
        
        
        df.out = df.out[match(df[1:nrow(df.out),]$Coeff,df.out$Coeff),
                        c('Coeff','True.Est','True.ESE','True.SE','True.Cov', 
                          'Naive.Est','Naive.ESE','Naive.SE','Naive.Cov',
                          'Raking.Est','Raking.ESE','Raking.SE','Raking.Cov')]
        df.out[,2:ncol(df.out)] %<>% round(digits)
        return(df.out)
    }
    
    if(plots){
        # This needs to be revised
        # fig.CalibEst =  ggplot(data=df.calib.est)+geom_histogram(position='identity',aes(x=Calib.Est,y=..density..))+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+theme_bw()+
        #     geom_vline(data=subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),aes(xintercept=Estimate))+
        #     xlab('Estimated coefficients from outcome-calibrated model')
        # ggsave(filename = paste('./Misc/Hist',model,paste0(method,'Est','_',func.avg,'_','.eps'),sep='_'),plot = fig.CalibEst)
        # 
        # fig.AvgEst =  ggplot(data=df.avg.est)+geom_histogram(position='identity',aes(x=Avg.Est,y=..density..))+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+theme_bw()+
        #     geom_vline(data=subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),aes(xintercept=Estimate))+
        #     xlab('Estimated coefficients from outcome-average model')
        # ggsave(filename = paste('./Misc/Hist',model,'AvgEst.eps',sep='_'),plot = fig.AvgEst)
        # 
        # fig.TrueEst =  ggplot(data=df.true.est)+geom_histogram(position='identity',aes(x=True.Est,y=..density..))+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+theme_bw()+
        #     geom_vline(data=subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),aes(xintercept=Estimate))+
        #     xlab('Estimated coefficients from outcome-true model')
        # ggsave(filename = paste('./Misc/Hist',model,'TrueEst.eps',sep='_'),plot = fig.TrueEst)
        # 
        # fig.BioEst =  ggplot(data=df.bio.est)+geom_histogram(position='identity',aes(x=Bio.Est,y=..density..))+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+theme_bw()+
        #     geom_vline(data=subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),aes(xintercept=Estimate))+
        #     xlab('Estimated coefficients from outcome-biomarker (phase 2) model')
        # ggsave(filename = paste('./Misc/Hist',model,'BioEst.eps',sep='_'),plot = fig.BioEst)
        # 
        # fig.vif1 = ggplot(data=df[!(Coeff=='Intercept'),])+geom_histogram(position='identity',aes(x=VIF,y=..density..))+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+theme_bw()+geom_vline(aes(xintercept=5),color='red')+
        #     xlab('Variance Inflation Factor')+ylab('Density')
        # ggsave(filename = paste('./Misc/VIF',model,'Hist.eps',sep='_'),plot = fig.vif1)
        # 
        # fig.vif2 = ggplot(data=df[!(Coeff=='Intercept'),],aes(x=VIF,y=Est))+
        #     geom_smooth(se=T)+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+geom_vline(aes(xintercept=5),color='red')+
        #     xlab('Variance Inflation Factor')+ylab('Parameter estimates from Outcome-calibrated model')
        # ggsave(filename = paste('./Misc/VIF',model,'Scatter.eps',sep='_'),plot = fig.vif2,dpi='screen')
    }
}

# Tables - Average/SD
m1.boot = summ(df=df.m1.boot,method='Bootstrap',model='Hypertension',raking=F,digits = 3,
               func.avg = 'mean',func.sd = 'sd',
               func.var.avg = 'median',func.var.sd = 'mad')
m1.mi = summ(df=df.m1.mi,method='MI',model='Hypertension',raking=F,digits = 3,
             func.avg = 'mean',func.sd = 'sd',
             func.var.avg = 'median',func.var.sd = 'mad')

m2.boot = summ(df=df.m2.boot,method='Bootstrap',model='SBP',raking=F,digits = 3,
               func.avg = 'mean',func.sd = 'sd',
               func.var.avg = 'median',func.var.sd = 'mad')
m2.mi = summ(df=df.m2.mi,method='MI',model='SBP',raking=F,digits = 3,
             func.avg = 'mean',func.sd = 'sd',
             func.var.avg = 'median',func.var.sd = 'mad')

### Organizing the output
key = Coeff

m1.boot = merge(subset(pop.coeff,Model=='Hypertension',select=c(Coeff,Estimate)),m1.boot,by='Coeff',all.x=T);m1.boot.out <- m1.boot[order(match(Coeff,key)),];m1.boot.out[,Estimate := round(Estimate,3)]
m1.mi = merge(subset(pop.coeff,Model=='Hypertension',select=c(Coeff,Estimate)),m1.mi,by='Coeff',all.x=T);m1.mi.out <- m1.mi[order(match(Coeff,key)),];m1.mi.out[,Estimate := round(Estimate,3)]

m2.boot = merge(subset(pop.coeff,Model=='SBP',select=c(Coeff,Estimate)),m2.boot,by='Coeff',all.x=T);m2.boot.out <- m2.boot[order(match(Coeff,key)),];m2.boot.out[,Estimate := round(Estimate,3)]
m2.mi = merge(subset(pop.coeff,Model=='SBP',select=c(Coeff,Estimate)),m2.mi,by='Coeff',all.x=T);m2.mi.out <- m2.mi[order(match(Coeff,key)),];m2.mi.out[,Estimate := round(Estimate,3)]

# Saving 
write.csv(m1.boot.out,file=paste0('./Output/summarize_bootstrap_logistic.csv'))
write.csv(m2.boot.out,file=paste0('./Output/summarize_bootstrap_linear.csv'))
write.csv(m1.mi.out,file=paste0('./Output/summarize_MI_logistic.csv'))
write.csv(m2.mi.out,file=paste0('./Output/summarize_MI_linear.csv'))

### Organizing the data to plot

m1.plot <- rbindlist(list(m1.boot.out[,c('Coeff','Estimate','True.Est','True.ESE','True.SE','True.Cov')][,Method := 'True (Phase1, Unobservable)'],
                          m1.boot.out[,c('Coeff','Estimate','Naive.Est','Naive.ESE','Naive.SE','Naive.Cov')][,Method := 'Naive 2-day Mean (Phase 1)'],
                          m1.boot.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Calib.SE','Calib.Cov')][,Method := 'Naive Calibrated (Phase 1)'],
                          m1.boot.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Corrected.SE','Corrected.Cov.Normal')][,Method := 'Naive Calibrated w/ Bootstrap (Phase 1)'],
                          m1.mi.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Corrected.SE','Corrected.Cov.Normal')][,Method := 'Naive Calibrated w/ MI (Phase 1)']),use.names=FALSE)
setnames(m1.plot,c('Coeff','True','Estimate','ESE','SE','Cov','Method'))

m1.plot$Coeff %<>% factor(levels = unique(.))
m1.plot$Method %<>% as.factor() %<>% factor(levels = unique(.))
m1.plot[,Lbl := min(Estimate-SE,Estimate-ESE)-max(SE,ESE),by='Coeff']

m2.plot <- rbindlist(list(m2.boot.out[,c('Coeff','Estimate','True.Est','True.ESE','True.SE','True.Cov')][,Method := 'True (Phase1, Unobservable)'],
                          m2.boot.out[,c('Coeff','Estimate','Naive.Est','Naive.ESE','Naive.SE','Naive.Cov')][,Method := 'Naive 2-day Mean (Phase 1)'],
                          m2.boot.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Calib.SE','Calib.Cov')][,Method := 'Naive Calibrated (Phase 1)'],
                          m2.boot.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Corrected.SE','Corrected.Cov.Normal')][,Method := 'Naive Calibrated w/ Bootstrap (Phase 1)'],
                          m2.mi.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Corrected.SE','Corrected.Cov.Normal')][,Method := 'Naive Calibrated w/ MI (Phase 1)']),use.names=FALSE)
setnames(m2.plot,c('Coeff','True','Estimate','ESE','SE','Cov','Method'))

m2.plot$Coeff %<>% factor(levels = unique(.))
m2.plot$Method %<>% as.factor() %<>% factor(levels = unique(.))
m2.plot[,Lbl := min(Estimate-SE,Estimate-ESE)-max(SE,ESE),by='Coeff']

### Subsetting to the relevant methods
m1.plot.sub <- m1.plot[Method %in% c('True (Phase1, Unobservable)','Naive 2-day Mean (Phase 1)','Naive Calibrated w/ Bootstrap (Phase 1)','Naive Calibrated w/ MI (Phase 1)'),]
m1.plot.sub$Method %<>% mapvalues(from = unique(.),to = c('True\n(Unobservable)','Naive','Bootstrap','MI'))
m1.plot.sub$Coeff %<>% factor(levels = c('Intercept','Log-Sodium',as.character(unique(m1.plot.sub$Coeff[!m1.plot.sub$Coeff%in%c('Intercept','Log-Sodium')]))))

m2.plot.sub <- m2.plot[Method %in% c('True (Phase1, Unobservable)','Naive 2-day Mean (Phase 1)','Naive Calibrated w/ Bootstrap (Phase 1)','Naive Calibrated w/ MI (Phase 1)'),]
m2.plot.sub$Method %<>% mapvalues(from = unique(.),to = c('True\n(Unobservable)','Naive','Bootstrap','MI'))
m2.plot.sub$Coeff %<>% factor(levels = c('Intercept','Log-Sodium',as.character(unique(m2.plot.sub$Coeff[!m2.plot.sub$Coeff%in%c('Intercept','Log-Sodium')]))))

#####################################################################
### Plotting Error bars
#####################################################################

sizeLetter.erroBars <- 30

fig.m1.ls <- lapply(sort(unique(m1.plot.sub[!Method=="True\n(Unobservable)",]$Coeff)),function(coef){
    s <- 0.25
    
    subdt <- m1.plot.sub[!Method=="True\n(Unobservable)" & Coeff==coef,]
    
    range <- c(subdt[,min(Estimate-ESE)],subdt[,max(Estimate+ESE)])
    range <- c(range[1]*ifelse(range[1]<0,1+s,1-s),range[2]*ifelse(range[2]<0,1-s,1+s))
    range <- c(plyr::round_any(range[1],ifelse(abs(diff(range))<0.1,0.01,0.1),floor),
               plyr::round_any(range[2],ifelse(abs(diff(range))<0.1,0.01,0.1),ceiling))
    
    subfig <- ggplot(data = subdt,aes(x = Method, y = Estimate,shape = Method)) +
        geom_hline(data = subdt, aes(yintercept = True),size=convertPt(1))+
        geom_pointrange(aes(ymin = Estimate-SE,ymax = Estimate+SE,color = Coeff),size=convertPt(1)) +
        geom_point(aes(y = Estimate+ESE,color = Coeff),shape = 5)+
        geom_point(aes(y = Estimate-ESE,color = Coeff),shape = 5)+
        theme_my(base_size = sizeLetter.erroBars)+
        scale_color_manual(values = coeffcol)+
        guides(shape = FALSE, color = FALSE)+
        theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.text = NULL,legend.title = element_blank())+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01,decimal.mark = '.'),limits = range,breaks = seq(range[1],range[2],length.out = 4))+
        theme(axis.text.x=element_text(angle=20,hjust=1),axis.title.x = element_blank(),plot.margin = unit(c(15, 5.5, 5.5, 15), "points"))
    
    return(subfig)
})

fig.m1 <- ggpubr::ggarrange(plotlist = fig.m1.ls,ncol = 3,nrow = 3,labels = as.list(paste0(LETTERS[1:9],')')),
                            hjust = 0,vjust = 1.25,font.label = list(size = sizeLetter.erroBars/.pt,face = 'plain'))

fig.m2.ls <- lapply(sort(unique(m2.plot.sub[!Method=="True\n(Unobservable)",]$Coeff)),function(coef){
    
    s <- ifelse(coef == 'Intercept',0.01,ifelse(coef == 'Background: Puerto Rico',0.4,0.25))
    
    subdt <- m2.plot.sub[!Method=="True\n(Unobservable)" & Coeff==coef,]
    
    range <- c(subdt[,min(Estimate-ESE)],subdt[,max(Estimate+ESE)])
    range <- c(range[1]*ifelse(range[1]<0,1+s,1-s),range[2]*ifelse(range[2]<0,1-s,1+s))
    range <- c(plyr::round_any(range[1],ifelse(abs(diff(range))<0.1,0.01,0.1),floor),
               plyr::round_any(range[2],ifelse(abs(diff(range))<0.1,0.01,0.1),ceiling))
    
    subfig <- ggplot(data = subdt,aes(x = Method, y = Estimate,shape = Method)) +
        geom_hline(data = subdt, aes(yintercept = True),size=convertPt(1))+
        geom_pointrange(aes(ymin = Estimate-SE,ymax = Estimate+SE,color = Coeff),size=convertPt(1)) +
        geom_point(aes(y = Estimate+ESE,color = Coeff),shape = 5)+
        geom_point(aes(y = Estimate-ESE,color = Coeff),shape = 5)+
        theme_my(base_size = sizeLetter.erroBars)+
        scale_color_manual(values = coeffcol)+
        guides(shape = FALSE, color = FALSE)+
        theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.text = NULL,legend.title = element_blank())+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01,decimal.mark = '.'),limits = range,breaks = seq(range[1],range[2],length.out = 4))+
        theme(axis.text.x=element_text(angle=20,hjust=1),plot.margin = unit(c(15, 5.5, 5.5, 15), "points"),axis.title.x = element_blank())
    
    return(subfig)
})

fig.m2 <- ggpubr::ggarrange(plotlist = fig.m2.ls,ncol = 3,nrow = 3,labels = as.list(paste0(LETTERS[10:18],')')),
                            hjust = 0,vjust = 1.25,font.label = list(size = sizeLetter.erroBars/.pt,face = 'plain'))

### Put them together

fig4 <- ggpubr::ggarrange(fig.m1,fig.m2,ncol = 1,nrow = 2)
fig4 <- annotate_figure(fig4,bottom = text_grob("Regression Model",just = 'center', size = sizeLetter.erroBars/.pt))
ggsave(fig4,filename = './Output/Figure4.pdf',dpi = 'retina',height = 9,width = 7)

### Saving the plots

for(i in 1:9){
    subfig <- ggpubr::ggarrange(fig.m1.ls[[i]],ncol = 1,nrow = 1,labels = as.list(paste0(LETTERS[i],')')),
                                hjust = 0,vjust = 1.25,font.label = list(size = sizeLetter.erroBars/.pt,face = 'plain'))
    ggsave(subfig,filename = paste0('./Output/Figure4',LETTERS[i],'.eps'),dpi = 'retina',height = 3,width = 3,device = grDevices::cairo_ps,fallback_resolution = 300)
}
for(i in 10:18){
    subfig <- ggpubr::ggarrange(fig.m2.ls[[i-9]],ncol = 1,nrow = 1,labels = as.list(paste0(LETTERS[i],')')),
                                hjust = 0,vjust = 1.25,font.label = list(size = sizeLetter.erroBars/.pt,face = 'plain'))
    ggsave(subfig,filename = paste0('./Output/Figure4',LETTERS[i],'.eps'),dpi = 'retina',height = 3,width = 3,device = grDevices::cairo_ps,fallback_resolution = 300)
}

# Computing time

df.time <- rbindlist(list(df.m1.boot.time[,Method := 'Resampling-based MI'][,Model := 'Logistic'],
                          df.m2.boot.time[,Method := 'Resampling-based MI'][,Model := 'Linear'],
                          df.m1.mi.time[,Method := 'Parametric MI'][,Model := 'Logistic'],
                          df.m2.mi.time[,Method := 'Parametric MI'][,Model := 'Linear']))
df.time[,Total.Time := Time.Calibration+Time.Model+Time.Resampling]

fig.time <- ggplot(df.time,aes(x = Model,y = log10(Total.Time)))+
    geom_boxplot(aes(fill = Method))+
    theme_bw()+
    ylab('Log-10 computing time (in seconds) across 1000 simulations')+
    geom_hline(yintercept = c(log10(15),log10(30)),linetype = 'dashed',color = 'grey',alpha = 0.75)+
    scale_fill_brewer(palette = "Set1")+
    annotate('text',x = c(2.5,2.5), y = 0.975*c(log10(15),log10(30)),label = c('15 sec.','30 sec.'))+
    theme(legend.position = c(0.25,0.90),axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),legend.background = element_rect(fill=alpha('white', 0)))

ggsave(fig.time,filename = './Output/SuppFigure2.eps',
       dpi = 'retina',height = 11/2,width = 8.5/1.25,device = grDevices::cairo_ps,fallback_resolution = 300)

### Saving the output

save(m1.boot.out,m2.boot.out,m1.mi.out,m2.mi.out,file = './Output/summarize.RData')

### Now, computing tables

t.m1.boot.out <- data.table(Parameter = names(m1.boot.out),data.table::transpose(m1.boot.out))
setnames(t.m1.boot.out,as.character(t.m1.boot.out[1,]))
t.m1.boot.out <- t.m1.boot.out[-1,]
t.m1.boot.out <- t.m1.boot.out[match(c('Estimate','True.Est','Naive.Est','Calib.Est',
                                       'True.ESE','Naive.ESE','Calib.ESE',
                                       'True.SE','Naive.SE','Calib.SE','Corrected.SE',
                                       'True.Cov','Naive.Cov','Calib.Cov','Corrected.Cov.Normal'),Coeff),]

t.m2.boot.out <- data.table(Parameter = names(m2.boot.out),data.table::transpose(m2.boot.out))
setnames(t.m2.boot.out,as.character(t.m2.boot.out[1,]))
t.m2.boot.out <- t.m2.boot.out[-1,]
t.m2.boot.out <- t.m2.boot.out[match(c('Estimate','True.Est','Naive.Est','Calib.Est',
                                       'True.ESE','Naive.ESE','Calib.ESE',
                                       'True.SE','Naive.SE','Calib.SE','Corrected.SE',
                                       'True.Cov','Naive.Cov','Calib.Cov','Corrected.Cov.Normal'),Coeff),]

t.m1.mi.out <- data.table(Parameter = names(m1.mi.out),data.table::transpose(m1.mi.out))
setnames(t.m1.mi.out,as.character(t.m1.mi.out[1,]))
t.m1.mi.out <- t.m1.mi.out[-1,]
t.m1.mi.out <- t.m1.mi.out[match(c('Estimate','True.Est','Naive.Est','Calib.Est',
                                   'True.ESE','Naive.ESE','Calib.ESE',
                                   'True.SE','Naive.SE','Calib.SE','Corrected.SE',
                                   'True.Cov','Naive.Cov','Calib.Cov','Corrected.Cov.Normal'),Coeff),]

t.m2.mi.out <- data.table(Parameter = names(m2.mi.out),data.table::transpose(m2.mi.out))
setnames(t.m2.mi.out,as.character(t.m2.mi.out[1,]))
t.m2.mi.out <- t.m2.mi.out[-1,]
t.m2.mi.out <- t.m2.mi.out[match(c('Estimate','True.Est','Naive.Est','Calib.Est',
                                   'True.ESE','Naive.ESE','Calib.ESE',
                                   'True.SE','Naive.SE','Calib.SE','Corrected.SE',
                                   'True.Cov','Naive.Cov','Calib.Cov','Corrected.Cov.Normal'),Coeff),]

t.m1.boot.out$Coeff %<>% mapvalues(from = c('Estimate','True.Est','Naive.Est','Calib.Est',
                                            'True.ESE','Naive.ESE','Calib.ESE',
                                            'True.SE','Naive.SE','Calib.SE','Corrected.SE',
                                            'True.Cov','Naive.Cov','Calib.Cov','Corrected.Cov.Normal'),
                                   to = c('','True','Naive','Calibrated',
                                          'True','Naive','Calibrated',
                                          'True','Naive','Calibrated','Corrected',
                                          'True','Naive','Calibrated','Corrected'))
t.m2.boot.out$Coeff %<>% mapvalues(from = c('Estimate','True.Est','Naive.Est','Calib.Est',
                                            'True.ESE','Naive.ESE','Calib.ESE',
                                            'True.SE','Naive.SE','Calib.SE','Corrected.SE',
                                            'True.Cov','Naive.Cov','Calib.Cov','Corrected.Cov.Normal'),
                                   to = c('','True','Naive','Calibrated',
                                          'True','Naive','Calibrated',
                                          'True','Naive','Calibrated','Corrected',
                                          'True','Naive','Calibrated','Corrected'))
t.m1.mi.out$Coeff %<>% mapvalues(from = c('Estimate','True.Est','Naive.Est','Calib.Est',
                                          'True.ESE','Naive.ESE','Calib.ESE',
                                          'True.SE','Naive.SE','Calib.SE','Corrected.SE',
                                          'True.Cov','Naive.Cov','Calib.Cov','Corrected.Cov.Normal'),
                                 to = c('','True','Naive','Calibrated',
                                        'True','Naive','Calibrated',
                                        'True','Naive','Calibrated','Corrected',
                                        'True','Naive','Calibrated','Corrected'))
t.m2.mi.out$Coeff %<>% mapvalues(from = c('Estimate','True.Est','Naive.Est','Calib.Est',
                                          'True.ESE','Naive.ESE','Calib.ESE',
                                          'True.SE','Naive.SE','Calib.SE','Corrected.SE',
                                          'True.Cov','Naive.Cov','Calib.Cov','Corrected.Cov.Normal'),
                                 to = c('','True','Naive','Calibrated',
                                        'True','Naive','Calibrated',
                                        'True','Naive','Calibrated','Corrected',
                                        'True','Naive','Calibrated','Corrected'))
t.m1.boot.out$Type <- c('',rep('Estimate',3),rep('ESE',3),rep('SE',4),rep('Coverage',4))
t.m2.boot.out$Type <- c('',rep('Estimate',3),rep('ESE',3),rep('SE',4),rep('Coverage',4))
t.m1.mi.out$Type <- c('',rep('Estimate',3),rep('ESE',3),rep('SE',4),rep('Coverage',4))
t.m2.mi.out$Type <- c('',rep('Estimate',3),rep('ESE',3),rep('SE',4),rep('Coverage',4))

t.m1.boot.out<-t.m1.boot.out[,c('Type', "Coeff","Intercept","Age","BMI","Log-Sodium","High Cholesterol","US Born","Sex: Female","Background: Puerto Rico","Background: Other")]
t.m2.boot.out<-t.m2.boot.out[,c('Type', "Coeff","Intercept","Age","BMI","Log-Sodium","High Cholesterol","US Born","Sex: Female","Background: Puerto Rico","Background: Other")]
t.m1.mi.out<-t.m1.mi.out[,c('Type', "Coeff","Intercept","Age","BMI","Log-Sodium","High Cholesterol","US Born","Sex: Female","Background: Puerto Rico","Background: Other")]
t.m2.mi.out<-t.m2.mi.out[,c('Type', "Coeff","Intercept","Age","BMI","Log-Sodium","High Cholesterol","US Born","Sex: Female","Background: Puerto Rico","Background: Other")]

t.m1.boot.out[,c(colnames(t.m1.boot.out)[3:ncol(t.m1.boot.out)]) := lapply(.SD, as.numeric), .SDcols = colnames(t.m1.boot.out)[3:ncol(t.m1.boot.out)]]
t.m2.boot.out[,c(colnames(t.m2.boot.out)[3:ncol(t.m2.boot.out)]) := lapply(.SD, as.numeric), .SDcols = colnames(t.m2.boot.out)[3:ncol(t.m2.boot.out)]]
t.m1.mi.out[,c(colnames(t.m1.mi.out)[3:ncol(t.m1.mi.out)]) := lapply(.SD, as.numeric), .SDcols = colnames(t.m1.mi.out)[3:ncol(t.m1.mi.out)]]
t.m2.mi.out[,c(colnames(t.m2.mi.out)[3:ncol(t.m2.mi.out)]) := lapply(.SD, as.numeric), .SDcols = colnames(t.m2.mi.out)[3:ncol(t.m2.mi.out)]]

sink('./Output/WebTable4.tex')
kable(t.m1.boot.out,format = 'latex',digits = 3,col.names = c("Metric",'Nutrient',"Intercept","Age","BMI","log-Sodium","High Chol","USBorn","Female","Bkg PRican","Bkg Other"),align = c('l','l','c','c','c','c','c','c','c','c','c'),booktabs = T,escape = F,caption = 'Simulation results. Logistic regression with resampling-based multiple imputation correction of standard errors.')%>%kable_styling(latex_options = 'scale_down')%>%collapse_rows(columns = 1,latex_hline = 'major',valign = 'top') %>%
    footnote(general = " For SE (standard error), 'Calibrated' indicates the average (across 1000 simulations) of the model-based uncorrected standard errors from the outcome model using biomarker calibrated nutrients, and  'Corrected' indicates the average (across 1000 simulations) of the corrected standard errors from the outcome model using biomarker calibrated nutrients.\n For Coverage, the table shows the proportion of 95 percent confidence intervals (Estimate $\\\\pm z_{0.975}$SE) covering the true parameter (first row), where $z_{0.975}$ is the $97.5$ percentile of the standard Gaussian distribution.",escape = F,threeparttable = T)
sink()

sink('./Output/WebTable5.tex')
kable(t.m2.boot.out,format = 'latex',digits = 3,col.names = c("Metric",'Nutrient',"Intercept","Age","BMI","log-Sodium","High Chol","USBorn","Female","Bkg PRican","Bkg Other"),align = c('l','l','c','c','c','c','c','c','c','c','c'),booktabs = T,escape = F,caption = 'Simulation results. Linear regression with resampling-based multiple imputation correction of standard errors.')%>%kable_styling(latex_options = 'scale_down')%>%collapse_rows(columns = 1,latex_hline = 'major',valign = 'top') %>%
    footnote(general = " For SE (standard error), 'Calibrated' indicates the average (across 1000 simulations) of the model-based uncorrected standard errors from the outcome model using biomarker calibrated nutrients, and  'Corrected' indicates the average (across 1000 simulations) of the corrected standard errors from the outcome model using biomarker calibrated nutrients.\n For Coverage, the table shows the proportion of 95 percent confidence intervals (Estimate $\\\\pm z_{0.975}$SE) covering the true parameter (first row), where $z_{0.975}$ is the $97.5$ percentile of the standard Gaussian distribution.",escape = F,threeparttable = T)
sink()

sink('./Output/WebTable6.tex')
kable(t.m1.mi.out,format = 'latex',digits = 3,col.names = c("Metric",'Nutrient',"Intercept","Age","BMI","log-Sodium","High Chol","USBorn","Female","Bkg PRican","Bkg Other"),align = c('l','l','c','c','c','c','c','c','c','c','c'),booktabs = T,escape = F,caption = 'Simulation results. Logistic regression with parametric multiple imputation correction of standard errors.')%>%kable_styling(latex_options = 'scale_down')%>%collapse_rows(columns = 1,latex_hline = 'major',valign = 'top') %>%
    footnote(general = " For SE (standard error), 'Calibrated' indicates the average (across 1000 simulations) of the model-based uncorrected standard errors from the outcome model using biomarker calibrated nutrients, and  'Corrected' indicates the average (across 1000 simulations) of the corrected standard errors from the outcome model using biomarker calibrated nutrients.\n For Coverage, the table shows the proportion of 95 percent confidence intervals (Estimate $\\\\pm z_{0.975}$SE) covering the true parameter (first row), where $z_{0.975}$ is the $97.5$ percentile of the standard Gaussian distribution.",escape = F,threeparttable = T)
sink()

sink('./Output/WebTable7.tex')
kable(t.m2.mi.out,format = 'latex',digits = 3,col.names = c("Metric",'Nutrient',"Intercept","Age","BMI","log-Sodium","High Chol","USBorn","Female","Bkg PRican","Bkg Other"),align = c('l','l','c','c','c','c','c','c','c','c','c'),booktabs = T,escape = F,caption = 'Simulation results. Linear regression with parametric multiple imputation correction of standard errors.')%>%kable_styling(latex_options = 'scale_down')%>%collapse_rows(columns = 1,latex_hline = 'major',valign = 'top') %>%
    footnote(general = " For SE (standard error), 'Calibrated' indicates the average (across 1000 simulations) of the model-based uncorrected standard errors from the outcome model using biomarker calibrated nutrients, and  'Corrected' indicates the average (across 1000 simulations) of the corrected standard errors from the outcome model using biomarker calibrated nutrients.\n For Coverage, the table shows the proportion of 95 percent confidence intervals (Estimate $\\\\pm z_{0.975}$SE) covering the true parameter (first row), where $z_{0.975}$ is the $97.5$ percentile of the standard Gaussian distribution.",escape = F,threeparttable = T)
sink()

### Checking bias & coverage of true regression calibration paramaters across 1000 simulated data

load('./RData/TargetPopulationData.RData')
pop <- pop[(v.num==1),-1]

# Fitting the model in the population

lm.regcalib = glm(c_ln_na_bio1 ~ c_age + c_bmi + c_ln_na_avg + high_chol + usborn + female + bkg_pr + bkg_o,data=pop)

df.lm.regcalib <- data.table(Coeff = names(lm.regcalib$coefficients),True.Est = as.numeric(lm.regcalib$coefficients))

df.lm.regcalib$Coeff %<>% plyr::mapvalues(from = c('(Intercept)', 'c_age', 'c_bmi', 'c_ln_na_avg', 'high_chol', 'usborn','female', 'bkg_pr', 'bkg_o'),
                                          to = c('Intercept', 'Age (Centered)', 'BMI (Centered)', 'Log-Sodium (Naive)', 'High Cholesterol', 'US Born','Sex: Female', 'Background: Puerto Rico', 'Background: Other'))

# Sumarizing the results across simulations

df.regcalib <- df.m1.boot[,list(RegCalib.Est = unique(RegCalib.Est),
                                RegCalib.SE = unique(RegCalib.SE)),by=c('Sim','Coeff')]

df.regcalib$Coeff %<>% plyr::mapvalues(from = c('(Intercept)', 'c_age', 'c_bmi', 'c_ln_na_calib.boot', 'high_chol', 'usborn','female', 'bkg_pr', 'bkg_o'),
                                       to = c('Intercept', 'Age (Centered)', 'BMI (Centered)', 'Log-Sodium (Naive)', 'High Cholesterol', 'US Born','Sex: Female', 'Background: Puerto Rico', 'Background: Other'))

# Merging

df.regcalib <- merge(df.regcalib,df.lm.regcalib,by = 'Coeff',all.x = T)

# Calculating metrics

df.regcalib <- df.regcalib[,Coverage := 1*((True.Est>=(RegCalib.Est-1.96*RegCalib.SE)) | (True.Est<=(RegCalib.Est+1.96*RegCalib.SE)))]
df.regcalib <- df.regcalib[,Bias := RegCalib.Est - True.Est]

df.regcalib <- df.regcalib[,list(True.Est = unique(True.Est),
                                 Average.Est = mean(RegCalib.Est),
                                 Relative.Bias = mean((RegCalib.Est-True.Est)/True.Est),
                                 Average.SE = mean(RegCalib.SE),
                                 Empirical.SE = sd(RegCalib.Est),
                                 Coverage = mean(1*((True.Est>RegCalib.Est-1.96*RegCalib.SE)&(True.Est<RegCalib.Est+1.96*RegCalib.SE)))),by = 'Coeff']

df.regcalib$Coeff %<>% factor(levels = c('Intercept', 'Log-Sodium (Naive)',  'Age (Centered)', 'BMI (Centered)', 'High Cholesterol', 'US Born','Sex: Female', 'Background: Puerto Rico', 'Background: Other'))
df.regcalib <- df.regcalib[order(Coeff),]

# Column names
setnames(df.regcalib,c('Coefficient','True Parameter','Estimate','Relative Bias','SE','Empirical SE','Coverage'))

# Exporting table

sink('./Output/WebTable8.tex')
df.regcalib %>%
    kable(format = 'latex',digits = 3,booktabs = T,escape = F,
          caption = 'Simulation Metrics of Regression Calibration Across 1000 Simulations.')%>%kable_styling(latex_options = 'scale_down')%>%
    footnote(general = "Average values of Estimate, Relative Bias, and SE (standard error) calculated across 1000 simulations of estimates from the regression calibration fitted in SOLNAS simulated phase 2 subsample. Empirical SE calculated as the standard deviation of estimated coefficients. Coverage indicates the proportion of 95 percent confidence intervals covering the true parameter.",escape = F,threeparttable = T)
sink()

###############################################################
############# Creating Plots for Population Data ##############
###############################################################

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
    guides(color = guide_legend(override.aes = list(size = rel(0.3))))+ylim(0,1.35)+xlim(5.5,NA)

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
ggsave(figure1,filename = './Output/WebFigure1.eps',dpi = 'retina',width = 8.5,height = 7,device = cairo_ps,fallback_resolution = 300)

## Sodium-only figure

sdnutr = round(apply(pop[pop$v.num==1,c('ln_na_true','ln_na_bio1','ln_na_avg')],2,sd),2)
pop.na <- as_tibble(pop[pop$v.num==1,c('ln_na_true','ln_na_bio1','ln_na_avg')]) %>%
    gather(Type,Sodium,ln_na_true:ln_na_avg) %>%
    mutate(Type = recode(Type,ln_na_true = 'True', ln_na_bio1 = 'Biomarker', ln_na_avg='Self-Reported'))
pop.na$Type %<>% factor(levels = c('True','Biomarker','Self-Reported'))

figure1.sodiumonly.naplot = ggplot(pop.na,aes(x=Sodium,linetype = Type))+
    stat_density(geom="line", position="identity",size = convertPt(1))+
    theme_my(base_size = 52)+ylab('Density')+xlab('Log(Sodium), mg')+
    ylim(0,1)+xlim(5.450926,11.067032)+
    labs(linetype = 'Simulated Intake')+
    scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
    theme(legend.title.align=0.5,legend.position = c(0.8,0.8),legend.background = element_rect(linetype = 'solid',colour = 'black',size = convertPt(1)))

ggsave(figure1.sodiumonly.naplot,filename = './Output/Figure1.pdf',dpi = 'retina',width = 7,height = 5)
ggsave(figure1.sodiumonly.naplot,filename = './Output/Figure1.eps',dpi = 'retina',width = 7,height = 5,device = cairo_ps,fallback_resolution = 300)

dt1 <- melt(pop[pop$v.num==1,c('subid','sex','bkg','age','bmi','ln_na_true','ln_na_bio1','ln_na_avg')],
            id.vars = c('subid','sex','bkg'),measure.vars = c('ln_na_true','ln_na_bio1','ln_na_avg'),variable.name = 'Nutrient',value.name = 'Value')
dt1$sex %<>% plyr::mapvalues(from = c('M','F'),to = c('Male','Female')) %<>% factor(levels = c('Female','Male'))
dt1$bkg %<>% plyr::mapvalues(from = c('D','PR','O'),to = c('Dominican','Puerto Rican','Other')) %<>% factor(levels = c('Dominican','Puerto Rican','Other'))
dt1$Nutrient %<>% plyr::mapvalues(from = c('ln_na_true','ln_na_bio1','ln_na_avg'),to = c('True','Biomarker','Self-Reported')) %<>% factor(levels = c('True','Biomarker','Self-Reported'))

figure1.sodiumonly.boxplot.ls <- lapply(unique(dt1[,paste0(bkg,'-',sex)]),function(x){
    splitx <- strsplit(x,'-')[[1]]
    rg <- range(dt1$Value)*c(0.925,1.05)
    rg[1] <- floor(rg[1])
    rg[2] <- ceiling(rg[2])
    
    ggplot(dt1[bkg == splitx[1] & sex == splitx[2],],aes(y = Value,x = Nutrient))+
        geom_boxplot(outlier.alpha = 0.25) + #outlier.size = 0.75
        ylab('Log(Sodium), mg') +
        theme_my(base_size = 36)+
        scale_y_continuous(breaks = c(5,7,9,11),limits = c(5,11))+
        theme(axis.text.x=element_text(angle=20,hjust=1),
              plot.margin = unit(c(15, 5.5, 5.5, 15), "points"),axis.title.x = element_blank())
})

figure1.sodiumonly.boxplot <- ggarrange(plotlist = figure1.sodiumonly.boxplot.ls,nrow = 2,ncol = 3,labels = as.list(paste0(LETTERS[1:6],')')),
                                        hjust = 0,vjust = 1.25,font.label = list(size = 36/.pt,face = 'plain'))


figure1.sodiumonly.boxplot <- annotate_figure(figure1.sodiumonly.boxplot,bottom = text_grob("Simulated Intake",just = 'center', size = 36/.pt))
ggsave(figure1.sodiumonly.boxplot,filename = './Output/Figure2.pdf',dpi = 'retina',width = 7,height = 6)
ggsave(figure1.sodiumonly.boxplot,filename = './Output/Figure2.eps',dpi = 'retina',width = 7,height = 6,device = grDevices::cairo_ps,fallback_resolution = 300)

### Saving the plots

for(i in 1:6){
    subfig <- ggpubr::ggarrange(figure1.sodiumonly.boxplot.ls[[i]],ncol = 1,nrow = 1,labels = as.list(paste0(LETTERS[i],')')),
                                hjust = 0,vjust = 1.25,font.label = list(size = 36/.pt,face = 'plain'))
    ggsave(subfig,filename = paste0('./Output/Figure2',LETTERS[i],'.eps'),dpi = 'retina',height = 3,width = 3,device = grDevices::cairo_ps,fallback_resolution = 300)
}

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

sink('./Output/WebTable1.tex')
kable(dt,booktab = T,caption = '\\label{svysamp}Characteristics of the simulated target population',
      col.names=c('','(N)','(N)','(N)','','','(N)','(\\%)'),escape=F,format = "latex",
      align = c('c','c','c','c','l','l','r','r')) %>%
    add_header_above(c('Stratum'=1,'Block Groups'=1,'Households'=1,'Household Type'=1,'Background'=1,'Sex'=1,'Individuals'=1,'45+ years old'=1), line = F) %>%
    collapse_rows(columns = 1:5,latex_hline = 'major',valign = 'top') %>%
    row_spec(0,align = 'c') %>%
    kable_styling(latex_options = "scale_down")
sink()

# Creating WebTable3

load("./Output/htn_parameters.RData")
load("./Output/sbp_parameters.RData")

out <- data.table(Variable = sbp_parameters$coefficients$Variable)
out <- merge(out,sbp_parameters$coefficients,by = 'Variable',all.x = T)
out <- merge(out,htn_parameters,by = 'Variable',all.x = T)
out <- rbind(out,data.table(Variable = 'Sigma2',Estimate.x = sbp_parameters$variance,Estimate.y = NA))

out$Variable %<>% plyr::mapvalues(from = c("BKGRD_OTHER","BKGRD_PR","C_AGE","C_BMI","C_LN_NA_CALIBR","FEMALE","HIGH_TOTAL_CHOL","Intercept","US_BORN","Sigma2"),
                                  to = c('Background: Other','Background: Puerto Rican','Age (centered)','BMI (centered)','Log-sodium (centered)','Female','Hypercholesterolemia','Intercept','Nativity','Sigma2'))
out$Variable %<>% factor(levels = c('Intercept','Log-sodium (centered)','Age (centered)','BMI (centered)','Female','Background: Puerto Rican','Background: Other','Hypercholesterolemia','Nativity','Sigma2'))
out <- out[order(Variable),]
out[Variable=='Sigma2',Variable := '$\\sigma^{2}$']

sink('./Output/WebTable3.tex')
kable(out,format='latex',caption='\\label{simoutcoeff}Model coefficients from linear and logistic regression models for the outcome simulation.',
      digits=3,booktabs=T,escape = F,col.names = c('Coefficient','SBP','Hypertension Status')) %>%
    row_spec(9,hline_after=T)
sink()

# Creating WebTable2

load("./data/calibrCoeff.RData")

sodicoeff <- as.data.table(sodicoeff)[Coeff == 'LSODI_2DMEAN_C',Coeff := 'Nutrient'][,Estimate_NA := Estimate][,Estimate := NULL]
potacoeff <- as.data.table(potacoeff)[Coeff == 'LPOTA_2DMEAN_C',Coeff := 'Nutrient'][,Estimate_K := Estimate][,Estimate := NULL]
kcalcoeff <- as.data.table(kcalcoeff)[Coeff == 'LKCAL_2DMEAN_C',Coeff := 'Nutrient'][,Estimate_KCAL := Estimate][,Estimate := NULL]
protcoeff <- as.data.table(protcoeff)[Coeff == 'LPROT_2DMEAN_C',Coeff := 'Nutrient'][,Estimate_PROT := Estimate][,Estimate := NULL]

out <- data.table(Coeff = sodicoeff$Coeff)
out <- merge(out,sodicoeff,by = 'Coeff',all.x = T)
out <- merge(out,potacoeff,by = 'Coeff',all.x = T)
out <- merge(out,kcalcoeff,by = 'Coeff',all.x = T)
out <- merge(out,protcoeff,by = 'Coeff',all.x = T)

out$Coeff %<>% plyr::mapvalues(from = c("BKGRD_OTHER","BKGRD_PR","AGE_C","BMI_C","Nutrient","FEMALE","(Intercept)","US_BORN","HIGH_TOTAL_CHOL"),
                               to = c('Background: Other','Background: Puerto Rican','Age (centered)','BMI (centered)','Naive Intake (centered)','Female','Intercept','Nativity','Hypercholesterolemia'))
out$Coeff %<>% factor(levels = c('Intercept','Naive Intake (centered)','Age (centered)','BMI (centered)','Female','Background: Puerto Rican','Background: Other','Hypercholesterolemia','Nativity'))
out <- out[order(Coeff),]

sink('./Output/WebTable2.tex')
out %>%
    bind_rows(tibble(Coeff='$\\sigma^{2}$',Estimate_NA=var.df[1,2],Estimate_K=var.df[2,2],Estimate_KCAL=var.df[3,2],Estimate_PROT=var.df[4,2])) %>%
    bind_rows(tibble(Coeff='$\\sigma^{2}_{X}$',Estimate_NA=var.df[1,4],Estimate_K=var.df[2,4],Estimate_KCAL=var.df[3,4],Estimate_PROT=var.df[4,4])) %>%
    bind_rows(tibble(Coeff='$\\sigma^{2}_{\\epsilon}$',Estimate_NA=var.df[1,3],Estimate_K=var.df[2,3],Estimate_KCAL=var.df[3,3],Estimate_PROT=var.df[4,3])) %>%
    kable(booktab=T,format='latex',caption='\\label{tab:calibreg}Gaussian linear model coefficients and variance parameters used to simulate the log-transformed true dietary exposures and their respective biomarkers.',
          escape = F,digits=3,col.names = c('Coefficient','Log-sodium','Log-potassium','Log-energy','Log-protein')) %>%
    row_spec(7,hline_after=T)
sink()
