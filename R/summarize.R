library(plyr)
library(magrittr)
library(ggplot2)
library(haven)
library(data.table)
library(readr)
library(ggpubr)
library(kableExtra)

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

# load('./Output/raking_BetaVersion.RData')
# 
# df.m1.raking <- m1.list.raking %<>% rbindlist()
# df.m2.raking <- m2.list.raking %<>% rbindlist()

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
        # ggsave(filename = paste('./Misc/Hist',model,paste0(method,'Est','_',func.avg,'_','.pdf'),sep='_'),plot = fig.CalibEst)
        # 
        # fig.AvgEst =  ggplot(data=df.avg.est)+geom_histogram(position='identity',aes(x=Avg.Est,y=..density..))+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+theme_bw()+
        #     geom_vline(data=subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),aes(xintercept=Estimate))+
        #     xlab('Estimated coefficients from outcome-average model')
        # ggsave(filename = paste('./Misc/Hist',model,'AvgEst.pdf',sep='_'),plot = fig.AvgEst)
        # 
        # fig.TrueEst =  ggplot(data=df.true.est)+geom_histogram(position='identity',aes(x=True.Est,y=..density..))+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+theme_bw()+
        #     geom_vline(data=subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),aes(xintercept=Estimate))+
        #     xlab('Estimated coefficients from outcome-true model')
        # ggsave(filename = paste('./Misc/Hist',model,'TrueEst.pdf',sep='_'),plot = fig.TrueEst)
        # 
        # fig.BioEst =  ggplot(data=df.bio.est)+geom_histogram(position='identity',aes(x=Bio.Est,y=..density..))+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+theme_bw()+
        #     geom_vline(data=subset(pop.coeff,Model==model,select=c(Coeff,Estimate)),aes(xintercept=Estimate))+
        #     xlab('Estimated coefficients from outcome-biomarker (phase 2) model')
        # ggsave(filename = paste('./Misc/Hist',model,'BioEst.pdf',sep='_'),plot = fig.BioEst)
        # 
        # fig.vif1 = ggplot(data=df[!(Coeff=='Intercept'),])+geom_histogram(position='identity',aes(x=VIF,y=..density..))+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+theme_bw()+geom_vline(aes(xintercept=5),color='red')+
        #     xlab('Variance Inflation Factor')+ylab('Density')
        # ggsave(filename = paste('./Misc/VIF',model,'Hist.pdf',sep='_'),plot = fig.vif1)
        # 
        # fig.vif2 = ggplot(data=df[!(Coeff=='Intercept'),],aes(x=VIF,y=Est))+
        #     geom_smooth(se=T)+
        #     facet_wrap(~Coeff,ncol=3,scales='free')+geom_vline(aes(xintercept=5),color='red')+
        #     xlab('Variance Inflation Factor')+ylab('Parameter estimates from Outcome-calibrated model')
        # ggsave(filename = paste('./Misc/VIF',model,'Scatter.pdf',sep='_'),plot = fig.vif2,dpi='screen')
    }
}

# Tables - Average/SD
m1.boot = summ(df=df.m1.boot,method='Bootstrap',model='Hypertension',raking=F,digits = 3,
               func.avg = 'mean',func.sd = 'sd',
               func.var.avg = 'median',func.var.sd = 'mad')
m1.mi = summ(df=df.m1.mi,method='MI',model='Hypertension',raking=F,digits = 3,
             func.avg = 'mean',func.sd = 'sd',
             func.var.avg = 'median',func.var.sd = 'mad')
# m1.raking = summ(df=df.m1.raking,method='Raking',model='Hypertension',raking=T,digits = 3,
#                  func.avg = 'mean',func.sd = 'sd',
#                  func.var.avg = 'median',func.var.sd = 'mad')

m2.boot = summ(df=df.m2.boot,method='Bootstrap',model='SBP',raking=F,digits = 3,
               func.avg = 'mean',func.sd = 'sd',
               func.var.avg = 'median',func.var.sd = 'mad')
m2.mi = summ(df=df.m2.mi,method='MI',model='SBP',raking=F,digits = 3,
             func.avg = 'mean',func.sd = 'sd',
             func.var.avg = 'median',func.var.sd = 'mad')
# m2.raking = summ(df=df.m2.raking,method='Raking',model='SBP',raking=T,digits = 3,
#                  func.avg = 'mean',func.sd = 'sd',
#                  func.var.avg = 'median',func.var.sd = 'mad')

### Organizing the output
key = Coeff

m1.boot = merge(subset(pop.coeff,Model=='Hypertension',select=c(Coeff,Estimate)),m1.boot,by='Coeff',all.x=T);m1.boot.out <- m1.boot[order(match(Coeff,key)),];m1.boot.out[,Estimate := round(Estimate,3)]
m1.mi = merge(subset(pop.coeff,Model=='Hypertension',select=c(Coeff,Estimate)),m1.mi,by='Coeff',all.x=T);m1.mi.out <- m1.mi[order(match(Coeff,key)),];m1.mi.out[,Estimate := round(Estimate,3)]
# m1.raking = merge(subset(pop.coeff,Model=='Hypertension',select=c(Coeff,Estimate)),m1.raking,by='Coeff',all.x=T);m1.raking.out <- m1.raking[order(match(Coeff,key)),];m1.raking.out[,Estimate := round(Estimate,3)]

m2.boot = merge(subset(pop.coeff,Model=='SBP',select=c(Coeff,Estimate)),m2.boot,by='Coeff',all.x=T);m2.boot.out <- m2.boot[order(match(Coeff,key)),];m2.boot.out[,Estimate := round(Estimate,3)]
m2.mi = merge(subset(pop.coeff,Model=='SBP',select=c(Coeff,Estimate)),m2.mi,by='Coeff',all.x=T);m2.mi.out <- m2.mi[order(match(Coeff,key)),];m2.mi.out[,Estimate := round(Estimate,3)]
# m2.raking = merge(subset(pop.coeff,Model=='SBP',select=c(Coeff,Estimate)),m2.raking,by='Coeff',all.x=T);m2.raking.out <- m2.raking[order(match(Coeff,key)),];m2.raking.out[,Estimate := round(Estimate,3)]

# Saving 
write.csv(m1.boot.out,file=paste0('./Output/summarize_bootstrap_logistic.csv'))
write.csv(m2.boot.out,file=paste0('./Output/summarize_bootstrap_linear.csv'))
write.csv(m1.mi.out,file=paste0('./Output/summarize_MI_logistic.csv'))
write.csv(m2.mi.out,file=paste0('./Output/summarize_MI_linear.csv'))
# write.csv(m1.raking.out,file=paste0('./Output/summarize_raking_logistic.csv'))
# write.csv(m2.raking.out,file=paste0('./Output/summarize_raking_linear.csv'))

### Organizing the data to plot

m1.plot <- rbindlist(list(m1.boot.out[,c('Coeff','Estimate','True.Est','True.ESE','True.SE','True.Cov')][,Method := 'True (Phase1, Unobservable)'],
                          # m1.raking.out[,c('Coeff','Estimate','Bio.Moments.Est','Bio.Moments.ESE','Bio.Moments.SE','Bio.Moments.Cov')][,Method := 'Biomarker (Phase 2)'],
                          # m1.raking.out[,c('Coeff','Estimate','Bio.Moments.Est','Bio.Moments.ESE','Bio.Moments.Corrected.SE','Bio.Moments.Corrected.Cov')][,Method := 'Biomarker w/ Bootstrap (Phase 2)'],
                          m1.boot.out[,c('Coeff','Estimate','Naive.Est','Naive.ESE','Naive.SE','Naive.Cov')][,Method := 'Naive 2-day Mean (Phase 1)'],
                          m1.boot.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Calib.SE','Calib.Cov')][,Method := 'Naive Calibrated (Phase 1)'],
                          m1.boot.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Corrected.SE','Corrected.Cov.Normal')][,Method := 'Naive Calibrated w/ Bootstrap (Phase 1)'],
                          m1.mi.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Corrected.SE','Corrected.Cov.Normal')][,Method := 'Naive Calibrated w/ MI (Phase 1)']),use.names=FALSE)
setnames(m1.plot,c('Coeff','True','Estimate','ESE','SE','Cov','Method'))
# m1.plot$Coeff %<>% mapvalues(from = unique(.),to = c('Intercept','Age','BMI','Log-Sodium','High Chol.','US Born','Female','Bkg: Puerto Rican','Bkg: Other'))
m1.plot$Coeff %<>% factor(levels = unique(.))
m1.plot$Method %<>% as.factor() %<>% factor(levels = unique(.))
m1.plot[,Lbl := min(Estimate-SE,Estimate-ESE)-max(SE,ESE),by='Coeff']

m2.plot <- rbindlist(list(m2.boot.out[,c('Coeff','Estimate','True.Est','True.ESE','True.SE','True.Cov')][,Method := 'True (Phase1, Unobservable)'],
                          # m2.raking.out[,c('Coeff','Estimate','Bio.Moments.Est','Bio.Moments.ESE','Bio.Moments.SE','Bio.Moments.Cov')][,Method := 'Biomarker (Phase 2)'],
                          # m2.raking.out[,c('Coeff','Estimate','Bio.Moments.Est','Bio.Moments.ESE','Bio.Moments.Corrected.SE','Bio.Moments.Corrected.Cov')][,Method := 'Biomarker w/ Bootstrap (Phase 2)'],
                          m2.boot.out[,c('Coeff','Estimate','Naive.Est','Naive.ESE','Naive.SE','Naive.Cov')][,Method := 'Naive 2-day Mean (Phase 1)'],
                          m2.boot.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Calib.SE','Calib.Cov')][,Method := 'Naive Calibrated (Phase 1)'],
                          m2.boot.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Corrected.SE','Corrected.Cov.Normal')][,Method := 'Naive Calibrated w/ Bootstrap (Phase 1)'],
                          m2.mi.out[,c('Coeff','Estimate','Calib.Est','Calib.ESE','Corrected.SE','Corrected.Cov.Normal')][,Method := 'Naive Calibrated w/ MI (Phase 1)']),use.names=FALSE)
setnames(m2.plot,c('Coeff','True','Estimate','ESE','SE','Cov','Method'))
# m2.plot$Coeff %<>% mapvalues(from = unique(.),to = c('Intercept','Age','BMI','Log-Sodium','High Chol.','US Born','Female','Bkg: Puerto Rican','Bkg: Other'))
m2.plot$Coeff %<>% factor(levels = unique(.))
m2.plot$Method %<>% as.factor() %<>% factor(levels = unique(.))
m2.plot[,Lbl := min(Estimate-SE,Estimate-ESE)-max(SE,ESE),by='Coeff']

### Subsetting to the relevant methods
m1.plot.sub <- m1.plot[Method %in% c('True (Phase1, Unobservable)','Naive 2-day Mean (Phase 1)','Naive Calibrated w/ Bootstrap (Phase 1)','Naive Calibrated w/ MI (Phase 1)'),]
m1.plot.sub$Method %<>% mapvalues(from = unique(.),to = c('True\n(Unobservable)','Naive','Calibrated\n(Bootstrap)','Calibrated\n(MI)'))
m1.plot.sub$Coeff %<>% factor(levels = c('Intercept','Log-Sodium',as.character(unique(m1.plot.sub$Coeff[!m1.plot.sub$Coeff%in%c('Intercept','Log-Sodium')]))))

m2.plot.sub <- m2.plot[Method %in% c('True (Phase1, Unobservable)','Naive 2-day Mean (Phase 1)','Naive Calibrated w/ Bootstrap (Phase 1)','Naive Calibrated w/ MI (Phase 1)'),]
m2.plot.sub$Method %<>% mapvalues(from = unique(.),to = c('True\n(Unobservable)','Naive','Calibrated\n(Bootstrap)','Calibrated\n(MI)'))
m2.plot.sub$Coeff %<>% factor(levels = c('Intercept','Log-Sodium',as.character(unique(m2.plot.sub$Coeff[!m2.plot.sub$Coeff%in%c('Intercept','Log-Sodium')]))))

### Plotting Error bars
# names(methcol) = c('True\n(Unobservable)','Naive','Calibrated\n(Bootstrap)','Calibrated\n(MI)')

fig.m1 = ggplot(data = m1.plot.sub[!Method=="True\n(Unobservable)",],aes(x = Method, y = Estimate,shape = Method)) +
    facet_wrap(~Coeff,scales = 'free_y') +
    geom_hline(data = m1.plot.sub, aes(yintercept = True))+
    geom_pointrange(aes(ymin = Estimate-SE,ymax = Estimate+SE,color = Coeff),size=0.5) +
    geom_point(aes(y = Estimate+ESE,color = Coeff),shape = 5)+
    geom_point(aes(y = Estimate-ESE,color = Coeff),shape = 5)+
    geom_text(aes(label = paste0(formatC(round(100*Cov,1),format='f',digits=1),'%'),y = Lbl),size = 3.5,position = position_dodge(0.9),vjust = 0,fontface = "bold") +
    theme_bw()+
    scale_color_manual(values = coeffcol)+
    guides(shape = FALSE, color = FALSE)+
    theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.text = NULL,axis.title.x = element_blank(),legend.title = element_blank())

fig.m2 =  ggplot(data = m2.plot.sub[!Method=="True\n(Unobservable)",],aes(x = Method, y = Estimate,shape = Method)) +
    facet_wrap(~Coeff,scales = 'free_y') +
    geom_hline(data = m2.plot.sub, aes(yintercept = True))+
    geom_pointrange(aes(ymin = Estimate-SE,ymax = Estimate+SE,color = Coeff),size=0.5) +
    geom_point(aes(y = Estimate+ESE,color = Coeff),shape = 5)+
    geom_point(aes(y = Estimate-ESE,color = Coeff),shape = 5)+
    geom_text(aes(label = paste0(formatC(round(100*Cov,1),format='f',digits=1),'%'),y = Lbl),size = 3.5,position = position_dodge(0.9),vjust = 0,fontface = "bold") +
    theme_bw()+
    scale_color_manual(values = coeffcol)+
    guides(shape = FALSE, color = FALSE)+
    theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.text = NULL,axis.title.x = element_blank(),legend.title = element_blank())

### Saving the plots

ggsave(fig.m1,filename = './Output/Figure3A.eps',
       dpi = 'retina',height = 8.5,width = 11,device = cairo_ps,fallback_resolution = 300)

ggsave(fig.m2,filename = './Output/Figure3B.pdf',
       dpi = 'retina',height = 8.5,width = 11,device = cairo_ps,fallback_resolution = 300)

fig = ggarrange(fig.m1+theme(plot.margin = unit(c(5.5, 5.5, 7.5, 5.5), "points"))+scale_y_continuous(labels = function(x) sprintf("%.2f", x))+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                fig.m2+theme(plot.margin = unit(c(7.5, 5.5, 5.5, 5.5), "points"))+scale_y_continuous(labels = function(x) sprintf("%.2f", x))+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                nrow=2,ncol=1,legend = 'bottom',common.legend = T,labels = list('A','B'))

ggsave(fig,filename = './Output/Figure3.pdf',
       dpi = 'retina',height = 11,width = 8.5,device = cairo_ps,fallback_resolution = 300)

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
    geom_hline(yintercept = c(log10(10),log10(15)),linetype = 'dashed',color = 'grey',alpha = 0.75)+
    scale_fill_brewer(palette = "Set1")+
    annotate('text',x = c(2.5,2.5), y = 0.95*c(log10(10),log10(15)),label = c('10 sec.','15 sec.'))+
    theme(legend.position = c(0.75,0.75),axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))

ggsave(fig.time,filename = './Output/SuppFigure1.pdf',
       dpi = 'retina',height = 11/2,width = 8.5/1.25,device = cairo_ps,fallback_resolution = 300)

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

sink('./Output/SuppTable1.tex')
kable(t.m1.boot.out,format = 'latex',digits = 3,col.names = c("Metric",'Nutrient',"Intercept","Age","BMI","log-Sodium","High Chol","USBorn","Female","Bkg PRican","Bkg Other"),align = c('l','l','c','c','c','c','c','c','c','c','c'),booktabs = T,escape = F,caption = 'Simulation results. Logistic regression with resampling-based multiple imputation correction of standard errors.')%>%kable_styling(latex_options = 'scale_down')%>%collapse_rows(columns = 1,latex_hline = 'major',valign = 'top') %>%
    footnote(general = " For SE (standard error), 'Calibrated' indicates the average (across 1000 simulations) of the model-based uncorrected standard errors from the outcome model using biomarker calibrated nutrients, and  'Corrected' indicates the average (across 1000 simulations) of the corrected standard errors from the outcome model using biomarker calibrated nutrients.\n For Coverage, the table shows the proportion of 95 percent confidence intervals (Estimate $\\\\pm z_{0.975}$SE) covering the true parameter (first row), where $z_{0.975}$ is the $97.5$ percentile of the standard Gaussian distribution.",escape = F,threeparttable = T)
sink()

sink('./Output/SuppTable2.tex')
kable(t.m2.boot.out,format = 'latex',digits = 3,col.names = c("Metric",'Nutrient',"Intercept","Age","BMI","log-Sodium","High Chol","USBorn","Female","Bkg PRican","Bkg Other"),align = c('l','l','c','c','c','c','c','c','c','c','c'),booktabs = T,escape = F,caption = 'Simulation results. Linear regression with resampling-based multiple imputation correction of standard errors.')%>%kable_styling(latex_options = 'scale_down')%>%collapse_rows(columns = 1,latex_hline = 'major',valign = 'top') %>%
    footnote(general = " For SE (standard error), 'Calibrated' indicates the average (across 1000 simulations) of the model-based uncorrected standard errors from the outcome model using biomarker calibrated nutrients, and  'Corrected' indicates the average (across 1000 simulations) of the corrected standard errors from the outcome model using biomarker calibrated nutrients.\n For Coverage, the table shows the proportion of 95 percent confidence intervals (Estimate $\\\\pm z_{0.975}$SE) covering the true parameter (first row), where $z_{0.975}$ is the $97.5$ percentile of the standard Gaussian distribution.",escape = F,threeparttable = T)
sink()

sink('./Output/SuppTable3.tex')
kable(t.m1.mi.out,format = 'latex',digits = 3,col.names = c("Metric",'Nutrient',"Intercept","Age","BMI","log-Sodium","High Chol","USBorn","Female","Bkg PRican","Bkg Other"),align = c('l','l','c','c','c','c','c','c','c','c','c'),booktabs = T,escape = F,caption = 'Simulation results. Logistic regression with parametric multiple imputation correction of standard errors.')%>%kable_styling(latex_options = 'scale_down')%>%collapse_rows(columns = 1,latex_hline = 'major',valign = 'top') %>%
    footnote(general = " For SE (standard error), 'Calibrated' indicates the average (across 1000 simulations) of the model-based uncorrected standard errors from the outcome model using biomarker calibrated nutrients, and  'Corrected' indicates the average (across 1000 simulations) of the corrected standard errors from the outcome model using biomarker calibrated nutrients.\n For Coverage, the table shows the proportion of 95 percent confidence intervals (Estimate $\\\\pm z_{0.975}$SE) covering the true parameter (first row), where $z_{0.975}$ is the $97.5$ percentile of the standard Gaussian distribution.",escape = F,threeparttable = T)
sink()

sink('./Output/SuppTable4.tex')
kable(t.m2.mi.out,format = 'latex',digits = 3,col.names = c("Metric",'Nutrient',"Intercept","Age","BMI","log-Sodium","High Chol","USBorn","Female","Bkg PRican","Bkg Other"),align = c('l','l','c','c','c','c','c','c','c','c','c'),booktabs = T,escape = F,caption = 'Simulation results. Linear regression with parametric multiple imputation correction of standard errors.')%>%kable_styling(latex_options = 'scale_down')%>%collapse_rows(columns = 1,latex_hline = 'major',valign = 'top') %>%
    footnote(general = " For SE (standard error), 'Calibrated' indicates the average (across 1000 simulations) of the model-based uncorrected standard errors from the outcome model using biomarker calibrated nutrients, and  'Corrected' indicates the average (across 1000 simulations) of the corrected standard errors from the outcome model using biomarker calibrated nutrients.\n For Coverage, the table shows the proportion of 95 percent confidence intervals (Estimate $\\\\pm z_{0.975}$SE) covering the true parameter (first row), where $z_{0.975}$ is the $97.5$ percentile of the standard Gaussian distribution.",escape = F,threeparttable = T)
sink()