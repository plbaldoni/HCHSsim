library(plyr)
library(magrittr)
library(ggplot2)
library(haven)
library(data.table)
library(readr)

load('./Output/bootstrap.RData')

df.m1.boot <- m1.list.boot %<>% rbindlist()
df.m2.boot <- m2.list.boot %<>% rbindlist()

load('./Output/MI.RData')

df.m1.mi <- m1.list.mi %<>% rbindlist()
df.m2.mi <- m2.list.mi %<>% rbindlist()

load('./Output/raking_BetaVersion.RData')

df.m1.raking <- m1.list.raking %<>% rbindlist()
df.m2.raking <- m2.list.raking %<>% rbindlist()

# Color of plots
Coeff = c('Intercept','Age','BMI','Log-Sodium','High Chol.','US Born','Female','Background: PR','Background: Other')
coeffcol = c('#000000','#000000','#000000',"#E41A1C",'#000000','#000000','#000000','#000000','#000000')
names(coeffcol) = Coeff
# methcol = c('#000000',gray.colors(12)[c(1,2,4)])
# names(methcol) <- c('True (Unobservable)','Naive 2-day Mean','Calibrated w/ Bootstrap','Calibrated w/ MI')

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
m1.raking = summ(df=df.m1.raking,method='Raking',model='Hypertension',raking=T,digits = 3,
                 func.avg = 'mean',func.sd = 'sd',
                 func.var.avg = 'median',func.var.sd = 'mad')

m2.boot = summ(df=df.m2.boot,method='Bootstrap',model='SBP',raking=F,digits = 3,
               func.avg = 'mean',func.sd = 'sd',
               func.var.avg = 'median',func.var.sd = 'mad')
m2.mi = summ(df=df.m2.mi,method='MI',model='SBP',raking=F,digits = 3,
             func.avg = 'mean',func.sd = 'sd',
             func.var.avg = 'median',func.var.sd = 'mad')
m2.raking = summ(df=df.m2.raking,method='Raking',model='SBP',raking=T,digits = 3,
                 func.avg = 'mean',func.sd = 'sd',
                 func.var.avg = 'median',func.var.sd = 'mad')

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

pdf(file = './Output/summarize_logistic.pdf',height = 8.5,width = 11)
fig.m1
dev.off()

pdf(file = './Output/summarize_linear.pdf',height = 8.5,width = 11)
fig.m2
dev.off()

pt.title = 10
pt.text = 8

fig = ggarrange(fig.m1+theme(plot.margin = unit(c(5.5, 5.5, 7.5, 5.5), "points"))+scale_y_continuous(labels = function(x) sprintf("%.2f", x))+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                fig.m2+theme(plot.margin = unit(c(7.5, 5.5, 5.5, 5.5), "points"))+scale_y_continuous(labels = function(x) sprintf("%.2f", x))+theme(axis.text = element_text(size = pt.text),axis.title = element_text(size = pt.title)),
                nrow=2,ncol=1,legend = 'bottom',common.legend = T,labels = list('A','B'))

ggsave(filename = './Output/Figure3.pdf',plot = fig,height = 11,width = 8.5)

### Saving the output

save(m1.boot.out,m2.boot.out,m1.mi.out,m2.mi.out,file = './Output/summarize.RData')

# ### Checking the appropriate number of MI/Bootstrap iterations
# 
# traceplot <- function(df,method,x = c(2,3,4,5,seq(10,500,10)),
#                       func.avg = 'mean',func.sd='sd',
#                       func.var.avg = 'mean',func.var.sd = 'sd'){
#     
#     df$Coeff %<>% mapvalues(from=unique(df$Coeff),to=Coeff) %<>% as.character()
#     
#     funclist = list(func.avg,func.sd)
#     f.avg = match.fun(funclist[[1]]);f.sd = match.fun(funclist[[2]])
#     
#     funcvarlist = list(func.var.avg,func.var.sd)
#     f.var.avg = match.fun(funcvarlist[[1]]);f.var.sd = match.fun(funcvarlist[[2]])
#     
#     # Calculating the empirical SE
#     df.calib.est.ese <- df[,.(Calib.Est = unique(Calib.Est)),by=c('Coeff','Sim')][,.(Calib.Est = f.sd(Calib.Est)),by=c('Coeff')]
#     df.calib.est.ese$Coeff %<>% factor(levels = c('Intercept','Log-Sodium',as.character(unique(df.calib.est.ese$Coeff[!df.calib.est.ese$Coeff%in%c('Intercept','Log-Sodium')]))))
#     
#     # Calculating the empirical ESE
#     df.avg.est.ese <- df[,.(Avg.Est = unique(Avg.Est)),by=c('Coeff','Sim')][,.(Avg.Est = f.sd(Avg.Est)),by=c('Coeff')]
#     df.avg.est.ese$Coeff %<>% factor(levels = c('Intercept','Log-Sodium',as.character(unique(df.avg.est.ese$Coeff[!df.avg.est.ese$Coeff%in%c('Intercept','Log-Sodium')]))))
#     
#     # Bootstrap/MI correction for SE
#     df.calib.se <- rbindlist(lapply(x,FUN = function(y){
#         if(method=='Bootstrap'){sub.df <- df[Boot<=y,]}
#         if(method=='MI'){sub.df <- df[MI<=y,]}
#         
#         df.var.within = sub.df[,.(withinVar = f.var.avg(SE^2)),by=c('Coeff','Sim')]
#         df.var.between = sub.df[,.(betweenSE = f.var.sd(Est)),by=c('Coeff','Sim')]
#         df.var = merge(x=df.var.within,y=df.var.between,by=c('Coeff','Sim'))
#         df.var[,Corrected.SE := sqrt(withinVar + betweenSE^2)]
#         df.out <- df.var[,.(Median = median(Corrected.SE),Pct_025 = quantile(Corrected.SE,0.025),Pct_975 = quantile(Corrected.SE,0.975)),by = 'Coeff']
#         df.out[,Cutoff := y]
#     }))
#     df.calib.se$Coeff %<>% factor(levels = c('Intercept','Log-Sodium',as.character(unique(df.calib.se$Coeff[!df.calib.se$Coeff%in%c('Intercept','Log-Sodium')]))))
#     
#     df.ese <- merge(df.calib.est.ese,df.avg.est.ese,by = 'Coeff')
#     df.ese <- melt(df.ese,id.vars = 'Coeff',measure.vars = c('Calib.Est','Avg.Est'),variable.name = 'Empirical',value.name = 'SE')
#     df.ese$Empirical %<>% mapvalues(from = c('Calib.Est','Avg.Est'), to = c('Calibrated','Naive'))
#     
#     return(list('trace' = df.calib.se, 'ese' = df.ese))
# }
# 
# trace.m1.boot <- traceplot(df = df.m1.boot,method = 'Bootstrap')
# trace.m2.boot <- traceplot(df = df.m2.boot,method = 'Bootstrap')
# trace.m1.mi <- traceplot(df = df.m1.mi,method = 'MI')
# trace.m2.mi <- traceplot(df = df.m2.mi,method = 'MI')
# 
# fig.trace.m1.boot <- ggplot(data = trace.m1.boot$trace,aes(x = Cutoff))+
#     facet_wrap(~Coeff,scales = 'free_y') +
#     geom_hline(data = trace.m1.boot$ese, aes(yintercept = SE,linetype = Empirical,color = Empirical))+
#     geom_ribbon(aes(ymin = Pct_025, ymax = Pct_975, fill = "2.5 and 97.5 percentiles"),alpha = 0.5) +
#     guides(fill=guide_legend(title=element_blank(),order = 2),color = guide_legend(order = 1),linetype = guide_legend(order = 1))+
#     scale_fill_manual(values = c("2.5 and 97.5 percentiles" = 'grey'))+
#     geom_line(aes(y = Median)) +
#     theme_bw() + labs(x = 'Number of Imputations',y = 'Corrected SE (median estimate across 1000 simulations)', color = 'Empirical SE',linetype = 'Empirical SE')+
#     theme(legend.position = 'top',legend.direction = 'horizontal',panel.grid.major = element_blank(),panel.grid.minor = element_blank())
# 
# fig.trace.m2.boot <- ggplot(data = trace.m2.boot$trace,aes(x = Cutoff))+
#     facet_wrap(~Coeff,scales = 'free_y') +
#     geom_hline(data = trace.m2.boot$ese, aes(yintercept = SE,linetype = Empirical,color = Empirical))+
#     geom_ribbon(aes(ymin = Pct_025, ymax = Pct_975, fill = "2.5 and 97.5 percentiles"),alpha = 0.5) +
#     guides(fill=guide_legend(title=element_blank(),order = 2),color = guide_legend(order = 1),linetype = guide_legend(order = 1))+
#     scale_fill_manual(values = c("2.5 and 97.5 percentiles" = 'grey'))+
#     geom_line(aes(y = Median)) +
#     theme_bw() + labs(x = 'Number of Imputations',y = 'Corrected SE (median estimate across 1000 simulations)', color = 'Empirical SE',linetype = 'Empirical SE')+
#     theme(legend.position = 'top',legend.direction = 'horizontal',panel.grid.major = element_blank(),panel.grid.minor = element_blank())
# 
# fig.trace.m1.mi <- ggplot(data = trace.m1.mi$trace,aes(x = Cutoff))+
#     facet_wrap(~Coeff,scales = 'free_y') +
#     geom_hline(data = trace.m1.mi$ese, aes(yintercept = SE,linetype = Empirical,color = Empirical))+
#     geom_ribbon(aes(ymin = Pct_025, ymax = Pct_975, fill = "2.5 and 97.5 percentiles"),alpha = 0.5) +
#     guides(fill=guide_legend(title=element_blank(),order = 2),color = guide_legend(order = 1),linetype = guide_legend(order = 1))+
#     scale_fill_manual(values = c("2.5 and 97.5 percentiles" = 'grey'))+
#     geom_line(aes(y = Median)) +
#     theme_bw() + labs(x = 'Number of Imputations',y = 'Corrected SE (median estimate across 1000 simulations)', color = 'Empirical SE',linetype = 'Empirical SE')+
#     theme(legend.position = 'top',legend.direction = 'horizontal',panel.grid.major = element_blank(),panel.grid.minor = element_blank())
# 
# fig.trace.m2.mi <- ggplot(data = trace.m2.mi$trace,aes(x = Cutoff))+
#     facet_wrap(~Coeff,scales = 'free_y') +
#     geom_hline(data = trace.m2.mi$ese, aes(yintercept = SE,linetype = Empirical,color = Empirical))+
#     geom_ribbon(aes(ymin = Pct_025, ymax = Pct_975, fill = "2.5 and 97.5 percentiles"),alpha = 0.5) +
#     guides(fill=guide_legend(title=element_blank(),order = 2),color = guide_legend(order = 1),linetype = guide_legend(order = 1))+
#     scale_fill_manual(values = c("2.5 and 97.5 percentiles" = 'grey'))+
#     geom_line(aes(y = Median)) +
#     theme_bw() + labs(x = 'Number of Imputations',y = 'Corrected SE (median estimate across 1000 simulations)', color = 'Empirical SE',linetype = 'Empirical SE')+
#     theme(legend.position = 'top',legend.direction = 'horizontal',panel.grid.major = element_blank(),panel.grid.minor = element_blank())
# 
# ### Saving plots
# 
# ggsave(filename = './Output/traceplot.logistic.bootstrap.pdf',plot = fig.trace.m1.boot,width = 6,height = 5)
# ggsave(filename = './Output/traceplot.linear.bootstrap.pdf',plot = fig.trace.m2.boot,width = 6,height = 5)
# ggsave(filename = './Output/traceplot.logistic.MI.pdf',plot = fig.trace.m1.mi,width = 6,height = 5)
# ggsave(filename = './Output/traceplot.linear.MI.pdf',plot = fig.trace.m2.mi,width = 6,height = 5)