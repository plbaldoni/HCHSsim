library(stringr)

S=1000
set.seed(20190414)

### Generate Samples for Simulation Study (based on plan from 10/14/16) ###

popname = 'TargetPopulationData'
sampname = 'SampleData'

####### read-in population datasets #######
load(paste0('./RData/',popname,".RData"))

####### function to generate sample #######
samp.hh.bystrat=function(s,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other){
  hh.select.s=rep(FALSE,dim(pop.unq[pop.unq$strat==s,])[1])         #initially set hh.select.s=FALSE for all subjects in stratum s, so that unselected BG's will be set to FALSE
  for (j in bg.select[bg.select.s==s]){
    hh.list.hisp=hh.list$hhid[hh.list$BGid==j & hisp.strat==TRUE]   #list of unique HHs in BGid j, w/ Hisp surname
    hh.list.other=hh.list$hhid[hh.list$BGid==j & hisp.strat==FALSE] #list of unique HHs in BGid j, w/ other surname
    hh.select.hisp=sample(hh.list.hisp,round(prob.hh.hisp[s]*length(hh.list.hisp)))      #list of randomly sampled HHs in BGid w/ Hispanic surname
    hh.select.other=sample(hh.list.other,round(prob.hh.other[s]*length(hh.list.other)))  #list of randomly sampled HHs in BGid w/ other surname
    hh.select.s[pop.unq$hhid[pop.unq$strat==s] %in% c(hh.select.hisp,hh.select.other)]=TRUE  #select subjects from randomly sampled HHs in BGid j (one entry per subject in BGid j)
  }
  return(hh.select.s)
}

samp.gen = function(pop,prob.bg,num.bg,hisp.prop,other.prop,prob.hh.hisp,prob.hh.other,prob.age){
  bgstrat = cumsum(num.bg)
  pop.unq=pop[pop$v.num==1,]   #this is the population only including visit 1 records (to simplify sampling)  
  
  ### re-create hh.size ###
  hh.list=unique(pop.unq[,c("BGid","hhid")])  # list of BGid & hhid for each unique hhid
  bg.size=table(hh.list$BGid)   # number of HHs per BG, 
  hh.size=table(pop.unq$hhid)   # number of subjects per HH
  num.hisp.strat=round(c(bg.size[as.numeric(names(bg.size))<=bgstrat[1]]*hisp.prop[1]/(hisp.prop[1]+other.prop[1]),bg.size[as.numeric(names(bg.size))>=(bgstrat[1]+1) & as.numeric(names(bg.size))<=bgstrat[2]]*hisp.prop[2]/(hisp.prop[2]+other.prop[2]),bg.size[as.numeric(names(bg.size))>=(bgstrat[2]+1) & as.numeric(names(bg.size))<=bgstrat[3]]*hisp.prop[3]/(hisp.prop[3]+other.prop[3]),bg.size[as.numeric(names(bg.size))>=(bgstrat[3]+1)]*hisp.prop[4]/(hisp.prop[4]+other.prop[4])))
  num.other.strat=bg.size-num.hisp.strat
  A=matrix(rep(NA,times=max(bg.size)*length(bg.size)),nrow=max(bg.size),ncol=length(bg.size))
  for (i in 1:length(bg.size)){
    A[,i]=c(rep(TRUE,times=num.hisp.strat[i]),rep(FALSE,times=num.other.strat[i]),rep(NA,times=max(bg.size)-bg.size[i]))    # create matrix A with each column corresponding to a BG (containing 1's for each Hispanic HH followed by 0's for each other HH)
  }
  hisp.strat=na.omit(c(A))   #indicator for Hispanic surname (one entry per HH)
  
  ### generate raw weights (these do not depend on the sample (only depend on sampling probabilities)) ###
  pop$W.bg=1/prob.bg[pop$strat]   # BG stage raw weight (based on BG sampling fraction; 1/sampling fraction for that stratum)
  pop$W.hh=ifelse(pop$hhid %in% hh.list$hhid[hisp.strat],1/prob.hh.hisp[pop$strat],1/prob.hh.other[pop$strat])   #HH stage raw weight
  pop$W.sub=ifelse(pop$age.strat,1/prob.age[2],1/prob.age[1])    # subject stage raw weight
  
  ### select random sample from population ###
  #select stratified random sample of BGs from pop & save list of BGs
  bg.select=c(sample(1:bgstrat[1],round(prob.bg[1]*num.bg[1])),sample((bgstrat[1]+1):bgstrat[2],round(prob.bg[2]*num.bg[2])),sample((bgstrat[2]+1):bgstrat[3],round(prob.bg[3]*num.bg[3])),sample((bgstrat[3]+1):bgstrat[4],round(prob.bg[4]*num.bg[4])))
  bg.select.s=1*(bg.select<=bgstrat[1])+2*(bg.select>=(bgstrat[1]+1) & bg.select<=bgstrat[2])+3*(bg.select>=(bgstrat[2]+1) & bg.select<=bgstrat[3])+4*(bg.select>=(bgstrat[3]+1) & bg.select<=bgstrat[4])  #stratum for each BG
  
  #select stratified random sample of HHs from each selected BG & save indicator of HH selection for all subjects within selected HHs
  hh.select=c(samp.hh.bystrat(1,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other),samp.hh.bystrat(2,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other),
              samp.hh.bystrat(3,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other),samp.hh.bystrat(4,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other))  #indicator of HH selection
  
  #select random sample of subjects from each selected HH & save indicator of subject selection
  sub.select=rep(FALSE,dim(pop.unq)[1])           # initially set sub.select=FALSE for all subjects, so that unselected HH's will be set to FALSE
  sub.select[pop.unq$subid %in% sample(pop.unq$subid[!pop.unq$age.strat & hh.select],round(prob.age[1]*dim(pop.unq[!pop.unq$age.strat & hh.select,])[1]))]=TRUE   # randomly sample younger subjects among sampled HH's
  sub.select[pop.unq$subid %in% sample(pop.unq$subid[pop.unq$age.strat & hh.select],round(prob.age[2]*dim(pop.unq[pop.unq$age.strat & hh.select,])[1]))]=TRUE     # randomly sample older subjects among sampled HH's
  
  samp=pop[rep(sub.select,each=2),]   # create sample by restricting pop to selected subjects (replicate sub.select vector to select both visits from selected subjects)
  
  # generate normalized weight
  samp$W.bghhsub=samp$W.bg*samp$W.hh*samp$W.sub     # raw combined weights
  samp$bghhsub_s2=samp$W.bghhsub/mean(samp$W.bghhsub[samp$v.num==1])    # normalized weight (raw combined weight/mean combined weight)
  
  #samp=samp[, colnames(samp)!="hhid"]
  
  return(samp)
}

####### generate S samples from same population for each scenario #######
for(i in 1:S){cat(i,'\n')
  samp=samp.gen(pop,
                prob.bg=c(.25,.25,.6,.6),
                num.bg=c(58,21,130,167),
                hisp.prop=c(.2,.2125,.2975,.2975),
                other.prop=c(.15,.225,.26,.2925),
                prob.hh.hisp=c(.18,.225,.14,.14),
                prob.hh.other=c(.025,.0175,.035,.04),
                prob.age=c(.3575,.55))
  samp=subset(samp,v.num==1)
  samp$dat.num=i
  samp$solnas = F; samp$solnas[sample(x=nrow(samp),size=450,replace=F)] = T
  
  ### Sorting
  samp = samp[order(samp$strat,samp$BGid,samp$hhid,samp$subid),]
  
  save(samp,file=paste0('./RData/',sampname,"_",str_pad(i,nchar(S),pad=0),".RData"),compress = 'xz')
}



