####### generate HCHS population #######
library(data.table)

# Set the seed
set.seed(123)

# Number of BG in the population
Nbg=376

# Cumulative sum of BG in the population across all 4 strata
bgstrat = cumsum(c(58,21,130,167))

# Number of HH per BG, at least 1 HH per BG, mean number of HHs/BG=451 (this is the number of ALL HHs in each BG including non-eligible HHs)
bg.size.temp=1+round(rexp(Nbg,1/450))    

# Number of eligible HHs with Hispanic surname in each BG
# The proportion of eligible HHs with Hispanic surname per stratum is 0.2, 0.2125, 0.2975, and 0.2975
num.hisp.strat=round(c(bg.size.temp[1:bgstrat[1]]*.2,bg.size.temp[(bgstrat[1]+1):bgstrat[2]]*.2125,bg.size.temp[(bgstrat[2]+1):bgstrat[3]]*.2975,bg.size.temp[(bgstrat[3]+1):bgstrat[4]]*.2975)) 

# Number of eligible HHs with Other surname in each BG
# The proportion of eligible HHs with Hispanic surname per stratum is 0.15, 0.225, 0.26, and 0.2925
num.other.strat=round(c(bg.size.temp[1:bgstrat[1]]*.15,bg.size.temp[(bgstrat[1]+1):bgstrat[2]]*.225,bg.size.temp[(bgstrat[2]+1):bgstrat[3]]*.26,bg.size.temp[(bgstrat[3]+1):bgstrat[4]]*.2925))  

# Number of eligible HHs in each BG
# Note that bg.size < bg.size.temp, as the latter contains non-eligible HHs
bg.size=num.hisp.strat+num.other.strat

# Number of HHs in target population
Nhh=sum(bg.size)                  

# HH size, at least 1 subject per HH, mean number of subjects/HH=2
hh.size=1+rpois(Nhh,1)   

# Number of subjects in target population
N=sum(hh.size)                          

# all ID's unique (e.g., subid=k only for one subject within one HH)
BGid=rep(rep(rep(1:Nbg, times=bg.size), times=hh.size), times=2)  
hhid=rep(rep(1:Nhh, times=hh.size), times=2)
subid=rep(1:N, times=2)
v.num=rep(c(1,2),each=N)
strat=rep(NA,times=Nbg); strat[BGid<=bgstrat[1]]=1; strat[BGid>bgstrat[1] & BGid<=bgstrat[2]]=2; strat[BGid>bgstrat[2] & BGid<=bgstrat[3]]=3; strat[BGid>bgstrat[3] & BGid<=bgstrat[4]]=4


# Now creating hips.strat variable
# It creates a matrix with Nbg columns (BGs) and it fills up with TRUE (if hispanic and in target population),
# FALSE (if other surname and in target population), and NA (if otherwise)
A=matrix(NA,nrow=max(bg.size),ncol=length(bg.size))
for (i in 1:length(bg.size)){
    A[,i]=c(rep(TRUE,times=num.hisp.strat[i]),rep(FALSE,times=num.other.strat[i]),rep(NA,times=max(bg.size)-bg.size[i]))    # create matrix A with each column corresponding to a BG (containing 1's for each Hispanic HH followed by 0's for each other HH)
}
# Indicator for Hispanic surname (na.omit(c(A)) gives the household-level hispanic surname)
hisp.strat=rep(rep(c(na.omit(c(A))),times=hh.size),times=2)

pop = data.table(strat,BGid,hhid,hisp.strat,subid,v.num)
pop = pop[order(strat,BGid,hhid,hisp.strat,subid,v.num),]

# Sanity check: Check if there is any household with multiple distinct hispanic surname (should be 0 or 1)
pop[,.(Check = mean(hisp.strat)%in%c(0,1)),by='hhid'][,all(Check)]
# Check the proportion of hispanic surname per Strat (should be somewhat close to 0.57, 0.48, 0.53, 0.50)
pop[,mean(hisp.strat),by='strat']
# Check how many BGid per stratum we have (it should be 58, 21, 130, and 167)
pop[,unique(BGid),by='strat'][,.N,by='strat']

# Saving output
if(!dir.exists('./RData')){system('mkdir RData')}
save(pop,file='./RData/TargetPopulation.RData',compress = 'xz')


