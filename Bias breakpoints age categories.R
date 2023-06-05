#### generate random dataset resembling Million et al. as regards age, death, and HCQ-AZ ####

dbr = data.frame(
  HCQAZ = c(rep(TRUE,23172),rep(FALSE,7030)),
  age_num = c(rnorm(n=23172,mean=47.0,sd=16.1),rnorm(n=7030,mean=54.6,sd=19.0))
)
#make sure everyone is between 18 and 100 y.o.
while(any(dbr$age_num<18 | dbr$age_num>100)){
  for(i in which(dbr$age_num<18 | dbr$age_num>100)) {
    if (dbr$HCQAZ[i]) {
      dbr$age_num[i] <- rnorm(n=1,mean=47.0,sd=16.1)
    } else {
      dbr$age_num[i] <- rnorm(n=1,mean=54.6,sd=19.0)
    }
  }
}
#cut in 2 types of age categories, either broad, like in Million et al. with >20y range
dbr$AGE = cut(dbr$age_num,breaks=c(18,50,70,90,100),include.lowest = T,right = F,labels=c(1:4))
prop.table(with(dbr,table(AGE,HCQAZ)),2)*100 #remarkably similar to Million et al.
#or in 5-y range
dbr$AGE_group5y = cut(dbr$age_num,breaks=c(18,seq(50,90,5),100),include.lowest = T,right = F)

# generate mortality depending on age from Supplementary Figure 2. Mortality rate according to age (n = 30,423)
dbr$DEATH = ifelse(dbr$age_num>90,rbinom(n=sum(dbr$age_num>90), size=1, prob=0.3278),0)
agegp = dbr$age_num<=90 & dbr$age_num>85
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.2404),dbr$DEATH)
agegp = dbr$age_num<=85 & dbr$age_num>80
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.1664),dbr$DEATH)
agegp = dbr$age_num<=80 & dbr$age_num>75
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.0842),dbr$DEATH)
agegp = dbr$age_num<=75 & dbr$age_num>70
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.0385),dbr$DEATH)
agegp = dbr$age_num<=70 & dbr$age_num>65
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.0200),dbr$DEATH)
agegp = dbr$age_num<=65 & dbr$age_num>60
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.0109),dbr$DEATH)
agegp = dbr$age_num<=60 & dbr$age_num>55
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.0041),dbr$DEATH)
agegp = dbr$age_num<=55 & dbr$age_num>50
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.0009),dbr$DEATH)
agegp = dbr$age_num<=50 & dbr$age_num>45
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.0022),dbr$DEATH)
agegp = dbr$age_num<=25
dbr$DEATH = ifelse(agegp,rbinom(n=sum(agegp), size=1, prob=0.0004),dbr$DEATH)
prop.table(with(dbr,table(DEATH,HCQAZ)),2)*100 #differ from Million, of course, because of the other biases (or the effect of treatment, if any)!!

# Effectiveness of treatment
model <- glm(DEATH ~ factor(AGE) + HCQAZ, data=dbr, family="binomial")
summary(model)
library("questionr")
odds.ratio(model) #significant OR 0.74 [0.60- 0.91], p = 0.005

# while model with age in 5-years age categories is [of course] non-significant
model <- glm(DEATH ~ AGE_group5y + HCQAZ, data=dbr, family="binomial")
odds.ratio(model) #non-significant OR 0.90 [0.73- 1.12], p = 0.356

#but also with age in 10-y categ
model <- glm(DEATH ~ AGE_group5y + HCQAZ, data=dbr, family="binomial")
odds.ratio(model) #non-significant OR 0.93 [0.76-1.14], p=	0.502

# and the same with age in numeric 
model <- glm(DEATH ~ age_num + HCQAZ, data=dbr, family="binomial")
odds.ratio(model) #non-significant OR 0.92 [0.74-1.14], p = 0.442
