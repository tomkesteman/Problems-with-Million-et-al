#### R script describing analyse of dataset IHU (Million et al, 2023) using hospitalization as an outcome
#### 05 Jun 2023 - Thomas Kesteman 

#### load database ####
# download database from Science Data Bank https://doi.org/10.57760/sciencedb.07803
db <- read.delim("Database_Cohort_30423_COVID-19_IHU_032923.txt")
colnames(db)
#### create variables treatment ####
#same as HCQ in logical
db$HCQAZ = ifelse(is.na(db$HCQ) | is.na(db$AZ),NA,db$HCQ %in% 1 & db$AZ %in% 1) 
# compare HCQ+AZ against standard of care
db$HCQAZ_vs_SOC = ifelse(db$HCQ %in% 0 & db$AZ %in% 0 & db$IVM %in% 0,F,
                         ifelse(db$HCQ %in% 1 & db$AZ %in% 1,T,NA ))
# compare 'other' treatments against standard of care
db$Other_vs_SOC = ifelse(db$HCQ %in% 0 & db$AZ %in% 0 & db$IVM %in% 0,F,
                         ifelse((db$HCQ %in% 1 | db$AZ %in% 1 | db$IVM %in% 1) & !(db$HCQ %in% 1 & db$AZ %in% 1),T,NA ))
# compare HCQ+AZ against 'other' treatments and against standard of care
db$Treatment = ifelse(db$HCQ %in% 0 & db$AZ %in% 0 & db$IVM %in% 0,"SOC",
                      ifelse(db$HCQ %in% 1 & db$AZ %in% 1,"HCQ-AZ","Other" ))
db$Treatment = ifelse(is.na(db$HCQ) | is.na(db$AZ),NA,db$Treatment)
db$Treatment = relevel(factor(db$Treatment),"SOC")

#### Debugging table 1 ####
# variable OPD/IPD
with(db,table(OUTPATIENT,INPATIENT,useNA = "ifany")) #explains why counts on Table 1 don't add up to 30423
db$OPD = ifelse(db$OUTPATIENT %in% 1,"OPD",NA)
db$OPD = ifelse(db$INPATIENT %in% 1, "IPD",db$OPD)
db$OPD = ifelse(db$OUTPATIENT %in% 1 & db$INPATIENT %in% 1,"OPD-IPD",db$OPD)
db$OPD = relevel(factor(db$OPD),"OPD")
with(db,table(OPD,Treatment,useNA = "ifany")) 

#### strange pval with B.1.7.7
with(db,table(VARIANT,HCQAZ,useNA = "ifany"))
with(db,table(VARIANT %in% "B.1.7.7",HCQAZ))
chisq.test(with(db,table(VARIANT %in% "B.1.7.7",HCQAZ)))

#### Multivariate model in Table 2 ####
model <- glm(DEATH ~ SEX + factor(AGE) + factor(PERIOD) + OUTPATIENT + HCQAZ, data=db, family="binomial")
library("questionr")
odds.ratio(model) 
#HCQ: OR 0.50 [0.41 - 0.60] --> reproduces roughly Million's finding (0.55 [0.45-0.68]) 
#Differences explained by use of other software???

#### Effect of treatment on Hospitalization of OPD patients ####
# 1. No adjustment on comorbidity
# 1.1 All
model <- glm(INPATIENT ~ SEX + factor(AGE)  + factor(PERIOD) + Treatment    , 
             data=db[db$OUTPATIENT %in% 1,], family="binomial")
odds.ratio(model)

# 1.2 By age groups
do.call(rbind,by(db,db$AGE,
                 FUN = function(x){
                   model <- glm(INPATIENT ~ Treatment + SEX + factor(PERIOD), data=x[x$OUTPATIENT %in% 1,], family="binomial")
                   odds.ratio(model)[2,] #uncomment for HCQ+AZ 
                   # odds.ratio(model)[3,] #uncomment for 'HCQ+AZ 'other' treatments
                   }))

# 1.3 By sex
do.call(rbind,by(db,db$SEX,
                 FUN = function(x){
                   model <- glm(INPATIENT ~ Treatment + factor(AGE) + factor(PERIOD), data=x[x$OUTPATIENT %in% 1,], family="binomial")
                   odds.ratio(model)[2,] #uncomment for HCQ+AZ 
                   # odds.ratio(model)[3,] #uncomment for 'HCQ+AZ 'other' treatments
                 }))

# 1.4 By period
do.call(rbind,by(db,db$PERIOD,
                 FUN = function(x){
                   model <- glm(INPATIENT ~ Treatment + factor(AGE) + SEX, data=x[x$OUTPATIENT %in% 1,], family="binomial")
                   odds.ratio(model)[2,] #uncomment for HCQ+AZ 
                   # odds.ratio(model)[3,] #uncomment for 'HCQ+AZ 'other' treatments
                 }))


# 2. Adjustment on comorbidity
# + exclude ChronicCardiacDiseases (indication bias)
# 2.1 All
model <- glm(INPATIENT ~ SEX + factor(AGE)  + factor(PERIOD) +
               VACCINATION + OBESITY + HBP + DIABETE + ASTHMA + COPD + CANCER + IMMUNODEFICIENCY + AutoImmuneDiseases +
               Treatment , 
             data=db[db$ChronicCardiacDiseases %in% 0 & db$OUTPATIENT %in% 1,], family="binomial")
odds.ratio(model)


# 2.2 By age groups
do.call(rbind,by(db,db$AGE,
                 FUN = function(x){
                   model <- glm(INPATIENT ~ Treatment + SEX + factor(PERIOD) +
                                  VACCINATION + OBESITY + HBP + DIABETE + ASTHMA + COPD + CANCER + IMMUNODEFICIENCY + AutoImmuneDiseases, 
                                data=x[x$ChronicCardiacDiseases %in% 0 & x$OUTPATIENT %in% 1,], family="binomial")
                   odds.ratio(model)[2,] #uncomment for HCQ+AZ 
                   # odds.ratio(model)[3,] #uncomment for 'HCQ+AZ 'other' treatments
                 }))

# 2.3 By sex
do.call(rbind,by(db,db$SEX,
                 FUN = function(x){
                   model <- glm(INPATIENT ~ Treatment + factor(AGE) + factor(PERIOD) +
                                  VACCINATION + OBESITY + HBP + DIABETE + ASTHMA + COPD + CANCER + IMMUNODEFICIENCY + AutoImmuneDiseases, 
                                data=x[x$ChronicCardiacDiseases %in% 0 & x$OUTPATIENT %in% 1,], family="binomial")
                   odds.ratio(model)[2,] #uncomment for HCQ+AZ 
                   # odds.ratio(model)[3,] #uncomment for 'HCQ+AZ 'other' treatments
                 }))

# 2.4 By period
# data not available in periods 1-3
prop.table(with(db[db$OUTPATIENT %in% 1,],table(PERIOD,ChronicCardiacDiseases)),1)*100
do.call(rbind,by(db[db$PERIOD>3,],db$PERIOD[db$PERIOD>3],
                 FUN = function(x){
                   model <- glm(INPATIENT ~ Treatment + factor(AGE) + SEX +
                                  VACCINATION + OBESITY + HBP + DIABETE + ASTHMA + COPD + CANCER + IMMUNODEFICIENCY + AutoImmuneDiseases, 
                                data=x[x$ChronicCardiacDiseases %in% 0 & x$OUTPATIENT %in% 1,], family="binomial")
                   odds.ratio(model)[2,] #uncomment for HCQ+AZ 
                   # odds.ratio(model)[3,] #uncomment for 'HCQ+AZ 'other' treatments
                 }))
