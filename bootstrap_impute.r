##   Copyright 2022:
##      Neil Ferguson, Imperial College London
##      Tommy Nyberg, University of Cambridge
##
##   Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

##############################
### Bootstrap if required ####
##############################

if(do_bootstrap) { # bootstrap data
  n=nrow(d2all)
  d2b=d2all[sample.int(n, n, replace = TRUE),]

} else {
  d2b=d2all
}


d2b$agebf=as.factor(d2b$ageband1) # 10 year bands
d2b$age_ib=d2b$age-d2b$ageband1 # age within ageband (subtract off lowest age)
d2b$Tgf=as.factor(d2b$Tg)

if(do_impute) { # impute missing VOC
  d2s=d2b %>% filter(omicron==1 | delta==1)

  # voc.model=glm(omicron ~ agebf + agebf:age_ib + vaccinated:reinfection_flag + 
  #                 ethnicity_group2 + imd_decile + sex + dose_vacc + 
  #                 hosp + hosp4 + hospA + dead + Tgf + NHSER_name,
  #               data = d2s, family = binomial)
  
  voc.model=parglm(omicron ~ agebf + agebf:age_ib + vaccinated:reinfection_flag + 
                  ethnicity_group2 + imd_decile + sex + dose_vacc + 
                  hosp + hosp4 + hospA + dead + Tgf + NHSER_name,
                data = d2s, family = binomial,
                control=parglm.control(epsilon = 1e-08, maxit = 25, trace = FALSE,
                 nthreads = 20L, block_size = NULL, method = "FAST"))
  

  # d2s$om_pred = predict(voc.model,type = "response")
  # voc.model.roc = roc(d2s$omicron ~ d2s$om_pred)
  # print(paste0("AUC=",voc.model.roc$auc))
  
  d2b$om_pred = predict(voc.model,newdata=d2b,type = "response")
  d2b$om_pred[d2b$omicron==1]=1
  d2b$om_pred[d2b$delta==1]=0
  d2b$om_rand=runif(nrow(d2b))
  d2b$omicron=ifelse(d2b$om_rand<=d2b$om_pred,1,0)
  d2b$delta=1-d2b$omicron
}

## retain cases which are known to be Omicron or Delta
d2=d2b %>% filter(omicron==1 | delta==1)
print(paste0("Dropping non Omicron/Delta:",nrow(d2)))

#d2all$voc=ifelse(d2all$omicron==1,1,ifelse(d2all$delta==1,0,0.5))

# d2all2=d2all %>% filter(last_specimen_date >= as.Date("2021-12-14") & last_specimen_date <= as.Date("2021-12-30"))
# table(d2all2$num_doses,d2all2$voc,d2all2$reinfection_flag)
# 
# table(d2all2$num_doses[d2all2$hospA==1],d2all2$voc[d2all2$hospA==1],d2all2$reinfection_flag[d2all2$hospA==1])

## more processing of selected data to add strata, agebands etc

if(TRUE) {

  d2$dose_type=droplevels(d2$dose_type)
  d2$id9=droplevels(d2$id9)
  d2$id10=droplevels(d2$id10)
  d2$id11=droplevels(d2$id11)
  d2$id12=droplevels(d2$id12)
  d2$id13=droplevels(d2$id13)
  d2$id14=droplevels(d2$id14)
  d2$id15=droplevels(d2$id15)
  
  d2$om=as.numeric(d2$omicron)
  d2$dlt=as.numeric(d2$delta)
  d2$hosp.om=d2$hosp*d2$om
  d2$hosp.dlt=d2$hosp*d2$dlt
  
  # compute proportion of data retained for different stratifications

  if(compute_retained) {
    var.list=c("id9","id10","id11","id12","id13","id14","id15")
    
    # all ages
    prop.used=NULL
    for(i in var.list) {
      sum.d2=d2 %>% group_by(!!as.name(i)) %>% summarise(hosp.om=sum(hosp.om),hosp.dlt=sum(hosp.dlt),om=sum(om),dlt=sum(dlt)) %>%
        filter(om >0 & dlt >0)
      prop=c(sum(sum.d2$om)/sum(d2$om),
                  sum(sum.d2$dlt)/sum(d2$dlt),
                  sum(sum.d2$hosp.om)/sum(d2$hosp.om),
                  sum(sum.d2$hosp.dlt)/sum(d2$hosp.dlt))
      prop.used=rbind(prop.used,prop)
    }
    rownames(prop.used)=c("10y","<1,1-4, 5-9, 10y","10y+vacc","<1,1-4, 5-9, 10y+vacc","10y, UTLA","day only","10y no eth")
    colnames(prop.used)=c("Omicron cases","Delta cases","Omicron hosp","Delta hosp")
    print(prop.used)
    
    # under 10s
    prop.used.lt10=NULL
    for(i in var.list) {
      sum.d2=d2 %>% group_by(!!as.name(i), ageband1) %>% summarise(hosp.om=sum(hosp.om),hosp.dlt=sum(hosp.dlt),om=sum(om),dlt=sum(dlt)) %>%
        filter(om >0 & dlt >0 & ageband1==0)
      prop=c(sum(sum.d2$om)/sum(d2$om[d2$ageband1==0]),
             sum(sum.d2$dlt)/sum(d2$dlt[d2$ageband1==0]),
             sum(sum.d2$hosp.om)/sum(d2$hosp.om[d2$ageband1==0]),
             sum(sum.d2$hosp.dlt)/sum(d2$hosp.dlt[d2$ageband1==0]))
      prop.used.lt10=rbind(prop.used.lt10,prop)
    }
    rownames(prop.used.lt10)=c("10y","<1,1-4, 5-9, 10y","10y+vacc","<1,1-4, 5-9, 10y+vacc","10y, UTLA","day only","10y no eth")
    colnames(prop.used.lt10)=c("Omicron cases","Delta cases","Omicron hosp","Delta hosp")
    print(prop.used.lt10)
  }
  
  # define vaccination-VOC interactions 
  
  d2$om_vacc3=interaction(d2$omicron,d2$dose_type)
  d2$om_vacc3=fct_collapse(d2$om_vacc3,"Unvaccinated"=c("0.None","1.None"))
  d2$om_vacc3=relevel(d2$om_vacc3,"Unvaccinated")
  
  d2$om_vacc4=interaction(d2$omicron,d2$dose_type)
  d2$om_vacc4=relevel(d2$om_vacc4,"0.None")
  
  d2$om_age=interaction(d2$omicron,d2$ageband)
  d2$om_age=fct_collapse(d2$om_age,"Delta"=c("0.[0,1)","0.[1,5)","0.[5,10)","0.[10,20)","0.[20,30)","0.[30,40)","0.[40,50)","0.[50,60)","0.[60,70)","0.[70,80)","0.[80,Inf]"))
  d2$om_age=relevel(d2$om_age,"Delta")
  d2$om_age3=d2$om_age
  
  d2$om_age1=fct_collapse(d2$om_age,"1.[0,10)"=c("1.[0,1)","1.[1,5)","1.[5,10)"))
  
  d2$om_age2=interaction(d2$omicron,d2$ageband2) # 0-4, 5-9 then 10 year bands
  d2$om_age2=fct_collapse(d2$om_age2,"Delta"=c("0.0","0.5","0.10","0.20","0.30","0.40","0.50","0.60","0.70","0.80"))
  d2$om_age2=relevel(d2$om_age2,"Delta")

  d2$om_dose=interaction(d2$omicron,d2$dose_vacc,d2$lv2s)
  d2$om_dose=droplevels(d2$om_dose)
  d2$om_dose=fct_collapse(d2$om_dose,"Unvaccinated"=c("0.None.8","1.None.8"))
  d2$om_dose=relevel(d2$om_dose,"Unvaccinated")
  
  d2$om_dose2=interaction(d2$omicron,d2$dose_vacc)
  d2$om_dose2=droplevels(d2$om_dose2)
  d2$om_dose2=fct_collapse(d2$om_dose2,"Unvaccinated"=c("0.None","1.None"))
  d2$om_dose2=relevel(d2$om_dose2,"Unvaccinated")

}

## processing to adjust for under ascertainment of reinfection
if(TRUE) {
  
  # estimates of prop of first infections reported from seb funk - <1, 1-4, 5-9, then 10-year
  # these use total first cases up to 28/11/2021 and his estimates of total numbers with 1+ infections to that point
  
  prop_rep=c(0.155, 0.155, 0.185, 0.366, 0.418, 0.437, 0.411, 0.386, 0.366, 0.343, 0.343)
  prop_rep_sd=c(0.0037, 0.0037, 0.0042, 0.0074, 0.0104, 0.0132, 0.0127, 0.0117, 0.0112, 0.0114, 0.0114)
  
  # From UKHSA Roche N seroprevalence data
  # prop_rep=c(0.282, 0.282, 0.313, 0.529, 0.675, 0.732, 0.751, 0.708, 0.741, 0.907, 0.907)
  # prop_rep_sd=c(0.0037, 0.0037, 0.0042, 0.0074, 0.0104, 0.0132, 0.0127, 0.0117, 0.0112, 0.0114, 0.0114)
  
  prop_rep_rand=rnorm(length(prop_rep),prop_rep,prop_rep_sd)  

  d2$prop_rep=prop_rep_rand[d2$ageband]
  d2$reinf_resamp_prob=1/d2$prop_rep-1
  if(do_model_reinf_impute) # use logistic regression to impute reinfection status [not used for paper]
  {
    reinf.model=parglm(reinfection_flag ~ omicron + dose_vacc + agebf + agebf:age_ib + 
                         ethnicity_group2 + imd_decile + sex +  
                      hosp + hosp4 + hospA + dead + Tgf + NHSER_name,
                    data = d2, family = binomial,
                    control=parglm.control(epsilon = 1e-08, maxit = 25, trace = FALSE,
                    nthreads = 20L, block_size = NULL, method = "FAST"))
    
    d2$reinf_pred = predict(reinf.model,type = "response")
    # reinf.model.roc = roc(d2$reinfection_flag ~ d2$reinf_pred)
    # print(paste0("Reinfection model AUC=",reinf.model.roc$auc))
  
    d2$reinf_pred=d2$reinf_pred/(1-d2$reinf_pred)*d2$reinf_resamp_prob
    d2$reinf_pred[d2$reinf_pred>=1|d2$reinfection_flag==1]=1
    
    d2$reinf_rand=runif(nrow(d2))
    d2$reinf_rand=ifelse(d2$reinf_rand<=d2$reinf_pred,1,0)
    d2$reinf_rand[d2$reinfection_flag==1]=2
    d2$reinf_rand=as.factor(d2$reinf_rand)
    d2$reinf_rand=relevel(d2$reinf_rand,"0")
    d2$vacc_reinf=interaction(d2$vaccinated,d2$reinf_rand)
    d2$vacc_reinf=fct_collapse(d2$vacc_reinf,"first"=c("Unvacc.0","Vacc.0"))
    d2$vacc_reinf=relevel(d2$vacc_reinf,"first")
    
  } else { # stratified approach to reinfection under ascertainment, as used in paper
    
    d2$rand=runif(nrow(d2)) # uniformly distributed random numbers between 0 and 1
    d2s=d2
    
    ## need to do imputation separately for each endpoint, given we stratify by endpoint in the imputation
    
    # hosp
    if("hosp" %in% endpoints) {
      # stratify by age, specimen date, simple vaccination group, variant and endpoint
      d2s$idr=interaction(d2s$ageband,d2s$Tg,d2s$dose_vacc,d2s$omicron,d2s$hosp) # <1, 1-4, 5-9, then 10 year bands + vacc
      d2s$idr=droplevels(d2s$idr)
      # calculate numbers of reinfections and total case numbers in each stratum
      d2.sum = d2s %>% group_by(idr) %>% summarise(idr_hosp_n=n(),idr_hosp_reinf=sum(reinfection_flag))
      # link back to main datatset
      d2s = left_join(d2s,d2.sum)
      # calculate probability each case is a reinfection. Always 1 for observed reinfections
      d2$hosp_reinf_var=ifelse(d2$reinfection_flag==1,1,d2$reinf_resamp_prob*d2s$idr_hosp_reinf/(d2s$idr_hosp_n-d2s$idr_hosp_reinf))
      # check for NA probabilities (none observed)
      print(nrow(d2[is.infinite(d2$hosp_reinf_var) | is.na(d2$hosp_reinf_var),]))
      print(nrow(d2[d2$hosp_reinf_var>1,]))
      # handle NA probabilities (but none observed)
      d2$hosp_reinf_var[is.infinite(d2$hosp_reinf_var) | is.na(d2$hosp_reinf_var)]=0
      # this line not really needed
      d2$hosp_reinf_var[d2$hosp_reinf_var>1]=1
      # sample 0/1 reinfection status using the estimated probability
      d2$hosp_reinf_rand=ifelse(d2$rand<=d2$hosp_reinf_var,1,0)
    }
    # hosp4
    if("hosp4" %in% endpoints) {
      d2s$idr=interaction(d2s$ageband,d2s$Tg,d2s$dose_vacc,d2s$omicron,d2s$hosp4) # <1, 1-4, 5-9, then 10 year bands + vacc
      d2s$idr=droplevels(d2s$idr)
      d2.sum = d2s %>% group_by(idr) %>% summarise(idr_hosp4_n=n(),idr_hosp4_reinf=sum(reinfection_flag))
      d2s = left_join(d2s,d2.sum) 
      d2$hosp4_reinf_var=ifelse(d2$reinfection_flag==1,1,d2$reinf_resamp_prob*d2s$idr_hosp4_reinf/(d2s$idr_hosp4_n-d2s$idr_hosp4_reinf))
      print(nrow(d2[is.infinite(d2$hosp4_reinf_var) | is.na(d2$hosp4_reinf_var),]))
      print(nrow(d2[d2$hosp4_reinf_var>1,]))
      d2$hosp4_reinf_var[is.infinite(d2$hosp4_reinf_var) | is.na(d2$hosp4_reinf_var)]=0
      d2$hosp4_reinf_var[d2$hosp4_reinf_var>1]=1
      d2$hosp4_reinf_rand=ifelse(d2$rand<=d2$hosp4_reinf_var,1,0)
    }
    # hospA
    if("hospA" %in% endpoints) {
      d2s$idr=interaction(d2s$ageband,d2s$Tg,d2s$dose_vacc,d2s$omicron,d2s$hospA) # <1, 1-4, 5-9, then 10 year bands + vacc
      d2s$idr=droplevels(d2s$idr)
      d2.sum = d2s %>% group_by(idr) %>% summarise(idr_hospA_n=n(),idr_hospA_reinf=sum(reinfection_flag))
      d2s = left_join(d2s,d2.sum)
      d2$hospA_reinf_var=ifelse(d2$reinfection_flag==1,1,d2$reinf_resamp_prob*d2s$idr_hospA_reinf/(d2s$idr_hospA_n-d2s$idr_hospA_reinf))
      print(nrow(d2[is.infinite(d2$hospA_reinf_var) | is.na(d2$hospA_reinf_var),]))
      print(nrow(d2[d2$hospA_reinf_var>1,]))
      d2$hospA_reinf_var[is.infinite(d2$hospA_reinf_var) | is.na(d2$hospA_reinf_var)]=0
      d2$hospA_reinf_var[d2$hospA_reinf_var>1]=1
      d2$hospA_reinf_rand=ifelse(d2$rand<=d2$hospA_reinf_var,1,0)
    }
    # dead
    if("dead" %in% endpoints) {
      d2s$idr=interaction(d2s$ageband,d2s$Tg,d2s$dose_vacc,d2s$omicron,d2s$dead) # <1, 1-4, 5-9, then 10 year bands + vacc
      d2s$idr=droplevels(d2s$idr)
      d2.sum = d2s %>% group_by(idr) %>% summarise(idr_dead_n=n(),idr_dead_reinf=sum(reinfection_flag))
      d2s = left_join(d2s,d2.sum)
      d2$dead_reinf_var=ifelse(d2$reinfection_flag==1,1,d2$reinf_resamp_prob*d2s$idr_dead_reinf/(d2s$idr_dead_n-d2s$idr_dead_reinf))
      print(nrow(d2[is.infinite(d2$dead_reinf_var) | is.na(d2$dead_reinf_var),]))
      print(nrow(d2[d2$dead_reinf_var>1,]))
      d2$dead_reinf_var[is.infinite(d2$dead_reinf_var) | is.na(d2$dead_reinf_var)]=0
      d2$dead_reinf_var[d2$dead_reinf_var>1]=1
      d2$dead_reinf_rand=ifelse(d2$rand<=d2$dead_reinf_var,1,0)
    }
  }
  
}

