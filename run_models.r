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

######################
### Run models     ###
######################

d2c=d2
datacut_c=datacut
### select only unvaccinated if required
if(do_unvacc) {
  d2=d2c %>% filter(vaccinated=="Unvacc")
  datacut=paste0(datacut_c,"_unvacc")
  d2$dose_lv=droplevels(d2$dose_lv)
  d2$dose_vacc=droplevels(d2$dose_vacc)
  d2$vaccinated=droplevels(d2$vaccinated)
}

if(do_impute) datacut=paste0(datacut,"_imp")

#### delays by age

if(do_delays) {

### onset to specimen date

  d2$omicron_n=ifelse(d2$omicron==1,"Omicron","Delta")

  d2$onset2sp=as.numeric(d2$last_specimen_date-d2$Onsetdate)
  d2.on2sp=d2[!is.na(d2$onset2sp) & d2$onset2sp>=0 & d2$reinfection_flag==0 & d2$onset2sp<=14,]
  onset2sp.tb = d2.on2sp %>% 
    group_by(ageband1,omicron_n,hosp) %>% 
    summarise(av=mean(onset2sp),se=sd(onset2sp)/sqrt(n()-1),n=n()) %>%
    pivot_wider(names_from = c(omicron_n,hosp), names_sep = ".",values_from = c(av,se,n))
  onset2sp.tb$endpoint="hosp"
  
  onset2sp.tbA = d2.on2sp %>% 
    group_by(ageband1,omicron_n,hospA) %>% 
    summarise(av=mean(onset2sp),se=sd(onset2sp)/sqrt(n()-1),n=n()) %>%
    pivot_wider(names_from = c(omicron_n,hospA), names_sep = ".",values_from = c(av,se,n))
  onset2sp.tbA$endpoint="hospA"
  
  onset2sp.tb=bind_rows(onset2sp.tb,onset2sp.tbA)
  write.csv(onset2sp.tb,paste0("Omicron/Onset2Test_10y_",datacut,"_",date_name,".csv"))
  
  ### onset to hosp
  
  d2$onset2hosp=as.numeric(d2$hosp_date-d2$Onsetdate)
  ## hosp
  d2.on2hosp=d2[d2$hosp==1 & !is.na(d2$onset2hosp) & d2$onset2sp>=0 & d2$onset2sp<=14 & d2$reinfection_flag==0 ,]
  onset2hosp.tb = d2.on2hosp %>% 
    group_by(ageband1,omicron_n) %>% 
    summarise(av=mean(onset2hosp),se=sd(onset2hosp)/sqrt(n()-1),n=n()) %>%
    pivot_wider(names_from = c(omicron_n), names_sep = ".",values_from = c(av,se,n))
  onset2hosp.tb$endpoint="hosp"
  
  d2$onset2hosp=as.numeric(d2$hospA_date-d2$Onsetdate)
  d2.on2hosp=d2[d2$hospA==1 & !is.na(d2$onset2hosp) & d2$onset2sp>=0 & d2$onset2sp<=14 & d2$reinfection_flag==0 ,]
  onset2hosp.tbA = d2.on2hosp %>% 
    group_by(ageband1,omicron_n) %>% 
    summarise(av=mean(onset2hosp),se=sd(onset2hosp)/sqrt(n()-1),n=n()) %>%
    pivot_wider(names_from = c(omicron_n), names_sep = ".",values_from = c(av,se,n))
  onset2hosp.tbA$endpoint="hospA"
  
  onset2sp.tb=bind_rows(onset2hosp.tb,onset2hosp.tbA)
  write.csv(onset2sp.tb,paste0("Omicron/Onset2Hosp_",datacut,"_",date_name,".csv"))
  
  
  ### LoS
  d2$los=d2$hosp_out_date - d2$hosp_date
  d2.los=d2[d2$hosp==1 & !is.na(d2$los) & d2$los>=0 & d2$los<=7,]
  los.tb = d2.los %>% 
    group_by(ageband1,omicron_n) %>% 
    summarise(av=mean(los),se=sd(los)/sqrt(n()-1),n=n()) %>%
    pivot_wider(names_from = c(omicron_n), names_sep = ".",values_from = c(av,se,n))
  los.tb$endpoint="hosp"
  
  d2$los=d2$hospA_out_date - d2$hospA_date
  d2.los=d2[d2$hospA==1 & !is.na(d2$los) & d2$los>=0 & d2$los<=7,]
  los.tbA = d2.los %>% 
    group_by(ageband1,omicron_n) %>% 
    summarise(av=mean(los),se=sd(los)/sqrt(n()-1),n=n()) %>%
    pivot_wider(names_from = c(omicron_n), names_sep = ".",values_from = c(av,se,n))
  los.tbA$endpoint="hospA"
  
  table(d2$los,d2$hospE)
  
  los.tb=bind_rows(los.tb,los.tbA)
  write.csv(los.tb,paste0("Omicron/LoS_lt7_",datacut,"_",date_name,".csv"))
  
  ### specimen to hosp
  d2$sp2hosp=as.numeric(d2$hosp_date-d2$last_specimen_date)
  d2.sp2hosp=d2[d2$hosp==1 & !is.na(d2$sp2hosp),]
  sp2hosp.tb = d2.sp2hosp %>% 
    group_by(ageband1,omicron_n) %>% 
    summarise(av=mean(sp2hosp),se=sd(sp2hosp)/sqrt(n()-1),n=n()) %>%
    pivot_wider(names_from = c(omicron_n), names_sep = ".",values_from = c(av,se,n))
  sp2hosp.tb$endpoint="hosp"
  
  d2$sp2hosp=as.numeric(d2$hospA_date-d2$last_specimen_date)
  d2.sp2hosp=d2[d2$hospA==1 & !is.na(d2$sp2hosp),]
  sp2hosp.tbA = d2.sp2hosp %>% 
    group_by(ageband1,omicron_n) %>% 
    summarise(av=mean(sp2hosp),se=sd(sp2hosp)/sqrt(n()-1),n=n()) %>%
    pivot_wider(names_from = c(omicron_n), names_sep = ".",values_from = c(av,se,n))
  sp2hosp.tbA$endpoint="hospA"
  
  sp2hosp.tb=bind_rows(sp2hosp.tb,sp2hosp.tbA)
  write.csv(sp2hosp.tb,paste0("Omicron/Sp2Hosp_",datacut,"_",date_name,".csv"))
}

### numbers of Delta and Omicron cases and hosp by age (including unadjusted HRs)

if(do_numbers) {
  
  tb1=as.data.frame.matrix(table(d2all$last_specimen_date[d2all$hosp==1],d2all$ageband[d2all$hosp==1]))
  tb1$date=rownames(tb1)
  tb1$voc="all"
  tb1$hosp="hosp4"
  tb2=as.data.frame.matrix(table(d2all$last_specimen_date,d2all$ageband))
  tb2$date=rownames(tb2)
  tb2$voc="all"
  tb2$hosp="case"
  
  tb3=as.data.frame.matrix(table(d2$last_specimen_date[d2$hosp==1 & d2$omicron==1],d2$ageband[d2$hosp==1 & d2$omicron==1]))
  tb3$date=rownames(tb3)
  tb3$voc="omicron"
  tb3$hosp="hosp4"
  tb4=as.data.frame.matrix(table(d2$last_specimen_date[d2$omicron==1],d2$ageband[d2$omicron==1]))
  tb4$date=rownames(tb4)
  tb4$voc="omicron"
  tb4$hosp="case"
  
  tb5=as.data.frame.matrix(table(d2$last_specimen_date[d2$hosp==1 & d2$delta==1],d2$ageband[d2$hosp==1 & d2$delta==1]))
  tb5$date=rownames(tb5)
  tb5$voc="delta"
  tb5$hosp="hosp4"
  tb6=as.data.frame.matrix(table(d2$last_specimen_date[d2$delta==1],d2$ageband[d2$delta==1]))
  tb6$date=rownames(tb6)
  tb6$voc="delta"
  tb6$hosp="case"
  num.age.date=rbind(tb1,tb2,tb3,tb4,tb5,tb6)
  num.age.date=cbind(num.age.date[,c("date","voc","hosp")],num.age.date[,1:(ncol(num.age.date)-3)])
  rownames(num.age.date)=c()
  
  write.csv(num.age.date,paste0("Omicron/num_age_date_",datacut,"_",date_name,".csv"))
  
  col.names=c("case.delta","case.omicron",
              "hosp.delta","hosp.omicron",
              "hosp1.delta","hosp1.omicron",
              "hosp4.delta","hosp4.omicron",
              "hosp2.delta","hosp2.omicron",
              "hospE.delta","hospE.omicron",
              "hospF.delta","hospF.omicron",
              "hospA.delta","hospA.omicron",
              "dead.delta","dead.omicron")
  ## all cases
  num.all=as.data.frame(cbind(t(table(d2$omicron)),
                              t(table(d2$omicron[d2$hosp==1])),
                              t(table(d2$omicron[d2$hosp1==1])),
                              t(table(d2$omicron[d2$hosp4==1])),
                              t(table(d2$omicron[d2$hosp2==1])),
                              t(table(d2$omicron[d2$hospE==1])),
                              t(table(d2$omicron[d2$hospF==1])),
                              t(table(d2$omicron[d2$hospA==1])),
                              t(table(d2$omicron[d2$dead==1]))))
  colnames(num.all)=col.names
  rownames(num.all)[1]="All"
  
  ## by age
  
  num.age=as.data.frame(cbind(table(d2$ageband,d2$omicron),
                              table(d2$ageband[d2$hosp==1],d2$omicron[d2$hosp==1]),
                              table(d2$ageband[d2$hosp1==1],d2$omicron[d2$hosp1==1]),
                              table(d2$ageband[d2$hosp4==1],d2$omicron[d2$hosp4==1]),
                              table(d2$ageband[d2$hosp2==1],d2$omicron[d2$hosp2==1]),
                              table(d2$ageband[d2$hospE==1],d2$omicron[d2$hospE==1]),
                              table(d2$ageband[d2$hospF==1],d2$omicron[d2$hospF==1]),
                              table(d2$ageband[d2$hospA==1],d2$omicron[d2$hospA==1]),
                              table(d2$ageband[d2$dead==1],d2$omicron[d2$dead==1])))
  colnames(num.age)=col.names
  
  
  write.csv(num.age,paste0("Omicron/num_age_",datacut,"_",date_name,".csv"))
  
  ### numbers of Delta and Omicron cases and hosp by vacc status
  
  num.vacc.simp=as.data.frame(cbind(table(d2$dose_vacc,d2$omicron),
                                    table(d2$dose_vacc[d2$hosp==1],d2$omicron[d2$hosp==1]),
                                    table(d2$dose_vacc[d2$hosp1==1],d2$omicron[d2$hosp1==1]),
                                    table(d2$dose_vacc[d2$hosp4==1],d2$omicron[d2$hosp4==1]),
                                    table(d2$dose_vacc[d2$hosp2==1],d2$omicron[d2$hosp2==1]),
                                    table(d2$dose_vacc[d2$hospE==1],d2$omicron[d2$hospE==1]),
                                    table(d2$dose_vacc[d2$hospF==1],d2$omicron[d2$hospF==1]),
                                    table(d2$dose_vacc[d2$hospA==1],d2$omicron[d2$hospA==1]),
                                    table(d2$dose_vacc[d2$dead==1],d2$omicron[d2$dead==1])))
  colnames(num.vacc.simp)=col.names
  
  
  num.vacc=as.data.frame(cbind(table(d2$dose_lv,d2$omicron),
                               table(d2$dose_lv[d2$hosp==1],d2$omicron[d2$hosp==1]),
                               table(d2$dose_lv[d2$hosp1==1],d2$omicron[d2$hosp1==1]),
                               table(d2$dose_lv[d2$hosp4==1],d2$omicron[d2$hosp4==1]),
                               table(d2$dose_lv[d2$hosp2==1],d2$omicron[d2$hosp2==1]),
                               table(d2$dose_lv[d2$hospE==1],d2$omicron[d2$hospE==1]),
                               table(d2$dose_lv[d2$hospF==1],d2$omicron[d2$hospF==1]),
                               table(d2$dose_lv[d2$hospA==1],d2$omicron[d2$hospA==1]),
                               table(d2$dose_lv[d2$dead==1],d2$omicron[d2$dead==1])))
  colnames(num.vacc)=col.names
  
  
  ### numbers of Delta and Omicron cases and hosp by reinfection_flag
  
  d2$reinf_f=as.factor(d2$reinfection_flag)
  num.reinf=as.data.frame(cbind(table(d2$reinf_f,d2$omicron),
                                table(d2$reinf_f[d2$hosp==1],d2$omicron[d2$hosp==1]),
                                table(d2$reinf_f[d2$hosp1==1],d2$omicron[d2$hosp1==1]),
                                table(d2$reinf_f[d2$hosp4==1],d2$omicron[d2$hosp4==1]),
                                table(d2$reinf_f[d2$hosp2==1],d2$omicron[d2$hosp2==1]),
                                table(d2$reinf_f[d2$hospE==1],d2$omicron[d2$hospE==1]),
                                table(d2$reinf_f[d2$hospF==1],d2$omicron[d2$hospF==1]),
                                table(d2$reinf_f[d2$hospA==1],d2$omicron[d2$hospA==1]),
                                table(d2$reinf_f[d2$dead==1],d2$omicron[d2$dead==1])))
  colnames(num.reinf)=col.names
  
  rownames(num.reinf)=paste0("Reinfection:",rownames(num.reinf))
  
  
  ## by sex
  
  d2$sexf=as.factor(d2$sex)
  num.sex=as.data.frame(cbind(table(d2$sexf,d2$omicron),
                        table(d2$sexf[d2$hosp==1],d2$omicron[d2$hosp==1]),
                        table(d2$sexf[d2$hosp1==1],d2$omicron[d2$hosp1==1]),
                        table(d2$sexf[d2$hosp4==1],d2$omicron[d2$hosp4==1]),
                        table(d2$sexf[d2$hosp2==1],d2$omicron[d2$hosp2==1]),
                        table(d2$sexf[d2$hospE==1],d2$omicron[d2$hospE==1]),
                        table(d2$sexf[d2$hospF==1],d2$omicron[d2$hospF==1]),
                        table(d2$sexf[d2$hospA==1],d2$omicron[d2$hospA==1]),
                        table(d2$sexf[d2$dead==1],d2$omicron[d2$dead==1])))
  colnames(num.sex)=col.names
  

  ## by IMD decile
  
  num.imd=as.data.frame(cbind(table(d2$imd_decile,d2$omicron),
                              table(d2$imd_decile[d2$hosp==1],d2$omicron[d2$hosp==1]),
                              table(d2$imd_decile[d2$hosp1==1],d2$omicron[d2$hosp1==1]),
                              table(d2$imd_decile[d2$hosp4==1],d2$omicron[d2$hosp4==1]),
                              table(d2$imd_decile[d2$hosp2==1],d2$omicron[d2$hosp2==1]),
                              table(d2$imd_decile[d2$hospE==1],d2$omicron[d2$hospE==1]),
                              table(d2$imd_decile[d2$hospF==1],d2$omicron[d2$hospF==1]),
                              table(d2$imd_decile[d2$hospA==1],d2$omicron[d2$hospA==1]),
                              table(d2$imd_decile[d2$dead==1],d2$omicron[d2$dead==1])))
  colnames(num.imd)=col.names
  rownames(num.imd)=paste0("IMD:",rownames(num.imd))
  
  
  ### numbers of Delta and Omicron cases and hosp by vacc_reinf
  
  d2$vacc_reinf=interaction(d2$vaccinated,d2$reinf_f)
  
  num.vacc_reinf=as.data.frame(cbind(table(d2$vacc_reinf,d2$omicron),
                                     table(d2$vacc_reinf[d2$hosp==1],d2$omicron[d2$hosp==1]),
                                     table(d2$vacc_reinf[d2$hosp1==1],d2$omicron[d2$hosp1==1]),
                                     table(d2$vacc_reinf[d2$hosp4==1],d2$omicron[d2$hosp4==1]),
                                     table(d2$vacc_reinf[d2$hosp2==1],d2$omicron[d2$hosp2==1]),
                                     table(d2$vacc_reinf[d2$hospE==1],d2$omicron[d2$hospE==1]),
                                     table(d2$vacc_reinf[d2$hospF==1],d2$omicron[d2$hospF==1]),
                                     table(d2$vacc_reinf[d2$hospA==1],d2$omicron[d2$hospA==1]),
                                     table(d2$vacc_reinf[d2$dead==1],d2$omicron[d2$dead==1])))
  colnames(num.vacc_reinf)=col.names
  
  ## unadjusted HRs
  
  d2$all="all"
  d2$all=as.factor(d2$all)
  
  var.list=c("all","ageband","sexf","reinf_f","imd_decile","dose_vacc","vacc_reinf","dose_lv")
  endpoints=c("hosp","hosp4","hospA","dead")
  res.unadjust=NULL
  for(k in var.list) {
    lv=levels(d2[[eval(k)]])
    for(j in lv) {
      d2s=d2[d2[eval(k)]==j,]
      res.cox=NULL  
      for(i in endpoints) {
        d2s[,"surv_status"]=d2s[,eval(i)]+1
        if(i=="dead") {
          d2s$surv_date=d2s$last_specimen_date+28
        } else {
          d2s$surv_date=d2s$last_specimen_date+14
        }
        d2s$surv_date[d2s$surv_date>as.Date(hosp_data_date)]=as.Date(hosp_data_date)
        d2s[d2s$surv_status==2,"surv_date"]=d2s[d2s$surv_status==2,eval(paste0(i,"_date"))]
        d2s$surv_time=d2s$surv_date-d2s$last_specimen_date
        d2s$surv_time[d2s$surv_time<=0]=0.5 # small positive value
        cox.mod = coxph(Surv(surv_time, surv_status) ~ omicron, data = d2s)
        cox.mod.s=summary(cox.mod)
        res.cox1=as.data.frame(cox.mod.s$coefficients)
        
        res.cox=cbind(res.cox,paste0(round(exp(res.cox1$coef[1]),2),
                                     " (", round(exp(res.cox1$coef[1]-1.96*res.cox1$'se(coef)'[1]),2),
                      "-",round(exp(res.cox1$coef[1]+1.96*res.cox1$'se(coef)'[1]),2),")"))
        rownames(res.cox)[1]=paste0(k,":",j)
    
      }
      res.unadjust=rbind(res.unadjust,res.cox)
    }
  }
  
  colnames(res.unadjust)=c("HR(hosp)","HR(hosp4)","HR(hospA)","HR(dead)")

  res.unadjust
  
  # combine counts with unadjusted HRs
  num.tab=cbind(rbind(num.all,num.age,num.sex,num.reinf,num.imd,num.vacc.simp,num.vacc_reinf,num.vacc),
                res.unadjust)
  
  write.csv(num.tab,paste0("Omicron/num_tab_",datacut,"_",date_name,".csv"))
  
}

### run cox models 

for(k in 2:2) { # option to run in <10s only (not used)
  for(j in model_set) {  # 5 is the "central" model
    if(k==2) {
      d2s=d2
      agedesc="allage"
    } else {
      d2s=d2[d2$age<10,]
      agedesc="lt10"
    }
    d2s$om_age=d2s$om_age3 # <1, 1-4, 5-9, then 10 year bands
    form.str=" ~ om_age + vaccinated:reinfection_flag + strata(idc)"
    epb=0.0 # epidemic phase bias adjustment
    if (j==1) {
      stratadesc="sev_vacc_noage"
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ omicron + vaccinated:reinfection_flag + imd_decile + sex + om_dose + agebf:age_ib + strata(idc)"
    } else if (j==2) {
      stratadesc="sev_noage"
      d2s$idc=d2s$id11 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ omicron + vaccinated:reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    } else if (j==3) {
      stratadesc="sev_vacc"
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + om_dose + agebf:age_ib + strata(idc)"
    } else if (j==4) {
      stratadesc="sev_vacc_fs"
      d2s$idc=d2s$id10 # <1, 1-9, then 10 year bands
      d2s$om_age=d2s$om_age3 # <1, 1-9, then 10 year bands
      d2s$agebf=as.factor(d2s$ageband3) # <1, 1-4, 5-9, then 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband3 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + om_dose + agebf:age_ib + strata(idc)"
    } else if (j==5) {
      stratadesc="sev"
      d2s$idc=d2s$id11 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    } else if (j==6) {
      stratadesc="sev_fs"
      d2s$idc=d2s$id12 # <1, 1-9, then 10 year bands
      d2s$om_age=d2s$om_age3 # <1, 1-9, then 10 year bands
      d2s$agebf=as.factor(d2s$ageband3) #  <1, 1-4, 5-9, then 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband3 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    }  else if (j==7) {
      stratadesc="sev_utla" # UTLA stratification
      d2s$idc=d2s$id13 # stratification with 10 year bands 
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    }  else if (j==8) {
      stratadesc="sev_ms" # minimal stratification
      d2s$idc=d2s$id14 # stratification with 10 year bands 
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf + agebf:age_ib + NHSER_name + dose_vacc + strata(idc)"
    } else if (j==9) {
      stratadesc="sev_area" # area not in strata
      d2s$idc=interaction(d2s$ageband1, d2s$ethnicity_group2,d2s$Tg,d2s$dose_vacc) # 10 year bands + vacc
      d2s$idc=droplevels(d2s$idc)
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf:age_ib + NHSER_name + strata(idc)"
    } else if (j==10) {
      stratadesc="sev_age" # age not in strata
      d2s$idc=interaction(d2s$NHSER_name, d2s$ethnicity_group2,d2s$Tg,d2s$dose_vacc) # 10 year bands + vacc
      d2s$idc=droplevels(d2s$idc)
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf + agebf:age_ib + strata(idc)"
    } else if (j==11) {
      stratadesc="sev_eth" # ethnicity not in strata
      d2s$idc=interaction(d2s$NHSER_name, d2s$ageband1, d2s$Tg,d2s$dose_vacc) # 10 year bands + vacc
      d2s$idc=droplevels(d2s$idc)
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + ethnicity_group2 + agebf:age_ib + strata(idc)"
    } else if (j==12) {
      stratadesc="sev_vc" # vacc not in strata
      d2s$idc=interaction(d2s$ageband1,d2s$NHSER_name, d2s$ethnicity_group2,d2s$Tg) # 10 year bands + vacc
      d2s$idc=droplevels(d2s$idc)
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + dose_vacc + agebf:age_ib + strata(idc)"
    } else if (j==13) {
      stratadesc="sev_sex" # sex in strata
      d2s$idc=interaction(d2s$ageband1,d2s$NHSER_name, d2s$ethnicity_group2,d2s$Tg,d2s$dose_vacc,d2s$sex)
      d2s$idc=droplevels(d2s$idc)
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + dose_vacc + agebf:age_ib + strata(idc)"
    } else if (j==14) {
      stratadesc="sev_imd" # imd in strata
      d2s$idc=interaction(d2s$ageband1,d2s$NHSER_name, d2s$ethnicity_group2,d2s$Tg,d2s$dose_vacc,d2s$imd_decile)
      d2s$idc=droplevels(d2s$idc)
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + sex + dose_vacc + agebf:age_ib + strata(idc)"
    } else if (j==15) {
      stratadesc="sev_reinf_voc"  # reinf by variant
      d2s$idc=d2s$id11 # stratification with 10 year bands
      d2s$reinf_var=interaction(d2s$reinfection_flag,d2s$omicron)
      d2s$reinf_var=fct_collapse(d2s$reinf_var,"No past infection"=c("0.0","0.1"),"Past infection:Delta"="1.0","Past infection:Omicron"="1.1")
      d2s$reinf_var=relevel(d2s$reinf_var,"No past infection")
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + reinf_var + imd_decile + sex + agebf:age_ib + strata(idc)"
    } else if (j==16) {
      stratadesc="sev_epb1" # with epidemic phase bias
      epb=1.0
      d2s$idc=d2s$id11 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    } else if (j==17) {
      stratadesc="sev_epb2" # with epidemic phase bias
      epb=2.0
      d2s$idc=d2s$id11 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    } else if (j==18) {
      stratadesc="sev_vam" # VAM only
      d2s$omicron=d2s$omicron3
      d2s$delta=d2s$delta3
      d2s=d2s %>% filter(omicron==1|delta==1)
      d2s$idc=d2s$id11 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    } else if (j==19) {
      stratadesc="sev_seq" # sequencing only
      d2s$omicron=d2s$B11529_conf
      d2s$delta=d2s$B16172_conf
      d2s=d2s %>% filter(omicron==1|delta==1)
      d2s$idc=d2s$id11 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    } else if (j==20) {
      stratadesc="sev_vacc_reinf_npR"
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinfvar_rand + imd_decile + sex + om_dose + agebf:age_ib + strata(idc)" # updated below
    }  else if (j==21) {
      stratadesc="sev_vacc2"
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + num_doses:reinfection_flag + imd_decile + sex + om_dose2 + agebf:age_ib + strata(idc)"
    } else if (j==22) {
      stratadesc="sev_vacc_dose_reinf"
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + num_doses:reinfection_flag + imd_decile + sex + om_dose + agebf:age_ib + strata(idc)"
    } else if (j==23) {
      stratadesc="sev_unvacc_noage"
      d2s=d2s %>% filter(num_doses==0)
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ omicron + reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    }  else if (j==24) {
      stratadesc="sev_unvacc"
      d2s=d2s %>% filter(num_doses==0)
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    }  else if (j==25) {
      stratadesc="sev_unvacc_fs"
      d2s=d2s %>% filter(num_doses==0)
      d2s$idc=d2s$id10 #  <1, 1-4, 5-9, then 10 year bands
      d2s$om_age=d2s$om_age3 #  <1, 1-4, 5-9, then 10 year bands
      d2s$agebf=as.factor(d2s$ageband3) # <1, 1-9, then 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband3 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    } else if (j==26) {
      stratadesc="sev_reinfdose" # finer stratification of reinfection interaction with vaccination in main analysis
      d2s$idc=d2s$id11 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + num_doses:reinfection_flag + imd_decile + sex + agebf:age_ib + strata(idc)"
    } else if (j==27) {
      stratadesc="sev_vacc_reinf_order" # accounting for order of first infection and vacc
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + reinf_vacc_order + imd_decile + sex + om_dose + agebf:age_ib + strata(idc)"
    } else if (j==28) {
      stratadesc="sev_vacc_reinfR"
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vacc_reinf + imd_decile + sex + om_dose + agebf:age_ib + strata(idc)" # updated below
    } else if (j==29) {
      stratadesc="sev_vacc_reinfP"
      d2s$idc=d2s$id9 # stratification with 10 year bands
      d2s$om_age=d2s$om_age1 # 10 year bands
      d2s$agebf=as.factor(d2s$ageband1) # 10 year bands
      d2s$age_ib=d2s$age-d2s$ageband1 # age within ageband (subtract off lowest age)
      d2s$age_ib2=d2s$age_ib*d2s$age_ib # quadratic term
      form.str=" ~ om_age + vaccinated:reinf_pred + imd_decile + sex + om_dose + agebf:age_ib + strata(idc)" # updated below
    }
    form.str2=paste0("Surv(surv_time, surv_status)",form.str)
    d2s$idc=droplevels(d2s$idc)
    d2s$om_age=droplevels(d2s$om_age)
    f=TRUE
    res.cox.age=NULL
    for(i in endpoints) {
      print(i)
      if(j==20) { # fit to real valued reinf var
        form.str2=paste0("Surv(surv_time, surv_status)",sub("reinfvar",paste0(i,"_reinf"),form.str))
      }
      d2s[,"surv_status"]=d2s[,eval(i)]+1
      if(epb>0) { # shift specimen dates for endpoint cases and re-stratify
        d2s$lesd_epb=d2s$last_episode_specimen_date
        d2s$lesd_epb[d2s$surv_status==2]=d2s$lesd_epb[d2s$surv_status==2]+epb
        d2s$idc=interaction(d2s$ageband1,d2s$NHSER_name, d2s$ethnicity_group2,d2s$lesd_epb,d2s$dose_vacc) # 10 year bands + vacc
        d2s$idc=droplevels(d2s$idc)
      }
      if(i=="dead") {
        d2s$surv_date=d2s$last_specimen_date+28
      } else {
        d2s$surv_date=d2s$last_specimen_date+14
      }
      d2s$surv_date[d2s$surv_date>as.Date(hosp_data_date)]=as.Date(hosp_data_date)
      d2s[d2s$surv_status==2,"surv_date"]=d2s[d2s$surv_status==2,eval(paste0(i,"_date"))]
      d2s$surv_time=d2s$surv_date-d2s$last_specimen_date
      d2s$surv_time[d2s$surv_time<=0]=0.5 # small positive value
      cox.mod = coxph(as.formula(form.str2), data = d2s, robust=FALSE, y=FALSE)
      if(do_se) {
        cox.mod.s=summary(cox.mod)
        # gather estimates into a nice table
        res.cox.age1=as.data.frame(cox.mod.s$coefficients)
        names(res.cox.age1)[names(res.cox.age1) == "coef"] <- "Estimate"
        names(res.cox.age1)[names(res.cox.age1) == "se(coef)"] <- "Std. Error"
        res.cox.age1$var=sub("om_age1.","",sub("om_dose","",rownames(res.cox.age1)))
        if(f) {
          res.cox.age=res.cox.age1[,c("var","Estimate")]
          res.cox.age.num=res.cox.age
        }
        res.cox.age.num[,eval(paste0(i,".mean"))]=exp(res.cox.age1$Estimate)
        res.cox.age.num[,eval(paste0(i,".se"))]=exp(res.cox.age1$'Std. Error')
        res.cox.age.num[,eval(paste0(i,".l95"))]=exp(res.cox.age1$Estimate-1.96*res.cox.age1$'Std. Error')
        res.cox.age.num[,eval(paste0(i,".u95"))]=exp(res.cox.age1$Estimate+1.96*res.cox.age1$'Std. Error')
        res.cox.age[,eval(i)]=paste0(round(res.cox.age.num[,eval(paste0(i,".mean"))],2),
                                     " (", round(res.cox.age.num[,eval(paste0(i,".l95"))],2),
                                     "-",round(res.cox.age.num[,eval(paste0(i,".u95"))],2),")")
        if(f) { # table headers
          res.cox.age.num=res.cox.age.num[,-which(names(res.cox.age.num)=="Estimate")]
          res.cox.age=res.cox.age[,-which(names(res.cox.age)=="Estimate")]
          rownames(res.cox.age)=c()
          rownames(res.cox.age.num)=c()
          f=FALSE
        } 
      } else {
        if(f) {
          f=FALSE
          res.cox.age=as.data.frame(cox.mod$coefficients)
        } else {
          res.cox.age=cbind(res.cox.age,as.data.frame(cox.mod$coefficients))
        }
      }
    }
    

    if(do_se) {
      write.csv(res.cox.age,paste0("Omicron/res_cox_age_",datacut,"_",date_name,"_",agedesc,"_",stratadesc,".csv"))
      write.csv(res.cox.age.num,paste0("Omicron/res_cox_age_num_",datacut,"_",date_name,"_",agedesc,"_",stratadesc,".csv"))
    } else {
      names(res.cox.age)=endpoints
      res.cox.age$var=rownames(res.cox.age)
    }
  }
}

# restore full dataset
datacut=datacut_c
if(do_unvacc) d2=d2c

