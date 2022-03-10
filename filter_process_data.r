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

####################
### Filter data  ###
####################


if(TRUE) { # if statement just to allow collapse of block in editor
  
  # last_specimen_date is key date variable in later analyses - set to the first test date in the last infection episode
  
  ll$last_specimen_date=ll$last_episode_specimen_date
  
  print(paste0("All cases:", nrow(ll))) # these print statements indicate effect of data filtering
  
  # Date range
  ll=ll %>% filter(last_specimen_date>=as.Date(date_start) & last_specimen_date<=as.Date(date_stop))
  print(paste0("All cases in date range:",nrow(ll)))
  
  ## filtering on variant - usually done later, but can uncomment these lines to generate numbers in Appendix fig 1
  #ll=ll %>% filter(omicron==1 | delta==1)
  #print(paste0("Omicron or Delta only:",nrow(ll)))
  
  # Missingness
  
  ll=ll %>% filter(!is.na(age) & !is.na(sex) & sex!="Unknown" & !is.na(imd_decile) &
                        !is.na(ethnicity_final))
  print(paste0("All cases without missing sex, age, imd_decile, ethnicity:",nrow(ll)))
  # Traveller status
  
  ll=ll %>% filter(traveller==0)
  print(paste0("All cases excl travellers:",nrow(ll)))
  
  ## remove cases who died before test
  ll=ll %>% filter(is.na(dod)|(!is.na(dod) & dod>=last_episode_specimen_date))
  print(paste0("Dropping people with DoD before test:",nrow(ll)))
        
  # remove people with 4+ doses
  ll = ll %>% filter(is.na(date_dose_4) & is.na(date_dose_5) & is.na(date_dose_6))
  print(paste0("Dropping people with 4+ doses:",nrow(ll)))
  
  # remove people with third dose <80 days after second dose - these are CXVs with more than 2 doses in primary course
  ll = ll %>% filter(is.na(date_dose_3) | (!is.na(date_dose_3) & !is.na(date_dose_2) & date_dose_3-date_dose_2>=80))
  print(paste0("Dropping people with <80 days between doses 2 and 3:",nrow(ll)))
  
  # use dose 2 vacc type if dose 1 type is missing
  ll$vacc_dose_1[is.na(ll$vacc_dose_1)]=ll$vacc_dose_2[is.na(ll$vacc_dose_1)]
  
  # only people with AZ, MD or PF vaccines (only MD and PF for dose 3)
  ll = ll %>% filter(is.na(vacc_dose_1) | vacc_dose_1=="AZ" | vacc_dose_1=="PF" | vacc_dose_1=="MD")
  ll = ll %>% filter(is.na(vacc_dose_2) | vacc_dose_2=="AZ" | vacc_dose_2=="PF" | vacc_dose_2=="MD")
  ll = ll %>% filter(is.na(vacc_dose_3) | vacc_dose_3=="PF" | vacc_dose_3=="MD")
  print(paste0("Dropping non AZ/PF/MD vaccines:",nrow(ll)))
  
}


####################################################
### Define required derived variables/strata     ###
####################################################

if(TRUE) {
  
  ########################################
  ### define vaccine status related vars #
  ########################################
  
  # delays & vacc status vars
  
  ll$vacc2spec1=ll$last_episode_specimen_date-ll$date_dose_1
  ll$vacc2spec2=ll$last_episode_specimen_date-ll$date_dose_2
  ll$vacc2spec3=ll$last_episode_specimen_date-ll$date_dose_3
  
  # assume that vaccine doses reach peak effectiveness 14 days after dose
  vacc2caseT=14
  vacc2caseD2T=14
  vacc2caseD3T=14
  
  ll$dose_done=ifelse(is.na(ll$vacc2spec1)|ll$vacc2spec1<0,"None",ifelse(ll$vacc2spec1<vacc2caseT,"D1:<14","D1:21+"))
  ll$dose_done=ifelse(is.na(ll$vacc2spec2)|ll$vacc2spec2<0,ll$dose_done,ifelse(ll$vacc2spec2<vacc2caseD2T,"D2:<14","D2:14+"))
  ll$dose_done=ifelse(is.na(ll$vacc2spec3)|ll$vacc2spec3<0,ll$dose_done,ifelse(ll$vacc2spec3<vacc2caseD3T,"D3:<14","D3:14+"))
  
  ll$num_doses=ifelse(is.na(ll$vacc2spec1)|ll$vacc2spec1<0,0,1)
  ll$num_doses=ifelse(is.na(ll$vacc2spec2)|ll$vacc2spec2<0,ll$num_doses,2)
  ll$num_doses=ifelse(is.na(ll$vacc2spec3)|ll$vacc2spec3<0,ll$num_doses,3)
  ll$num_doses=as.factor(ll$num_doses)
  
  ll$last_vacc2spec=ifelse(is.na(ll$vacc2spec1)|ll$vacc2spec1<0,0,ifelse(ll$vacc2spec1<vacc2caseT,-1,ll$vacc2spec1-vacc2caseT))
  ll$last_vacc2spec=ifelse(is.na(ll$vacc2spec2)|ll$vacc2spec2<0,ll$last_vacc2spec,ifelse(ll$vacc2spec2<vacc2caseD2T,-1,ll$vacc2spec2-vacc2caseD2T))
  ll$last_vacc2spec=ifelse(is.na(ll$vacc2spec3)|ll$vacc2spec3<0,ll$last_vacc2spec,ifelse(ll$vacc2spec3<vacc2caseD3T,-1,ll$vacc2spec3-vacc2caseD3T))
  
  ll$vacc_type=ifelse(is.na(ll$vacc_dose_1)| ll$dose_done=="None","None",ifelse(ll$vacc_dose_1=="AZ","AZ","PF/MD"))
  ll$dose_type=as.factor(ifelse(ll$dose_done=="None"|ll$vacc_type=="None","None",paste0(ll$vacc_type,":",ll$dose_done)))
  ll$dose_type=relevel(ll$dose_type,"None")
  
  ll$vacc_type[!is.na(ll$num_doses) & ll$num_doses==0]="None"
  ll$num_doses[ll$vacc_type=="None"]=0
  
  ll$dose_vacc=as.factor(ifelse(ll$num_doses==0|ll$vacc_type=="None","None",paste0(ll$vacc_type,":",ll$num_doses)))
  ll$dose_vacc=relevel(ll$dose_vacc,"None")
  
  ll$last_vacc2spec[ll$dose_vacc=="None"]=0
  
  # vaccination by follow-up time
  ll$lv2s=ifelse(ll$last_vacc2spec<0,0,ifelse(ll$last_vacc2spec<42,8,12+4*floor((ll$last_vacc2spec-42)/28)))
  ll$lv2s[ll$lv2s>=24]=24
  ll$lv2s[ll$lv2s>16 & ll$num_doses==3]=16
  ll$lv2s[ll$lv2s>16 & ll$num_doses==1]=16
  ll$lv2s[ll$num_doses==1 & ll$vacc_type=="AZ"]=0
  ll$lv2s=as.factor(ll$lv2s)
  
  # not stratified by VOC - for tables
  ll$dose_lv=interaction(ll$dose_vacc,ll$lv2s)
  ll$dose_lv=droplevels(ll$dose_lv)
  ll$dose_lv=fct_collapse(ll$dose_lv,"Unvaccinated"=c("None.8"))
  ll$dose_lv=relevel(ll$dose_lv,"Unvaccinated")
  
  # overall vacc status
  ll$vaccinated="Vacc"
  ll$vaccinated[ll$num_doses==0]="Unvacc"
  ll$vaccinated=as.factor(ll$vaccinated)
  
  ## other vars
  
  ll$imd_decile=as.factor(ll$imd_decile)
  
  ll$ageband=cut(ll$age, breaks = c(0,1,5,10,20,30,40,50,60,70,80,Inf),
                 right = FALSE, include.lowest = TRUE, ordered.result = TRUE)
  
  ll$ethnicity_final=as.factor(ll$ethnicity_final)
  ll$ethnicity_final=relevel(ll$ethnicity_final,"British (White)")
  
  ll$ethnicity_group2 <- fct_collapse(ll$ethnicity_final,
                                      white = c("Any other White background",
                                                "Irish (White)",
                                                "British (White)"),
                                      black=c("African (Black or Black British)",
                                              "Any other Black background",
                                              "White and Black African (Mixed)",
                                              "White and Black Caribbean (Mixed)",
                                              "Caribbean (Black or Black British)"),
                                      all_other=c("White and Asian (Mixed)",
                                                  "Any other Asian background",
                                                  "Any other ethnic group",
                                                  "Any other Mixed background",
                                                  "Chinese (other ethnic group)"),
                                      pakastani_or_bangladeshi =   c("Bangladeshi (Asian or Asian British)",
                                                                     "Pakistani (Asian or Asian British)"),
                                      indian =         c("Indian (Asian or Asian British)"),
                                      unknown =        c("Unknown"))
  
  ll$ethnicity_group2=relevel(ll$ethnicity_group2,"white")
  
  ## death endpoint
  ll$dead=ifelse(!is.na(ll$dod) & (ll$dod-ll$last_episode_specimen_date)<=MaxSp2Death & (ll$dod-ll$last_episode_specimen_date)>=MinSp2Hosp,1,0)
  
  # death date
  ll$dead_date=ll$dod
  
  # Tg =  days since earliest specimen date in study interval
  min_sp_date=min(ll$last_specimen_date)
  ll$Tg=as.numeric(ll$last_specimen_date-min_sp_date)
  
  # reinfection allowing for order of first infection and first dose
  
  ll$reinf_vacc_order=as.factor(ifelse(ll$reinfection_flag==0, "First infection",
                                       ifelse(ll$reinfection_flag==1 & ll$vaccinated=="Unvacc","Prior infection, unvaccinated",
                                              ifelse(!is.na(ll$date_dose_1) & ll$specimen_date<ll$date_dose_1,"Prior infection before dose 1","Prior infection after dose 1"))))
  
  # ageband definitions
  ll$ageband1=floor(ll$age/10)*10  # 10 year bands
  ll$ageband1[ll$ageband1>=80]=80
  ll$ageband2=ifelse(ll$age<5,0,ifelse(ll$age<10,5,floor(ll$age/10)*10)) # 0-4, 5-9 then 10 year bands
  ll$ageband2[ll$ageband2>=80]=80
  ll$ageband3=ifelse(ll$age<1,0,ifelse(ll$age<5,1,ifelse(ll$age<10,5,floor(ll$age/10)*10))) # <1, 1-4, 5-9, then 10 year bands
  ll$ageband3[ll$ageband3>=80]=80
  ll$ageband4=ifelse(ll$age<1,0,ifelse(ll$age<10,1,floor(ll$age/10)*10)) # <1, 1-9, then 10 year bands
  ll$ageband4[ll$ageband4>=80]=80
  

  # stratifications
  
  ll$id9=interaction(ll$ageband1,ll$NHSER_name, ll$ethnicity_group2, ll$Tg) # 10 year bands
  ll$id10=interaction(ll$ageband3,ll$NHSER_name, ll$ethnicity_group2,ll$Tg) # <1, 1-4, 5-9, then 10 year bands
  ll$id11=interaction(ll$ageband1,ll$NHSER_name, ll$ethnicity_group2,ll$Tg,ll$dose_vacc) # 10 year bands + vacc
  ll$id12=interaction(ll$ageband3,ll$NHSER_name, ll$ethnicity_group2,ll$Tg,ll$dose_vacc) # <1, 1-4, 5-9, then 10 year bands + vacc
  ll$id13=interaction(ll$ageband1,ll$UTLA_name, ll$ethnicity_group2,ll$Tg,ll$dose_vacc) # 10 year band, UTLA +vacc
  ll$id14=interaction(ll$Tg) # spec date only
  ll$id15=interaction(ll$ageband1,ll$NHSER_name, ll$Tg) # 10 year bands, no ethnicity

}



############################
### Data selection       ###
############################

# filter by Pillar if required

if(pillar2.filter) {
	d2all=ll %>% filter(pillar=="Pillar 2")
	d2all$omicron=d2all$omicron4
	d2all$delta=d2all$delta4
	datacut=paste0("pillar2",link_desc)
} else {
	d2all=ll 
	datacut=paste0("pillar1+2",link_desc)
	}

# remove cases with invalid NHS numbers
if(NHSnv.filter) {
	d2all=d2all %>% filter(NHS_nv==0)
	print(paste0("Dropping invalid NHS number:",nrow(d2all)))
} else {
	datacut=paste0(datacut,"nv")
}

