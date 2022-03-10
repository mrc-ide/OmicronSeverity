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


#################################
#### Load individual datasets ###
#################################

## name of linked linelist file
ll.name=paste0("data/full_linked_",date_name,link_desc,".fst")


## read linked dataset 
if(LoadPreviousLinkedList & file.exists(ll.name)) {
    ll=read.fst(ll.name)

} else {
  
  ### load line list
  
  if(file.exists(paste0("data/",ll_file,".fst"))) {
    llcopy <- read.fst(paste0("data/",ll_file,".fst"))
  } else {
    llcopy<-read_csv(paste0("data/",ll_file,".csv.xz"),guess_max = 1e6, progress = TRUE)
    if(SaveFST) write_fst(llcopy, paste0("data/",ll_file,".fst"), compress=90)
  }
  nrow(llcopy)
  llcopy = llcopy %>% filter(!is.na(finalid)) %>%
                      mutate(specimen_date=as.Date(specimen_date,"%d/%m/%Y"),
                        Onsetdate=as.Date(Onsetdate,"%d/%m/%Y"),
                        lab_report_date=as.Date(lab_report_date,"%d/%m/%Y"))
  
  ### load death file
  
  if(file.exists(paste0("data/",death_file,".fst"))) {
    dth <- read.fst(paste0("data/",death_file,".fst"))
  } else {
    dth<-read_csv(paste0("data/",death_file,".csv.xz"),guess_max = 1e6, progress = TRUE)
    if(SaveFST) write_fst(dth, paste0("data/",death_file,".fst"), compress=90)
  }
  nrow(dth)
  dth=dth %>% filter(!is.na(finalid)) # drop any without NHS number
  dth$dod=as.Date(dth$dod,"%d/%m/%Y")
  dth$dateadmission=as.Date(dth$dateadmission,"%d/%m/%Y")
  dth = dth %>% select(finalid,death_type28,dod,covidcod)
  
  
  ## load VAM file - no repeated finalids in VAM file
  
  if(file.exists(paste0("data/",vam_file,".fst"))) {
    vam <- read.fst(paste0("data/",vam_file,".fst"))
  } else {
    vam<-read_csv(paste0("data/",vam_file,".csv.xz"),guess_max = 1e6, progress = TRUE)
    if(SaveFST) write_fst(vam, paste0("data/",vam_file,".fst"), compress=90)
  }
  nrow(vam)
  vam=vam %>% filter(!is.na(finalid)) # drop any without NHS number
  vam$vam_specimen_date=as.Date(vam$specimen_date_sk)
  vam=vam %>% select(finalid,vam_specimen_date,variant,variant_alt_name,seq_result,vam_pillar=pillar)
  
  ## load SGTF file - some repeated finalids
  
  if(file.exists(paste0("data/",sgtf_file,".fst"))) {
    sgtf <- read.fst(paste0("data/",sgtf_file,".fst"))
  } else {
    sgtf=read_csv(paste0("data/",sgtf_file,".csv.xz"),guess_max = 1e6, progress = TRUE)
    if(SaveFST) write_fst(sgtf, paste0("data/",sgtf_file,".fst"), compress=90)
  }
  nrow(sgtf)
  sgtf=sgtf %>% rename(finalid=FINALID) %>% filter(!is.na(finalid)) # drop any without NHS number
  sgtf$sgtf_specimen_date=as.Date(sgtf$specimen_date,"%Y-%m-%d")
  
  # S-neg/pos only called if other genes have Ct<30
  sgtf$sneg=0
  sgtf$sneg[sgtf$P2CH3CQ==0 & sgtf$P2CH1CQ<30 & sgtf$P2CH1CQ>0 & sgtf$P2CH2CQ<30 & sgtf$P2CH2CQ>0]=1
  sgtf$spos=0
  sgtf$spos[sgtf$P2CH3CQ>0 & sgtf$P2CH3CQ<30 & sgtf$P2CH1CQ<30 & sgtf$P2CH1CQ>0 & sgtf$P2CH2CQ<30 & sgtf$P2CH2CQ>0]=1
  # only use SGTF data up to and including last_sgtf_date
  sgtf = sgtf %>% filter(specimen_date<=last_sgtf_date)
  # can be multiple SGTF test results per person, so create list with 1 row per person
  sgtf.sum = sgtf %>% group_by(finalid) %>% summarise(n_sgtf=n(),sgtf_last_specimen_date=max(sgtf_specimen_date))
  # join to main SGTF list
  sgtf2 = left_join(sgtf,sgtf.sum)
  # only keep definitive results
  sgtf2 = sgtf2 %>% filter((sneg==1 & spos==0)|(sneg==0 & spos==1))
  # only keep results within 28 days of last test
  sgtf2 = sgtf2 %>% filter(sgtf_last_specimen_date-sgtf_specimen_date<=28)
  # summarise all such tests for each individual. 
  sgtf.sum = sgtf2 %>% group_by(finalid) %>% 
    summarise(spos=sum(spos),sneg=sum(sneg),
              sgtf_first_specimen_date=min(sgtf_specimen_date),sgtf_last_specimen_date=max(sgtf_specimen_date))  
  # transform sums of spos and sneg into categorical 0/1 variables again
  sgtf.sum$sneg[sgtf.sum$sneg>0]=1
  sgtf.sum$spos[sgtf.sum$spos>0]=1
  # and drop results for people who get contradictory results within 28 days of last test
  sgtf.sum = sgtf.sum %>% filter((sneg==1 & spos==0)|(sneg==0 & spos==1))
  rm(sgtf2)
  
  ### load reinfection file
  
  if(file.exists(paste0("data/",reinf_file,".fst"))) {
    reinf <- read.fst(paste0("data/",reinf_file,".fst"))
  } else {
    reinf<-read_csv(paste0("data/",reinf_file,".csv.xz"),guess_max = 1e6, progress = TRUE)
    if(SaveFST) write_fst(reinf, paste0("data/",reinf_file,".fst"), compress=90)
  }
  nrow(reinf)
  reinf=reinf %>% filter(!is.na(finalid_primary_case))  # drop any without NHS number
  reinf$specimen_date=as.Date(reinf$SPECIMEN_DATE)
  reinf$first_specimen_date=as.Date(reinf$primary_specimen_date,"%d/%m/%Y")
  reinf$finalid=reinf$finalid_primary_case
  reinf$pillar=ifelse(reinf$pillar=="PILLAR 2 TESTING","Pillar 2","Pillar 1")
  
  ## load hospitalisation file
  
  if(file.exists(paste0("data/",hosp_file,".fst"))) {
    hosp <- read.fst(paste0("data/",hosp_file,".fst"))
  } else {
    hosp<-read_csv(paste0("data/",hosp_file,".csv.xz"),guess_max = 1e6, progress = TRUE)
    if(SaveFST) write_fst(hosp, paste0("data/",hosp_file,".fst"), compress=90)
  }
  nrow(hosp)
  hosp$specimen_date=as.Date(hosp$specimen_date,"%d%b%Y")
  hosp$spell_end_date=as.Date(hosp$spell_end_date,"%d%b%Y")
  hosp$spell_start_date=as.Date(hosp$spell_start_date,"%d%b%Y")
  hosp$arrival_date=as.Date(hosp$arrival_date,"%d%b%Y")
  hosp$departure_date=as.Date(hosp$departure_date,"%d%b%Y")
  hosp$hospital_in=as.Date(hosp$hospital_in,"%d%b%Y")
  hosp$hospital_out=as.Date(hosp$hospital_out,"%d%b%Y")
  
  hosp$hospital_in_orig=hosp$hospital_in
  hosp$hospital_out_orig=hosp$hospital_out
  
  # use raw SUS or ECDS hospitalisation dates
  hosp$hospital_in=as.Date(NA)
  hosp$hospital_out=as.Date(NA)
  hosp$hospital_in[!is.na(hosp$spell_start_date)]=hosp$spell_start_date[!is.na(hosp$spell_start_date)]
  hosp$hospital_out[!is.na(hosp$spell_end_date)]=hosp$spell_end_date[!is.na(hosp$spell_end_date)]
  hosp$hospital_in[!is.na(hosp$arrival_date)]=hosp$arrival_date[!is.na(hosp$arrival_date)]
  hosp$hospital_out[!is.na(hosp$departure_date)]=hosp$departure_date[!is.na(hosp$departure_date)]
  
  #remove rows with no hospital_in date and those prior to 2021 (to speed up processing)
  hosp=hosp[!is.na(hosp$hospital_in) & hosp$hospital_in>=as.Date("2021-01-04"),]
  
  #remove episodes where hospital_out<hospital_in (data errors)
  sel=(!is.na(hosp$hospital_in) & !is.na(hosp$hospital_out) & hosp$hospital_out<hosp$hospital_in)
  hosp$hospital_in[sel]=NA
  hosp$hospital_out[sel]=NA
  
  
  ## load immunisation file
  
  if(file.exists(paste0("data/",vacc_file,".fst"))) {
    vacc <- read.fst(paste0("data/",vacc_file,".fst"))
  } else {
    vacc<-read_csv(paste0("data/",vacc_file,".csv.xz"),guess_max = 1e6, progress = TRUE)
    if(SaveFST) write_fst(vacc, paste0("data/",vacc_file,".fst"), compress=90)
  }
  nrow(vacc)
  vacc = vacc %>% filter(!is.na(finalid))
  vacc$date_dose = as.Date(vacc$vaccination_date,"%d%b%Y")
  vacc = vacc %>% select(finalid,vacc_dose=product_display_type,date_dose,dose_number)
  # raw file is a list of all doses given. Transform to a list of people, with doses in separate columns
  vacc <- vacc %>% pivot_wider(id_cols=c(finalid),
                              names_from=dose_number, 
                              values_from=c(vacc_dose,date_dose)) 
  
  ## load travel file
  
  if(file.exists(paste0("data/",travel_file,".fst"))) {
    travel <- read.fst(paste0("data/",travel_file,".fst"))
  } else {
    travel<-read_csv(paste0("data/",travel_file,".csv.xz"),guess_max = 1e6, progress = TRUE)
    
  }
  nrow(travel)
  travel = travel %>% filter(!is.na(finalid) & overall_travel=="Yes") # only keep known travellers with NHS number
  travel$travel_specimen_date=ymd(travel$specimen_date_sk)
  # only keep last known travel date
  travel = travel %>% group_by(finalid) %>% summarise(travel_specimen_date=max(travel_specimen_date))
  
  
  ###############################
  ############ LINKAGE ##########
  ###############################
  

  ## reset linelist dataframe
  
  ll=llcopy
  
  ## reinfection_linkage - done first to determine most recent episode date
  
  if(pillar2.filter) {
    reinf.sum = reinf %>% filter(pillar=="Pillar 2")
  } else {
    reinf.sum=reinf
  }
  
  reinf.sum = reinf.sum %>% group_by(finalid) %>% 
    arrange(specimen_date,by_group=TRUE) %>%
    summarise(num_episodes=n(),last_episode_specimen_date=last(specimen_date),reinf_pillar=last(pillar),reinf_asymptomatic_indicator=last(Asymptomatic_Indicator))
  
  ll=left_join(ll,reinf.sum)
  rm(reinf.sum)

  # pillar 2 cases in reinfection list which are pillar 1 in main line list
  if(pillar2.filter) {
    sel=!is.na(ll$num_episodes) & ll$pillar=="Pillar 1"
    ll$num_episodes[sel]=ll$num_episodes[sel]-1
  }
  # update pillar for reinfections
  ll$pillar[!is.na(ll$last_episode_specimen_date)]=ll$reinf_pillar[!is.na(ll$last_episode_specimen_date)]
  # pillar 2 filter for linelist 
  if(pillar2.filter) {
    ll = ll %>% filter(pillar=="Pillar 2")
  }  
  # and asymp status
  ll$asymptomatic_indicator[!is.na(ll$last_episode_specimen_date)]=ll$reinf_asymptomatic_indicator[!is.na(ll$last_episode_specimen_date)]
  # if case not in reinfection list, use specimen date from main linelist as last_episode_specimen_date
  ll$last_episode_specimen_date[is.na(ll$last_episode_specimen_date)]=ll$specimen_date[is.na(ll$last_episode_specimen_date)]
  # set episode number for first infections
  ll$num_episodes[is.na(ll$num_episodes)]=0
  ll$num_episodes=ll$num_episodes+1
  # reinfection_flag = 1 only if num_episodes>0
  ll$reinfection_flag=ifelse(ll$num_episodes>1,1,0)
  # just in case any reinfections with <90 day interval slip into reinfection list
  ll$reinfection_flag[ll$last_episode_specimen_date-ll$specimen_date<=90]=0
  # save orig last episode date
  ll$orig_last_episode_specimen_date=ll$last_episode_specimen_date
  
  ## death linkage
  
  ll=left_join(ll,dth)
  ll$death_type28[is.na(ll$death_type28)]=0
  
  ## Add new reinfections from VAM and SGTF if not present in reinfection list?
  DoAddExtraReinf=TRUE
  
  ## VAM linkage

  if(pillar2.filter) {
    vam2=vam %>% filter(vam_pillar=="Pillar 2")
  } else {
    vam2=vam
  }
  ll=left_join(ll,vam2)
  rm(vam2)
  
  # only call VOC if VAM test date close to last_episode_date
  ll$last2vam=ll$vam_specimen_date-ll$last_episode_specimen_date
  
  # check if vam is >90 days after current last_episode_specimen_date derived from linelist and reinfection list
  # add as unrecognisied reinfection if so
  ll$vam_sgtf_reinf=0
  if(DoAddExtraReinf) {
    sel=(!is.na(ll$last2vam) & ll$last2vam > 90)
    ll$last_episode_specimen_date[sel]=ll$vam_specimen_date[sel]
    ll$last2vam[sel]=0
    ll$reinfection_flag[sel]=1
    ll$pillar[sel]=ll$vam_pillar[sel]
    ll$vam_sgtf_reinf[sel]=1
  }
  
  # only call VOC if VAM test date within 14 days of last_episode_date
  sel=(!is.na(ll$last2vam) & ll$last2vam <= 14 & ll$last2vam >= 0)
  ll$pillar[sel]=ll$vam_pillar[sel]
  # remove VAM info otherwise
  sel=(!is.na(ll$last2vam) & !(ll$last2vam <= 14 & ll$last2vam >= 0))
  ll$variant[sel]=NA
  ll$variant_alt_name[sel]=NA
  ll$seq_result[sel]=NA
  ll$vam_specimen_date[sel]=NA
  ll$vam_pillar[sel]=NA
  
  ## variant defs using genotyping
  ll$B16172=ifelse(!is.na(ll$variant) & !is.na(ll$seq_result) & 
                     (ll$seq_result=="Confirmed" | ll$seq_result=="Provisional genotyping") &
                     (ll$variant=="VOC-21APR-02"|ll$variant=="VUI-21OCT-01"),1,0)
  ll$B11529=ifelse(!is.na(ll$variant) & !is.na(ll$seq_result) &
                     (ll$seq_result=="Confirmed" | ll$seq_result=="Provisional genotyping") &
                     ll$variant=="VOC-21NOV-01",1,0)
  ll$other_var=ifelse(!is.na(ll$variant) & !is.na(ll$seq_result) &
                     (ll$seq_result=="Confirmed" | ll$seq_result=="Provisional genotyping") &
                     ll$variant!="VOC-21NOV-01" & ll$variant!="VOC-21APR-02" & ll$variant!="VUI-21OCT-01",1,0)
  ll$B16172_conf=ifelse(!is.na(ll$variant) & !is.na(ll$seq_result) & 
                     ll$seq_result=="Confirmed"  &
                     (ll$variant=="VOC-21APR-02"|ll$variant=="VUI-21OCT-01"),1,0)
  ll$B11529_conf=ifelse(!is.na(ll$variant) & !is.na(ll$seq_result) &
                     ll$seq_result=="Confirmed"  &
                     ll$variant=="VOC-21NOV-01",1,0)
  
  ## SGTF linkage
  
  if(!vam.filter) {

    ll=left_join(ll,sgtf.sum)
    
    # use both first and last SGTF dates within 28 days of the last SGTF test to match to most recent episode from linelist/reinfection list
    ll$lastsp2lastsgtf=ll$sgtf_last_specimen_date-ll$last_episode_specimen_date
    ll$lastsp2firstsgtf=ll$sgtf_first_specimen_date-ll$last_episode_specimen_date
    
    # sgtf is >90 days after current last_episode_specimen_date derived from linelist, reinfection and vam lists
    if(DoAddExtraReinf) {
      sel=(!is.na(ll$lastsp2firstsgtf) & ll$lastsp2firstsgtf > 90)
      # such instances get sgtf_first_specimen_date as the last_episode_specimen_date
      ll$last_episode_specimen_date[sel]=ll$sgtf_first_specimen_date[sel]
      ll$lastsp2firstsgtf[sel]=0
      ll$lastsp2lastsgtf[sel]=ll$sgtf_last_specimen_date[sel]-ll$sgtf_first_specimen_date[sel]
      ll$reinfection_flag[sel]=1
      ll$pillar[sel]="Pillar 2"
      ll$vam_sgtf_reinf[sel]=1
      # clear VAM fields in such cases
      ll$B11529[sel]=0
      ll$B16172[sel]=0
      ll$B11529_conf[sel]=0
      ll$B16172_conf[sel]=0
      ll$variant[sel]=NA
      ll$variant_alt_name[sel]=NA
      ll$seq_result[sel]=NA
      ll$vam_specimen_date[sel]=NA
      ll$vam_pillar[sel]=NA
    }
    # if SGTF test within 14 days of last_episode_specimen_date and there is no VAM result, class cases as pillar 2 (since all SGTF is pillar 2)
    sel=(!is.na(ll$lastsp2firstsgtf) & is.na(ll$vam_pillar) & 
           ((ll$lastsp2firstsgtf <= 14 & ll$lastsp2firstsgtf>=0)|(ll$lastsp2lastsgtf <= 14 & ll$lastsp2lastsgtf>=0)))
    ll$pillar[sel]="Pillar 2"
    
    # tests outside the 0-14 day window are ignored
    sel=(!is.na(ll$lastsp2firstsgtf) & 
           !((ll$lastsp2firstsgtf <= 14 & ll$lastsp2firstsgtf>=0)|(ll$lastsp2lastsgtf <= 14 & ll$lastsp2lastsgtf>=0)))
    ll$spos[sel]=NA
    ll$sneg[sel]=NA
    ll$sgtf_last_specimen_date[sel]=NA
    
  } else { # if vam.filter==TRUE
    
    ll$spos=NA
    ll$sneg=NA
    ll$sgtf_last_specimen_date=NA
    ll$sgtf_first_specimen_date=NA
    ll$sgtf2last=NA
  }
  
  ## adjust last episode date to earliest of VAM, SGTF and specimen date from main & reinfection lists
  
  ##ll$last_episode_specimen_date=pmin(ll$last_episode_specimen_date,ll$sgtf_first_specimen_date,ll$vam_specimen_date,na.rm=TRUE)
  
  ## Call Omicron and Delta - multiple alternative definitions
  
  # combine SGTF and VAM
  ll$omicron1=0
  ll$omicron1[(!is.na(ll$sneg) & ll$sneg==1 & ll$spos==0 & ll$B11529==0 & ll$B16172==0 & ll$other_var==0 ) | ll$B11529==1 ]=1
  
  ll$delta1=0
  ll$delta1[(!is.na(ll$spos) & ll$spos==1 & ll$sneg==0 & ll$B16172==0 & ll$B11529==0 & ll$other_var==0) | ll$B16172==1]=1
  
  # only SGTF
  ll$omicron2=0
  ll$omicron2[!is.na(ll$sneg) & ll$sneg==1]=1
  
  ll$delta2=0
  ll$delta2[!is.na(ll$spos) & ll$spos==1]=1
  
  #only VAM  
  ll$omicron3=0
  ll$omicron3[!is.na(ll$B11529) & ll$B11529==1]=1
  
  ll$delta3=0
  ll$delta3[!is.na(ll$B16172) & ll$B16172==1]=1
  
  #pillar 2 tests
  ll$omicron4=ifelse(!is.na(ll$vam_pillar) & ll$vam_pillar=="Pillar 2",ll$omicron1,ll$omicron2)
  ll$delta4=ifelse(!is.na(ll$vam_pillar) & ll$vam_pillar=="Pillar 2",ll$delta1,ll$delta2)
  
  # default to VAM+SGTF
  ll$omicron=ll$omicron1
  ll$delta=ll$delta1
  
  ## hospital linkage
  
  # analyses always use variable last_specimen_date - this now set to last_episode_specimen_date, which is the first specimen date for the most recent infection episode
  dr=ll %>% select(finalid, last_episode_specimen_date)
  # left_join cases with hosp
  dr=left_join(dr,hosp,by=c("finalid"="final_id"))
  dr=dr %>% filter(!is.na(hospital_in))
  dr$los=dr$hospital_out-dr$hospital_in
  dr$sp2adm=dr$hospital_in-dr$last_episode_specimen_date
  dr$sp2dis=dr$hospital_out-dr$last_episode_specimen_date
  # any attendance, 14 day threshold (dropped is.na(los)|los>=0 check, as los<0 already eliminated above)
  dr$hosp=0
  dr$hosp[!is.na(dr$sp2adm) & dr$sp2adm<=MaxSp2Hosp & dr$sp2adm>=MinSp2Hosp]=1
  # any attendance, 7 day threshold
  dr$hosp1=0
  dr$hosp1[!is.na(dr$sp2adm) & dr$sp2adm<=7 & dr$sp2adm>=MinSp2Hosp]=1
  # any attendance,LoS>=1 days
  dr$hosp2=0
  dr$hosp2[dr$hosp==1 & !is.na(dr$sp2dis) & dr$sp2dis>0 & !is.na(dr$los) & dr$los>0]=1
  # dashboard definition
  dr$hosp4=0
  dr$hosp4[!is.na(dr$sp2adm) & dr$sp2adm<=MaxSp2Hosp & ((!is.na(dr$sp2dis) & (dr$sp2dis>0)|(dr$sp2dis==0 & dr$sp2adm==0)) | is.na(dr$sp2dis))]=1
  # ECDS admission, 14 day threshold
  dr$hospE=0
  dr$hospE[dr$hosp==1 & !is.na(dr$ecds_discharge) & (dr$ecds_discharge=="Admitted"|dr$ecds_discharge=="Transfer"|dr$ecds_discharge=="Died")]=1
  # ECDS admission, 7 day threshold
  dr$hospF=0
  dr$hospF[dr$hosp1==1 & !is.na(dr$ecds_discharge) & (dr$ecds_discharge=="Admitted"|dr$ecds_discharge=="Transfer"|dr$ecds_discharge=="Died")]=1
  # SUS death
  dr$hospG=0
  dr$hospG[dr$hosp==1 & !is.na(dr$discharge_destination) & dr$discharge_destination==79]=1
  # ECDS admission or 1+ night stay
  dr$hospA=0
  dr$hospA[dr$hospE==1|dr$hosp2==1|dr$hospG==1 ]=1
  # past hosp (1-90 days prior to latest episode)
  dr$past_hosp=0
  dr$past_hosp[!is.na(dr$sp2adm) & dr$sp2adm>=-90 & dr$sp2adm<0 & (is.na(dr$sp2dis) | (!is.na(dr$sp2dis) & dr$sp2dis<0)) ]=1
  
  # track dates separately for each outcome (use date in future for those not meeting def)
  dr$hosp_date=dr$hospital_in
  dr$hosp1_date=dr$hospital_in
  dr$hosp2_date=dr$hospital_in
  dr$hosp4_date=dr$hospital_in
  dr$hospE_date=dr$hospital_in
  dr$hospF_date=dr$hospital_in
  dr$hospG_date=dr$hospital_in
  dr$hospA_date=dr$hospital_in
  dr$hosp_out_date=dr$hospital_out
  dr$hosp1_out_date=dr$hospital_out
  dr$hosp2_out_date=dr$hospital_out
  dr$hosp4_out_date=dr$hospital_out
  dr$hospE_out_date=dr$hospital_out
  dr$hospF_out_date=dr$hospital_out
  dr$hospG_out_date=dr$hospital_out
  dr$hospA_out_date=dr$hospital_out
  
  # set dates to NA for those who don't meet criteria
  dr$hosp_date[dr$hosp==0]=NA
  dr$hosp1_date[dr$hosp1==0]=NA
  dr$hosp2_date[dr$hosp2==0]=NA
  dr$hosp4_date[dr$hosp4==0]=NA
  dr$hospE_date[dr$hospE==0]=NA
  dr$hospF_date[dr$hospF==0]=NA
  dr$hospG_date[dr$hospG==0]=NA
  dr$hospA_date[dr$hospA==0]=NA
  dr$hosp_out_date[dr$hosp==0]=NA
  dr$hosp1_out_date[dr$hosp1==0]=NA
  dr$hosp2_out_date[dr$hosp2==0]=NA
  dr$hosp4_out_date[dr$hosp4==0]=NA
  dr$hospE_out_date[dr$hospE==0]=NA
  dr$hospF_out_date[dr$hospF==0]=NA
  dr$hospG_out_date[dr$hospG==0]=NA
  dr$hospA_out_date[dr$hospA==0]=NA
  
  # only keep hospitalisations meeting at least 1 definition
  dr=dr %>% filter(hosp==1|hosp1==1|hosp2==1|hosp4==1|hospE==1|hospF==1|hospG==1|hospA==1|past_hosp==1)
  
  # now group_by finalid and sum over numbers of each type of hosp for each finalid
  dr=dr %>% group_by(finalid) %>% summarise(hosp=sum(hosp),hosp1=sum(hosp1),
                                            hosp4=sum(hosp4),hosp2=sum(hosp2),
                                            hospE=sum(hospE),hospF=sum(hospF),
                                            hospG=sum(hospG),hospA=sum(hospA),
                                            past_hosp=sum(past_hosp),
                                            hosp_date=min(hosp_date,na.rm=TRUE),
                                            hosp1_date=min(hosp1_date,na.rm=TRUE),
                                            hosp2_date=min(hosp2_date,na.rm=TRUE),
                                            hosp4_date=min(hosp4_date,na.rm=TRUE),
                                            hospE_date=min(hospE_date,na.rm=TRUE),
                                            hospF_date=min(hospF_date,na.rm=TRUE),
                                            hospG_date=min(hospG_date,na.rm=TRUE),
                                            hospA_date=min(hospA_date,na.rm=TRUE),
                                            hosp_out_date=min(hosp_out_date,na.rm=TRUE),
                                            hosp1_out_date=min(hosp1_out_date,na.rm=TRUE),
                                            hosp2_out_date=min(hosp2_out_date,na.rm=TRUE),
                                            hosp4_out_date=min(hosp4_out_date,na.rm=TRUE),
                                            hospE_out_date=min(hospE_out_date,na.rm=TRUE),
                                            hospF_out_date=min(hospF_out_date,na.rm=TRUE),
                                            hospG_out_date=min(hospG_out_date,na.rm=TRUE),
                                            hospA_out_date=min(hospA_out_date,na.rm=TRUE),
                                            hosp_rank=mean(hospital_event_rank))
  
  # now join aggregated hosp list to main case list
  ll=left_join(ll,dr)
  # and finalise hosp status vars
  ll$hosp[is.na(ll$hosp)]=0
  ll$hosp1[is.na(ll$hosp1)]=0
  ll$hosp4[is.na(ll$hosp4)]=0
  ll$hosp2[is.na(ll$hosp2)]=0
  ll$hospE[is.na(ll$hospE)]=0
  ll$hospF[is.na(ll$hospF)]=0
  ll$hospG[is.na(ll$hospG)]=0
  ll$hospA[is.na(ll$hospA)]=0
  ll$past_hosp[is.na(ll$past_hosp)]=0
  ll$hosp[ll$hosp>0]=1
  ll$hosp1[ll$hosp1>0]=1
  ll$hosp4[ll$hosp4>0]=1
  ll$hosp2[ll$hosp2>0]=1
  ll$hospE[ll$hospE>0]=1
  ll$hospF[ll$hospF>0]=1
  ll$hospG[ll$hospG>0]=1
  ll$hospA[ll$hospA>0]=1

  # filter out any with past hosp, if required
  if(RemovePastHosp) ll=ll %>% filter(past_hosp==0)
  

  ## link to vacc data
  
  ll = left_join(ll,vacc)
  
  # remove doses after last_episode_specimen_date
  sel=(!is.na(ll$date_dose_1) & ll$last_episode_specimen_date<ll$date_dose_1)
  ll$vacc_dose_1[sel]=NA
  ll$date_dose_1[sel]=NA
  
  sel=(!is.na(ll$date_dose_2) & ll$last_episode_specimen_date<ll$date_dose_2)
  ll$vacc_dose_2[sel]=NA
  ll$date_dose_2[sel]=NA
  
  sel=(!is.na(ll$date_dose_3) & ll$last_episode_specimen_date<ll$date_dose_3)
  ll$vacc_dose_3[sel]=NA
  ll$date_dose_3[sel]=NA
  
  sel=(!is.na(ll$date_dose_4) & ll$last_episode_specimen_date<ll$date_dose_4)
  ll$vacc_dose_4[sel]=NA
  ll$date_dose_4[sel]=NA
  
  sel=(!is.na(ll$date_dose_5) & ll$last_episode_specimen_date<ll$date_dose_5)
  ll$vacc_dose_5[sel]=NA
  ll$date_dose_5[sel]=NA
  
  sel=(!is.na(ll$date_dose_6) & ll$last_episode_specimen_date<ll$date_dose_6)
  ll$vacc_dose_6[sel]=NA
  ll$date_dose_6[sel]=NA
  
  
  ## travel linkage
  
  ll=left_join(ll,travel)
  ll$traveller=0
  ll$traveller[!is.na(ll$travel_specimen_date) & ll$last_episode_specimen_date-ll$travel_specimen_date<14]=1
  
  
  ## cleanup
  rm(llcopy, vacc,dth,reinf,hosp,vam,sgtf,dr,travel,sgtf.sum,sel)
  gc()
  file.remove(list.files(tempdir(), full.names = TRUE))
  
  ## save linked linelist

  if(SaveLinkedList) write.fst(ll,ll.name,compress = 90)
}


