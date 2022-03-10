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

library(tidyverse)
library(forcats)
library(fst)
library(lubridate)
library(survival)
library(pROC)
library(parglm)

# working directory
setwd("D:/severity")

# raw data file names
ll_file="Anonymised Combined Line List 20220124"
hosp_file="20220124 SUS ECDS linked data"
reinf_file="20220124 reinfection SPIM"
death_file="20220124 COVID19 Deaths"
vacc_file="20220124 immunisations SPIM"
vam_file="20220124 VAM line list"
sgtf_file="SGTF_linelist_20220124"
travel_file="Modellers_rapid_travel_20220124"

# used in output file names
date_name="2022-01-24"
# assumed censoring data
hosp_data_date="2022-01-24"

# study period
date_start="2021-11-29"
date_stop="2022-01-09"

# last date when SGTF data is reliable to call delta/omicron
last_sgtf_date=as.Date("2021-12-30") 

# various flags
do_delays=FALSE # don't output delays
vam.filter=FALSE # default value - use both VAM and SGTF
do_unvacc=FALSE # default value - include vacc and unvacc
compute_retained=FALSE # compute data retained after stratification
do_bootstrap=FALSE # boostrap data?
do_impute=FALSE # impute VOC? [not used]
do_se=TRUE # calculate SE of Cox estimates?
do_model_reinf_impute=FALSE # impute reinfection status for "first" infections

# endpoints to run
# 	hosp = attendance or admission, hosp4 = as hosp, but incl diagnosis in hosp
# 	hospA =  admission only, dead =  death

endpoints=c("hosp","hosp4","hospA","dead")


################################################################
## main analysis - pillars 1 and 2, dropping invalid NHS numbers
################################################################

# set to TRUE once initial (slow) linkage done
LoadPreviousLinkedList=TRUE
# flag to save data files in FST format
SaveFST=FALSE
# save final linked dataset?
SaveLinkedList=TRUE
# minimum delay from test to hospital attendance - 0 by default, but also look at 1
MinSp2Hosp=0
# max delay from test to hosp attendance - 14 days by default
MaxSp2Hosp=14
# max delay from test to death - 28 days by default
MaxSp2Death=28

# Flag to remove instances of past hospitalisation
RemovePastHosp=FALSE
# set to "sp1" for MinSp2Hosp==1
link_desc=""
# Don't filter on pillar 2
pillar2.filter=FALSE
source("Omicron/R/load_link_data.r")
# Do filter on valid NHS number
NHSnv.filter=TRUE
source("Omicron/R/filter_process_data.r")
source("Omicron/R/bootstrap_impute.r")

# run main models
do_numbers=TRUE # output counts
model_set = c(1:3,5) # 5=primary analysis, 2=primary, no age for Om:Dlt HR, 3=secondary, 5=2ndary, no age for Om:Dlt HR
source("Omicron/R/run_models.r")

# most sensitivity analyses
do_numbers=FALSE # output counts
model_set = c(6:17,22,27)
source("Omicron/R/run_models.r")


# sensitivity analysis - unvaccinated only
do_unvacc=TRUE
do_numbers=TRUE # output counts
model_set = c(23,24)
source("Omicron/R/run_models.r")
do_unvacc=FALSE # reset unvacc flag


#############################################################################
## secondary analysis - pillars 1+2, including cases with invalid NHS numbers
#############################################################################
# no need to relink if run after main analysis
NHSnv.filter=FALSE
source("Omicron/R/filter_process_data.r")
source("Omicron/R/bootstrap_impute.r")
model_set=c(5) # just main model ("adjusted")
do_numbers=TRUE # output counts
source("Omicron/R/run_models.r")
NHSnv.filter=TRUE # reset flag

#######################################################
## secondary analysis - reinfection under-ascertainment
######################################################

# no need to relink if run after main analysis
NHSnv.filter=TRUE
source("Omicron/R/filter_process_data.r")

do_bootstrap=TRUE # enable bootstrapping
do_se=FALSE # don't bother to output SEs
model_set=c(20) # Secondary analysis, using imputed reinfections
do_numbers=FALSE # output counts

# random number seed. Results in paper split computation across 4 instances, using these seeds:
set.seed(10287918) # iter=1-50
#set.seed(9873021) # iter=51-100
#set.seed(51037871) # iter=101-150
#set.seed(29910943) # iter=151-200

boot.res=NULL # data frame for results
# run iterations
for(iter in 1:200) {
  tm=system.time({
    source("Omicron/R/bootstrap_impute.r")
    source("Omicron/R/run_models.r")
    res.cox.age$iter=iter
    boot.res=rbind(boot.res,res.cox.age) # res.cox.age stores results from run_models.r
  })
  print(iter)
  print(tm)
}
# exponentiate HRs
boot.res$hosp=exp(boot.res$hosp)
boot.res$hosp4=exp(boot.res$hosp4)
boot.res$hospA=exp(boot.res$hospA)
boot.res$dead=exp(boot.res$dead)

# quick summaries
boot.res.sum = boot.res %>% group_by(var) %>% 
  summarise(hosp.median=quantile(hosp,0.5),hosp.sd=sd(hosp),
            hosp.l95=quantile(hosp,0.025),hosp.u95=quantile(hosp,0.975),
            hosp4.median=quantile(hosp4,0.5),hosp4.sd=sd(hosp4),
            hosp4.l95=quantile(hosp4,0.025),hosp4.u95=quantile(hosp4,0.975),
            hospA.median=quantile(hospA,0.5),hospA.sd=sd(hospA),
            hospA.l95=quantile(hospA,0.025),hospA.u95=quantile(hospA,0.975),
            dead.median=quantile(dead,0.5,na.rm=TRUE),dead.sd=sd(dead, na.rm=TRUE),
            dead.l95=quantile(dead,0.025, na.rm=TRUE),dead.u95=quantile(dead,0.975, na.rm=TRUE),
            dead.na= sum(is.na(dead)))

# save raw results and summary
write.csv(boot.res,paste0("Omicron/boot_res_",datacut,"_",date_name,"_",agedesc,"_",stratadesc,".csv"))
write.csv(boot.res.sum,paste0("Omicron/boot_res_sum_",datacut,"_",date_name,"_",agedesc,"_",stratadesc,".csv"))

# reset flags
do_bootstrap=FALSE
do_se=TRUE

#############################################################
## secondary analysis - pillars 1+2, 1-14 day hosp definition
#############################################################

# need to relink
pillar2.filter=FALSE
MinSp2Hosp=1
RemovePastHosp=FALSE
# set link_desc to "sp1" for MinSp2Hosp==1
link_desc="_sp1"
source("Omicron/R/load_link_data.r")
NHSnv.filter=TRUE
source("Omicron/R/filter_process_data.r")
source("Omicron/R/bootstrap_impute.r")
model_set=c(5)  # just main model ("adjusted")
do_numbers=TRUE # output counts
source("Omicron/R/run_models.r")


#############################################################
## secondary analysis - VAM only (don't use SGTF)
#############################################################

# need to relink
pillar2.filter=FALSE
MinSp2Hosp=0
RemovePastHosp=FALSE
# Only use VAM data
vam.filter=TRUE
# set link_desc 
link_desc="_vam"
source("Omicron/R/load_link_data.r")
NHSnv.filter=TRUE
source("Omicron/R/filter_process_data.r")
source("Omicron/R/bootstrap_impute.r")
model_set=c(18,19)  # main model, 19=sequence only, 18=sequence+genotype
do_numbers=TRUE # output counts
source("Omicron/R/run_models.r")
# reset VAM filter flag
vam.filter=FALSE

###################################################################
## secondary analysis - pillar 2 only, dropping invalid NHS numbers
###################################################################

# need to relink
pillar2.filter=TRUE
MinSp2Hosp=0
RemovePastHosp=FALSE
link_desc="_p2"
source("Omicron/R/load_link_data.r")
NHSnv.filter=TRUE
source("Omicron/R/filter_process_data.r")
source("Omicron/R/bootstrap_impute.r")
model_set = c(5) # just main model ("adjusted")
do_numbers=TRUE # output counts
source("Omicron/R/run_models.r")



