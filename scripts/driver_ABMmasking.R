# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# this is the driver script to load specific settings (main.runs, maskathome.runs, nonphywt.runs, or R0x2.runs) 
# and run the simulations for different masking rates and timing

# note in the study, we run the simulations on a cluster
# here we use the main analysis as an exmple
# change the settings per your study need

dummy_dredn_no = 1 # example runs
dummy_dctc_no = 1 # example runs
dummy_wt_no = 1 # example runs

num_pop = 100000 # population size
num_ens = 1 # number of ensemble runs
num.yr.sim = 5 # number of years to run the model

n.age.grp.t = 8
date.start = as.Date('2020/3/1')

agg.by.cohort = F # whether to compute and return tallies by birth cohort (no need if run for a short time)
agg.by.sn = T # whether to compute and return tallies by season  (no need if run for a short time)
wt_nonphys.main = .2 # main setting for wt of nonphysical contact

# TO USE - FIRST SET THE WORKING DIRECTORY TO SOURCE FILE LOCATION
dir_data = '../data/'
dir_code = '../scripts/'
dir_res = '../results/test/' # set the result sub folder


source(paste0(dir_code, 'loadPackages.R'))

batch.tag = 'main' # diff tests


if(!file.exists(dir_res))
  dir.create(dir_res, recursive = T)

ens_no = 1; set.seed(ens_no)

# get specific settings for this run
if(grepl('main', batch.tag)){
  source(paste0(dir_code, 'get_main.runs.settings.R'))
} else if(grepl('maskathome', batch.tag)){
  source(paste0(dir_code, 'get_maskathome.runs.settings.R'))
} else if(grepl('nonphywt', batch.tag)){
  source(paste0(dir_code, 'get_nonphywt.runs.settings.R'))
  wt_nonphys.t = wt_nonphys_vec[dummy_wt_no] # .2  # assume non-phyical contact has 20% of the chance resulting infection
  # !!! change dummy_wt_no to specific numbers
} else if(grepl('R0x2', batch.tag)){
  source(paste0(dir_code, 'get_R0x2.runs.settings.R'))
} 

inclVax = inclVax.t # whether to include vaccination
UseContactDuration = UseContactDuration.t #  whether to use contact duration directly or use number of contacts
UseVxData = UseVxData.t # whether to use vax data for the vaccination module
UseAgeSpecSuscept = UseAgeSpecSuscept.t # whether to assume diff susceptibility by age. set it to FALSE to reduce uncertainty

MASKATHOME = MASKATHOME.t # set this to TRUE only when testing a very specific scenario that ppl also mask at home

prob.tx.per.ctc0 = prob.tx.per.ctc0.t # 12/29/23 this affects R0 (default = .5 -> R0 = ~2.6)

diff.ctc.redn.bymask = diff.ctc.redn.bymask_vec[dummy_dredn_no]
diff.ctc.bymask = diff.ctc.bymask_vec[dummy_dctc_no]


if(inclVax){
  status_vec = c('S','E','I','R','V')  # all possible status here
} else {
  status_vec = c('S','E','I','R')  # all possible status here
}

source(paste0(dir_code, 'loadPackages.R'))
source(paste0(dir_code, 'Fns_util_AMBmasking.R'))
source(paste0(dir_code, 'Fn_ABMmasking.R'))
options(dplyr.summarise.inform = FALSE)
options(datatable.verbose = F)



# read in age group setting
age.grps = read.csv(paste0(dir_code,'age.grps',n.age.grp.t,'.csv')) %>% data.table()
str_age_recode = age.grps[,c('code','value')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')

birth.rate = 1/age.grps[nrow(age.grps)]$age.end/365

mask.E.cut = c(.5, .9) # cut point for stratification of mask effectiveness
method.calHHrisk = 'like4like' # set method to calculate HH risk (like4like is probably more accurate than adjHHmask)

tm1 = Sys.time()

# load parameters for the simulation
source(paste0(dir_code, 'load_parms.R'))
# save a copy of parms used for this run
write.csv(parm.bounds, paste0(dir_res, 'parm.bounds.csv'), row.names = F)

# frac.masking_vec = seq(0, 1, by = .1)
# cut.mask.start_vec = c(0, .01/100, .02/100, .05/100, .1/100, .2/100, .5/100, 1/100, 2/100) # seq(0, 1, by = .1)
# cut.mask.start = seq(0, 100, by = 10) # diff starting point of masking

# test run, using the following settings:
dummy_fm_no = 2
dummy_cm_no = 2
ave.frac.masking = frac.masking_vec[dummy_fm_no];  
cut.mask.start = cut.mask.start_vec[dummy_cm_no]

set.seed(ens_no)  # use the same seed for the same ensemble run

{
  print(c(ave.frac.masking, cut.mask.start))
  
  model_inputs = list(tab_dates = tab_dates, # dates of simulation
                      da.vx.daily = da.vx.daily, # vaccination data
                      da.vx.weekly = da.vx.weekly,
                      da.p.ctc.redn = da.p.ctc.redn, # non-mask risk reduction (social distancing etc)
                      ctc.fit = ctc.fit # prob distributions for the contact patterns
  )
  model_settings = list(
    num_pop = num_pop, # population size
    num_ens = num_ens, # number of ensemble runs
    ave.frac.masking = ave.frac.masking, # fraction with population 'assigned' to wear mask
    cut.mask.start = cut.mask.start,  # threshold of population prevalence when mask wearing is enacted
    inclVax = inclVax, # whether to include vaccination
    tt.vax.start = tt.vax.start, # time when vaccination starts
    UseVxData = UseVxData, # whether to use vax data for the vaccination module
    UseContactDuration = UseContactDuration, #  whether to use contact duration directly or use number of contacts
    diff.ctc.redn.bymask = diff.ctc.redn.bymask, # whether assume diff level of redn in contact rate by masking status
    adj.ctc.redn.bymask = adj.ctc.redn.bymask, # specific settings for diff level of risk redn by masking status
    diff.ctc.bymask = diff.ctc.bymask, # whether assume diff level of contact rate for ppl have diff masking preference
    UseAgeSpecSuscept = UseAgeSpecSuscept, # whether to assume diff susceptibility by age. set it to FALSE to reduce uncertainty
    suscept.by.age = suscept.by.age, # specific settings
    age.grps = age.grps, # age group settings
    str_age_recode =str_age_recode,
    status_vec = status_vec,  # all possible status here
    mask.E.cut = mask.E.cut, # cut point for stratification of mask effectiveness
    method.calHHrisk = method.calHHrisk, # set method to calculate HH risk (like4like is probably more accurate than adjHHmask)
    frac.ini.inf0 = frac.ini.inf0, # initial prevalence
    seed0 = seed0, # background seeding
    days.epi.seeding = days.epi.seeding,
    ave.p.redn.mask.byage = ave.p.redn.mask.byage,
    sd.p.redn.mask.byage = sd.p.redn.mask.byage,
    ave.num.nonHHctc.byage = ave.num.nonHHctc.byage, 
    sd.num.nonHHctc.byage = sd.num.nonHHctc.byage,
    ave.num.HHctc.byage = ave.num.HHctc.byage,
    sd.num.HHctc.byage = sd.num.HHctc.byage,
    birth.rate = birth.rate,
    parm.immLossInf = parm.immLossInf,
    parm.immLossVax = parm.immLossVax,
    parm.immLossVseverity = parm.immLossVseverity,
    parm.getVaxProb = parm.getVaxProb,
    MASKATHOME = MASKATHOME,
    cont_freq = cont_freq.t,
    wt_nonphys = wt_nonphys.t,
    parm.type = parm.type.t
  )
  
  fname = paste0(dir_res, 'res_',
                 ifelse(inclVax,'vx','novx'),'_',
                 ifelse(UseAgeSpecSuscept,'Sbyage','So'),'_',
                 ifelse(diff.ctc.bymask,'CTCbymask','CTCo'),'_',
                 ifelse(diff.ctc.redn.bymask,'MOBbymask','MOBo'),
                 ifelse(wt_nonphys.t == wt_nonphys.main, '', paste0('_wt', dummy_wt_no)),
                 '_fm', dummy_fm_no, '_cm', dummy_cm_no, 
                 '_ens', ens_no, '.RData'
  )
  tmp = try(load(fname), silent = T)
  if(class(tmp) == 'try-error'){
    res = Fn_ABMmasking(model_inputs, model_settings, agg.by.cohort = agg.by.cohort, agg.by.sn = agg.by.sn) # 
    
    save(res, file = fname)
    
  }
}

tm2 = Sys.time()
print(tm1)
print(tm2)
print(paste('It took', tm2 - tm1))
