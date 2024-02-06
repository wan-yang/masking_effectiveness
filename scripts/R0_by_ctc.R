# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# Get R0 for different contact settings
# Get example ctc pattern from the function: fn_initalizeABM_ens

dir_data = '../data/'
dir_code = '../scripts/'
dir_res = '../res/test/'

num_pop = 100000 # population size
num_ens = 1 # 100 # number of ensemble runs
num.yr.sim = 1 # number of years to run the model

n.age.grp.t = 8
date.start = as.Date('2020/3/1')

agg.by.cohort = F # whether to compute and return tallies by birth cohort (no need if run for a short time)
agg.by.sn = T # whether to compute and return tallies by season  (no need if run for a short time)
wt_nonphys.main = .2 # main setting for wt of nonphysical contact
diff.ctc.redn.bymask = F # whether assume diff level of redn in contact rate by masking status
diff.ctc.bymask = F # whether to assume diff (e.g. lower) contact rate for those prefer masking

ens_no = 1; set.seed(ens_no)
source(paste0(dir_code, 'loadPackages.R'))
source(paste0(dir_code, 'get_main.runs.settings.R'))
# source(paste0(dir_code, 'get_nonphywt.runs.settings.R'))
# source(paste0(dir_code, 'R0x2.runs.settings.R'))

prob.tx.per.ctc0 = prob.tx.per.ctc0.t # 12/29/23 this affects R0 (default = .5 -> R0 = ~2.6)

inclVax = inclVax.t # whether to include vaccination
UseContactDuration = UseContactDuration.t #  whether to use contact duration directly or use number of contacts
UseVxData = UseVxData.t # whether to use vax data for the vaccination module
UseAgeSpecSuscept = UseAgeSpecSuscept.t # whether to assume diff susceptibility by age. set it to FALSE to reduce uncertainty
# diff.ctc.redn.bymask = F # whether assume diff level of redn in contact rate by masking status
# diff.ctc.bymask = F # whether to assume diff (e.g. lower) contact rate for those prefer masking
MASKATHOME = MASKATHOME.t # set this to TRUE only when testing a very specific scenario that ppl also mask at home

if(inclVax){
  status_vec = c('S','E','I','R','V')  # all possible status here
} else {
  status_vec = c('S','E','I','R')  # all possible status here
}

source(paste0(dir_code, 'Fns_util_AMBmasking.R'))
source(paste0(dir_code, 'Fn_ABMmasking.R'))
options(dplyr.summarise.inform = FALSE)
options(datatable.verbose = F)



# read in age group setting
age.grps = read.csv(paste0(dir_code,'age.grps',n.age.grp.t,'.csv')) %>% data.table()
# age.grps$label = age.grps[,c('age.start','age.end')] %>% apply(1, paste, collapse = '-') %>%
#   gsub(pattern = '0-0', x = ., replacement = '<1')
# age.grps[nrow(age.grps)]$label = paste0(age.grps[nrow(age.grps)]$age.start,'+')
# age.grps$value = age.grps[,c('age.start','age.end')] %>% apply(1, paste, collapse = '-') %>% sQuote(q='') %>%
#   gsub(pattern = '0-0', x = ., replacement = '<1') %>%
#   gsub(pattern = '65-80', x = ., replacement = '65+')
# age.grps$lower = paste("`age` >=", age.grps$age.start)
# age.grps$upper = paste("`age` <=", age.grps$age.end)
# age.grps$upper[nrow(age.grps)] = "`age` <= 150" # increase the upper bound
# age.grps$code = age.grps[, c('lower', 'upper')] %>% apply(1, paste, collapse = ' & ')
str_age_recode = age.grps[,c('code','value')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')

birth.rate = 1/age.grps[nrow(age.grps)]$age.end/365

mask.E.cut = c(.5, .9) # cut point for stratification of mask effectiveness
method.calHHrisk = 'like4like' # set method to calculate HH risk (like4like is probably more accurate than adjHHmask)

tm1 = Sys.time()

# load parameters for the simulation
source(paste0(dir_code, 'load_parms.R'))

ave.frac.masking = 0;  cut.mask.start = 0


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


# get the model inputs from model_inputs
with(model_inputs,{
  args_names = names(model_inputs)
  for(arg in args_names){
    arg.t = get(arg)
    # assign(paste0(arg,'.t'), arg)
    assign(arg, arg.t)
  }
})
# get the model settings from model_settings
with(model_settings,{
  args_names = names(model_settings)
  for(arg in args_names){
    arg.t = get(arg)
    assign(arg, arg.t)
  }
})


# set parms accordingly
parms.t = list(num_ens = num_ens,
               num_pop = num_pop, # population size
               year.cur = tab_dates[1]$year, # year simulation starts
               UseContactDuration = UseContactDuration, #  whether to use contact duration directly or use number of contacts
               frac.ini.inf = frac.ini.inf0, # initial fraction of population being infectious 
               frac.masking = c(0, rep(ave.frac.masking, nrow(age.grps)-1)), # initial fraction of population adopt masking
               diff.ctc.bymask = diff.ctc.bymask,
               ave.p.redn.mask = ave.p.redn.mask.byage$value, # .7, # average effectiveness
               sd.p.redn.mask = sd.p.redn.mask.byage$value, # .25, 
               age.grps = age.grps,
               ave.num.nonHHctc = ave.num.nonHHctc.byage$value, # average number of non hh contact for diff age groups
               sd.num.nonHHctc = sd.num.nonHHctc.byage$value,  # std
               ave.num.HHctc = ave.num.HHctc.byage$value, # average number of hh contact (exclude self); # average HH size in NYC is 2.55
               sd.num.HHctc = sd.num.HHctc.byage$value,
               birth.rate =birth.rate)  # std

ndays = tab_dates %>% nrow()


# initialize the system for each ens member and combine
tmp = fn_initalizeABM_ens(parms.t)
rec = tmp$rec %>% data.table()
n.byage_ens.t = tmp$n.byage_ens

R0ini.t = tmp$R0ini %>% round(2)
R0ini.no.mask.t = tmp$R0ini.no.mask  %>% round(2)

summary(R0ini.no.mask.t)
quantile(R0ini.no.mask.t, probs = c(.5, .025, .975))  # main: 2.63  2.59  2.67 

# for wt_nonphys.t = 1: 3.17, 3.12, 3.22, mean = 3.17
# for wt_nonphys.t = .5: 2.63, 2.59, 2.66, mean = 2.63
# for wt_nonphys.t = .1: 2.54, 2.5, 2.58, mean = 2.54

